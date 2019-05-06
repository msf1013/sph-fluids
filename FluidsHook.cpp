#include "FluidsHook.h"
#include "RigidBodyTemplate.h"
#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include "igl/opengl/glfw/imgui/ImGuiHelpers.h"
#include "CollisionDetection.h"
#include "Particle.h"
#include <Eigen/Geometry>

using namespace Eigen;

#define PI 3.14159265358979323846

FluidsHook::FluidsHook() : PhysicsHook(), sceneFile_("box.scn")
{
    birdTemplate_ = NULL;
    launch_ = false;
}

void FluidsHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Scene", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputText("Filename", sceneFile_);
    }
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputFloat("Timestep", &params_.timeStep, 0, 0, 3);
        ImGui::DragFloat("Newton Tolerance", &params_.NewtonTolerance, 0.01, 1e-16, 1e-1, "%.3e", 10);
        ImGui::InputInt("Newton Max Iters", &params_.NewtonMaxIters);
    }
    if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Gravity Enabled", &params_.gravityEnabled);
        ImGui::InputFloat("Gravity G", &params_.gravityG, 0, 0, 3);
        ImGui::Checkbox("Penalty Forces Enabled", &params_.penaltyEnabled);
        ImGui::InputFloat("Penalty Stiffness", &params_.penaltyStiffness, 0, 0, 3);
        ImGui::Checkbox("Impulses Enabled", &params_.impulsesEnabled);
        ImGui::InputFloat("CoR", &params_.CoR, 0, 0, 3);
    }    
}

void FluidsHook::updateRenderGeometry()
{
    int totverts = 0;
    int totfaces = 0;

    totverts = sphereTemplate_->getVerts().rows() * particles_.size();
    totfaces = sphereTemplate_->getFaces().rows() * particles_.size();

    renderQ.resize(totverts, 3);
    renderF.resize(totfaces, 3);
    int voffset = 0;
    int foffset = 0;

    for (int x = 0; x < particles_.size(); x ++)
    {
        int nverts = sphereTemplate_->getVerts().rows();
        for (int i = 0; i < nverts; i++)
            renderQ.row(voffset + i) = (particles_[x]->position + sphereTemplate_->getVerts().row(i).transpose()).transpose();
        int nfaces = sphereTemplate_->getFaces().rows();
        for (int i = 0; i < nfaces; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                renderF(foffset + i, j) = sphereTemplate_->getFaces()(i, j) + voffset;
            }
        }
        voffset += nverts;
        foffset += nfaces;
    }
}


void FluidsHook::initSimulation()
{
    time_ = 0;    
    loadScene();
    updateRenderGeometry();
}

void FluidsHook::tick()
{    
    launchMutex_.lock();

    double launchVel = 100.0;

    if (launch_)
    {
        Eigen::Vector3d cvel(0, 0, 0);
        Eigen::Vector3d w(0, 0, 0);
        RigidBodyInstance *rbi = new RigidBodyInstance(*birdTemplate_, launchPos_, cvel, launchVel * launchDir_, w, 1.0);
        bodies_.push_back(rbi);
    
        launch_ = false;
    }

    launchMutex_.unlock();
}

void FluidsHook::computeForces(VectorXd &Fc, VectorXd &Ftheta)
{
    Fc.resize(3*bodies_.size());
    Ftheta.resize(3*bodies_.size());
    Fc.setZero();
    Ftheta.setZero();    

    if(params_.gravityEnabled)
    {
        for(int i=0; i<bodies_.size(); i++)
        {
            double m = bodies_[i]->density * bodies_[i]->getTemplate().getVolume();
            Fc[3 * i + 1] -= params_.gravityG*m;
        }
    }
}

bool FluidsHook::mouseClicked(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d dir, int button)
{
    if (button != 2)
        return false;

    launchMutex_.lock();
    launch_ = true;
    Eigen::Matrix4f view = viewer.core.view;
    Eigen::Vector4f eye = view.inverse() * Eigen::Vector4f(0, 0, 0, 1.0f);
    for (int i = 0; i < 3; i++)
    {

        launchPos_[i] = eye[i] + dir[i];
        launchDir_[i] = dir[i];
    }

    launchMutex_.unlock();
    return true;
}

bool FluidsHook::simulateOneStep()
{   
    time_ += params_.timeStep;
    int nbodies = (int)bodies_.size();

    std::vector<Vector3d> oldthetas;
    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {
        RigidBodyInstance &body = *bodies_[bodyidx];
        body.c += params_.timeStep*body.cvel;
        Matrix3d Rhw = VectorMath::rotationMatrix(params_.timeStep*body.w);
        Matrix3d Rtheta = VectorMath::rotationMatrix(body.theta);

        Vector3d oldtheta = body.theta;
        body.theta = VectorMath::axisAngle(Rtheta*Rhw);  
        if (body.theta.dot(oldtheta) < 0 && oldtheta.norm() > M_PI/2.0)
        {
            double oldnorm = oldtheta.norm();
            oldtheta = (oldnorm - 2.0*M_PI)*oldtheta/oldnorm;
        }
        oldthetas.push_back(oldtheta);
    }

    std::set<Collision> collisions;
    collisionDetection(bodies_, collisions);

    Eigen::VectorXd cForce(3 * nbodies);
    Eigen::VectorXd thetaForce(3 * nbodies);
    computeForces(cForce, thetaForce);
    
    // TODO compute and add penalty forces
    // TODO apply collision impulses
    
    for(int bodyidx=0; bodyidx < (int)bodies_.size(); bodyidx++)
    {        
        RigidBodyInstance &body = *bodies_[bodyidx];
        Matrix3d Mi = body.getTemplate().getInertiaTensor();

        body.cvel += params_.timeStep*cForce.segment<3>(3 * bodyidx) / body.density / body.getTemplate().getVolume();

        Vector3d newwguess(body.w);

        int iter = 0;
        for (iter = 0; iter < params_.NewtonMaxIters; iter++)
        {
            Vector3d term1 = (-VectorMath::TMatrix(-params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * newwguess;
            Vector3d term2 = (VectorMath::TMatrix(params_.timeStep*body.w).inverse()*VectorMath::TMatrix(oldthetas[bodyidx])).transpose() * Mi * body.density * body.w;
            Vector3d term3 = -params_.timeStep * thetaForce.segment<3>(3 * bodyidx);
            Vector3d fval = term1 + term2 + term3;
            if (fval.norm() / body.density / Mi.trace() <= params_.NewtonTolerance)
                break;

            Matrix3d Df = (-VectorMath::TMatrix(-params_.timeStep*newwguess).inverse() * VectorMath::TMatrix(body.theta)).transpose() * Mi * body.density;

            Vector3d deltaw = Df.inverse() * (-fval);
            newwguess += deltaw;
        }
        std::cout << "Converged in " << iter << " Newton iterations" << std::endl;
        body.w = newwguess;

    }

    return false;
}

void FluidsHook::loadScene()
{
    for (Particle *p : particles_)
        delete p;
    particles_.clear();

    double width = 2.0, height = 1.0, depth = 1.0;
    int num_w = 3, num_h = 3, num_d = 3;

    for (int i = 0; i < num_w; i ++) {
        double x = -width/2.0 + i * width/(num_w - 1.0);
        for (int j = 0; j < num_h; j ++) {
            double y = -height/2.0 + j * height/(num_h - 1.0) + 1.0;
            for (int k = 0; k < num_d; k ++) {
                double z = -depth/2.0 + k * depth/(num_d - 1.0);
                particles_.push_back(new Particle(Eigen::Vector3d(x, y, z), Eigen::Vector3d(0.0, 0.0, 0.0), 1.0));
            }
        }
    }

    tankV.resize(8,3);
    tankV << -2.0, -1.0, -1.0,
             -2.0, -1.0,  1.0,
              2.0, -1.0, -1.0,
              2.0, -1.0,  1.0,
             -2.0,  1.0, -1.0,
             -2.0,  1.0,  1.0,
              2.0,  1.0, -1.0,
              2.0,  1.0,  1.0;

    tankE.resize(12,2);
    tankE <<
          0, 1,
          0, 2,
          1, 3,
          2, 3,
          1, 5,
          3, 7,
          2, 6,
          0, 4,
          4, 6,
          4, 5,
          5, 7,
          6, 7;

    sphereTemplate_ = new RigidBodyTemplate("../meshes/sphere.obj", 0.1);
}

double FluidsHook::viscosityKernelLaplacian(double distance, double h) {
    return 45.0 / (PI * h * h * h * h * h * h) * (h - distance);
}

double FluidsHook::pressureKernelGradient(double distance, double h) {
    return -90.0 * (h - distance) * (h - distance) * (h - distance) / (PI *h*h*h*h*h*h*h);
}
