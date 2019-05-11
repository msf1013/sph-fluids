#include "FluidsHook.h"
#include "RigidBodyTemplate.h"
#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include "igl/opengl/glfw/imgui/ImGuiHelpers.h"
#include "CollisionDetection.h"
#include "Particle.h"
#include <Eigen/Geometry>
#include <stdlib.h>

using namespace Eigen;
using Eigen::Vector3d;
using std::vector;

#define PI 3.14159265358979323846

FluidsHook::FluidsHook() : PhysicsHook(), sceneFile_("box.scn")
{
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

}


void FluidsHook::initSimulation()
{
    time_ = 0;    
    loadScene();
}

void FluidsHook::tick()
{    
    applyForceMutex_.lock();

    if (applyForce)
    {
        std::cout << "Start vector\n";
        std::cout << startV;
        std::cout << "End vector\n";
        std::cout << endV;

        applyForce = false;
    }

    applyForceMutex_.unlock();
}

void FluidsHook::computeAcc(vector<Vector3d> &Acc)
{
    for (auto &acc : Acc) acc.setZero();

    computeGravityAcc(Acc);
    computeFloorWallAcc(Acc);
    computePressureAcc(Acc);
    computeViscosityAcc(Acc);
    computeSurfaceTensionAcc(Acc);

}


void FluidsHook::computeGravityAcc(vector<Vector3d> &Acc) {
    // TODO. Where is rho?
    for(int i = 0; i < particles_.size(); i ++)
    {
        Acc[i][1] -= params_.gravityG;
    }
}

void FluidsHook::computeFloorWallAcc(vector<Vector3d> &Acc) {

    // Floor force
    // TODO. Should this be in params_?
    double basestiffness = 10000;
    double basedrag = 1000.0;

    for(int i = 0; i < particles_.size(); i ++)
    {
        if(particles_[i]->position[1] < -t_height/2.0)
        {
            double vel = (particles_[i]->position[1] - particles_[i]->prev_position[1])/params_.timeStep;
            double dist = -t_height/2.0 - particles_[i]->position[1];

            Acc[i][1] += basestiffness*dist - basedrag*dist*vel;
        }
    }

    // Wall forces
    for(int i = 0; i < particles_.size(); i ++)
    {
        if(particles_[i]->position[0] < -t_width/2.0)
        {
            double vel = (particles_[i]->position[0] - particles_[i]->prev_position[0])/params_.timeStep;
            double dist = -t_width/2.0 - particles_[i]->position[0];

            Acc[i][0] += basestiffness*dist - basedrag*dist*vel;
        }
        if(particles_[i]->position[0] > t_width/2.0)
        {
            double vel = (particles_[i]->position[0] - particles_[i]->prev_position[0])/params_.timeStep;
            double dist = particles_[i]->position[0] - t_width/2.0;

            Acc[i][0] += basedrag*dist*vel - basestiffness*dist;
        }
        if(particles_[i]->position[2] < -t_depth/2.0)
        {
            double vel = (particles_[i]->position[2] - particles_[i]->prev_position[2])/params_.timeStep;
            double dist = -t_depth/2.0 - particles_[i]->position[2];

            Acc[i][2] += basestiffness*dist - basedrag*dist*vel;
        }
        if(particles_[i]->position[2] > t_depth/2.0)
        {
            double vel = (particles_[i]->position[2] - particles_[i]->prev_position[2])/params_.timeStep;
            double dist = particles_[i]->position[2] - t_depth/2.0;

            Acc[i][2] +=  basedrag*dist*vel - basestiffness*dist;
        }
    }
}

void FluidsHook::computePressureAcc(vector<Vector3d> &Acc) {
    
}

void FluidsHook::computeViscosityAcc(vector<Vector3d> &Acc) {

}

void FluidsHook::computeSurfaceTensionAcc(vector<Vector3d> &Acc) {

}

bool FluidsHook::mouseClicked(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d dir)
{
    if (pressed) true;

    std::cout << "CLICKED\n";

    pressed = true;
    startV = dir;

    return true;
}

bool FluidsHook::mouseReleased(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d dir)
{
    pressed = false;

    applyForceMutex_.lock();

    std::cout << "RELEASED\n";
    applyForce = true;
    endV = dir;

    applyForceMutex_.unlock();
    
    return true;
}

bool FluidsHook::simulateOneStep()
{   
    // Implicit Euler
    time_ += params_.timeStep;

    // Calculate position_(i+1)
    for (int i = 0; i < particles_.size(); i ++) {
        particles_[i]->prev_position = particles_[i]->position;
        particles_[i]->position += params_.timeStep*particles_[i]->velocity;
    }

    // TODO. This is not useful.. making 3*size eigen, rather use vector of 3D eigen.
    // Calculate velocity_(i+1) with position_(i+1)
    vector<Vector3d> Acc(particles_.size());
    computeAcc(Acc);

    for (int i = 0; i < particles_.size(); i ++) {
        particles_[i]->velocity += params_.timeStep*Acc[i];
    }

    return false;
}

void FluidsHook::loadScene()
{
    for (Particle *p : particles_)
        delete p;
    particles_.clear();

    double width = 2.0, height = 1.0, depth = 1.0;
    int num_w = 10, num_h = 10, num_d = 10;

    for (int i = 0; i < num_w; i ++) {
        double x = -width/2.0 + i * width/(num_w - 1.0);
        for (int j = 0; j < num_h; j ++) {
            double y = -height/2.0 + j * height/(num_h - 1.0) + 1.0;
            for (int k = 0; k < num_d; k ++) {
                double z = -depth/2.0 + k * depth/(num_d - 1.0);
                particles_.push_back(new Particle(Eigen::Vector3d(x, y, z), Eigen::Vector3d(((double) rand() / (RAND_MAX)), ((double) rand() / (RAND_MAX)), ((double) rand() / (RAND_MAX))) * 2.0 - Eigen::Vector3d(1,1,1), 1.0));
            }
        }
    }

    tankV.resize(8,3);
    tankV << -t_width/2.0, -t_height/2.0, -t_depth/2.0,
             -t_width/2.0, -t_height/2.0,  t_depth/2.0,
              t_width/2.0, -t_height/2.0, -t_depth/2.0,
              t_width/2.0, -t_height/2.0,  t_depth/2.0,
             -t_width/2.0,  t_height/2.0, -t_depth/2.0,
             -t_width/2.0,  t_height/2.0,  t_depth/2.0,
              t_width/2.0,  t_height/2.0, -t_depth/2.0,
              t_width/2.0,  t_height/2.0,  t_depth/2.0;

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
}

double FluidsHook::viscosityKernelLaplacian(double distance, double h) {
    return 45.0 / (PI * h * h * h * h * h * h) * (h - distance);
}

double FluidsHook::pressureKernelGradient(double distance, double h) {
    return -90.0 * (h - distance) * (h - distance) * (h - distance) / (PI *h*h*h*h*h*h*h);
}
