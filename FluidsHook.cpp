#include "FluidsHook.h"
#include "igl/opengl/glfw/imgui/ImGuiHelpers.h"
#include "Particle.h"
#include <Eigen/Geometry>
#include <stdlib.h>
#include <cmath>

using namespace Eigen;

#define PI 3.14159265358979323846
#define GAS_CONSTANT 8.3144598

FluidsHook::FluidsHook() : PhysicsHook(), sceneFile_("box.scn")
{
}

void FluidsHook::drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
    if (ImGui::CollapsingHeader("Simulation Options", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Combo("Number of particles", &params_.numParticles_index, params_.numParticles, sizeof(params_.numParticles) / sizeof(params_.numParticles[0]));
        ImGui::InputFloat("Particle mass", &params_.particleMass, 0, 0, 4);
        ImGui::InputFloat("Timestep", &params_.timeStep, 0, 0, 4);
    }
    if (ImGui::CollapsingHeader("Forces", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::Checkbox("Pressure Enabled", &params_.pressureEnabled);
        ImGui::Checkbox("Viscosity Enabled", &params_.viscosityEnabled);
        ImGui::Checkbox("Surface Tension Enabled", &params_.surfaceTensionEnabled);
        ImGui::Checkbox("Gravity Enabled", &params_.gravityEnabled);
    }
    if (ImGui::CollapsingHeader("Forces Parameters", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::InputFloat("Gravity G", &params_.gravityG, 0, 0, 4);
        ImGui::InputFloat("Smoothing Length", &params_.smoothingLength, 0, 0, 4);
        ImGui::InputFloat("Rest Density", &params_.restDensity, 0, 0, 4);
        ImGui::InputFloat("Viscosity Coefficient", &params_.viscosityCoefficient, 0, 0, 4);
        ImGui::InputFloat("Tension Coefficient", &params_.tensionCoefficient, 0, 0, 4);
        ImGui::InputFloat("Epsilon Color Normal", &params_.epsColorNormal, 0, 0, 4);
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
        applyExternalForce = true;
        applyForce = false;
    }

    applyForceMutex_.unlock();
}

void FluidsHook::computeAcc(vector<Vector3d> &Acc)
{
    for (auto &acc : Acc) acc.setZero();

    // Recompute Grid.
    grid.clear();
    computeGrid();

    if (applyExternalForce) {
        computeExternalAcc(Acc);
    }
    if (params_.gravityEnabled) {
        computeGravityAcc(Acc);
    }
    computeFloorWallAcc(Acc);

    vector<double> Density(particles_.size());
    computeDensity(Density);

    if (params_.pressureEnabled) {
        computePressureAcc(Density, Acc);
    }
    if (params_.viscosityEnabled) {
        computeViscosityAcc(Density, Acc);
    }
    if (params_.surfaceTensionEnabled) {
        computeSurfaceTensionAcc(Density, Acc);
    }
}

void FluidsHook::computeGrid() {
    for (int i = 0; i < particles_.size(); ++i) {
        grid[point_to_coord(particles_[i]->position)].push_back(i);
    }
}

vector<int> FluidsHook::grid_neighbors(Vector3d R) {
    vector<int> ret = grid[point_to_coord(R)];
    double step[2] = {-params_.smoothingLength, params_.smoothingLength};

    for (double x : step) {
        for (double y : step) {
            for (double z : step) {
                Vector3d dR(x,y,z);
                coord c = point_to_coord(R+dR);
                ret.insert(ret.end(), grid[c].begin(), grid[c].end());
            }
        }
    }
    return ret;
}

FluidsHook::coord FluidsHook::point_to_coord(Vector3d R) {
    assert (params_.smoothingLength >= 0);
    coord a;
    a.x = R[0] / params_.smoothingLength;
    a.y = R[1] / params_.smoothingLength;
    a.z = R[2] / params_.smoothingLength;
    return a;
}

void FluidsHook::computeDensity(vector<double> &Density) {
    Density.resize(particles_.size());
    for (int i = 0; i < particles_.size(); ++i) {
        double density = 0;
        for (int j : grid_neighbors(particles_[i]->position)) {
            density += params_.particleMass * kernelPoly6( 
                (particles_[i]->position - particles_[j]->position).norm(), params_.smoothingLength );
        }
        Density[i] = density;
    }
}

void FluidsHook::computeExternalAcc(vector<Vector3d> &Acc)
{
    Eigen::Vector3d acc = accAlongPlane(startV, endV) * 3200.0;
    for (int i = 0; i < particles_.size(); i ++) {
        double distance = pointToPlaneDistance(particles_[i]->position, startV, endV); 

        if (distance < 0.1) {
            Acc[i] += acc * (1.0 - distance) * (1.0 - distance);
        }
    }
    applyExternalForce = false;
}

Eigen::Vector3d FluidsHook::accAlongPlane(Eigen::Vector3d startV, Eigen::Vector3d endV)
{
    Eigen::Vector3d startU = startV / startV.norm();
    Eigen::Vector3d endU = endV / endV.norm();

    Eigen::Vector3d acc = endU - startU;

    return acc / acc.norm();
}

double FluidsHook::pointToPlaneDistance(Eigen::Vector3d p, Eigen::Vector3d v1, Eigen::Vector3d v2) {
    Eigen::Vector3d n = v1.cross(v2);

    double a = n[0],
           b = n[1],
           c = n[2],
           d = -n[0]*eye_[0] -n[1]*eye_[1] -n[2]*eye_[2];

    return abs(a*p[0] + b*p[1] + c*p[2] + d) / sqrt(a*a + b*b + c*c);
}

void FluidsHook::computeGravityAcc(vector<Vector3d> &Acc) {
    // TODO. Where is rho?
    for(int i = 0; i < particles_.size(); i ++)
    {
        Acc[i][1] -= params_.gravityG;
    }
}

void FluidsHook::computeFloorWallAcc(vector<Vector3d> &Acc) {
    // Compute wall/floor bouncing forces
    // TODO. Should this be in params_?
    double basestiffness = 10000;
    double basedrag = 1000.0;
    double basedragLat = 1000.0;

    // Floor force
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
        // Left wall force
        if(particles_[i]->position[0] < -t_width/2.0)
        {
            double vel = (particles_[i]->position[0] - particles_[i]->prev_position[0])/params_.timeStep;
            double dist = -t_width/2.0 - particles_[i]->position[0];

            Acc[i][0] += basestiffness*dist - basedragLat*dist*vel;
        }
        // Right wall force
        if(particles_[i]->position[0] > t_width/2.0)
        {
            double vel = (particles_[i]->position[0] - particles_[i]->prev_position[0])/params_.timeStep;
            double dist = particles_[i]->position[0] - t_width/2.0;

            Acc[i][0] += basedragLat*dist*vel - basestiffness*dist;
        }
        // Back wall force
        if(particles_[i]->position[2] < -t_depth/2.0)
        {
            double vel = (particles_[i]->position[2] - particles_[i]->prev_position[2])/params_.timeStep;
            double dist = -t_depth/2.0 - particles_[i]->position[2];

            Acc[i][2] += basestiffness*dist - basedragLat*dist*vel;
        }
        // Front wall force
        if(particles_[i]->position[2] > t_depth/2.0)
        {
            double vel = (particles_[i]->position[2] - particles_[i]->prev_position[2])/params_.timeStep;
            double dist = particles_[i]->position[2] - t_depth/2.0;

            Acc[i][2] +=  basedragLat*dist*vel - basestiffness*dist;
        }
        // Roof/top wall force
        if(particles_[i]->position[1] > t_height/2.0)
        {
            double vel = (particles_[i]->position[1] - particles_[i]->prev_position[1])/params_.timeStep;
            double dist = particles_[i]->position[1] - t_height/2.0;

            Acc[i][1] += basedragLat*dist*vel - basestiffness*dist;
        }
    }
}

void FluidsHook::computePressureAcc(const vector<double> &Density, vector<Vector3d> &Acc) {

    // Computing Pressure.
    vector<double> Pressure(particles_.size());
    for (int i = 0; i < particles_.size(); ++i) {
        Pressure[i] = GAS_CONSTANT * (Density[i] - params_.restDensity);
    }

    for (int i = 0; i < particles_.size(); ++i) {
        Vector3d force = Vector3d::Zero();
        for (int j : grid_neighbors(particles_[i]->position)) {
            force += -params_.particleMass * ( (Pressure[i] + Pressure[j]) / (2 * Density[j]) ) *
                kernelSpikyGradient(particles_[i]->position - particles_[j]->position, params_.smoothingLength);
        }
        Acc[i] += force / Density[i];
    }
}

void FluidsHook::computeViscosityAcc(const vector<double> &Density, vector<Vector3d> &Acc) {

    for (int i = 0; i < particles_.size(); ++i) {
        Vector3d force = Vector3d::Zero();
        for (int j : grid_neighbors(particles_[i]->position)) {
            force += params_.particleMass * ( (particles_[j]->velocity - particles_[i]->velocity) /
                Density[j] ) * kernelViscosityLaplacian( (particles_[i]->position 
                    - particles_[j]->position).norm(), params_.smoothingLength );
        }
        Acc[i] += params_.viscosityCoefficient * force / Density[i];
    }
}

void FluidsHook::computeSurfaceTensionAcc(const vector<double> &Density, vector<Vector3d> &Acc) {

    vector<Vector3d> GradientColorField(particles_.size());

    for (int i = 0; i < particles_.size(); ++i) {
        Vector3d grad = Vector3d::Zero();
        for (int j : grid_neighbors(particles_[i]->position)) {
            grad += (params_.particleMass / Density[j]) * kernelPoly6Gradient( particles_[i]->position 
                    - particles_[j]->position, params_.smoothingLength );
        }
        GradientColorField[i] = grad;
    }

    for (int i = 0; i < particles_.size(); ++i) {
        double force = 0;
        for (int j : grid_neighbors(particles_[i]->position)) {
            if (GradientColorField[i].norm() > params_.epsColorNormal) break;
            force += kernelPoly6Laplacian( (particles_[i]->position 
                    - particles_[j]->position).norm(), params_.smoothingLength );
        }
        Acc[i] += -params_.tensionCoefficient * force * GradientColorField[i] / (
            Density[i] / GradientColorField[i].norm());
    }
}


double FluidsHook::kernelPoly6(double r, double h) {
    assert (r >= 0);
    if (r <= h) {
        return ( 315.0 / (64 * PI * pow(h,9)) ) * pow(pow(h,2) - pow(r,2), 3);
    }
    return 0;
}

Vector3d FluidsHook::kernelPoly6Gradient(Vector3d R, double h) {
    
    if (R.norm() <= h) {
        return ( 315.0 / (64 * PI * pow(h,9)) ) * pow(pow(h,2) - pow(R.norm(),2), 2) *
            3 * ( -2 * R);
    }
    return Vector3d::Zero();
}

double FluidsHook::kernelPoly6Laplacian(double r, double h) {
    assert (r >= 0);
    if (r <= h) {
        return 6 * ( 315.0 / (64 * PI * pow(h,9)) ) * (pow(h,2) - pow(r,2)) *
            ( 7 * pow(r,2) - 3 * pow(h,2) );
    }
    return 0;
}

Vector3d FluidsHook::kernelSpikyGradient(Vector3d R, double h) {
    if (R.norm() == 0) return Vector3d::Zero();

    if (R.norm() <= h) {
        return ( -45.0 / (PI * pow(h,6)) ) * pow(h - R.norm(), 2) * R / R.norm();
    }

    return Vector3d::Zero();
}

double FluidsHook::kernelViscosityLaplacian(double r, double h) {
    assert (r >= 0);
    if (r <= h) {
        return ( 45.0 / (PI * pow(h,6)) ) * (h - r);
    }
    return 0;
}

bool FluidsHook::mouseClicked(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d dir)
{
    if (pressed) return true;

    pressed = true;

    // Calculate startV
    Eigen::Matrix4f view = viewer.core.view;
    Eigen::Vector4f eye = view.inverse() * Eigen::Vector4f(0, 0, 0, 1.0f);
    for (int i = 0; i < 3; i++)
    {

        eye_[i] = eye[i] + dir[i];
        startV[i] = dir[i];
    }

    return true;
}

bool FluidsHook::mouseReleased(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d dir)
{
    pressed = false;

    applyForceMutex_.lock();

    applyForce = true;
    
    // Calculate endV
    Eigen::Matrix4f view = viewer.core.view;
    Eigen::Vector4f eye = view.inverse() * Eigen::Vector4f(0, 0, 0, 1.0f);
    for (int i = 0; i < 3; i++)
    {

        endV[i] = dir[i];
    }

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

    // Calculate velocity_(i+1) with position_(i+1)
    vector<Vector3d> Acc(particles_.size());
    computeAcc(Acc);

    for (int i = 0; i < particles_.size(); i ++) {
        particles_[i]->velocity += params_.timeStep*Acc[i];
        // Bound the amount of kinetic energy of particles.
        // Otherwise, external forces (mouse drag) can indefinitely add energy to the system.
        if (particles_[i]->velocity.norm() > 500.0) {
            particles_[i]->velocity = particles_[i]->velocity / particles_[i]->velocity.norm() * 500.0;
        }
    }

    return false;
}

void FluidsHook::loadScene()
{
    for (Particle *p : particles_)
        delete p;
    particles_.clear();

    // Dimensions of particles 'box'.
    double width = 1.3, height = 1.3, depth = 1.3;

    // Number of particles per dimension of the box.
    int num_p = atoi(params_.numParticles[params_.numParticles_index]);
    int num_w = cbrt(num_p), num_h = cbrt(num_p), num_d = cbrt(num_p);

    for (int i = 0; i < num_w; i ++) {
        double x = -width/2.0 + i * width/(num_w - 1.0) - 0.4;
        for (int j = 0; j < num_h; j ++) {
            double y = -height/2.0 + j * height/(num_h - 1.0);
            for (int k = 0; k < num_d; k ++) {
                double z = -depth/2.0 + k * depth/(num_d - 1.0);
                particles_.push_back(new Particle(Eigen::Vector3d(x, y, z),
                                                  Eigen::Vector3d(((double) rand() / (RAND_MAX)),
                                                                  ((double) rand() / (RAND_MAX)),
                                                                  ((double) rand() / (RAND_MAX))) * 2.0 - Eigen::Vector3d(1,1,1)));
            }
        }
    }

    // Define position of vertices of tank.
    tankV.resize(8,3);
    tankV << -t_width/2.0-0.08, -t_height/2.0-0.08, -t_depth/2.0-0.08,
             -t_width/2.0-0.08, -t_height/2.0-0.08,  t_depth/2.0+0.08,
              t_width/2.0+0.08, -t_height/2.0-0.08, -t_depth/2.0-0.08,
              t_width/2.0+0.08, -t_height/2.0-0.08,  t_depth/2.0+0.08,
             -t_width/2.0-0.08,  t_height/2.0+0.08, -t_depth/2.0-0.08,
             -t_width/2.0-0.08,  t_height/2.0+0.08,  t_depth/2.0+0.08,
              t_width/2.0+0.08,  t_height/2.0+0.08, -t_depth/2.0-0.08,
              t_width/2.0+0.08,  t_height/2.0+0.08,  t_depth/2.0+0.08;

    // Define edges of tank.
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
