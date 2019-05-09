#include "PhysicsHook.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/ViewerData.h>
#include <deque>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include "SimParameters.h"
#include <set>
#include "CollisionDetection.h"
#include "Particle.h"

class RigidBodyTemplate;
class RigidBodyInstance;

class FluidsHook : public PhysicsHook
{
public:
    FluidsHook();    

    virtual void drawGUI(igl::opengl::glfw::imgui::ImGuiMenu &menu);

    virtual void initSimulation();

    virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d dir);

    virtual bool mouseReleased(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d dir);

    virtual void updateRenderGeometry();

    virtual void tick();

    virtual bool simulateOneStep();

    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
    {
        viewer.data().clear();
        //viewer.data().set_mesh(renderQ, renderF);

        //viewer.data().add_points(tankV,Eigen::RowVector3d(1.0, 0.0, 0.0));

        Eigen::MatrixXd particlesPos(particles_.size(), 3);

        for (int i = 0; i < particles_.size(); i ++) {
            particlesPos.row(i) = particles_[i]->position;
        }

        viewer.data().add_points(particlesPos,Eigen::RowVector3d(0.0, 0.0, 1.0));
        
        for (unsigned i=0;i<tankE.rows(); ++i)
            viewer.data().add_edges
            (
              tankV.row(tankE(i,0)),
              tankV.row(tankE(i,1)),
              Eigen::RowVector3d(1,0,0)
        );
    }

private:
    void loadScene();
    void computeForces(Eigen::VectorXd &F);    

    std::mutex applyForceMutex_;
    Eigen::Vector3d launchPos_;
    Eigen::Vector3d launchDir_;

    double time_;
    SimParameters params_;
    std::string sceneFile_;

    std::vector<RigidBodyTemplate *> templates_;
    std::vector<RigidBodyInstance *> bodies_;

    RigidBodyTemplate *birdTemplate_;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;

    std::vector<Particle *> particles_;
    Eigen::MatrixXd tankV;
    Eigen::MatrixXi tankE;

    double t_width=4.0, t_height=2.0, t_depth=2.0;

    bool pressed = false;
    bool applyForce = false;
    Eigen::Vector3d startV;
    Eigen::Vector3d endV;

    double viscosityKernelLaplacian(double distance, double h);
    double pressureKernelGradient(double distance, double h);
};