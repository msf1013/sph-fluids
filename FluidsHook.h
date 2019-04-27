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

    virtual bool mouseClicked(igl::opengl::glfw::Viewer &viewer, Eigen::Vector3d dir, int button);

    virtual void updateRenderGeometry();

    virtual void tick();

    virtual bool simulateOneStep();

    virtual void renderRenderGeometry(igl::opengl::glfw::Viewer &viewer)
    {
        viewer.data().clear();
        viewer.data().set_mesh(renderQ, renderF);
    }

private:
    void loadScene();
    void computeForces(Eigen::VectorXd &Fc, Eigen::VectorXd &Ftheta);    

    std::mutex launchMutex_;
    bool launch_;
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

    double viscosityKernelLaplacian(double distance, double h);
    double pressureKernelGradient(double distance, double h);
};