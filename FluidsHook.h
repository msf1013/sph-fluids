#include "PhysicsHook.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/ViewerData.h>
#include <deque>
#include <Eigen/Sparse>
#include <Eigen/StdVector>
#include "SimParameters.h"
#include <set>
#include "Particle.h"

class RigidBodyTemplate;
class RigidBodyInstance;

using Eigen::Vector3d;
using std::vector;

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

        viewer.data().add_points(particlesPos,Eigen::RowVector3d(35.0/255.0,137.0/255.0,218/255.0));

        // Render outer tank edges
        for (unsigned i=0;i<oTankE.rows(); ++i)
            viewer.data().add_edges
            (
              oTankV.row(oTankE(i,0)),
              oTankV.row(oTankE(i,1)),
              Eigen::RowVector3d(1,1,1)
        );
        
        // Render inner tank edges
        for (unsigned i=0;i<iTankE.rows(); ++i)
            viewer.data().add_edges
            (
              iTankV.row(iTankE(i,0)),
              iTankV.row(iTankE(i,1)),
              Eigen::RowVector3d(1,1,1)
        );
    }

private:
    void loadScene();

    std::mutex applyForceMutex_;
    Eigen::Vector3d launchPos_;
    Eigen::Vector3d launchDir_;

    double time_;
    SimParameters params_;
    std::string sceneFile_;

    Eigen::MatrixXd renderQ;
    Eigen::MatrixXi renderF;

    std::vector<Particle *> particles_;

    // Dimensions of outer tank.
    double ot_width=4.5, ot_height=3.2, ot_depth=1.5;
    Eigen::MatrixXd oTankV;
    Eigen::MatrixXi oTankE;

    // Dimensions of inner tank.
    double it_width=0.6, it_height=0.15, it_depth=1.5;
    Eigen::MatrixXd iTankV;
    Eigen::MatrixXi iTankE;

    bool pressed = false;
    bool applyForce = false;
    bool applyExternalForce = false;
    Eigen::Vector3d startV;
    Eigen::Vector3d endV;
    Eigen::Vector3d eye_;

    Eigen::Vector3d accAlongPlane(Eigen::Vector3d startV, Eigen::Vector3d endV);
    double pointToPlaneDistance(Eigen::Vector3d p, Eigen::Vector3d v1, Eigen::Vector3d v2);


    void computeAcc(vector<Vector3d> &Acc);
    void computeGravityAcc(vector<Vector3d> &Acc);
    void computeFloorWallAcc(vector<Vector3d> &Acc);
    void computePressureAcc(const vector<double> &Density, vector<Vector3d> &Acc);
    void computeViscosityAcc(const vector<double> &Density, vector<Vector3d> &Acc);
    void computeSurfaceTensionAcc(const vector<double> &Density, vector<Vector3d> &Acc);
    void computeExternalAcc(vector<Vector3d> &Acc);

    void computeDensity(vector<double> &Density);

    double kernelPoly6(double r, double h);
    Vector3d kernelPoly6Gradient(Vector3d R, double h);
    double kernelPoly6Laplacian(double r, double h);
    Vector3d kernelSpikyGradient(Vector3d R, double h);
    double kernelViscosityLaplacian(double r, double h);

    struct coord { 
        int x, y, z; 

        bool operator<(const coord &o) const {
            return x < o.x || (x == o.x && y < o.y) || (x == o.x && y == o.y && z < o.z);
        }
    };

    std::map<coord, vector<int>> grid;
    void computeGrid();
    coord point_to_coord(Vector3d);
    vector<int> grid_neighbors(Vector3d);
};