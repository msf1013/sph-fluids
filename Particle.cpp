#include "RigidBodyInstance.h"
#include "VectorMath.h"
#include "RigidBodyTemplate.h"
#include <Eigen/Geometry>
#include <iostream>
#include "CollisionDetection.h"
#include "Particle.h"

using namespace Eigen;
using namespace std;

Particle::Particle(Eigen::Vector3d position, Eigen::Vector3d velocity, double mass)
    : position(position), velocity(velocity), mass(mass)
{
}

Particle::~Particle()
{    
}