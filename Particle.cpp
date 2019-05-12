#include <Eigen/Geometry>
#include <iostream>
#include "Particle.h"

using namespace Eigen;
using namespace std;

Particle::Particle(Eigen::Vector3d position, Eigen::Vector3d velocity)
    : position(position), velocity(velocity)
{
	prev_position = Eigen::Vector3d(0.0, 0.0, 0.0);
}

Particle::~Particle()
{    
}