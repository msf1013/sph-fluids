#ifndef PARTICLE_H
#define PARTICLE_H

#include <Eigen/Core>
#include <list>
#include <vector>

class Particle
{
public:
    Particle(Eigen::Vector3d position, Eigen::Vector3d velocity, double mass);
    ~Particle();

    Eigen::Vector3d position;
    Eigen::Vector3d velocity;

    double mass;
};

#endif // PARTICLE_H
