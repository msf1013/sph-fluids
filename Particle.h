#ifndef PARTICLE_H
#define PARTICLE_H

#include <Eigen/Core>
#include <list>
#include <vector>

class Particle
{
public:
    Particle(Eigen::Vector3d position, Eigen::Vector3d velocity);
    ~Particle();

    Eigen::Vector3d position;
    Eigen::Vector3d prev_position;
    Eigen::Vector3d velocity;
};

#endif // PARTICLE_H
