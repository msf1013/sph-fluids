#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

struct SimParameters
{
    SimParameters()
    {
        numParticles_index = 6;
        particleMass = 1;
        timeStep = 0.0015;

        pressureEnabled = true;
        viscosityEnabled = true;
        surfaceTensionEnabled = true;
        gravityEnabled = true;

        gravityG = 9.8;
        smoothingLength = 0.2;
        restDensity = 1.0;
        viscosityCoefficient = 0.1;
        tensionCoefficient = 1.0;
        epsColorNormal = 0.001;

    }

    int numParticles_index;
    const char* numParticles[9] = {"8", "27", "64", "125", "216", "343", "512", "729", "1000"};
    float particleMass;
    float timeStep = 0.0015;

    bool pressureEnabled;
    bool viscosityEnabled;
    bool surfaceTensionEnabled;
    bool gravityEnabled;

    float gravityG;
    float smoothingLength;
    float restDensity;
    float viscosityCoefficient;
    float tensionCoefficient;
    float epsColorNormal;
};

#endif