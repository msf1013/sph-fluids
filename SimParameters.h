#ifndef SIMPARAMETERS_H
#define SIMPARAMETERS_H

struct SimParameters
{
    SimParameters()
    {
        timeStep = 0.001;
        NewtonMaxIters = 50;
        NewtonTolerance = 1e-8;
        
        gravityEnabled = true;
        gravityG = 9.8;                
        penaltyEnabled = true;
        penaltyStiffness = 1.0;
        impulsesEnabled = true;
        CoR = 0.5;
    }

    float timeStep;
    float NewtonTolerance;
    int NewtonMaxIters;
    
    bool gravityEnabled;
    float gravityG;    
    bool penaltyEnabled;
    float penaltyStiffness;
    bool impulsesEnabled;
    float CoR;
};

#endif