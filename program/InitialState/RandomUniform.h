#pragma once
#include "InitialState.h"

class RandomUniform : public InitialState
{

  public:
    RandomUniform ( System* system,
		    int nDimensions, 
		    int nParticles);

    void setupInitialState ();

  protected:
    double my_scaling = 0; 
    unsigned  seed;
    typedef   std::chrono::high_resolution_clock clock;
    clock::time_point my_start = clock::now();

};
