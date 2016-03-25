#pragma once
#include "Hamiltonian.h"

class HarmonicOscillator : public Hamiltonian
{
  public:
    
    HarmonicOscillator( System* system);

    double  computeLocalEnergy	();
};

