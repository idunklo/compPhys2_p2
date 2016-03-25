#pragma once
#include "WaveFunction.h"

class TrialWaveFunction : public WaveFunction
{
  public:
    TrialWaveFunction		        (class System* system);
    double evaluate		          ();
    double computeQuantumForce  (int p);
};
