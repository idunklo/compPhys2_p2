#pragma once
#include "WaveFunction.h"

class TrialWaveFunction : public WaveFunction
{
  public:
    TrialWaveFunction		        (class System* system);
    double evaluate		          ();
    double computeGradient      (int p, int d);
    double computeQuantumForce  (int p, int d);
    void   computePsiBars       (double &psiBar_alpha, double &psiBar_beta);
};
