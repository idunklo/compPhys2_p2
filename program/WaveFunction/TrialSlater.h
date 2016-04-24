#pragma once
#include "WaveFunction.h"

class TrialSlater: public WaveFunction
{
  public:
    TrialSlater                   (class System* system);
    double evaluate               (int p, int nx, int ny);
    double computeGradient        (int row, int level);
    double computeLaplacian       (int row, int level);
    double computeQuantumForce    (int p, int d);
    void   computePsiBars         (double &psiBar_alpha, double &psiBar_beta);
};
