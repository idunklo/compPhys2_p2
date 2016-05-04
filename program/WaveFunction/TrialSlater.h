#pragma once
#include "WaveFunction.h"

class TrialSlater: public WaveFunction
{
  public:
    TrialSlater                   (class System* system);
    double Phi                    (int p, int nx, int ny);
    double GradPhi                (int pos, int d);
    double GradJas                (int k,int d);
    double LapPhi                 (int pos, int orbital);
    double LapJas                 ();
    double computeJastrow         ();
    double computeQuantumForce    (int p, int d);
    void   computePsiBars         (double &psiBar_alpha, double &psiBar_beta);
};
