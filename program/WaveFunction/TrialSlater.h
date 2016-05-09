#pragma once
#include "WaveFunction.h"

class TrialSlater: public WaveFunction
{
  public:
    TrialSlater                   (class System* system, double jf);
    double Phi                    (int p, int nx, int ny);
    double GradPhi                (int pos, int d);
    double GradJas                (int k,int d);
    double LapPhi                 (int pos, int nx, int ny);
    double LapJas                 ();
    double computeJastrow         ();
    double H                      (int state, double x);
    double dH                     (int state, double x);
    double ddH                    (int state, double x);
    double computeQuantumForce    (int p, int d);
    void   computePsiBars         (double &psiBar_alpha, double &psiBar_beta);
  protected:
    double my_jf = 0;
};
