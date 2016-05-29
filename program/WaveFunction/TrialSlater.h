#pragma once
#include "WaveFunction.h"

class TrialSlater: public WaveFunction
{
  public:
    TrialSlater                   (class System* system);
    double Phi                    (int p, int nx, int ny);
    double GradPhi                (int pos, int d);
    double GradJas                (int k,int d);
    double LapPhi                 (int pos, int nx, int ny);
    double LapJas                 ();
    double computeJastrow         ();
    double H                      (int state, double x);
    double dH                     (int state, double x);
    double ddH                    (int state, double x);
    double dPhi_alpha             (int p, int nx, int ny);
    double dlnjast_beta           ();
    Eigen::Vector2d computeQuantumForce    (int p, int d);
    void   computePsiBars         (double &psiBar_alpha, double &psiBar_beta);

  protected:
};
