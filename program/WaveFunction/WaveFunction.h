#pragma once
#include "../System/System.h"

class WaveFunction
{
  public:
  WaveFunction (class System* system);

    virtual double Phi                  (int p, int nx, int ny) = 0;
    virtual double computeQuantumForce  (int p, int d) = 0;
    virtual double GradPhi              (int pos, int d) = 0;
    virtual double LapPhi               (int pos, int orbital) = 0;
    virtual double GradJas              (int k, int d) = 0;
    virtual double LapJas               () = 0;
    virtual double computeJastrow       () = 0;
    virtual void   computePsiBars       (double &psiBar_alpha, double &psiBar_beta)=0;
    double computeDerivative            (int p, int d, double waveFunctionCurrent);
    double computeDoubleDerivative      (int p, int d, double waveFunctionCurrent);
    

  protected:
    class System* my_system = nullptr;
};
