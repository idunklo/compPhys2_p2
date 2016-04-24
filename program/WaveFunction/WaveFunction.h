#pragma once
#include "../System/System.h"

class WaveFunction
{
  public:
  WaveFunction (class System* system);

    virtual double evaluate		          (int p, int nx, int ny) = 0;
    virtual double computeQuantumForce  (int p, int d) = 0;
    virtual double computeGradient      (int row, int level) = 0;
    virtual double computeLaplacian     (int row, int level) = 0;
    virtual void   computePsiBars       (double &psiBar_alpha, double &psiBar_beta)=0;
    double computeDerivative            (int p, int d, double waveFunctionCurrent);
    double computeDoubleDerivative      (int p, int d, double waveFunctionCurrent);
    

  protected:
    class System* my_system = nullptr;
};
