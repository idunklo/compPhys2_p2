#pragma once
#include "../System/System.h"

class WaveFunction
{
  public:
  WaveFunction (class System* system);

    virtual double evaluate		          () = 0;
    virtual double computeQuantumForce  (int p, int d) = 0;
    virtual double computeGradient      (int p, int d) = 0;
    double computeDerivative            (int p, int d, double waveFunctionCurrent);
    double computeDoubleDerivative      (int p, int d, double waveFunctionCurrent);

  protected:
    class System* my_system = nullptr;
};
