#include "WaveFunction.h"

WaveFunction::WaveFunction (System* system)
{
  my_system = system;
}


double WaveFunction::computeDoubleDerivative(int p, int d, double waveFunctionCurrent)
{
  double doubleDerivative     = 0;
  double waveFunctionPlus     = 0;
  double waveFunctionMinus    = 0;
 

  my_system->get_particle()[p]->changePosition(d, my_system->get_derivativeStep());
  waveFunctionPlus = evaluate();

  my_system->get_particle()[p]->changePosition(d, -2*my_system->get_derivativeStep());
  waveFunctionMinus = evaluate();

  my_system->get_particle()[p]->changePosition(d, my_system->get_derivativeStep());
  
  doubleDerivative += waveFunctionPlus + waveFunctionMinus - 2*waveFunctionCurrent;
  
  doubleDerivative = doubleDerivative/my_system->get_derivativeStep2();
  return doubleDerivative;
}

double WaveFunction::computeDerivative(int p, int d, double waveFunctionCurrent)
{
  double derivative	      = 0;
  double waveFunctionPlus = 0;

  my_system->get_particle()[p]->changePosition(d, my_system->get_derivativeStep());

  waveFunctionPlus = evaluate();

  my_system->get_particle()[p]->changePosition(d,-my_system->get_derivativeStep());

  derivative = (waveFunctionPlus - waveFunctionCurrent)/my_system->get_derivativeStep();

  return derivative;
}

