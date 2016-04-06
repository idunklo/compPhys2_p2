#include "HarmonicOscillator.h"

HarmonicOscillator::HarmonicOscillator (System* system):
  Hamiltonian(system)
{}

double HarmonicOscillator::computeLocalEnergy ()
{
  double  HOLap         = 0;
  double  HOext         = 0;
  double  Hrep          = 0;
  double  ddPsi         = 0;
  double  waveFuncNow   = 0;
  double  omega         = my_system->get_parameters()[2];
  double  omega2        = omega*omega;
  int     nP            = my_system->get_nParticles();
  int     nD            = my_system->get_nDimensions();

  waveFuncNow = my_system->get_waveFunction()->evaluate();

  for (int p1 = 0 ; p1 < nP ; p1++){
    ddPsi += my_system->get_waveFunction()->
      computeDoubleDerivative(p1,0,waveFuncNow);
    ddPsi += my_system->get_waveFunction()->
      computeDoubleDerivative(p1,1,waveFuncNow);
    const double x = my_system->get_particle()[p1]->get_position()[0];
    const double y = my_system->get_particle()[p1]->get_position()[1];
    const double r2 = x*x + y*y; 

    HOext += r2;
    
    for (int p2 = p1 + 1 ; p2 < nP ; p2++){
      const double deltaX = x - my_system->get_particle()[p2]->get_position()[0];
      const double deltaY = y - my_system->get_particle()[p2]->get_position()[1];
      Hrep += 1/sqrt(deltaX*deltaX + deltaY*deltaY);
    }
  }
  HOLap = -0.5*ddPsi/waveFuncNow;
  HOext = 0.5*HOext*omega2;
  
  return HOLap + HOext + Hrep;
}
