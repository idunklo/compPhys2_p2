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
  double  rep           = 0;
  double  omega         = my_system->get_parameters()[2];
  double  omega2        = omega*omega;
  int     nP            = my_system->get_nParticles();
  int     nD            = my_system->get_nDimensions();

  waveFuncNow = my_system->get_waveFunction()->evaluate();

  for (int p1 = 0 ; p1 < nP ; p1++){
    for (int d = 0 ; d < nD ; d++){
      ddPsi += my_system->get_waveFunction()->
               computeDoubleDerivative(p1,d,waveFuncNow);
      const double x = my_system->get_particle()[p1]->get_position()[d];
      HOext += x*x;
    }
    
    for (int p2 = p1 + 1 ; p2 < nP ; p2++){
      rep = 0;
      for (int d = 0 ; d < nD ; d++){
        const double deltaX = my_system->get_particle()[p1]->get_position()[d]- 
                              my_system->get_particle()[p2]->get_position()[d];
        rep += deltaX*deltaX;
      }
      Hrep += 1/sqrt(rep);
    }
  }
  HOLap = -0.5*ddPsi/waveFuncNow;
  HOext = 0.5*HOext*omega2;
  Hrep = 0; 
  return HOLap + HOext + Hrep;
}
