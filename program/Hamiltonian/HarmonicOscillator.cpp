#include "HarmonicOscillator.h"

HarmonicOscillator::HarmonicOscillator (System* system):
  Hamiltonian(system)
{}

double HarmonicOscillator::computeLocalEnergy ()
{
  double  HOpsi         = 0;
  double  HOterm1       = 0;
  double  HOterm2       = 0;
  double  Hrep          = 0;
  double  alpha         = my_system->get_parameters()[0];
  double  omega         = my_system->get_parameters()[2];
  double  alom          = alpha*omega;
  double  omega2        = omega*omega;
  int     nP            = my_system->get_nParticles();
  int     nD            = my_system->get_nDimensions();


  for (int p1 = 0 ; p1 < nP ; p1++){
    const double x = my_system->get_particle()[p1]->get_position()[0]; 
    const double y = my_system->get_particle()[p1]->get_position()[1];
    const double r = sqrt(x*x + y*y);

    HOterm1  += -r*(alom*r - 1);
    HOterm2  += omega*r*r;
    for (int p2 = p1 + 1 ; p2 < nP ; p2++){
      const double deltaX = x - my_system->get_particle()[p2]->get_position()[0];
      const double deltaY = y - my_system->get_particle()[p2]->get_position()[1];
      Hrep += 1/sqrt(deltaX*deltaX + deltaY*deltaY);
    }
  }

  HOpsi = 0.5*omega*(alpha*HOterm1 + HOterm2);

  return HOpsi + Hrep;
}
