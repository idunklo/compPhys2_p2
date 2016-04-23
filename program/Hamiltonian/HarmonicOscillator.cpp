#include "HarmonicOscillator.h"

HarmonicOscillator::HarmonicOscillator (System* system):
  Hamiltonian(system)
{}

double HarmonicOscillator::computeNumLocalEnergy ()
{
/*
  double  HOLap         = 0;
  double  HOExt         = 0;
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
      HOExt += x*x;
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
  HOExt = 0.5*HOExt*omega2;

  return HOLap + HOExt + Hrep;
  */
}

double HarmonicOscillator::computeAnaLocalEnergy()
{
  /*
  const int nD          = my_system->get_nDimensions();
  const double alpha    = my_system->get_parameters()[0];
  const double beta     = my_system->get_parameters()[1];
  const double omega    = my_system->get_parameters()[2];
  const double a        = my_system->get_parameters()[3];
  const double alom     = alpha*omega;
  double       HOLap    = 0;
  double       HOExt    = 0;
  double       Hrep     = 0;
  double       r2       = 0;
  double       sep2     = 0;
  
  for (int d = 0 ; d < nD ; d++){
    const double x1   = my_system->get_particle()[0]->get_position()[d];
    const double x2   = my_system->get_particle()[1]->get_position()[d];
    const double sep  = x1-x2;
    r2   += x1*x1+x2*x2;
    sep2 += sep*sep;
  }

  const double r12        = sqrt(sep2);
  const double opbr12     = 1 + beta*r12;
  const double opbr122    = opbr12*opbr12;
  const double a_opbr122  = a/opbr122;

  HOLap = alom*alom*r2 
          -4*alom 
          -2*alom*r12*a_opbr122
          +2*(a_opbr122)*(a_opbr122 + (1-beta*r12)/(r12*opbr12));

  HOLap = -0.5*HOLap;
  HOExt = 0.5*omega*omega*(r2);
  Hrep  = 1/r12;
  return (HOLap + HOExt + Hrep);
  */
}
