#include "HarmonicOscillator.h"

HarmonicOscillator::HarmonicOscillator (System* system):
  Hamiltonian(system)
{}

double HarmonicOscillator::computeAnaLocalEnergy()
{
  const int orbitals = my_system->get_orbitals();
  const int spin     = my_system->get_spin();
  const int nP       = my_system->get_nParticles();
  const int nD       = my_system->get_nDimensions();
  double Grads  = 0.0;
  double SD_LAP = 0.0;
  double SD_LAP_new = 0.0;
  int startPos  = 0;
  int orbital   = 0;
  if(spin){
    SD_LAP = my_system->get_SDLap_old_dn();
    startPos=nP/2;
  }
  else{
    SD_LAP = my_system->get_SDLap_old_up();
  }
  for (int shell = 0 ; shell <= orbitals ; shell++){
    int nx = shell; int ny = 0;
    for (int state = 0 ; state <= shell ; state++){
      for (int pos = startPos ; pos<(2*startPos) ; pos++){
        for (int d = 0 ; d < nD ; d++){
          if(d==0){orbital=nx;}
          else{orbital=ny;}
          SD_LAP_new += my_system->get_waveFunction()->LapPhi(pos,orbital);
        }
      }
      nx -= 1; ny += 1;
    }
  }
  const double JastrowLap= my_system->get_waveFunction()->LapJas();
  double GradSum = 0.0;
  for (int i = 0 ; i<nP ; i++){
    for (int d = 0 ; d < nD ; d++){
      const double phi = my_system->get_waveFunction()->GradPhi(i,d);
      const double jas = my_system->get_waveFunction()->GradJas(i,d);
      GradSum += phi*jas; 
    }
  }

  if(spin){
    my_system->set_SDLap_up(SD_LAP_new);
  }
  else{
    my_system->set_SDLap_dn(SD_LAP_new);
  }
  SD_LAP += SD_LAP_new; 
  return SD_LAP + JastrowLap + 2*GradSum;
}

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

