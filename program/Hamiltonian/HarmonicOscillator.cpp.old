#include "HarmonicOscillator.h"


HarmonicOscillator::HarmonicOscillator (System* system) :
  Hamiltonian(system)
{
}

double HarmonicOscillator::computeLocalEnergy ()
{
  bool	 analytical	      = my_system->get_analytical();
  int	 nD		      = my_system->get_nDimensions();
  int	 nP		      = my_system->get_nParticles();
  double potEnergy	      = 0;
  double kinEnergy	      = 0;
  double doubleDerivative     = 0;
  double alpha		      = my_system->get_parameters()[0];
  double omega2		      = my_system->get_parameters()[1]*
				my_system->get_parameters()[1];
  double gamma		      = my_system->get_parameters()[2];
  double gamma2		      = gamma*gamma;
  double waveFunctionCurrent  = my_system->get_waveFunction()->evaluate();
  
  
 
  if (analytical){
    for (int p = 0 ; p < nP ; p++){
      for (int d = 0 ; d < nD ; d++){
        const double x = my_system->get_particle()[p]->get_position()[d];

        if (d==2){
          potEnergy += x*x*gamma2;
	  kinEnergy -= (2*alpha*gamma*x*x - 1)*gamma;
        }
        else{
          potEnergy += x*x*omega2;
	  kinEnergy -= (2*alpha*x*x - 1);
        }
      }
    }
    kinEnergy = kinEnergy*alpha;
    potEnergy = 0.5*potEnergy;
  }		  
  else{
    for (int p = 0 ; p < nP ; p++){
      for (int d = 0 ; d < nD ; d++){
        const double x = my_system->get_particle()[p]->get_position()[d];
        if (d==2){
          potEnergy += x*x*gamma2;
        }
        else{
          potEnergy += x*x*omega2;
        }
        doubleDerivative += my_system->get_waveFunction()->
	  computeDoubleDerivative(p,d,waveFunctionCurrent);
      }
    }
    kinEnergy = -0.5*doubleDerivative/waveFunctionCurrent;
    potEnergy = 0.5*potEnergy;
  }		  

  return kinEnergy + potEnergy;
}
 

