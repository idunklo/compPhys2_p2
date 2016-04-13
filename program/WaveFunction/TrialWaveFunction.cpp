#include "TrialWaveFunction.h"

using std::cout;
using std::endl;

TrialWaveFunction::TrialWaveFunction (System* system) :
  WaveFunction(system){
  }

double TrialWaveFunction::evaluate ()
{
  int nP            = my_system->get_nParticles();
  int nD            = my_system->get_nDimensions();
  double alpha      = my_system->get_parameters()[0];
  double beta       = my_system->get_parameters()[1];
  double omega      = my_system->get_parameters()[2];
  double a          = my_system->get_parameters()[3]; 

  double argument1  = 0;
  double argument2  = 0;
  double sep        = 0;


  for (int p1 = 0 ; p1 < nP ; p1++){
    for (int d = 0 ; d < nD ; d++){
      const double x = my_system->get_particle()[p1]->get_position()[d];
      argument1 += x*x;
    }

    for (int p2 = p1+1 ; p2 < nP ; p2++){
      sep = 0;
      for (int d = 0 ; d < nD ; d++){
        const double deltaX = my_system->get_particle()[p1]->get_position()[d]- 
                              my_system->get_particle()[p2]->get_position()[d];
        sep += deltaX*deltaX;
      }
      sep = sqrt(sep);
      argument2 += sep/(1+beta*sep);
    }
  }
  argument1 = argument1*(-0.5*alpha*omega);
  argument2 *= a;

  return exp(argument1 + argument2);
}

double TrialWaveFunction::computeGradient(int p, int d)
{
  int          sign   = (2*(p==0)-1);
  double       sep2   = 0;
  double       x      = 0;
  const double alpha  = my_system->get_parameters()[0];
  const double beta   = my_system->get_parameters()[1];
  const double omega  = my_system->get_parameters()[2];
  const double a      = my_system->get_parameters()[3];
  const double xi     = my_system->get_particle()[p]->get_position()[d];
  for (int di = 0 ; di < my_system->get_nDimensions() ; di++){
    const double x1  = my_system->get_particle()[0]->get_position()[di];
    const double x2  = my_system->get_particle()[1]->get_position()[di];
    x    += x1;
    sep  += (x1-x2);
  }
  const double r12 = sep*sep;
  gradient = -alpha*omega*xi + sign*a*sep/(r12*(1+beta*r12)*(1+beta*r12));
  return gradient;
}

double TrialWaveFunction::computeQuantumForce(int p, int d)
{
  //This is wrong, must change!
  /*
  double quantumForce	= 0;
  double alom         = my_system->get_parameters()[0]*
                        my_system->get_parameters()[2];

  const double x = my_system->get_particle()[p]->get_position()[0]; 
  const double y = my_system->get_particle()[p]->get_position()[1];
  const double r = sqrt(x*x + y*y);
  quantumForce = -r*alom;
  return quantumForce;
  */
}

