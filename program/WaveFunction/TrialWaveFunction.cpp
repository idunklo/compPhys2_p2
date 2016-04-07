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
  argument2 = 0;
  return exp(argument1 + argument2);
}

double TrialWaveFunction::computeGradient(int p, int d)
{
  int sign      = (1-2*(p%2==0));
  double d1     = 0;
  double d2     = 0;
  double expression = 0; 
  double alpha  = my_system->get_parameters()[0];
  double beta   = my_system->get_parameters()[1];
  double a      = my_system->get_parameters()[3];
  double x1     = my_system->get_particle()[0]->get_position()[0];
  double x2     = my_system->get_particle()[1]->get_position()[0];
  double y1     = my_system->get_particle()[0]->get_position()[1];
  double y2     = my_system->get_particle()[1]->get_position()[1];
  double r12    = sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
  if (d==0){
    d1 = x1;
    d2 = x2;
  }
  else{
    d1 = y1;
    d2 = y2;
  }
  expression = -alpha*2*d + sign*(a*(d1-d2))/((beta*r12+1)*(beta*r12+1)*r12);
  return expression;
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

