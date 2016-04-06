#include "TrialWaveFunction.h"

using std::cout;
using std::endl;

TrialWaveFunction::TrialWaveFunction (System* system) :
  WaveFunction(system){
  }

double TrialWaveFunction::evaluate ()
{
  int nP            = my_system->get_nParticles();
  double a          = 1; 
  double alpha      = my_system->get_parameters()[0];
  double beta       = my_system->get_parameters()[1];
  double omega      = my_system->get_parameters()[2];

  double argument1  = 0;
  double argument2  = 0;


  for (int p1 = 0 ; p1 < nP ; p1++){
    const double x = my_system->get_particle()[p1]->get_position()[0];
    const double y = my_system->get_particle()[p1]->get_position()[1];
    argument1 += x*x + y*y;

    for (int p2 = p1+1 ; p2 < nP ; p2++){
      const double deltaX = x - my_system->get_particle()[p2]->get_position()[0];
      const double deltaY = y - my_system->get_particle()[p2]->get_position()[1];
      const double sep    = sqrt(deltaX*deltaX + deltaY*deltaY);
      argument2 += sep/(1+beta*sep);
    }
  }
  argument1 = argument1*(-0.5*alpha*omega);
  argument2 *= a;

  return exp(argument1 + argument2);
}

double TrialWaveFunction::computeQuantumForce(int p)
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

