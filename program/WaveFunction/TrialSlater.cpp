#include "TrialSlater.h"

using std::cout;
using std::endl;

TrialSlater::TrialSlater (System* system) :
  WaveFunction(system){
  }

double TrialSlater::evaluate (int p, int nx, int ny)
{
  const double x = my_system->get_particle().at(p)->get_position().at(0);
  const double y = my_system->get_particle().at(p)->get_position().at(1);
  const double omega = my_system->get_parameters().at(2);
  const double omegasq = sqrt(omega);
  const double Hnx = my_system->Hermite_n(nx,omegasq*x);
  const double Hny = my_system->Hermite_n(ny,omegasq*y);

  return Hnx*Hny*exp(-omega*(x*x+y*y)*0.5);
}


double TrialSlater::computeGradient(int row, int level)
{
  int nP = my_system->get_nParticles();
  const double omega = my_system->get_parameters()[2];
  double grad = 0.0;
  for (int i=0; i<nP ; i++){
    const double x = my_system->get_particle()[i]->get_position()[0];
    const double y = my_system->get_particle()[i]->get_position()[1];
    switch (level)
    {
      case 0:
        grad += -omega*x;
        break;
      case 1:
        grad += (1-omega*x*x)*2;
        break;
      case 2:
        grad += (4 + omega*x - 2*omega*x*x)*2;
        break;
      case 3:
        grad += (6*x*x + 3*omega*x*x - 2*omega*x*x*x*x - 3)*4;
        break;
    }
    grad *= exp(-omega*(x*x + y*y)*0.5);
  }
  return grad;
}

double TrialSlater::computeLaplacian(int row, int level)
{
  int nP = my_system->get_nParticles();
  const double omega = my_system->get_parameters()[2];
  double lap = 0.0;
  for (int i=0; i<nP ;i++){
    const double x = my_system->get_particle()[i]->get_position()[0];
    const double y = my_system->get_particle()[i]->get_position()[1];
    switch (level)
    {
      case 0:
        lap += (-omega + omega*x*x);
        break;
      case 1:
        lap += (omega*omega*x*x*x-3*omega*x)*2;
        break;
      case 2:
        lap += (2*omega*omega*x*x*x + omega - 6*omega*x);
        break;
      case 3:
        lap += (48*x + 36*omega*x - 56*omega*x*x*x - 12*omega*omega*x*x*x 
                + 8*omega*omega*x*x*x*x*x + 12*omega*x);
        break;
    }
    lap *= exp(-omega*(x*x + y*y)*0.5);
  }
  return lap;

}
double TrialSlater::computeQuantumForce(int p, int d)
{
  return 2*computeGradient(p,d);
}

void TrialSlater::computePsiBars(double &psiBar_alpha,
                                       double &psiBar_beta)
{
  const double beta  = my_system->get_parameters()[1];
  const double omega = my_system->get_parameters()[2];
  const double a     = my_system->get_parameters()[3];
  const int nD   = my_system->get_nDimensions();
  double    r2   = 0;
  double    sep2 = 0;
  for (int d = 0 ; d < nD ; d++){
    const double x1  = my_system->get_particle()[0]->get_position()[d];
    const double x2  = my_system->get_particle()[1]->get_position()[d];
    sep2 += (x1-x2)*(x1-x2);
    r2  += x1*x1+x2*x2;
  }
  const double r12 = sqrt(sep2);
  psiBar_alpha = -0.5*omega*r2;
  psiBar_beta  = -a*r12*r12/((1+beta*r12)*(1+beta*r12));
}

