#include "TrialSlater.h"

using std::cout;
using std::endl;

TrialSlater::TrialSlater (System* system) :
  WaveFunction(system)
{
}


double TrialSlater::Phi(int p, int nx, int ny)
{
  const double x = my_system->get_particles()(p,0);
  const double y = my_system->get_particles()(p,1);
  const double omega = my_system->get_parameters()[0];
  const double alpha = my_system->get_parameters()[2];
  const double omegasq = sqrt(omega*alpha);
  const double Hnx = H(nx,omegasq*x);
  const double Hny = H(ny,omegasq*y);

  return Hnx*Hny*exp(-omega*alpha*(x*x+y*y)*0.5);
}

double TrialSlater::computeJastrow ()
{
  const int nP        = my_system->get_nParticles();
  const int nP_2      = nP/2;
  const double beta   = my_system->get_parameters()[1];
  const Eigen::MatrixXd rij = my_system->get_r_ij();
  const Eigen::MatrixXd aij = my_system->get_a();
  double argument     = 0.0;
  for (int i=0 ; i<nP; i++){
    for (int j=i+1 ; j<nP; j++){
      const double sep = rij(i,j);
      const double a   = aij(i,j);
      argument += a*sep/(1+beta*sep); 
    }
  }
  return exp(argument);
}

double TrialSlater::GradPhi(int k, int d)
{
  const int nP       = my_system->get_nParticles();
  const int orbitals = my_system->get_orbitals();
  const double omega = my_system->get_parameters()[0];
  const double alpha = my_system->get_parameters()[2];
  const double xi    = my_system->get_particles()(k,d);
  const double yi    = my_system->get_particles()(k,(d==0));
  const double e     = exp(-omega*alpha*(xi*xi+yi*yi)*0.5);
  const double de    = -omega*alpha*xi;
  const double x     = sqrt(omega*alpha)*xi;
  const double y     = sqrt(omega*alpha)*yi;
  int col     = 0;
  double Grad = 0.0;
  Eigen::VectorXd d_inv(nP);
  if(k<nP/2)
    d_inv = my_system->get_DMatrix_up_inv().col(k);
  else
    d_inv = my_system->get_DMatrix_dn_inv().col(k-nP/2);
  for (int shell = 0 ; shell <= orbitals ; shell++){
    int nx = shell; int ny = 0;
    for (int state = 0 ; state <= shell ; state++){
      const double grad = dH(nx,x)*H(ny,y) + H(nx,x)*H(ny,y)*de;
      Grad += grad*d_inv(col); 
      col++; nx--; ny++;
    }
  }
  return Grad*e;
}

double TrialSlater::GradJas(int k, int d)
{
  const int nP         = my_system->get_nParticles();
  const int nP_2       = nP/2;
  const double beta    = my_system->get_parameters()[1];
  const Eigen::MatrixXd r_ij = my_system->get_r_ij();
  const Eigen::MatrixXd a_ij = my_system->get_a();
  const double xk = my_system->get_particles()(k,d);
  double Grad = 0.0;
  for (int j = 0 ; j < k ; j++){
    const double xj = my_system->get_particles()(j,d);
    const double onePbetarkj = (1 + beta*r_ij(j,k));
    Grad += a_ij(j,k)*(xk-xj)/(r_ij(j,k)*onePbetarkj*onePbetarkj);
  }
  for (int j = k+1 ; j<nP ; j++){
    const double xj = my_system->get_particles()(j,d);
    const double onePbetarkj = (1 + beta*r_ij(k,j));
    Grad += a_ij(k,j)*(xk-xj)/(r_ij(k,j)*onePbetarkj*onePbetarkj);
  }
  return Grad;
}

double TrialSlater::LapPhi(int pos,int nx, int ny)
{
  double Lap = 0.0;
  const double omega = my_system->get_parameters()[0];
  const double alpha = my_system->get_parameters()[2];
  const double alom  = alpha*omega;
  const double xi  = my_system->get_particles()(pos,0);
  const double yi  = my_system->get_particles()(pos,1);
  const double e   = exp(-alom*(xi*xi+yi*yi)*0.5);
  const double dex = -alom*xi;
  const double dey = -alom*yi;
  const double dde = alom*alom*(xi*xi+yi*yi)-2*alom;
  const double x = sqrt(alom)*xi;
  const double y = sqrt(alom)*yi;

  Lap = ddH(nx,x)*H(ny,y) + ddH(ny,y)*H(nx,x) + H(ny,y)*H(nx,x)*dde
        + 2*(dH(nx,x)*H(ny,y)*dex + dH(ny,y)*H(nx,x)*dey);

  return Lap*e;
}


double TrialSlater::LapJas()
{
  const int nP         = my_system->get_nParticles();
  const int nP_2       = nP/2;
  const double beta    = my_system->get_parameters()[1];
  const Eigen::MatrixXd r_ij = my_system->get_r_ij();
  const Eigen::MatrixXd a    = my_system->get_a();
  double Lap  = 0.0;

  for (int k=0 ; k<nP ; k++){
    const double xk = my_system->get_particles()(k,0);
    const double yk = my_system->get_particles()(k,1);
    double term1= 0.0;
    double term2= 0.0;

    for (int i=k+1 ; i<nP ; i++){
      const double xi = my_system->get_particles()(i,0);
      const double yi = my_system->get_particles()(i,1);
      const double aki= a(k,i);
      const double r_ki = r_ij(k,i);
      const double bri  = 1+beta*r_ki;
      
      const double jast_ki = aki/(r_ki*bri*bri);
        
      for(int j=k+1 ; j<nP ; j++){
        const double xj = my_system->get_particles()(j,0);
        const double yj = my_system->get_particles()(j,1);
        const double akj= a(k,j);
        const double r_kj = r_ij(k,j); 
        const double brj  = 1+beta*r_kj;

        const double jast_kj = akj/(r_kj*brj*brj);
        const double rkri_rkrj = (xk-xi)*(xk-xj)+(yk-yi)*(yk-yj);

        term1 += rkri_rkrj*jast_kj;
      }
      term1 *= jast_ki;
      term2 += aki*(1-beta*r_ki)/(r_ki*bri*bri*bri);
    }
    Lap += (term1+term2);
  }
  return Lap;
}

double TrialSlater::H(int state, double x)
{
  switch(state)
  {
    case 0:
      return 1;
      break;

    case 1:
      return 2*x;
      break;

    case 2:
      return 4*x*x - 2;
      break;

    case 3:
      return 8*x*x*x - 12*x;
      break;
  }
  return 0.0;
}

double TrialSlater::dH(int state, double x)
{
  switch(state)
  {
    case 0:
      return 0;
      break;

    case 1:
      return 2;
      break;

    case 2:
      return 8*x;
      break;

    case 3:
      return 24*x*x - 12;
      break;
  }
  return 0.0;
}

double TrialSlater::ddH(int state, double x)
{
  switch(state)
  {
    case 0:
      return 0;
      break;

    case 1:
      return 0;
      break;

    case 2:
      return 8;
      break;

    case 3:
      return 48*x;
      break;
  }
  return 0.0;
}

Eigen::Vector2d TrialSlater::computeQuantumForce(int p, int d)
{
  Eigen::Vector2d Qforce (GradPhi(p,d),GradJas(p,d));
  return 2*Qforce;
  //return 2*(GradPhi(p,d)+GradJas(p,d));
}

void TrialSlater::computePsiBars(double &psiBar_alpha,
                                       double &psiBar_beta)
{
  /*
  const double beta  = my_system->get_parameters()[1];
  const double omega = my_system->get_parameters()[0];
  const double a     = my_system->get_parameters()[3];
  const int nD   = my_system->get_nDimensions();
  double    r2   = 0;
  double    sep2 = 0;
  for (int d = 0 ; d < nD ; d++){
    const double x1  = my_system->get_particles()[0]->get_position()[d];
    const double x2  = my_system->get_particles()[1]->get_position()[d];
    sep2 += (x1-x2)*(x1-x2);
    r2  += x1*x1+x2*x2;
  }
  const double r12 = sqrt(sep2);
  psiBar_alpha = -0.5*omega*r2;
  psiBar_beta  = -a*r12*r12/((1+beta*r12)*(1+beta*r12));
  */
}

double TrialSlater::dPhi_alpha(int p, int nx, int ny)
{
  const double x = my_system->get_particles()(p,0);
  const double y = my_system->get_particles()(p,1);
  const double omega = my_system->get_parameters()[0];
  const double alpha = my_system->get_parameters()[2];
  const double alom  = omega*alpha;
  const double omegasq = sqrt(alom);
  const double Hnx = H(nx,omegasq*x);
  const double Hny = H(ny,omegasq*y);
 
  return Hnx*Hny*(-omega*(x*x+y*y)*0.5)*exp(-alom*(x*x+y*y)*0.5);
}

double TrialSlater::dlnjast_beta ()
{
  const int nP        = my_system->get_nParticles();
  const int nP_2      = nP/2;
  const double beta   = my_system->get_parameters()[1];
  const Eigen::MatrixXd rij = my_system->get_r_ij();
  const Eigen::MatrixXd aij = my_system->get_a();
  double argument     = 0.0;
  for (int i=0 ; i<nP; i++){
    for (int j=i+1 ; j<nP; j++){
      const double sep = rij(i,j);
      const double a   = aij(i,j);
      argument += a*sep*sep/((1+beta*sep)*(1+beta*sep)); 
    }
  }
  return argument;
}

