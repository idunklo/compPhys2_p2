#include "TrialSlater.h"

using std::cout;
using std::endl;

TrialSlater::TrialSlater (System* system, double jf) :
  WaveFunction(system)
{
  my_jf = jf;
}


double TrialSlater::Phi(int p, int nx, int ny)
{
  const double x = my_system->get_particle().at(p)->get_position().at(0);
  const double y = my_system->get_particle().at(p)->get_position().at(1);
  const double omega = my_system->get_parameters().at(0);
  const double omegasq = sqrt(omega);
  const double Hnx = my_system->Hermite_n(nx,omegasq*x);
  const double Hny = my_system->Hermite_n(ny,omegasq*y);
 
  //cout << Hnx*Hny*exp(-omega*(x*x+y*y)*0.5)<<endl;
  return Hnx*Hny*exp(-omega*(x*x+y*y)*0.5);
}

double TrialSlater::computeJastrow ()
{
  const int nP        = my_system->get_nParticles();
  const int nP_2      = nP/2;
  const double beta   = my_system->get_parameters().at(1);
  const Eigen::MatrixXd rij = my_system->get_r_ij();
  double argument     = 0.0;
  for (int i=0 ; i<nP; i++){
    for (int j=i+1 ; j<nP; j++){
      const double sep = rij(i,j);
      const double a   = 1-((i<nP_2)*(j>=nP_2)!=1)*(2.0/3.0);
      argument += a*sep/(1+beta*sep); 
    }
  }
  return exp(argument/my_jf);
}

double TrialSlater::GradPhi(int pos, int d)
{
  const int nP       = my_system->get_nParticles();
  const double omega = my_system->get_parameters()[0];
  const int orbitals = my_system->get_orbitals();
  const double x     = my_system->get_particle()[pos]->get_position()[d];
  const double xi    = my_system->get_particle()[pos]->get_position()[0];
  const double yi    = my_system->get_particle()[pos]->get_position()[1];
  int orbital = 0;
  int col     = 0;
  double Grad = 0.0;
  double grad = 0.0;
  double arg  = 0.0;
  Eigen::VectorXd d_inv(nP);
  if(pos<nP/2)
    d_inv = my_system->get_DMatrix_up_inv().col(pos);
  else
    d_inv = my_system->get_DMatrix_dn_inv().col(pos-nP/2);

  for (int shell = 0 ; shell <= orbitals ; shell++){
    int nx = shell; int ny = 0;
    for (int state = 0 ; state <= shell ; state++){
      if(d==0){orbital=nx;}
      else{orbital=ny;}

      switch (orbital)
      {
        case 0:
          grad = -omega*x;
          break;
        case 1:
          grad = (2-2*omega*x*x);
          break;
        case 2:
          grad = (8 + 2*omega*x - 4*omega*x*x);
          break;
        case 3:
          grad = (24*x*x + 12*omega*x*x - 8*omega*x*x*x*x - 12);
          break;
      }
      Grad += grad*d_inv(col); 
      col++; nx--; ny++;
    }
  }
  return Grad*exp(-omega*(xi*xi+yi*yi)*0.5);
}

double TrialSlater::GradJas(int k, int d)
{
  const int nP         = my_system->get_nParticles();
  const int nP_2       = nP/2;
  const double beta    = my_system->get_parameters().at(1);
  const Eigen::MatrixXd r_ij = my_system->get_r_ij();
  double Grad = 0.0;
  const double xk = my_system->get_particle().at(k)->get_position().at(d);

  for (int j = 0 ; j < k ; j++){
    const double a   = 1-((k<nP_2)*(j>=nP_2)!=1)*(2.0/3.0);
    const double xj = my_system->get_particle().at(j)->get_position().at(d);
    const double onePbetarkj = my_jf*(1 + beta*r_ij(k,j));

    Grad += a*(xk-xj)/(r_ij(j,k)*onePbetarkj*onePbetarkj);
  }
  for (int j = k+1 ; j<nP ; j++){
    const double a   = 1-((k<nP_2)*(j>=nP_2)!=1)*(2.0/3.0);
    const double xj = my_system->get_particle().at(j)->get_position().at(d);
    const double onePbetarkj = my_jf*(1 + beta*r_ij(k,j));

    Grad += a*(xk-xj)/(r_ij(k,j)*onePbetarkj*onePbetarkj);
  }
  return Grad;
}

double TrialSlater::LapPhi(int pos,int orbital)
{
  double Lap = 0.0;
  double arg = 0.0;
  const double omega = my_system->get_parameters().at(0);
  for (int xi = 0 ; xi < 2 ; xi++){
    const double x = my_system->get_particle()[pos]->get_position()[xi];
    switch (orbital)
    {
      case 0:
        Lap += (omega*x*x - 1);
        break;
      case 1:
        Lap += (2*omega*x*x*x - 6*x);
        break;
      case 2:
        Lap += (4*omega*x*x*x*x - 2*omega*x*x - 20*x*x +2 +8/omega);
        break;
      case 3:
        Lap += (48*x/omega + 36*x - 56*x*x*x - 12*omega*x*x*x 
                + 8*omega*x*x*x*x*x);
        break;
    }
    arg += x*x;
  }
  return Lap*omega*exp(-omega*arg*0.5);
}

double TrialSlater::LapJas()
{
  const int nP         = my_system->get_nParticles();
  const int nP_2       = nP/2;
  const double beta    = my_system->get_parameters().at(1);
  const Eigen::MatrixXd r_ij = my_system->get_r_ij();
  double Lap  = 0.0;

  for (int k=0 ; k<nP ; k++){
    const double xk = my_system->get_particle().at(k)->get_position().at(0);
    const double yk = my_system->get_particle().at(k)->get_position().at(1);
    double term1= 0.0;
    double term2= 0.0;

    for (int i=k+1 ; i<nP ; i++){
      const double xi = my_system->get_particle().at(i)->get_position().at(0);
      const double yi = my_system->get_particle().at(i)->get_position().at(1);
      const double aki= 1-((k<nP_2)*(i>=nP_2)!=1)*(2.0/3.0);
      const double r_ki = r_ij(k,i);
      const double bri  = my_jf*(1+beta*r_ki);

      const double jast_ki = aki/(r_ki*bri*bri);
        
      for(int j=k+1 ; j<nP ; j++){
        const double xj = my_system->get_particle().at(j)->get_position().at(0);
        const double yj = my_system->get_particle().at(j)->get_position().at(1);
        const double akj= 1-((k<nP_2)*(j>=nP_2)!=1)*(2.0/3.0);
        const double r_kj = r_ij(k,j); 
        const double brj  = my_jf*(1+beta*r_kj);

        const double jast_kj = akj/(r_kj*brj*brj);
        const double rkri_rkrj = (xk-xi)*(xk-xj)+(yk-yi)*(yk-yj);

        term1 += rkri_rkrj*jast_kj;
      }
      term1 *= jast_ki;
      term2 += aki/(r_ki*bri*bri*bri);
    }
    Lap += (term1+2*term2);
  }
  return Lap;
}

double TrialSlater::computeQuantumForce(int p, int d)
{
  //return 2*computeGradient(p,d);
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
    const double x1  = my_system->get_particle()[0]->get_position()[d];
    const double x2  = my_system->get_particle()[1]->get_position()[d];
    sep2 += (x1-x2)*(x1-x2);
    r2  += x1*x1+x2*x2;
  }
  const double r12 = sqrt(sep2);
  psiBar_alpha = -0.5*omega*r2;
  psiBar_beta  = -a*r12*r12/((1+beta*r12)*(1+beta*r12));
  */
}
