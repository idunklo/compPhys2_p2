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
  const double omega = my_system->get_parameters().at(0);
  const double omegasq = sqrt(omega);
  const double Hnx = my_system->Hermite_n(nx,omegasq*x);
  const double Hny = my_system->Hermite_n(ny,omegasq*y);
  
  //cout << Hnx*Hny*exp(-omega*(x*x+y*y)*0.5)<<endl;
  return Hnx*Hny*exp(-omega*(x*x+y*y)*0.5);
}


double TrialSlater::GradPhi(int pos, int d)
{
  const int nP       = my_system->get_nParticles();
  const double omega = my_system->get_parameters()[0];
  const int orbitals = my_system->get_orbitals();
  const double x = my_system->get_particle()[pos]->get_position()[d];
  int orbital = 0;
  double Grad = 0.0;
  for (int shell = 0 ; shell <= orbitals ; shell++){
    int nx = shell; int ny = 0;
    for (int state = 0 ; state <= shell ; state++){
      if(d==0){orbital=nx;}
      else{orbital=ny;}
      switch (orbital)
      {
        case 0:
          Grad += -omega*x;
          break;
        case 1:
          Grad += (1-omega*x*x)*2;
          break;
        case 2:
          Grad += (4 + omega*x - 2*omega*x*x)*2;
          break;
        case 3:
          Grad += (6*x*x + 3*omega*x*x - 2*omega*x*x*x*x - 3)*4;
          break;
      }
      nx -= 1; ny += 1;
    }
  }
  return Grad;
}

double TrialSlater::GradJas(int k, int d)
{
  const int nP         = my_system->get_nParticles();
  const double beta    = my_system->get_parameters().at(1);
  const double a       = my_system->get_parameters().at(2);
  Eigen::MatrixXd r_ij = my_system->get_r_ij();
  double Grad = 0.0;
  const double xk = my_system->get_particle().at(k)->get_position().at(d);
  for (int j = 0 ; j<nP ; j++){
    if(j!=k){
      const double xj = my_system->get_particle().at(j)->get_position().at(d);
      Grad += (xk-xj)/(r_ij(k,j)*(1+beta*r_ij(k,j))*(1+beta*r_ij(k,j)));
    }
  }
  return Grad*a;
}

double TrialSlater::LapPhi(int pos,int orbital)
{
  double Lap = 0.0;
  const double omega = my_system->get_parameters().at(0);
  for (int xi = 0 ; xi < 2 ; xi++){
    const double x = my_system->get_particle()[pos]->get_position()[xi];
    switch (orbital)
    {
      case 0:
        Lap += (x*x - 1);
        break;
      case 1:
        Lap += (omega*x*x*x-3*x)*2;
        break;
      case 2:
        Lap += (2*omega*x*x*x + 1 - 6*x);
        break;
      case 3:
        Lap += (48*x/omega + 36*x - 56*x*x*x - 12*omega*x*x*x 
                + 8*omega*x*x*x*x*x + 12*x);
        break;
    }
  }
  return Lap*omega;
}

double TrialSlater::LapJas()
{
  const int nP         = my_system->get_nParticles();
  const double beta    = my_system->get_parameters().at(1);
  const double a       = my_system->get_parameters().at(2);
  Eigen::MatrixXd r_ij = my_system->get_r_ij();
  double term1= 0.0;
  double term2= 0.0;
  double Lap  = 0.0;

  for (int k=0 ; k<nP ; k++){
    const double xk = my_system->get_particle().at(k)->get_position().at(0);
    const double yk = my_system->get_particle().at(k)->get_position().at(1);
    for (int i=k+1 ; i<nP ; i++){
      const double xi = my_system->get_particle().at(i)->get_position().at(0);
      const double yi = my_system->get_particle().at(i)->get_position().at(1);
      const double r_ki = r_ij(k,i);
      const double jast_ki = 1/((1+beta*r_ki)*(1+beta*r_ki));
      for(int j=k+1 ; j<nP ; j++){
        const double xj = my_system->get_particle().at(j)->get_position().at(0);
        const double yj = my_system->get_particle().at(j)->get_position().at(1);
        //const double r_kj = sqrt((xk-xj)*(xk-xj) + (yk-yj)*(yk-yj));
        const double r_kj = r_ij(k,j); 
        const double jast_kj = 1/((1+beta*r_kj)*(1+beta*r_kj));
        const double rkri_rkrj = (xk-xi)*(xk-xj)+(yk-yi)*(yk-yj);
        term1 += rkri_rkrj*jast_ki*jast_kj/(r_kj*r_ki);
      }
      term2 += (1+2*beta*r_ki)/(r_ki*(1+beta*r_ki)*(1+beta*r_ki)*(1+beta*r_ki));
    }
    term2 *= 2;
    Lap += a*(term1+term2);
  }
  return Lap;
}

double TrialSlater::computeJastrow ()
{
  const int nP      = my_system->get_nParticles();
  const double beta = my_system->get_parameters().at(1);
  const double a    = my_system->get_parameters().at(2);
  double sep = 0.0;
  for (int p1=0 ; p1<nP; p1++){
    const double x1 = my_system->get_particle().at(p1)->get_position().at(0);
    const double y1 = my_system->get_particle().at(p1)->get_position().at(1);
    for (int p2=p1+1 ; p2<nP; p2++){
      const double x2 = my_system->get_particle().at(p2)->get_position().at(0);
      const double y2 = my_system->get_particle().at(p2)->get_position().at(1);
      sep += sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
    }
  }
  return (a*sep)/(1+beta*sep);
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

