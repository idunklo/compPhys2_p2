#include "HarmonicOscillator.h"

HarmonicOscillator::HarmonicOscillator (System* system):
  Hamiltonian(system)
{
  int rank = my_system->get_rank();
  std::string name = "hamil_"+std::to_string(rank)+".out";
  //my_hamil.open(name,std::ios::out);
}

double HarmonicOscillator::computeLocalEnergy()
{
  //cout << HOLap()<<" + "<<HOExt()<<" = "<<HOLap() + HOExt() << endl;
  return HOLap() + HOExt() + HRep();
}

double HarmonicOscillator::HOLap()
{
  const int orbitals = my_system->get_orbitals();
  const int spin     = my_system->get_spin();
  const int nP       = my_system->get_nParticles();
  int startPos       = 0;
  int stopPos        = 0;
  int row            = 0;
  double SD_LAP      = 0.0;
  double SD_LAP_new  = 0.0;
  Eigen::MatrixXd D_inv(nP/2,nP/2);
  D_inv.setZero();
  if(spin){
    SD_LAP   = my_system->get_SDLap_old_dn();
    startPos = 0;
    stopPos  = nP/2;
    D_inv    = my_system->get_DMatrix_up_inv();
  }
  else{
    SD_LAP   = my_system->get_SDLap_old_up();
    D_inv    = my_system->get_DMatrix_dn_inv();
    startPos = nP/2;
    stopPos  = nP;
  }
  const double JastrowLap= my_system->get_waveFunction()->LapJas();
  for (int shell = 0 ; shell <= orbitals ; shell++){
    int nx = shell; int ny = 0;
    for (int state = 0 ; state <= shell ; state++){
      for (int r = startPos; r < stopPos; r++){
        SD_LAP_new += my_system->get_waveFunction()->LapPhi(r,nx,ny)
                    *  D_inv(row,r-startPos);
      }
      row++; nx--;ny++;
    }
  }
  
  double totGrad = 0.0;
  //my_hamil << SD_LAP << "  " << JastrowLap << "  " << totGrad << endl;
  for (int i = 0 ; i<nP ; i++){
    totGrad += my_system->get_waveFunction()->GradPhi(i,0)
            *  my_system->get_waveFunction()->GradJas(i,0)
            +  my_system->get_waveFunction()->GradPhi(i,1)
            *  my_system->get_waveFunction()->GradJas(i,1);
  }

  if(spin){
    my_system->set_SDLap_up(SD_LAP_new);
  }
  else{
    my_system->set_SDLap_dn(SD_LAP_new);
  }
  SD_LAP += SD_LAP_new; 
  //if (fabs(SD_LAP) > 200)
  //  cout << D_inv << endl<<endl;
  //my_hamil << SD_LAP << "  " << JastrowLap << "  " << totGrad << endl;
  //cout << SD_LAP << "  " << JastrowLap << "  " << 2*totGrad << endl;
  return -0.5*(SD_LAP + JastrowLap + 2*totGrad);
}

double HarmonicOscillator::HOExt()
{
  const int nP       = my_system->get_nParticles();
  const double omega = my_system->get_parameters()[0];
        double r_all = 0.0;
  for (int p = 0 ; p < nP ; p++){
    const double x = my_system->get_particles()(p,0);
    const double y = my_system->get_particles()(p,1);
    r_all += x*x + y*y; 
  }
  return 0.5*omega*omega*r_all;

}
double HarmonicOscillator::HRep()
{
  const int nP  = my_system->get_nParticles();
  const Eigen::MatrixXd r_ij = my_system->get_r_ij();
  double rep = 0.0;
  for (int p1 = 0 ; p1 < nP ; p1++){
    for (int p2 = p1+1 ; p2 < nP ; p2++){
      rep += 1.0/r_ij(p1,p2);
    }
  }
  return rep;
}
