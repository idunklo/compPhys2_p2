#include "HarmonicOscillator.h"

HarmonicOscillator::HarmonicOscillator (System* system):
  Hamiltonian(system)
{}

double HarmonicOscillator::computeLocalEnergy()
{
  //cout << HOLap() << " " << HOExt() << " " << HRep() << endl;
  return HOLap() + HOExt();// + HRep();
}

double HarmonicOscillator::HOLap()
{
  const int orbitals = my_system->get_orbitals();
  const int spin     = my_system->get_spin();
  const int nP       = my_system->get_nParticles();
  const int nD       = my_system->get_nDimensions();
  int startPos       = 0;
  int stopPos        = 0;
  int row            = 0;
  int orbital        = 0;
  double Grads       = 0.0;
  double SD_LAP      = 0.0;
  double SD_LAP_new  = 0.0;
  Eigen::MatrixXd D_inv(nP/2,nP/2);;
  D_inv.setZero();
  if(spin){
    SD_LAP = my_system->get_SDLap_old_dn();
    startPos = 0;
    stopPos  = nP/2;
    D_inv = my_system->get_DMatrix_up_inv();
  }
  else{
    SD_LAP = my_system->get_SDLap_old_up();
    D_inv = my_system->get_DMatrix_dn_inv();
    startPos = nP/2;
    stopPos  = nP;
  }
  //cout << D_inv << endl;
  for (int shell = 0 ; shell <= orbitals ; shell++){
    int nx = shell; int ny = 0;
    for (int state = 0 ; state <= shell ; state++){
      for (int r = startPos; r < stopPos; r++){
        //for (int part = startPos; part < stopPos ; part++){
          SD_LAP_new += my_system->get_waveFunction()->LapPhi(r,nx)
                       *D_inv(row,r-startPos);
          SD_LAP_new += my_system->get_waveFunction()->LapPhi(r,ny)
                       *D_inv(row,r-startPos);
          if(my_system->get_waveFunction()->LapPhi(r,nx) > 7){

          cout << "LaPh: " << my_system->get_waveFunction()->LapPhi(r,nx) << "  "
            << "D_inv: "<< D_inv(row,r-startPos) << endl;
          cout << "shell: " << shell;
          cout << "  state: " << state;
          cout << "  nx: " << nx;
          cout << "  ny: " << ny;
          cout << "  phi_"<< row <<"("<<r<<")"<<endl;
        }
      }
      row++; nx--;ny++;
    }
  }
  
  //cout << "Sd " << SD_LAP_new << endl;
  const double JastrowLap= my_system->get_waveFunction()->LapJas();
  //cout << "Jas " << JastrowLap << endl;
  double totGrad = 0.0;
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
  //if (SD_LAP_new > 400)
  //  cout << "IIIK " << SD_LAP << "  " << SD_LAP_new << endl;
  SD_LAP += SD_LAP_new; 
  //cout << SD_LAP << "  " << JastrowLap << "  " << 2*totGrad << endl;
  return -0.5*(SD_LAP);// + JastrowLap + 2*totGrad);
}

double HarmonicOscillator::HOExt()
{
  const int nP       = my_system->get_nParticles();
  const double omega = my_system->get_parameters()[0];
        double r_all = 0.0;
  for (int p = 0 ; p < nP ; p++){
    const double x = my_system->get_particle().at(p)->get_position().at(0);
    const double y = my_system->get_particle().at(p)->get_position().at(1);
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
      rep += 1/r_ij(p1,p2);
    }
  }
  return rep;
}
