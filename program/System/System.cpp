#include "System.h"

using std::cout;
using std::endl;

void System::runMetropolis ()
{
  unsigned  seed;
  bool accepted = false;
  my_sampler  = new Sampler(this, my_File);
  clock::duration d = clock::now() - my_start;
  seed = d.count();
  my_generator.seed(seed);

  for (int cycle = 0 ; cycle < my_nCycles ; cycle++){
    //accepted = metropolis();
    accepted = importanceSampling();
    if (cycle > my_equilibrationFraction * my_nCycles)
      my_sampler->sample(accepted);
  }
  my_sampler->printResults();
}


bool System::metropolis ()
{
  int detSize   = my_nParticles/2;
  int i         = 0;
  double R_C    = 0.0;
  double R_SD   = 0;
  double randX, randY;
  my_spin = 0;
  Eigen::Vector2d old_pos;
  Eigen::VectorXd d_inv (detSize);
  Eigen::VectorXd SD_row_i (detSize);
  Eigen::MatrixXd DMatrix_old_up (detSize,detSize);
  Eigen::MatrixXd DMatrix_old_dn (detSize,detSize);
  

  /* Choosing particle to move and saving old positions */
  std::uniform_int_distribution<int> particle  (0,my_nParticles-1);
  my_elected = particle(my_generator);
  if(my_elected<detSize){my_spin=1;}
  old_pos = my_particles.row(my_elected);
  randX = my_stepLength*(my_uniform(my_generator)-0.5);
  randY = my_stepLength*(my_uniform(my_generator)-0.5);
  
  /* Getting column from inverse SD */
  if (my_spin)
    d_inv = my_DMatrix_up_inv.col(my_elected);
  else
    d_inv = my_DMatrix_dn_inv.col(my_elected-detSize);
  /* Calculating correlation ratio */
  DMatrix_old_up = my_DMatrix_up;
  DMatrix_old_dn = my_DMatrix_dn;
  //double R_SD_old = my_DMatrix_up.determinant()*
  //                 my_DMatrix_dn.determinant();

  R_C = my_waveFunction->computeJastrow();
  my_particles(my_elected,0) += randX;
  my_particles(my_elected,1) += randY;
  update_r_ij(my_elected);
  R_C = my_waveFunction->computeJastrow()/R_C;

  /* Calculating new row for SD */
  for (int n=0;n<=my_orbitals;n++){
    int nx = n; int ny = 0;
    for (int state=0;state<=n;state++){
      const double phi = my_waveFunction->Phi(my_elected,nx,ny);
      R_SD += phi*d_inv(i);
      SD_row_i(i) = phi;
      i++; nx--; ny++;
    }
  }
  //if(my_spin)
  //  my_DMatrix_up.row(my_elected) = SD_row_i;
  //else
  //  my_DMatrix_dn.row(my_elected-detSize) = SD_row_i;
  

  //double R_SD_new = my_DMatrix_up.determinant()*
  //                  my_DMatrix_dn.determinant();

  //R_SD = (R_SD_new*R_SD_new)/(R_SD_old*R_SD_old);
  R_SD*= R_SD;
  R_C *= R_C;
  //R_C = 1;
  const double RATIO = fabs(R_SD*R_C);

  if (RATIO < my_uniform(my_generator)){
    rejects += 1;
    my_particles(my_elected,0) -= randX;
    my_particles(my_elected,1) -= randY;
    update_r_ij(my_elected);
    //my_DMatrix_up = DMatrix_old_up;
    //my_DMatrix_dn = DMatrix_old_dn;
    return false;
  }
  else{
    /* Updating SD matrices after accepting new step */
    //if(my_spin){
    //  update_inverse(my_DMatrix_up_inv,my_DMatrix_up,
    //                 SD_row_i, my_elected, R_SD);
    //  my_DMatrix_up.row(my_elected) = SD_row_i;
    //}
    //else{
    //  update_inverse(my_DMatrix_dn_inv,my_DMatrix_dn,
    //                 SD_row_i, my_elected-detSize, R_SD);
    //  my_DMatrix_dn.row(my_elected-detSize) = SD_row_i;
    //}
    //cout << my_DMatrix_up_inv << endl;
    my_DMatrix_up_inv = my_DMatrix_up.inverse();
    my_DMatrix_dn_inv = my_DMatrix_dn.inverse();
    if(my_spin)
      my_DMatrix_up.row(my_elected) = SD_row_i;
    else
      my_DMatrix_dn.row(my_elected-detSize) = SD_row_i;
    return true;
  }
}
void System::update_r_ij (int pos)
{
  const double xk = my_particles(pos,0);
  const double yk = my_particles(pos,1);
  for (int i = 0 ; i < pos; i++){
    const double xi = my_particles(i,0);
    const double yi = my_particles(i,1);
    const double sep= sqrt((xk-xi)*(xk-xi) + (yk-yi)*(yk-yi));
    my_r_ij(i,pos) = sep; 
  }
  for (int i = pos+1; i < my_nParticles ; i++){
    const double xi = my_particles(i,0);
    const double yi = my_particles(i,1);
    const double sep= sqrt((xk-xi)*(xk-xi) + (yk-yi)*(yk-yi));
    my_r_ij(pos,i) = sep;
  }
}
void System::update_inverse(Eigen::MatrixXd& matrix_inv, 
                            Eigen::MatrixXd& matrix,
                            Eigen::VectorXd& new_row,
                            int i, double RATIO)
{
  const int detSize    = my_nParticles/2;
  const Eigen::MatrixXd old_matrix_inv = matrix_inv;
  double element = 0.0;
  double tempTerm= 0.0;
  for (int k=0; k<detSize ; k++){
    const double factor = old_matrix_inv(k,i)/RATIO;

    for (int j=0 ; j<detSize ; j++){
      element = factor;
      if (i != j){
        tempTerm = 0.0;
        for (int l=0 ; l<detSize ; l++)
          tempTerm += new_row(l)*old_matrix_inv(l,j);

        element = old_matrix_inv(k,j) - element*tempTerm;
      }
      matrix_inv(k,j) = element;
      //cout << "mai " <<matrix_inv << "\nhh\n" << "old " <<old_matrix_inv << endl;
    }
  }
  //cout << endl;
}

bool System::importanceSampling()
{
  
  int     i                 = 0;
  int     detSize           = my_nParticles/2;
  double    greensArg       = 0;
  double    R_C             = 0;
  double    R_SD            = 0;
  Eigen::Vector2d RandMove;
  Eigen::Vector2d OldPos;
  Eigen::Vector2d QforceX;
  Eigen::Vector2d QforceY;
  Eigen::Vector2d QforceOld;
  Eigen::Vector2d QforceNew;
  Eigen::Vector2d StepDiff;
  Eigen::VectorXd d_inv (detSize);
  Eigen::VectorXd SD_row_i (detSize);

  std::uniform_int_distribution<int> particle (0,my_nParticles-1);
  std::normal_distribution<double>   my_normal (0,sqrt(my_stepLength));

  my_elected = particle(my_generator);
  my_spin = 0;
  if(my_elected<detSize){my_spin=1;}

  OldPos = my_particles.row(my_elected);

  if (my_spin)
    d_inv = my_DMatrix_up_inv.col(my_elected);
  else
    d_inv = my_DMatrix_dn_inv.col(my_elected-detSize);
  
  
  QforceX = my_waveFunction->computeQuantumForce(my_elected,0); 
  QforceY = my_waveFunction->computeQuantumForce(my_elected,1); 
  QforceOld << QforceX(0)+QforceX(1), QforceY(0)+QforceY(1);
  
  RandMove << (my_normal(my_generator) + 0.5*my_stepLength*QforceOld(0)),
              (my_normal(my_generator) + 0.5*my_stepLength*QforceOld(1));

  R_C = my_waveFunction->computeJastrow();

  my_particles.row(my_elected) += RandMove;
  update_r_ij(my_elected); 
  R_C = my_waveFunction->computeJastrow()/R_C;
  for (int n=0;n<=my_orbitals;n++){
    int nx = n; int ny = 0;
    for (int state=0;state<=n;state++){
      const double phi = my_waveFunction->Phi(my_elected,nx,ny);
      R_SD += phi*d_inv(i);
      SD_row_i(i) = phi;
      i++; nx--; ny++;
    }
  }
  //QforceNew << my_waveFunction->computeQuantumForce(my_elected, 0),
  //             my_waveFunction->computeQuantumForce(my_elected, 1);
  QforceX = my_waveFunction->computeQuantumForce(my_elected,0); 
  QforceY = my_waveFunction->computeQuantumForce(my_elected,1); 
  QforceNew << QforceX(0)/(R_SD) + QforceX(1),
               QforceY(0)/(R_SD) + QforceY(1);

  StepDiff << OldPos(0) - my_particles.row(my_elected)(0),
              OldPos(1) - my_particles.row(my_elected)(1);

  greensArg  = ((QforceOld+QforceNew).dot
               (StepDiff + (QforceOld-QforceNew)*0.25*my_stepLength))*0.5;
  
  const double compared = exp(greensArg) * R_SD*R_SD * R_C*R_C;

  if (compared < my_uniform(my_generator)){
    my_particles.row(my_elected) -= RandMove;
    update_r_ij(my_elected);
    return false;
  }
  else{
    //if(my_spin){
    //  update_inverse(my_DMatrix_up_inv,my_DMatrix_up,
    //                 SD_row_i, my_elected, R_SD);
    //  my_DMatrix_up.row(my_elected) = SD_row_i;
    //}
    //else{
    //  update_inverse(my_DMatrix_dn_inv,my_DMatrix_dn,
    //                 SD_row_i, my_elected-detSize, R_SD);
    //  my_DMatrix_dn.row(my_elected-detSize) = SD_row_i;
    //}
    //cout << my_DMatrix_up_inv << endl;
    my_DMatrix_up_inv = my_DMatrix_up.inverse();
    my_DMatrix_dn_inv = my_DMatrix_dn.inverse();
    if(my_spin)
      my_DMatrix_up.row(my_elected) = SD_row_i;
    else
      my_DMatrix_dn.row(my_elected-detSize) = SD_row_i;
    return true;
  }
  
}

void System::OPTIMIZE()
{
  
  int maxIters = 500;
  int optCycles= 1e5;
  int nP_2     = my_nParticles/2;
  double step  = 0.001;
  for (int i = 0 ; i < maxIters ; i++){
    double Ebar_alpha     = 0;
    double Ebar_beta      = 0;
    double dSD_alpha      = 0;
    double dJas_beta      = 0;
    double cumElocal	    = 0;
    double cumEdSD_alpha  = 0;
    double cumEdJas_beta  = 0;
    double cumdSD_alpha   = 0;
    double cumdJas_beta   = 0;

    my_initialState->setupInitialState(); 
    Eigen::VectorXd d_inv (nP_2);
    set_DMatrix();

    for (int i = 0 ; i < optCycles ; i++){
      dSD_alpha = 0;
      dJas_beta = 0;
      int row   = 0;
      double phi = 0.0;
      double jas = 0.0;
      if (importanceSampling()){
        double localEnergy = get_hamiltonian()->computeLocalEnergy();
        for (int p=0 ; p<nP_2 ; p++){
          row = 0;
          for (int n=0;n<=my_orbitals;n++){
            int nx = n; int ny = 0;
            for (int state=0;state<=n;state++){
              const double phi_up = my_waveFunction->Phi(p,nx,ny);
              const double phi_dn = my_waveFunction->Phi(my_nParticles-p-1,nx,ny);
              const double dphi_up = my_waveFunction->dPhi_alpha(p,nx,ny);
              const double dphi_dn = my_waveFunction->dPhi_alpha(my_nParticles-p-1,nx,ny);
              dSD_alpha += dphi_up*my_DMatrix_up_inv(row,p)+
                           dphi_dn*my_DMatrix_dn_inv(row,nP_2-p-1);
              phi += phi_up*phi_dn;
              row++; nx--; ny++;
            }
          }
        }

        dJas_beta = my_waveFunction->dlnjast_beta ();
        jas  = my_waveFunction->computeJastrow();
        cumElocal     += localEnergy;
        cumEdSD_alpha += localEnergy*dSD_alpha*jas;
        cumEdJas_beta += localEnergy*dJas_beta*phi;
        cumdSD_alpha  += dSD_alpha*jas;
        cumdJas_beta  += dJas_beta*phi;
      }
    }
    Ebar_alpha = 2*(cumEdSD_alpha - cumdSD_alpha*cumElocal/(double)optCycles)
                   /(double)optCycles;

    Ebar_beta  = 2*(cumEdJas_beta - cumdJas_beta*cumElocal/(double)optCycles)
                   /(double)optCycles;
    cout << Ebar_alpha << " | " << Ebar_beta << endl;
    if (fabs(Ebar_alpha) > 1e-3)//{
      my_parameters[2] = my_parameters[2] - step*Ebar_alpha;
      //cout<<Ebar_alpha << " | ";}
    if (fabs(Ebar_beta) > 1e-3)//{
      my_parameters[1] = my_parameters[1] - step*Ebar_beta;
      //cout << Ebar_beta;}
    //cout << endl;
  }
  std::cout << "alpha:  " << my_parameters[2];
  std::cout << "  beta:  " << my_parameters[1] << std::endl;
  
}

void System::set_DMatrix()
{
  int nP = my_nParticles;
  Eigen::MatrixXd DMatrix1 (nP/2,nP/2);
  Eigen::MatrixXd DMatrix2 (nP/2,nP/2);
  Eigen::MatrixXd r_ij     (nP,nP);
  DMatrix1.setZero();
  DMatrix2.setZero();
  r_ij.setZero();
  //r_ij = Eigen::MatrixXd::Zero(nP,nP);
  int row = 0;
  for (int n=0;n<=my_orbitals;n++){
    int nx = n; int ny = 0;
    for (int i=0;i<n+1;i++){
      for (int j=0;j<my_nParticles/2;j++){
        const double phi1 = my_waveFunction->Phi(j,nx,ny);
        const double phi2 = my_waveFunction->Phi(j+nP/2,nx,ny);
        DMatrix1(row,j) = phi1;
        DMatrix2(row,j) = phi2;
      }
      row+=1; nx -= 1; ny += 1;
    }
  }
  my_DMatrix_up     = DMatrix1;
  my_DMatrix_dn     = DMatrix2;
  my_DMatrix_up_inv = DMatrix1.inverse();
  my_DMatrix_dn_inv = DMatrix2.inverse();
  //cout << my_DMatrix_up_inv << endl;
  //cout << my_DMatrix_dn_inv << endl;
  //cout << my_DMatrix_up << endl;
  //cout << my_DMatrix_dn << endl;

  for (int i = 0 ; i < nP ; i++){
    const double xi = my_particles(i,0);
    const double yi = my_particles(i,1);
    for (int j = i+1 ; j < nP ; j++){
      const double xj = my_particles(j,0);
      const double yj = my_particles(j,1);
      r_ij(i,j) = sqrt((xi-xj)*(xi-xj) + (yi-yj)*(yi-yj));
    }
  }
  my_r_ij = r_ij;
}
