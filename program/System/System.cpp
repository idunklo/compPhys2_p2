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
    accepted = metropolis();
    //accepted = importanceSampling();
    if (cycle > my_equilibrationFraction * my_nCycles)
      my_sampler->sample(accepted);
  }
  my_sampler->printResults();
}


bool System::metropolis ()
{
  int elected   = 0;
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
  elected = particle(my_generator);
  if(elected<detSize){my_spin=1;}
  old_pos = my_particles.row(elected);
  randX = my_stepLength*(my_uniform(my_generator)-0.5);
  randY = my_stepLength*(my_uniform(my_generator)-0.5);
  
  /* Getting column from inverse SD */
  if (my_spin)
    d_inv = my_DMatrix_up_inv.col(elected);
  else
    d_inv = my_DMatrix_dn_inv.col(elected-detSize);
  /* Calculating correlation ratio */
  DMatrix_old_up = my_DMatrix_up;
  DMatrix_old_dn = my_DMatrix_dn;
  //double R_SD_old = my_DMatrix_up.determinant()*
  //                 my_DMatrix_dn.determinant();

  R_C = my_waveFunction->computeJastrow();
  my_particles(elected,0) += randX;
  my_particles(elected,1) += randY;
  update_r_ij(elected);
  R_C = my_waveFunction->computeJastrow()/R_C;

  /* Calculating new row for SD */
  for (int n=0;n<=my_orbitals;n++){
    int nx = n; int ny = 0;
    for (int state=0;state<=n;state++){
      const double phi = my_waveFunction->Phi(elected,nx,ny);
      R_SD += phi*d_inv(i);
      SD_row_i(i) = phi;
      i++; nx--; ny++;
    }
  }
  if(my_spin)
    my_DMatrix_up.row(elected) = SD_row_i;
  else
    my_DMatrix_dn.row(elected-detSize) = SD_row_i;
  

  //double R_SD_new = my_DMatrix_up.determinant()*
  //                  my_DMatrix_dn.determinant();

  //R_SD = (R_SD_new*R_SD_new)/(R_SD_old*R_SD_old);
  R_SD*= R_SD;
  R_C *= R_C;
  //R_C = 1;
  const double RATIO = fabs(R_SD*R_C);

  if (RATIO < my_uniform(my_generator)){
    rejects += 1;
    my_particles(elected,0) -= randX;
    my_particles(elected,1) -= randY;
    update_r_ij(elected);
    my_DMatrix_up = DMatrix_old_up;
    my_DMatrix_dn = DMatrix_old_dn;
    return false;
  }
  else{
    /* Updating SD matrices after accepting new step */
    //if(my_spin){
    //  update_inverse(my_DMatrix_up_inv,my_DMatrix_up,
    //                 SD_row_i, elected, RATIO);
    //  my_DMatrix_up.row(elected) = SD_row_i;
    //}
    //else{
    //  update_inverse(my_DMatrix_dn_inv,my_DMatrix_dn,
    //                 SD_row_i, elected-detSize, RATIO);
    //  my_DMatrix_dn.row(elected-detSize) = SD_row_i;
    //}
    //cout << my_DMatrix_up_inv << endl;
    my_DMatrix_up_inv = my_DMatrix_up.inverse();
    my_DMatrix_dn_inv = my_DMatrix_dn.inverse();
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
  for (int k=0; k<detSize ; k++){
    const double factor = old_matrix_inv(k,i)/RATIO;

    for (int j=0 ; j<detSize ; j++){
      element = 0.0;
      if (i != j){
        element = old_matrix_inv(k,j);

        for (int l=0 ; l<detSize ; l++){
          element -= factor*new_row(l)*old_matrix_inv(l,j);
        }
      }
      else
        element = factor;
        //element += matrix(i,l)*old_matrix_inv(l,j);
      matrix_inv(k,j) = element;
      }
    }
  
}

bool System::importanceSampling()
{
  /*
  int     elected           = 0;
  int     chosenDimension          = 0;
  double    stepDifference         = 0;
  double    randomMove             = 0;
  double    oldPosition            = 0;
  double    greensExp              = 0;
  double    waveFunctionOld        = 0;
  double    waveFunctionNew        = 0;
  double    quantumForceOld        = 0;
  double    quantumForceNew        = 0;
  double    waveFunctionsCompared  = 0;
  double    greensFunctionCompared = 0;
  double    compared               = 0;

  std::uniform_int_distribution<int> particle (0,my_nParticles-1);
  std::uniform_int_distribution<int> dimension  (0,my_nDimensions-1);

  elected  = particle(my_generator);
  chosenDimension = dimension(my_generator);

  oldPosition = my_particles[elected]->get_position()[chosenDimension];

  waveFunctionOld = my_waveFunction->evaluate();

  quantumForceOld = my_waveFunction->
           computeQuantumForce(elected, chosenDimension);

  randomMove  = (my_normal(my_generator)) * sqrt(my_stepLength)
              +
              quantumForceOld * my_stepLength * 0.5;

  my_particles[elected]->changePosition(chosenDimension, randomMove);

  waveFunctionNew = my_waveFunction->evaluate();

  quantumForceNew = my_waveFunction->
        computeQuantumForce(elected, chosenDimension);

  stepDifference  = my_particles[elected]->get_position()[chosenDimension] - oldPosition;

  greensExp  = (0.5*(quantumForceOld + quantumForceNew)*((quantumForceOld - quantumForceNew)*
    0.25*my_stepLength - stepDifference));

  greensFunctionCompared = exp(greensExp);

  waveFunctionsCompared  = (waveFunctionNew*waveFunctionNew)/
         (waveFunctionOld*waveFunctionOld);

  compared = greensFunctionCompared * waveFunctionsCompared;

  if (compared < my_uniform(my_generator)){
    my_particles[elected]->get_position()[chosenDimension] = oldPosition;
    return false;
  }
  return true;
  */
}

void System::OPTIMIZE()
{
  /*
  int maxIters = 1e2;
  int optCycles= 1e4;
  double step  = 0.01;
  for (int i = 0 ; i < maxIters ; i++){
    double Ebar_alpha       = 0;
    double Ebar_beta        = 0;
    double psiBar_alpha     = 0;
    double psiBar_beta      = 0;
    double cumElocal	      = 0;
    double cumEPsi_alpha    = 0;
    double cumEPsi_beta     = 0;
    double cumPsiBar_alpha  = 0;
    double cumPsiBar_beta   = 0;
    my_particles.clear();
    my_initialState->setupInitialState(); 
    //std::vector <double> temp1 = {1,1};
    //std::vector <double> temp2 = {2,2};
    //my_particles[0]->set_position(temp1);
    //my_particles[1]->set_position(temp2);

    //for (int p = 0 ; p < my_nParticles ; p++){
    //  for (int d = 0 ; d<my_nDimensions ; d++)
    //    std::cout << my_particles[p]->get_position()[d] << " ";
    //  std::cout << std::endl; 
    //}
    for (int i = 0 ; i < optCycles ; i++){
      psiBar_alpha     = 0;
      psiBar_beta      = 0;
      if (metropolis()){
        my_waveFunction->computePsiBars(psiBar_alpha,psiBar_beta);
        
        double localEnergy = get_hamiltonian()->computeLocalEnergy();
        cumElocal       += localEnergy;
        cumEPsi_alpha   += localEnergy*psiBar_alpha;
        cumEPsi_beta    += localEnergy*psiBar_beta;
        cumPsiBar_alpha += psiBar_alpha;
        cumPsiBar_beta  += psiBar_beta;
        //std::cout << psiBar_alpha << "  " << psiBar_beta << 
        //  "  " << localEnergy << "  " << localEnergy*psiBar_alpha <<
        //  "  " << localEnergy*psiBar_beta << std::endl; 
      }
    }
    //std::cout << "(" << cumEPsi_alpha << " - " << cumPsiBar_alpha*cumElocal/optCycles;
    //std::cout << ")/ " << optCycles << std::endl; 
    Ebar_alpha = 2*(cumEPsi_alpha - cumPsiBar_alpha*cumElocal/(double)optCycles)
                   /(double)optCycles;

    Ebar_beta  = 2*(cumEPsi_beta - cumPsiBar_beta*cumElocal/(double)optCycles)
                   /(double)optCycles;
    if (abs(Ebar_alpha) > 1e-10){
      my_parameters[0] = my_parameters[0] - step*Ebar_alpha;
      std::cout<<Ebar_alpha << " | ";}
    if (abs(Ebar_beta) > 1e-10){
      my_parameters[1] = my_parameters[1] - step*Ebar_beta;
      std::cout << Ebar_beta;}
    std::cout << std::endl;
  }
  std::cout << "alpha:  " << my_parameters[0];
  std::cout << "  beta:  " << my_parameters[1] << std::endl;
  */
}


double System::Hermite_n ( int n, double x)
{
  /*
  double Hi = 0.0;
  double Hb = 1.0;
  double Hc = 2*x;
  if (n==0)
    return Hb;
  else if (n == 1)
    return Hc; 
  else
  {
    for (int i = 2 ; i <= n ; i++)
    {
      Hi = 2*x*Hc - 2*(i-1)*Hb;
      Hb = Hc;
      Hc = Hi;
    }
    return Hi; 
  }
  */
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
