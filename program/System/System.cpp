#include "System.h"

using std::cout;
using std::endl;

void System::runMetropolis ()
{
  /*
  unsigned  seed;
  bool accepted = false;
  my_sampler  = new Sampler(this, my_File);
  clock::duration d = clock::now() - my_start;
  seed = d.count();
  my_generator.seed(seed);

  for (int cycle = 0 ; cycle < my_nCycles ; cycle++){
    accepted = metropolis();
    if (cycle > my_equilibrationFraction * my_nCycles)
    {
      my_sampler->sample(accepted);
    }
  }
  my_sampler->printResults();
  */
}

void System::runImportanceSampling ()
{
  /*
  unsigned  seed;
  bool accepted = false;
  my_sampler  = new Sampler(this,my_File);
  clock::duration d = clock::now() - my_start;
  seed = d.count();
  my_generator.seed(seed);

  for (int cycle = 0 ; cycle < my_nCycles ; cycle++){
    accepted = importanceSampling();
    if (cycle > my_equilibrationFraction * my_nCycles)
    {
      my_sampler->sample(accepted);
    }
  }
  my_sampler->printResults();
  */
}

bool System::metropolis ()
{
  /*
  int     chosenParticle        = 0;
  int     chosenDimension       = 0;
  double  waveFunctionOld       = 0;
  double  waveFunctionNew       = 0;
  double  waveFunctionsCompared = 0;
  double  randomMove            = 0;

  std::uniform_int_distribution<int> particle  (0,my_nParticles-1);
  std::uniform_int_distribution<int> dimension (0,my_nDimensions-1);

  randomMove      = my_stepLength*(my_uniform(my_generator)-0.5);
       
  chosenParticle  = particle(my_generator);
  chosenDimension = dimension(my_generator);
  
  waveFunctionOld = my_waveFunction->evaluate();

  my_particles[chosenParticle]->changePosition(chosenDimension, randomMove);

  waveFunctionNew = my_waveFunction->evaluate();

  waveFunctionsCompared = (waveFunctionNew*waveFunctionNew)/
                          (waveFunctionOld*waveFunctionOld);

  if (waveFunctionsCompared < my_uniform(my_generator)){
    my_particles[chosenParticle]->changePosition(chosenDimension, -randomMove);
    return false;
  }
  return true;
  */
}

bool System::importanceSampling()
{
  /*
  int     chosenParticle           = 0;
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

  chosenParticle  = particle(my_generator);
  chosenDimension = dimension(my_generator);

  oldPosition = my_particles[chosenParticle]->get_position()[chosenDimension];

  waveFunctionOld = my_waveFunction->evaluate();

  quantumForceOld = my_waveFunction->
           computeQuantumForce(chosenParticle, chosenDimension);

  randomMove  = (my_normal(my_generator)) * sqrt(my_stepLength)
              +
              quantumForceOld * my_stepLength * 0.5;

  my_particles[chosenParticle]->changePosition(chosenDimension, randomMove);

  waveFunctionNew = my_waveFunction->evaluate();

  quantumForceNew = my_waveFunction->
        computeQuantumForce(chosenParticle, chosenDimension);

  stepDifference  = my_particles[chosenParticle]->get_position()[chosenDimension] - oldPosition;

  greensExp  = (0.5*(quantumForceOld + quantumForceNew)*((quantumForceOld - quantumForceNew)*
    0.25*my_stepLength - stepDifference));

  greensFunctionCompared = exp(greensExp);

  waveFunctionsCompared  = (waveFunctionNew*waveFunctionNew)/
         (waveFunctionOld*waveFunctionOld);

  compared = greensFunctionCompared * waveFunctionsCompared;

  if (compared < my_uniform(my_generator)){
    my_particles[chosenParticle]->get_position()[chosenDimension] = oldPosition;
    return false;
  }
  return true;
  */
}

void System::OPTIMIZE()
{
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
        
        double localEnergy = get_hamiltonian()->computeNumLocalEnergy();
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
}

double System::Hermite_n ( int n, double x)
{
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
}

void System::set_DMatrix(Eigen::MatrixXd& DMatrix,
                         Eigen::MatrixXd& DMatrix_inv,
                         int levels)
{
  my_DMatrix     = DMatrix;
  my_DMatrix_inv = DMatrix_inv;
  int row = 0;
  for (int n=0;n<=levels;n++){
    int nx = n;
    int ny = 0;
    for (int i=0;i<n+1;i++){
      for (int j=0;j<my_nParticles;j++){
        const double phi = my_waveFunction->evaluate(j,nx,ny);
        my_DMatrix(row,j)   = phi;
        my_DMatrix(row+1,j) = phi;
      }
      row+=2;
      nx -= 1;
      ny += 1;
    }
  }
  my_DMatrix_inv = my_DMatrix.inverse();
}
