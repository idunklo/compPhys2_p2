#include "System.h"

using std::cout;
using std::endl;

System::System (bool File)
{
  my_File = File;
}
System::~System ()
{
  for (int i = 0 ; i < my_particles.size() ; i++)
    delete my_particles.at(i);
  delete my_sampler;
}
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
    if (cycle > my_equilibrationFraction * my_nCycles)
    {
      my_sampler->sample(accepted);
    }
  }
  my_sampler->printResults();
}

void System::runImportanceSampling ()
{
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
}

bool System::metropolis ()
{
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
}

bool System::importanceSampling()
{
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

void System::set_nDimensions (int nDimensions)
{ my_nDimensions = nDimensions; }

void System::set_nParticles (int nParticles)
{ my_nParticles = nParticles; }

void System::set_nCycles (int nCycles)
{ my_nCycles = nCycles; }

void System::set_rank (int rank)
{ my_rank = rank; }

void System::set_procs (int procs)
{ num_procs = procs; }

void System::set_stepLength (double stepLength)
{ my_stepLength = stepLength; }

void System::set_derivativeStep (double h)
{ my_derivativeStep = h; my_derivativeStep2 = h*h;}

void System::set_equilibrationFraction (double equilibraFraction)
{ my_equilibrationFraction = equilibraFraction; }

void System::set_Hamiltonian (Hamiltonian* hamiltonian)
{ my_hamiltonian = hamiltonian; }

void System::set_WaveFunction (WaveFunction* waveFunction)
{ my_waveFunction = waveFunction; }

void System::set_InitialState (InitialState* initialState)
{ my_initialState = initialState; }

void System::set_Timer (Timer* timer)
{ my_timer = timer; }

void System::set_parameters (std::vector<double> parameters)
{ my_parameters = parameters; }

void System::add_particle ()
{ my_particles.push_back(new Particle()); }
  /*
  int       nD                      = my_nDimensions;
  int       nP                      = my_nParticles;
  double    stepDifference          = 0;
  double    randomMove              = 0;
  double    greensExp               = 0;
  double    waveFunctionOld         = 0;
  double    waveFunctionNew         = 0;
  double    quantumForceNew         = 0;
  double    waveFunctionsCompared   = 0;
  double    greensFunctionCompared  = 0;
  double    compared                = 0;
  

  std::vector<class Particle*> oldPositions = std::vector<class Particle*>();
  std::vector<double> tempPos (nD);
  std::vector<double> quantumForceOld (nP);

  for (int p = 0 ; p < nP ; p++){
    for (int d = 0 ; d < nD ; d++){
      tempPos[d] = my_particles[p]->get_position()[d];
    }
    oldPositions.push_back(new Particle ());
    oldPositions.at(p)->set_nDimensions (nD);
    oldPositions.at(p)->set_position    (tempPos);
  }
    
  waveFunctionOld = my_waveFunction->evaluate();

  for (int p = 0 ; p < nP ; p++){
    for (int d = 0 ; d < nD ; d++){
      quantumForceOld[p] = my_waveFunction->computeQuantumForce(p, d);
      randomMove  = (my_normal(my_generator)) * sqrt(my_stepLength)
                   + quantumForceOld[p] * my_stepLength * 0.5;
      my_particles[p]->changePosition(d,randomMove);
    }
  }

  waveFunctionNew = my_waveFunction->evaluate();

  for (int p = 0 ; p < nP ; p++){
    for (int d = 0 ; d < nD ; d++){
      stepDifference  = my_particles[p]->get_position()[d] - oldPositions[p]->get_position()[d];

      quantumForceNew = my_waveFunction->computeQuantumForce(p, d);

      greensExp  += (0.5*(quantumForceOld[p] + quantumForceNew)*
                    ((quantumForceOld[p] - quantumForceNew)*
                    0.25*my_stepLength - stepDifference));
    }
  }
  greensFunctionCompared = exp(greensExp);

  waveFunctionsCompared  =  (waveFunctionNew*waveFunctionNew)/
                            (waveFunctionOld*waveFunctionOld);

  compared = greensFunctionCompared * waveFunctionsCompared;

  if (compared < 1.0){
    if (compared < my_uniform(my_generator)){
      my_particles = oldPositions;
      return false;
    }
  }

  return true;
  */
