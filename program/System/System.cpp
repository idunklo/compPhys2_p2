#include "System.h"

using std::cout;
using std::endl;

System::System (bool File)
{
  my_File = File;
}

bool System::metropolis ()
{
  int     chosenParticle    = 0;
  int     chosenDimension   = 0;
  double  waveFunctionOld   = 0;
  double  waveFunctionNew   = 0;
  double  waveFunctionsCompared = 0;
  double  randomMove      = 0;

  std::uniform_int_distribution<int> particle  (0,my_nParticles-1);
  std::uniform_int_distribution<int> dimension (0,my_nDimensions-1);

  randomMove    = my_normal(my_generator)*my_stepLength;//*(my_uniform(my_generator)-0.5);
  chosenParticle  = particle(my_generator);
  chosenDimension = dimension(my_generator);

  waveFunctionOld = my_waveFunction->evaluate();

  my_particles[chosenParticle]->changePosition(chosenDimension, randomMove);

  waveFunctionNew = my_waveFunction->evaluate();

  waveFunctionsCompared = (waveFunctionNew*waveFunctionNew)/
             (waveFunctionOld*waveFunctionOld);


  if (waveFunctionsCompared < 1.0){
    if (waveFunctionsCompared < my_uniform(my_generator)){
      my_particles[chosenParticle]->changePosition(chosenDimension, -randomMove);
      return false;
    }
  }
  return true;
}

bool System::importanceSampling()
{
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
      randomMove  = (my_normal(my_generator)) * sqrt(my_stepLength)
                   + quantumForceOld[p] * my_stepLength * 0.5;
      my_particles[p]->changePosition(d,randomMove);
    }
    quantumForceOld[p] = my_waveFunction->computeQuantumForce(p);
  }

  waveFunctionNew = my_waveFunction->evaluate();

  for (int p = 0 ; p < nP ; p++){
    for (int d = 0 ; d < nD ; d++){
      stepDifference  = my_particles[p]->get_position()[d] - oldPositions[p]->get_position()[d];

      quantumForceNew = my_waveFunction->computeQuantumForce(p);

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
}

void System::runMetropolis ()
{
  //ProgressBar *bar = new ProgressBar(my_nCycles);
  //bar->SetFrequencyUpdate(9);

  unsigned  seed;
  bool accepted = false;
  my_sampler  = new Sampler(this, my_File);
  clock::duration d = clock::now() - my_start;
  seed = d.count();
  my_generator.seed(seed);

  for (int cycle = 0 ; cycle < my_nCycles ; cycle++){
    accepted = metropolis();
    //bar->Progressed(cycle);
    if (cycle > my_equilibrationFraction * my_nCycles)
    {
      my_sampler->sample(accepted);
    }
  }
  my_sampler->printResults();
}


void System::runImportanceSampling ()
{
  //ProgressBar *bar = new ProgressBar(my_nCycles);
  //bar->SetFrequencyUpdate(50);
  //randomMove.reserve    (my_nDimensions);
  //quantumForceOld.reserve (my_nDimensions);
  //oldPosition.reserve   (my_nDimensions);
  unsigned  seed;
  bool accepted = false;
  my_sampler  = new Sampler(this,my_File);
  clock::duration d = clock::now() - my_start;
  seed = d.count();
  my_generator.seed(seed);

  for (int cycle = 0 ; cycle < my_nCycles ; cycle++){
    accepted = importanceSampling();
    //bar->Progressed(cycle);
    if (cycle > my_equilibrationFraction * my_nCycles)
    {
      my_sampler->sample(accepted);
    }
  }
  my_sampler->printResults();
}

void System::set_nDimensions (int nDimensions)
{ my_nDimensions = nDimensions; }

void System::set_nParticles (int nParticles)
{ my_nParticles = nParticles; }

void System::set_nCycles (int nCycles)
{ my_nCycles = nCycles; }

void System::set_Comm (MPI_Comm comm)
{ my_comm = comm; }

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

void System::add_particle (Particle* particle)
{ my_particles.push_back(particle); }
