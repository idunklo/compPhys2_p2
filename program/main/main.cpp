#include <stdlib.h>
#include "mpi.h"
#include "../System/System.h"
#include "../System/Particle.h"
#include "../Sampler/Sampler.h"
#include "../Hamiltonian/Hamiltonian.h"
#include "../WaveFunction/WaveFunction.h"
#include "../InitialState/InitialState.h"
#include "../Hamiltonian/HarmonicOscillator.h"
#include "../WaveFunction/TrialWaveFunction.h"
#include "../InitialState/RandomUniform.h"

using std::cout;

void importanceSampling(bool,int,int,int,int,int,double,double,double,std::vector<double>);
void metropolis(bool,int,int,int,int,int,double,double,double,std::vector<double>);
 
int main (int argc,char* argv[]){

  int my_rank = 0;
  int num_procs;// = 1;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size (MPI_COMM_WORLD,&num_procs);

  bool    File            = false;
  int     nDimensions     = 2;
  int     nParticles      = 2;
  int 	  nCycles	        = (int) 1e6;
  double  omega		        = 1.0;//0.5;
  double  alpha		        = 1.0;//0.95455;
  double  beta		        = 0.3;//0.50905;
  double  a               = 1.0;
  double  equilibration	  = 0.1;
  double  derivativeStep  = 0.001;

  std::vector<double> parameters {alpha, beta, omega, a};

  int chosenOne = 1;

  double  stepLength = 1.0;
  switch (chosenOne)
  {
    case 0:
      stepLength = 0.01;
      importanceSampling(File,nCycles,nParticles,nDimensions,
          my_rank, num_procs,
          stepLength,equilibration,derivativeStep,
          parameters);
      break;
    case 1:
      stepLength = 2.0;
      metropolis(File,nCycles,nParticles,nDimensions,
          my_rank, num_procs,
          stepLength,equilibration,derivativeStep,
          parameters);
      break;
  }

  MPI_Finalize();

  return 0;
}

void importanceSampling(bool File,
			int nCycles,
			int nParticles,
			int nDimensions,
      int my_rank,
      int num_procs,
			double stepLength,
			double equilibration,
			double derivativeStep,
			std::vector<double>parameters)
{
  
  if (my_rank>0)
    File=false;

  nCycles = (double) nCycles;
  System* system = new System(File);

  system->set_parameters		        (parameters);
  system->set_stepLength		        (stepLength);
  system->set_equilibrationFraction (equilibration);
  system->set_derivativeStep		    (derivativeStep);
  system->set_nCycles			          (nCycles);
  system->set_rank                  (my_rank);
  system->set_procs                 (num_procs);
  system->set_InitialState  (new RandomUniform		  (system, nDimensions, nParticles));
  system->set_Hamiltonian	  (new HarmonicOscillator (system));
  system->set_WaveFunction  (new TrialWaveFunction  (system));
  if (my_rank==-1){
    system->set_Timer		      (new Timer			        (system));
    cout << "Starting timer...\n";
    system->get_timer()->startTimer  ();
  }
  system->runImportanceSampling    ();
  if (my_rank==-1){
    system->get_timer()->stopTimer   ();
    cout << "----------------------------------\n";
    cout << "              Timers              \n";
    cout << "       Importance Sampling        \n";
    cout << "----------------------------------\n";

    printf("Seconds     : %i \n",system->get_timer()->elapsedTimeSeconds());
    printf("Milliseconds: %i \n",system->get_timer()->elapsedTimeMilli());
    printf("Microseconds: %i \n",system->get_timer()->elapsedTimeMicro());
  }
}

void metropolis(bool File,
		int nCycles,
		int nParticles,
		int nDimensions,
    int my_rank,
    int num_procs,
		double stepLength,
		double equilibration,
		double derivativeStep,
		std::vector<double>parameters)
{
// Brute forece Metropolis 

  System* system = new System(File);
   
  system->set_parameters		        (parameters);
  system->set_stepLength		        (stepLength);
  system->set_equilibrationFraction	(equilibration);
  system->set_derivativeStep		    (derivativeStep);
  system->set_nCycles               (nCycles);
  system->set_rank                  (my_rank);
  system->set_procs                 (num_procs);


  system->set_InitialState	(new RandomUniform	 (system, nDimensions, nParticles));
  system->set_Hamiltonian	(new HarmonicOscillator  (system));
  system->set_WaveFunction	(new TrialWaveFunction   (system));
  //system->OPTIMIZE();

  if (my_rank==-1){
    system->set_Timer		(new Timer               (system));
    cout << "Starting timer...\n";
    system->get_timer()->startTimer	();
  }
  system->runMetropolis		();
  if (my_rank==-1){
    system->get_timer()->stopTimer	();
    cout << "----------------------------------\n";
    cout << "              Timers              \n";
    cout << "            Brute force           \n";
    cout << "----------------------------------\n";
    printf("Seconds     : %i \n",system->get_timer()->elapsedTimeSeconds());
    printf("Milliseconds: %i \n",system->get_timer()->elapsedTimeMilli());
    printf("Microseconds: %i \n",system->get_timer()->elapsedTimeMicro());
    cout << "----------------------------------\n\n\n";
  }
}

