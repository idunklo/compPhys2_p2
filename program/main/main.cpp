#include <stdlib.h>
#include "mpi.h"
#include "Eigen/Core"
#include "../System/System.h"
#include "../System/Particle.h"
#include "../Sampler/Sampler.h"
#include "../Hamiltonian/Hamiltonian.h"
#include "../WaveFunction/WaveFunction.h"
#include "../InitialState/InitialState.h"
#include "../Hamiltonian/HarmonicOscillator.h"
#include "../WaveFunction/TrialSlater.h"
#include "../InitialState/RandomUniform.h"

using std::cout;

//void importanceSampling(bool,int,int,int,int,int,double,double,double,std::vector<double>);
void metropolis(bool,int,int,int,int,int,double,double,double,std::vector<double>);
 
int main (int argc,char* argv[]){

  int my_rank = 0;
  int num_procs;// = 1;
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD,&my_rank);
  MPI_Comm_size (MPI_COMM_WORLD,&num_procs);

  bool    File            = false;
  const int     nDimensions     = 2;
  const int 	  nCycles	        = (int) 1e5;
  const double  omega		        = 1.0;
  const double  alpha           = 0.95;
  const double  beta		        = 0.45;//0.50905;
  const double  equilibration	  = 0.1;
  const double  derivativeStep  = 0.001;
  const double  stepLength      = 0.005;

  std::vector<double> parameters {omega, beta, alpha};
  int nParticles = 20;
  int orbitals   = 0;
  switch (nParticles)
  {
    case 2:
      orbitals   = 0;
      break;
    case 6:
      orbitals   = 1;
      break;

    case 12:
      orbitals   = 2;
      break;

    case 20:
      orbitals   = 3;
      break;
    default:
      cout << "Only 2, 6, 12 or 20 particles allowed" << endl;
      exit(1);
  }

  System* system = new System(File);
   
  system->set_parameters		        (parameters);
  system->set_orbitals              (orbitals);
  system->set_stepLength		        (stepLength);
  system->set_equilibrationFraction	(equilibration);
  system->set_derivativeStep		    (derivativeStep);
  system->set_nCycles               (nCycles);
  system->set_rank                  (my_rank);
  system->set_procs                 (num_procs);
  system->set_Hamiltonian	          (new HarmonicOscillator  (system));
  system->set_WaveFunction	        (new TrialSlater (system));
  system->set_InitialState	        (new RandomUniform	 (system, nDimensions, nParticles));

  system->set_DMatrix();
  //system->OPTIMIZE();
  //system->get_waveFunction()->LapJas();
  system->runMetropolis		();
  MPI_Finalize();
  return 0;
}
/*
  int chosenOne = 1;
  double  stepLength = 1.0;
  switch (chosenOne)
  {
    //case 0:
    //  stepLength = 0.01;
    //  importanceSampling(File,nCycles,nParticles,nDimensions,
    //      my_rank, num_procs,
    //      stepLength,equilibration,derivativeStep,
    //      parameters);
    //  break;
    case 1:
      stepLength = 5.7;
      
      //metropolis(File,nCycles,nParticles,nDimensions,
      //    my_rank, num_procs,
      //    stepLength,equilibration,derivativeStep,
      //    parameters);
      break;
  }

}
*/
/*
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
*/
/*
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
  system->set_Hamiltonian	  (new HarmonicOscillator  (system));
  system->set_WaveFunction	(new TrialSlater (system));
  Eigen::Matrix<double, nParticles,nParticles> DMatrix;
  system->set_DMatrix(DMatrix,2);
  std::cout << system->get_DMatrix()<< std::endl;
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
}*/
