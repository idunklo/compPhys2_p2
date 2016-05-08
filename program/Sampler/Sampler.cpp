#include "Sampler.h"

using std::cout;
using std::endl;

Sampler::Sampler(System* system, bool File){
  my_system	= system;
  my_stepNumber = 0;
  if (File){
    my_oEnergies.open("energies.out", std::ios::out | std::ios::binary);
    my_oPositions.open("positions.out",std::ios::out);
    my_oPositions << my_system->get_nCycles()*(1-my_system->get_equilibration()) << " "
		  << my_system->get_nParticles()  << " "
		  << my_system->get_nDimensions() << "\n";
  }
    
}

void Sampler::sample (bool accepted)
{
  if (my_stepNumber == 0)
  {
    cumulativeEnergy	        = 0;
    cumulativeEnergy2	        = 0;
    cumulativeAcceptanceRate  = 0;
  }

  double localEnergy = my_system->get_hamiltonian()->computeLocalEnergy();
  //for (int p = 0 ; p < my_system->get_nParticles() ; p++){
  //  cout << my_system->get_particle()[p]->get_position()[0] << "  " ;
  //  cout << my_system->get_particle()[p]->get_position()[0] << endl;
  //} 
  /*
  if (my_oEnergies.is_open())
  {
    my_oEnergies.write(reinterpret_cast<char *>(&localEnergy),sizeof(double));
    for (int p = 0 ; p < my_system->get_nParticles() ; p++){
      for (int d = 0 ; d < my_system->get_nDimensions() ; d++){
	      my_oPositions << my_system->get_particle()[p]->get_position()[d];
	      if (d+1 < my_system->get_nDimensions())
	        my_oPositions << ",";
      }
      if (p+1 < my_system->get_nParticles())
	      my_oPositions << ";";
    }
    my_oPositions << "\n";
  }
  */
  //else
  //{
    cumulativeEnergy	        += localEnergy;
    cumulativeEnergy2	        += localEnergy*localEnergy;
    cumulativeAcceptanceRate  += accepted;
    my_stepNumber		          += 1;
    //std::cout << my_stepNumber << std::endl;
  //}
} 

void Sampler::printResults ()
{

  if (my_oEnergies.is_open())
    my_oEnergies.close();
  if (my_oPositions.is_open())
    my_oPositions.close();
  int rank = my_system->get_rank();
  //cout << cumulativeAcceptanceRate << endl; 
  for (int p = 0 ; p<my_system->get_procs() ; p++){
    MPI_Barrier(MPI_COMM_WORLD);
    if(rank==p){
      double expectationValue    = cumulativeEnergy/(double)my_stepNumber;
      double expectationValue2   = cumulativeEnergy2/(double)my_stepNumber;
      double variance	       = (expectationValue2 - expectationValue * expectationValue)/
                                (double)my_stepNumber;
      cout << "==================" << endl;
      //cout << "Rank:             " << my_system->get_rank() << endl;
      cout << "Accept Ratio:     " << cumulativeAcceptanceRate/(double)my_stepNumber << endl;
      cout << "Expect Value:     " << expectationValue<< endl;
      cout << "Varian Ceeee:     " << variance << endl;
      cout << endl;
    }
  }
  
  //if(my_system->get_rank()==0){
/*
    double totalExpect  = 0;
    double totalExpect2 = 0;
    double totalVar	    = 0;
    double totalAccept  = 0;
  //}
  double expectationValue   = cumulativeEnergy/(double)my_stepNumber;
  double expectationValue2  = cumulativeEnergy2/(double)my_stepNumber;
  double variance	          = (expectationValue2 - expectationValue * expectationValue)/
                              (double)my_stepNumber;
  double acceptanceRatio    = cumulativeAcceptanceRate/(double)my_stepNumber;

  MPI_Reduce (&expectationValue, &totalExpect, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //MPI_Reduce (&expectationValue2, &totalExpect2, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //MPI_Reduce (&variance, &totalVar, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  //MPI_Reduce (&acceptanceRatio, &totalAccept, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    
  if (my_system->get_rank()==0){
    int	 nParticles	        = my_system->get_nParticles();
    int	 nDimensions	      = my_system->get_nDimensions(); 
    int  num_procs          = my_system->get_procs();
    double nCycles	        = my_system->get_nCycles();
    double alpha	          = my_system->get_parameters()[0];
    double beta             = my_system->get_parameters()[1];
    double omega	          = my_system->get_parameters()[2];
    double a                = my_system->get_parameters()[3];
//    double timeStep	       = my_system->get_timeStep();
    double stepLength	      = my_system->get_stepLength();
    double derivativeStep    = my_system->get_derivativeStep();
    //double expectationValue    = cumulativeEnergy/(double)my_stepNumber;
    double expectationValue2   = cumulativeEnergy2/(double)my_stepNumber;
    double variance	       = (expectationValue2 - expectationValue * expectationValue);
    double acceptanceRatio     = cumulativeAcceptanceRate/(double)my_stepNumber;
    cout << "taut: " << (totalExpect/num_procs - 3.0/4.0)/2.0 << endl;
    printf("\n");
    printf("\033[1;44m====================  System Data ====================\033[1;m\n");
    printf("\033[0;93mNumber of processes:     %i\033[0;m\n",num_procs);
    printf("\033[0;93mNumber of particles:     %i\033[0;m\n",nParticles);
    printf("\033[0;93mNumber of dimensions:    %i\033[0;m\n",nDimensions);
    printf("\033[0;93mNumber of cycles:        %0.1e\033[0;m\n",nCycles);
    printf("\033[0;93mAlpha:                   %f\033[0;m\n",alpha);
    printf("\033[0;93mOmega:                   %f\033[0;m\n",omega);
    printf("\033[0;93mBeta:                    %f\033[0;m\n",beta);
    printf("\033[0;93ma:                       %f\033[0;m\n",a);
    printf("\033[0;93mStep length:             %f\033[0;m\n",stepLength);
    printf("\033[0;93mNumber of steps          %i\033[0;m\n",my_stepNumber);
//    printf("\033[0;93mTime step:               %f\033[0;m\n",timeStep);
    printf("\033[0;93mDerivative step:         %f\033[0;m\n",derivativeStep);
    printf("\033[1;105m~~~~~~~~~~~~~~~~~~~~~ Results ~~~~~~~~~~~~~~~~~~~~~~~~\033[1;m\n");
    printf("\033[0;91mExpectation Value:       %e\033[0;m\n",totalExpect/num_procs);
    printf("\033[0;91mVariance:                %e\033[0;m\n",variance/nCycles);
    printf("\033[0;91mAcceptance ratio:        %f\033[0;m\n",acceptanceRatio);
    //printf("\033[0;91mExpectation Value:       %e\033[0;m\n",expectationValue);
    //printf("\033[0;91mVariance:                %e\033[0;m\n",variance);
    //printf("\033[0;91mAcceptance ratio:        %f\033[0;m\n",acceptanceRatio);
  }
  */
}
