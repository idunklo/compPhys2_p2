#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
#include <cmath>
#include <string>
#include "mpi.h"
#include "Eigen/Core"
#include "Particle.h"
#include "../Sampler/Sampler.h"
#include "../Hamiltonian/Hamiltonian.h"
#include "../WaveFunction/WaveFunction.h"
#include "../InitialState/InitialState.h"

using std::cout;
using std::endl;
class System
{
 
  public:
    System			                    (bool File);
    ~System                         ();
    bool metropolis		              ();
    bool importanceSampling	        ();
    void runMetropolis		          ();
    void OPTIMIZE                   ();
    void runImportanceSampling	    ();
    void set_nDimensions	          (int nDimensions);
    void set_nParticles	  	        (int nParticles);
    void set_nCycles	  	          (int nCycles);
    void set_rank                   (int my_rank);
    void set_procs                  (int num_procs);
    void set_stepLength		          (double stepLength);
    void set_derivativeStep	        (double h);
    void set_equilibrationFraction  (double equilibractionFraction);
    void set_Hamiltonian	          (class Hamiltonian* hamiltonian);
    void set_WaveFunction	          (class WaveFunction* waveFunction);
    void set_InitialState	          (class InitialState* initalState);
    void set_Timer		              (class Timer* timer);
    void set_parameters		          (std::vector<double> parameters);
    void set_matrix                 ();
    
    int	    get_nDimensions		  (){return my_nDimensions;}
    int     get_nParticles	   	(){return my_nParticles;}
    int	    get_nCycles			    (){return my_nCycles;}
    int     get_rank            (){return my_rank;}
    int     get_procs           (){return num_procs;}
    double  get_stepLength	    (){return my_stepLength;}
    double  get_derivativeStep	(){return my_derivativeStep;}
    double  get_derivativeStep2	(){return my_derivativeStep2;}
    double  get_equilibration		(){return my_equilibrationFraction;}
    std::vector<double>& get_parameters	(){return my_parameters;}


    class Hamiltonian*    get_hamiltonian   (){return my_hamiltonian;}
    class WaveFunction*		get_waveFunction  (){return my_waveFunction;}
    class Sampler*		    get_sampler       (){return	my_sampler;}
    class Timer*		      get_timer         (){return my_timer;}
    std::vector<class Particle*>&  get_particle()	{return my_particles;}

    void  add_particle (); 
    
  protected:
    //std::ofstream my_oFile;
    bool    my_File		                = false;
    int     my_rank                   = 0;
    int     num_procs                 = 1;
    int     my_nDimensions		        = 0;
    int     my_nParticles		          = 0;
    int     my_nCycles		            = 0;
    double  my_stepLength	            = 0.0;
    double  my_derivativeStep	        = 0.0;
    double  my_derivativeStep2	      = 0.0;
    double  my_equilibrationFraction  = 0.0;
    
    Eigen::MatrixXd my_invSlater;
    std::vector<double> my_parameters	= std::vector<double>();

    class Hamiltonian*    my_hamiltonian	  = nullptr;
    class WaveFunction*	  my_waveFunction  	= nullptr;
    class InitialState*	  my_initialState   = nullptr;
    class Sampler*	      my_sampler	      = nullptr;
    class Timer*	        my_timer		      = nullptr;

    std::vector<class Particle*> my_particles = std::vector<class Particle*>();

    typedef std::chrono::high_resolution_clock clock;
    clock::time_point my_start = clock::now();
    std::uniform_real_distribution<double> my_uniform {std::uniform_real_distribution<double>(0.0,1.0)};
    std::normal_distribution<double> my_normal {std::normal_distribution<double>(0,1.0/sqrt(2))};
    std::mt19937 my_generator;
};


