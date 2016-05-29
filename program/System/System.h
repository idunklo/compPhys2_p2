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
#include "Eigen/LU"
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
    System			                    (bool File)
      {my_File = File;}
    //~System                         ()
    // {for (int i = 0 ; i < my_particles.size() ; i++)
    //    delete my_particles.at(i);
    //  delete my_sampler;
    // }

    bool   metropolis		      ();
    bool   importanceSampling	();
    void   runMetropolis		  ();
    void   OPTIMIZE           ();
    void   update_r_ij        (int pos);
    void   update_inverse     (Eigen::MatrixXd& matrix_inv,
                               Eigen::MatrixXd& matrix,
                               Eigen::VectorXd& new_row,
                               int i, double R_SD);

    /* Set functions */
    void set_nDimensions	          (int nDimensions)
      { my_nDimensions = nDimensions; }
    void set_nParticles	  	        (int nParticles)
      { my_nParticles = nParticles; }
    void set_nCycles	  	          (int nCycles)
      { my_nCycles = nCycles; }
    void set_orbitals               (int orbitals)
      { my_orbitals = orbitals; }
    void set_rank                   (int rank)
      { my_rank = rank; }
    void set_procs                  (int procs)
      { num_procs = procs; }
    void set_stepLength		          (double stepLength)
      { my_stepLength = stepLength;}
        //norm_dist my_normal(0,sqrt(stepLength));}
    void set_derivativeStep	        (double h)
      { my_derivativeStep = h; my_derivativeStep2 = h*h;}
    void set_equilibrationFraction  (double equilibraFraction)
      { my_equilibrationFraction = equilibraFraction; }
    void set_SDLap_up               (double SDLap_up)
      { my_SDLap_old_up = SDLap_up; }
    void set_SDLap_dn               (double SDLap_dn)
      { my_SDLap_old_dn = SDLap_dn; }
    void set_Hamiltonian	          (class Hamiltonian* hamiltonian)
      { my_hamiltonian = hamiltonian; }
    void set_WaveFunction	          (class WaveFunction* waveFunction)
      { my_waveFunction = waveFunction; }
    void set_InitialState	          (class InitialState* initialState)
      { my_initialState = initialState; }
    void set_Timer		              (class Timer* timer)
      { my_timer = timer; }
    void set_parameters		          (std::vector<double> parameters)
      { my_parameters = parameters; }
    //void add_particle (class Particle* particle)
    //  { my_particles.push_back(particle); }
    void set_DMatrix();
    void set_particles              (Eigen::MatrixXd& particles)
      { my_particles = particles;}
    
    /* Return functions */
    int	    get_nDimensions		  (){return my_nDimensions;}
    int     get_nParticles	   	(){return my_nParticles;}
    int	    get_nCycles			    (){return my_nCycles;}
    int     get_rank            (){return my_rank;}
    int     get_procs           (){return num_procs;}
    int     get_orbitals        (){return my_orbitals;}
    int     get_spin            (){return my_spin;}
    double  get_stepLength	    (){return my_stepLength;}
    double  get_derivativeStep	(){return my_derivativeStep;}
    double  get_derivativeStep2	(){return my_derivativeStep2;}
    double  get_equilibration		(){return my_equilibrationFraction;}
    double  get_SDLap_old_up    (){return my_SDLap_old_up;}
    double  get_SDLap_old_dn    (){return my_SDLap_old_dn;}
    std::vector<double>& get_parameters	(){return my_parameters;}
    Eigen::MatrixXd& get_DMatrix_up     (){return my_DMatrix_up;}
    Eigen::MatrixXd& get_DMatrix_up_inv (){return my_DMatrix_up_inv;}
    Eigen::MatrixXd& get_DMatrix_dn     (){return my_DMatrix_dn;}
    Eigen::MatrixXd& get_DMatrix_dn_inv (){return my_DMatrix_dn_inv;}
    Eigen::MatrixXd& get_r_ij           (){return my_r_ij;}
    Eigen::MatrixXd& get_a              (){return my_a;}
    Eigen::MatrixXd& get_particles      (){return my_particles;}

    class Hamiltonian*    get_hamiltonian   (){return my_hamiltonian;}
    class WaveFunction*		get_waveFunction  (){return my_waveFunction;}
    class Sampler*		    get_sampler       (){return	my_sampler;}
    class Timer*		      get_timer         (){return my_timer;}
    //std::vector<class Particle*>&  get_particle()	{return my_particles;}
    
  protected:
    std::ofstream my_oFile;
    bool   my_File                  = false;
    int    my_rank                  = 0;
    int    num_procs                = 1;
    int    my_elected               = 0;
    int    my_orbitals              = 0;
    int    my_nDimensions           = 0;
    int    my_nParticles            = 0;
    int    my_nCycles               = 0;
    int    my_spin                  = 0;
    double my_stepLength            = 0.0;
    double my_derivativeStep        = 0.0;
    double my_derivativeStep2       = 0.0;
    double my_equilibrationFraction = 0.0;
    double my_SDLap_old_up          = 0.0;
    double my_SDLap_old_dn          = 0.0;
    
    Eigen::MatrixXd my_DMatrix_up;
    Eigen::MatrixXd my_DMatrix_dn;
    Eigen::MatrixXd my_DMatrix_up_inv;
    Eigen::MatrixXd my_DMatrix_dn_inv;
    Eigen::MatrixXd my_r_ij;
    Eigen::MatrixXd my_a;
    Eigen::MatrixXd my_particles;
    std::vector<double> my_parameters	= std::vector<double>();

    class Hamiltonian*    my_hamiltonian    = nullptr;
    class WaveFunction*	  my_waveFunction  	= nullptr;
    class InitialState*	  my_initialState   = nullptr;
    class Sampler*	      my_sampler        = nullptr;
    class Timer*          my_timer          = nullptr;
    //std::vector<class Particle*> my_particles = std::vector<class Particle*>();


    typedef std::chrono::high_resolution_clock clock;
    clock::time_point my_start = clock::now();
    std::uniform_real_distribution<double> my_uniform {std::uniform_real_distribution<double>(0.0,1.0)};
    //std::normal_distribution<double> my_normal {std::normal_distribution<double>(0,1.0/sqrt(2))};
    std::mt19937 my_generator;
};


