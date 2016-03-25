#pragma once
#include "../System/System.h"
#include "Timer.h"
#include "ProgressBar.h"

class Sampler {
  public:
    Sampler(class System* system, bool File);

    void sample		(bool accepted);
    void printResults	();
    
    double get_cumuEn (){return cumulativeEnergy;}

  protected:
    std::ofstream my_oEnergies;
    std::ofstream my_oPositions;
    int	    my_nCycles		      = 0;
    int	    my_stepNumber	      = 0;
    double  cumulativeEnergy	      = 0;
    double  cumulativeEnergy2	      = 0;
    double  cumulativeAcceptanceRate  = 0;
    class   System* my_system	      = nullptr;
};

