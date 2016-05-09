#pragma once
#include "../System/System.h"

class InitialState
{
  public:
    InitialState (class System* system);
    
    virtual void setupInitialState () = 0;

    //std::vector<class Particle*> get_particles() {return my_particles;}
    
  protected:
    
    int my_nDimensions	= 0;
    int my_nParticles	= 0;

    class System* my_system = nullptr;

    //std::vector<class Particle*> my_particles;
};
