#pragma once
#include "../System/System.h"

class Hamiltonian
{
  public:

    Hamiltonian (class System* system);

    virtual double computeLocalEnergy () = 0;

  protected:
    
    class System* my_system = nullptr;

};
