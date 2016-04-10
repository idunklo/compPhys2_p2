#pragma once
#include "../System/System.h"

class Hamiltonian
{
  public:

    Hamiltonian (class System* system);

    virtual double computeNumLocalEnergy () = 0;
    virtual double computeAnaLocalEnergy () = 0;

  protected:
    
    class System* my_system = nullptr;

};
