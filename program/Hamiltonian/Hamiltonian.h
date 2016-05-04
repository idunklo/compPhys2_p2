#pragma once
#include "../System/System.h"

class Hamiltonian
{
  public:

    Hamiltonian (class System* system);

    virtual double computeLocalEnergy    () = 0;
    virtual double HOLap                 () = 0;
    virtual double HOExt                 () = 0;
    virtual double HRep                  () = 0;

  protected:
    
    class System* my_system = nullptr;

};
