#pragma once
#include "System.h"

class Particle 
{
  public:
    Particle();
    void set_nDimensions  (int nDimensions);
    void set_position	    (const std::vector<double> &position);
    void changePosition	  (int dimension, double adjustment);
    std::vector<double>& get_position() {return my_position;}

  private:
    int my_nDimensions	  = 0;
    std::vector<double>	my_position = std::vector<double>();

};


