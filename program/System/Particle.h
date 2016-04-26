#pragma once
#include "System.h"

class Particle 
{
  public:
    Particle();
    void set_nDimensions  (int nDimensions);
    void set_position	    (const std::vector<double> &position);
    void changePosition	  (double x, double y);
    std::vector<double>& get_position() {return my_position;}

  private:
    int my_nDimensions	  = 0;
    std::vector<double>	my_position = std::vector<double>();

};


