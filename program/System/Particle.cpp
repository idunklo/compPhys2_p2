#include "Particle.h"

Particle::Particle(){}

void Particle::set_nDimensions (int nDimensions)
{
  my_nDimensions = nDimensions;
}

void Particle::set_position (const std::vector<double> &position)
{
  my_position = position;
}

void Particle::changePosition (double x, double y)
{
  my_position.at(0) += x;
  my_position.at(1) += y;
}


