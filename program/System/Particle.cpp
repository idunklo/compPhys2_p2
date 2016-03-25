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

void Particle::changePosition (int dimension, double movement)
{
  my_position.at(dimension) += movement;
}


