#include "RandomUniform.h"

RandomUniform::RandomUniform( System* system,
			      int nDimensions,
			      int nParticles):
			      InitialState(system)
{
  
  my_nDimensions  = nDimensions;
  my_nParticles	  = nParticles;

  my_system->set_nDimensions  (nDimensions);
  my_system->set_nParticles   (nParticles);
  
  setupInitialState();
}

void RandomUniform::setupInitialState ()
{
  Eigen::MatrixXd particles (my_nParticles,my_nDimensions);  
  std::default_random_engine		  generator;
  std::uniform_real_distribution<double>  uniform(0.0,1.0);
  clock::duration dur = clock::now() - my_start;
  seed = dur.count();
  generator.seed(seed);

  //std::vector<double> position(my_nDimensions);
  for (int i = 0 ; i < my_nParticles ; i++){
    for(int j = 0 ; j < my_nDimensions ; j++){
      particles(i,j) = uniform(generator);
      //position[j] = (uniform(generator) - 0.5);
    }
    //my_system->add_particle(new Particle());
    //my_system->get_particle().at(i)->set_nDimensions    (my_nDimensions);
    //my_system->get_particle().at(i)->set_position       (position);
  }
  my_system->set_particles(particles); 
}
