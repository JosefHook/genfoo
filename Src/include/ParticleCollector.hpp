//
// Particle Collector class
//

#pragma once
#ifndef PARTICLE_COLLECTOR_H
#define PARTICLE_COLLECTOR_H

#include "Particle.hpp"



class ParticleCollector 
{


public:
  ParticleCollector(Particle *p) : _p(p) {}
  int collect(int nr);
  void knn(Particle *p, int knn);
  
private:
  Particle *_p;
    
    


};


#endif
