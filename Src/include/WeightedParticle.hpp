#pragma once
#ifndef WEIGHTEDPARTICLE_H
#define WEIGHTEDPARTICLE_H

#include "Particle.hpp"
#include "mArray.hpp"



class WeightedParticle 
{
  
public:
  WeightedParticle() { _p = NULL; _w = NULL; }
  WeightedParticle(Particle *p, mArray *w) { _p = p; _w = w; }
  mArray* getWeight() const { return _w; }
  Particle* getParticle() const { return _p; }
  void setParticle(Particle *p)  { _p = p;  }
  void setWeight(mArray *w)  { _w = w;  }
private: 
  Particle *_p;
  mArray *_w;
  
};




#endif
