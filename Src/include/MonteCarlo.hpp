#pragma once
#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "Particle.hpp"
#include "mArray.hpp"
#include "WeightedParticle.hpp"
#include "Operator.hpp"
#include "ReflectedBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"



class MonteCarlo 
{

 public:
  MonteCarlo() {}
  MonteCarlo(Operator *ops, double dt) : _ops(ops), _dt(dt) {  }
  MonteCarlo(Operator *ops) : _ops(ops) {
    _drift = new double[ops->getDim()];
    _diff = new double[(ops->getDim()*ops->getDim())];
  }
  MonteCarlo(Operator *ops, RandomNumberGenerator *rng) : _ops(ops), _dW(rng) {
    _drift = new double[ops->getDim()];
    _diff = new double[(ops->getDim()*ops->getDim())];
  }
  ~MonteCarlo() {
    delete _drift;
    delete _diff;
  }
  int assemble( Operator *ops, ReflectedBoundaryCondition *bc, 
		RandomNumberGenerator *rng);
  int init(Particle *p);
  int step(Particle *p, double dt);
  int adaptivestep(Particle *p, mArray *dt );
  int initsource(Particle *p, double dl);
  int initsource(Particle *p, double dl, int burnin);
  int initsource(WeightedParticle *wp, double dl, int burnin);
  int sample(Particle *p, Particle *s, double np);
  int sample(WeightedParticle *wp, WeightedParticle *ws, double np);
  
private: 
  double _dt;
  Operator *_ops;
  double *_drift, *_diff;
  RandomNumberGenerator *_dW; 
  ReflectedBoundaryCondition *_bc;
};




#endif
