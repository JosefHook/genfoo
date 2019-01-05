#pragma once
#ifndef DATIO_H
#define DATIO_H

#include "Particle.hpp"
#include "WeightedParticle.hpp"


class DatIO
{
public:
  DatIO(std::string filename);
  void setNumParticles(long long n, int D);
  void operator<<(const WeightedParticle &wp);
  void operator<<(const Particle &p);
  void setIterationStep(int t);
private:
  std::string _filename;
  int _time_step;
  long long _np;


};




#endif
