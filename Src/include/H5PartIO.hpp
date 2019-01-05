#pragma once
#ifndef H5PartIO_H
#define H5PartIO_H

#define PARALLEL_IO 1

#include <H5Part.h>
#include "Particle.hpp"
#include "WeightedParticle.hpp"
class H5PartIO 
{
  
public:
  

  
 
  H5PartIO(std::string);
  H5PartIO(std::string, long long n);
  H5PartIO(std::string, long long n, int D);
  int file(std::string);
  void setTimeSteps(int);
  void setNumParticles(long long n);
  void setNumParticles(long long n, int D);
  int close(void);
  void operator<<(const Particle &p);    
  void operator<<(const WeightedParticle &p);  

private:
  std::string _filename;
  H5PartFile *_file_h;
  int _time_steps;
  long long _np;

  


};



#endif
