//
// A Random number generator class
//
#pragma once
#ifndef RANDOMNUMBERGENERATOR_H
#define RANDOMNUMBERGENERATOR_H
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>

class RandomNumberGenerator
{

public:
  RandomNumberGenerator() {
    boost::normal_distribution<> _nd(0.0, 1.0);  
    _rng = new boost::variate_generator<boost::mt19937&, boost::normal_distribution<> >(_rng_seed, _nd);
  }
  double operator()(void) {
    return (*_rng)();
  }
  ~RandomNumberGenerator() {
    delete _rng;
  }
private:
  boost::mt19937 _rng_seed;
  boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > *_rng;
  

};
#endif
