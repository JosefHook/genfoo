#pragma once
#ifndef PARTICLE_H
#define PARTICLE_H

// Disable Boost::ublas debugging
#define NDEBUG


#include <boost/numeric/ublas/matrix.hpp>



//typedef boost::numeric::ublas::matrix<double> Particle; 

class Particle : public boost::numeric::ublas::matrix<double>
{
  
public:
  Particle():boost::numeric::ublas::matrix<double>() {}
  Particle(int row, int col) : boost::numeric::ublas::matrix<double>(row, col) ,  _n_pos(0) {}
  void set_pos(long long np) { _n_pos = np; }
  long long pos() const { return _n_pos; }
  
private: 
  long long _n_pos; //particle position in matrix
  
  
};




#endif
