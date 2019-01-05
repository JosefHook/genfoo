#pragma once
#ifndef BALLRESAMPLE_H
#define BALLRESAMPLE_H
#include <ANN/ANN.h>
#include <vector>
#include "WeightedParticle.hpp"
//#include "mArray.hpp"
#include "RandomNumberGenerator.hpp"

#define DEFAULT_NEIGHBOURS 5
#define DEFAULT_REDUCE 0.10 // 10 % reduction


class BallResample 
{
public:
  BallResample( RandomNumberGenerator *rng ) : _nneigh(DEFAULT_NEIGHBOURS),
					       _reduce(DEFAULT_REDUCE),
					       _rng(rng) {};
  BallResample( int nneigh, double reduce, RandomNumberGenerator *rng 
		) : _nneigh(nneigh), 
		    _reduce(reduce),  
		    _rng(rng) {};
  int reduce(WeightedParticle *owp, 
	     double radius, int maxneighbors, double reducefactor );
  ~BallResample() {};
  

private:

  int calcParticleWeights(ANNpointArray newp,
			  mArray *neww,
			  ANNpointArray oldp,
			  mArray *oldw, 
			  ANNidxArray nnIdx,
			  double radius,
			  ANNpoint qp,
			  double qw,
			  int n_ball, 
			  int n_reduced,
			  long gIdx,
			  int dim);

  int _nneigh; 
  double _reduce;
  RandomNumberGenerator *_rng;
  

  
  

};




#endif
