// $Id$ 
// Drift and diffusion coefficient file for 
// the Ornstein-Uhlenbeck process in 2D 
// With diffusion in X-direction
// By Josef Höök Copyright 2010 (C) all rights reserved
//
#include <cmath>
#include "Operator.hpp"
#define DIM 2



// The operator class 
class OrnUhl2D : public Operator
{

public:
  
  OrnUhl2D();
  void evalDrift(double* val, const double* x) const;
  void evalDiffusion(double* val, const double* x)const;
  void evalIC(double* val, const double* x) const ;
  void evalSource(double* val, const double* x) const;


};



// Default constructor
OrnUhl2D::OrnUhl2D()
{
  _dim = DIM;
  _info = "Ornstein-Uhlenbeck process constant diagional diffusion.";

}




  
void OrnUhl2D::evalIC(double* values, const double* x) const
{
  
  values[0] = 20.0*1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-( (x[0]-0.0)*(x[0]-0.0) + (x[1]-0.0)*(x[1]-0.0) + (x[2]-0.0)*(x[2]-0.0))/0.2 ) ; 
  
  
}


// Drift to (x,y,z)=(5,5,5)
//populate drift vector  
void OrnUhl2D::evalDrift(double* values, const double* x) const
  {

	values[0] = -1.0*(x[0]- 0.0);
	values[1] = -1.0*(x[1]- 0.0);
	values[2] = -1.0*(x[2]- 0.0);
    
  }

//
// Alignment of Matrix in array
//
// [1,1 1,2 1,3
//  2,1 2,2 2,3
//  3,1 3,2 3,3]
// => 
// [ 0   1   2   3   4   5   6   7   8 ] 
// [1,1 1,2 1,3 2,1 2,2 2,3 3,1 3,2 3,3]
// row_major ordering 
// M_(m x n) => M_ij = [i*n + j]
//



void OrnUhl2D::evalDiffusion(double* values, const double* x) const
  {
    values[0] = 1.0; // B11 X-diffusion
    values[1] = 1.0; // B12

    values[2] = 1.0; // B21 
    values[3] = 1.0; // B22 Y-diffusion  



  }


void OrnUhl2D::evalSource(double* val, const double* x) const {

  val[0] = 0.0;
  
};

// Code snippet needed for object initalization
// Function always needed!
extern "C" {
  Operator *load()  { return new OrnUhl2D;  }
}



