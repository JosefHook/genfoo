// $Id$ 
// Drift and diffusion coefficient file for 
// the Ornstein-Uhlenbeck process in 3D 
// With diffusion in X-direction
// By Josef Höök Copyright 2010 (C) all rights reserved
//
#include <cmath>
#include "Operator.hpp"
#define DIM 3



// The operator class 
class OrnUhl3D : public Operator
{

public:
  
  OrnUhl3D();
  void evalDrift(double* val, const double* x) const;
  void evalDiffusion(double* val, const double* x)const;
  void evalIC(double* val, const double* x) const ;
  void evalSource(double* val, const double* x) const;


};



// Default constructor
OrnUhl3D::OrnUhl3D()
{
  _dim = DIM;
  _info = "Ornstein-Uhlenbeck process constant diagional diffusion.";

}




  
void OrnUhl3D::evalIC(double* values, const double* x) const
{
  
  values[0] = 20.0*1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-( (x[0]-0.0)*(x[0]-0.0) + (x[1]-0.0)*(x[1]-0.0) + (x[2]-0.0)*(x[2]-0.0))/0.2 ) ; 
  
  
}


// Drift to (x,y,z)=(5,5,5)
//populate drift vector  
void OrnUhl3D::evalDrift(double* values, const double* x) const
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



void OrnUhl3D::evalDiffusion(double* values, const double* x) const
  {
    values[0] = 0.0; // B11 X-diffusion
    values[1] = 0.0; // B12
    values[2] = 0.0; // B13  

    values[3] = 0.0; // B21 
    values[4] = 0.0; // B22 Y-diffusion  
    values[5] = 0.0; // B23 ?

    values[6] = 0.0; // B31
    values[7] = 0.0; // B32 Y- diffusion (with z-component)
    values[8] = 100.0; // B33


  }


void OrnUhl3D::evalSource(double* val, const double* x) const {

  val[0] = 0.0;
  
};

// Code snippet needed for object initalization
// Function always needed!
extern "C" {
  Operator *load()  { return new OrnUhl3D;  }
}



