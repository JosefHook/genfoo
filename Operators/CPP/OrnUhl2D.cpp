//
//
//
// $Id$ 
// Drift and diffusion coefficient file for 
// the Ornsteing-Uhlenbeck process 2D version
// By Josef Höök
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
  _info = "2D Ornstein-Uhlenbeck process constant diagional diffusion.";

}

  
void OrnUhl2D::evalIC(double* values, const double* x) const
{
    
  values[0] = 20.0*1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-( (x[0]-0.0)
							  *(x[0]-0.0) + (x[1]-0.0)*(x[1]-0.0) )/0.2 ); 
  
  
}

// Drift to (x,y,z)=(5,5,5)
//populate drift vector  
void OrnUhl2D::evalDrift(double* values, const double* x) const
  {
  
  for(int i=0; i<2; i++)
      values[i] = -2.0*(x[i]-1.5);
  
    
  }

void OrnUhl2D::evalDiffusion(double* values, const double* x) const
  {
    values[0] = 2.0; // B11 X-diffusion
    values[1] = 0.0; // 
    
    values[2] = 0.0; // 
    values[3] = 2.0; // B22 Y-diffusion  
    
  }


void OrnUhl2D::evalSource(double* val, const double* x) const {

  val[0] = 0.0;
  
};

// Code snippet needed for object initalization
// Function always needed!
extern "C" {
  Operator *load()  { return new OrnUhl2D;  }
}

