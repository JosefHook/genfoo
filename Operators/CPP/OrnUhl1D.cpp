//
//
//
// $Id$ 
// Drift and diffusion coefficient file for 
// the Ornsteing-Uhlenbeck process 1D version
// By Josef Höök
//

#include <cmath>
#include "Operator.hpp"
#define DIM 1

// The operator class 
class OrnUhl1D : public Operator
{

public:
  
  OrnUhl1D();
  void evalDrift(double* val, const double* x) const;
  void evalDiffusion(double* val, const double* x)const;
  void evalIC(double* val, const double* x) const ;
  void evalSource(double* val, const double* x) const;


};


// Default constructor
OrnUhl1D::OrnUhl1D()
{


  _dim = DIM;
  _info = "1D Ornstein-Uhlenbeck process constant diagional diffusion.";

}

  
void OrnUhl1D::evalIC(double* values, const double* x) const
{
    
  values[0] = 20.0*1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-( (x[0]-0.0)*(x[0]-0.0)  )/0.2 ); 
   
  
}

// Drift to (x,y,z)=(5,5,5)
//populate drift vector  
void OrnUhl1D::evalDrift(double* values, const double* x) const
  {
      values[0] = -2.0*(x[0]-1.5);

    
  }

void OrnUhl1D::evalDiffusion(double* values, const double* x) const
  {
    values[0] = 2.0; // B11 X-diffusion

  }


void OrnUhl1D::evalSource(double* val, const double* x) const {

  val[0] = 0.0;
  
};

// Code snippet needed for object initalization
// Function always needed!
extern "C" {
  Operator *load()  { return new OrnUhl1D;  }
}

