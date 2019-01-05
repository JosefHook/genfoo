//  
// North-Source flow example 
//
// By Josef Höök Copyright 2010 (C) all rights reserved
//
#include <cmath>
#include "Operator.hpp"

#define DIM 1

#include <iostream>

// The operator class 
class NSFlow : public Operator
{

public:
  
  NSFlow();
  void evalDrift(double* val, const double* x) const;
  void evalDiffusion(double* val, const double* x)const;
  void evalIC(double* val, const double* x) const ;
  void evalSource(double* val, const double* x) const;


};



// Default constructor
NSFlow::NSFlow()
{
  _dim = DIM;
  _info = "North Source flow on a circle presented in Carvehill, Characteristic exponents for stochastic flows,Stochastic Processes ";

}


  
void NSFlow::evalIC(double* values, const double* x) const
{
  // Uniform distribution
  values[0] = 1.0;

}


//populate drift vector  
void NSFlow::evalDrift(double* values, const double* x) const
  {
    double K = 1.0;
    values[0]= -K*sin(x[0]);


    
  }

void NSFlow::evalDiffusion(double* values, const double* x) const
  {

    values[0] = 1.0; 


  }


void NSFlow::evalSource(double* val, const double* x) const {

	       

  
};

// Code snippet needed for object initalization
// Function always needed!
extern "C" {
  Operator *load()  { return new NSFlow;  }
}



