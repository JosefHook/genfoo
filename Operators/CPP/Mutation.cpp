//  
// Allele Mutation Diffusion process 
// from W, J, Ewens, Mathematical population genetics
// p. 157
// By Josef Höök Copyright 2010 (C) all rights reserved
//
#include <cmath>
#include "Operator.hpp"

#define DIM 1

#include <iostream>

// The operator class 
class MultiAllele : public Operator
{

public:
  
  MultiAllele();
  void evalDrift(double* val, const double* x) const;
  void evalDiffusion(double* val, const double* x)const;
  void evalIC(double* val, const double* x) const ;
  void evalSource(double* val, const double* x) const;


};



// Default constructor
MultiAllele::MultiAllele()
{
  _dim = DIM;
  _info = "Mutation Diffusion process from W, J, Ewens, Mathematicla population genetics";

}


  
void MultiAllele::evalIC(double* values, const double* x) const
{
  // Uniform distribution
  values[0] = 1.0;

}


//populate drift vector  
void MultiAllele::evalDrift(double* values, const double* x) const
  {
    double N = 1e3;
    double s =  0.01;
      double u = 0.02;  // One percent mutation A1->A2 in population     
      double v = 0.0; //01;  // mutation A2->A1 in population 
    double alpha = 2.0*N*s;
    double beta1 = 2.0*N*u;
    double beta2 = 2.0*N*v;
    double h = 0.5; // No dominance in Selection
    //    values[0]= -2.0*N*u;
    values[0] = alpha*x[0]*(1.0 -x[0])*( x[0] + h*(1.0 -2*x[0]))
      - beta1*x[0] + beta2*(1.0 - x[0]);
   

    
  }

void MultiAllele::evalDiffusion(double* values, const double* x) const
  {

    values[0] = x[0]*(1.0 -x[0]);


  }


void MultiAllele::evalSource(double* val, const double* x) const {

	       

  
};

// Code snippet needed for object initalization
// Function always needed!
extern "C" {
  Operator *load()  { return new MultiAllele;  }
}



