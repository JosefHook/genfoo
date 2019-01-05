//  
// Implementation of the Duffing-Van der Pol oscillator
// from P.E. Kloeden and E. Platen 
// Numerical Solution of Stochastic differential equations
// 
// By Josef Höök Copyright 2010 (C) all rights reserved
//
#include <cmath>
#include "Operator.hpp"

#define DIM 2

#include <iostream>

// The operator class 
class Duffing : public Operator
{

public:
  
  Duffing();
  void evalDrift(double* val, const double* x) const;
  void evalDiffusion(double* val, const double* x)const;
  void evalIC(double* val, const double* x) const ;
  void evalSource(double* val, const double* x) const;


};



// Default constructor
Duffing::Duffing()
{
  _dim = DIM;
  _info = "North Source flow on a circle presented in Carvehill, Characteristic exponents for stochastic flows,Stochastic Processes ";

}


  
void Duffing::evalIC(double* values, const double* x) const
{
  // (x1, x2) = (-k epsilon , 0) , k = 11, 12, ..., 20
  // epsilon = 0.2 thatis a 
  // Uniform distribution over -20*0.2..-11*0.2
  
  // Smooth approximation of the box using the logistic function approximation 
  // of the heaviside function. This is needed since 
  // the Metropolis-Hastings method used for sampling
  // cannot handle step-functions 

  double k = 10.0;
  values[0] = 1.0; 
  ///(1.0 + std::exp(-2*k*( x[0] + 4.0 ) )) - 1.0/(1.0 + std::exp(-2*k*( x[0]  + 2.2 )));
  
   
}


//populate drift vector  
void Duffing::evalDrift(double* values, const double* x) const
  {
    double alpha = 1.0;
    values[0]= x[1];
    values[1] = x[0]*(alpha - std::pow(x[0], 2)) - x[1];

    
  }

void Duffing::evalDiffusion(double* values, const double* x) const
  {
    double sigma = 0.2;
    values[0] = 0.0;
    values[0] = sigma*x[0]; 


  }


void Duffing::evalSource(double* val, const double* x) const {

	       

  
};

// Code snippet needed for object initalization
// Function always needed!
extern "C" {
  Operator *load()  { return new Duffing;  }
}



