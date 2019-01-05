// $Id$ 
// Drift and diffusion coefficient file for 
// the Ornstein-Uhlenbeck process in 3D 
// By Josef Höök Copyright 2010 (C) all rights reserved
//
#include <cmath>
#include "Operator.hpp"
#include "RF1DParams.hpp"

#define DIM 1

#include <iostream>

// The operator class 
class RF1D : public Operator
{

public:
  
  RF1D();
  void evalDrift(double* val, const double* x) const;
  void evalDiffusion(double* val, const double* x)const;
  void evalIC(double* val, const double* x) const ;
  void evalSource(double* val, const double* x) const;


};



// Default constructor
RF1D::RF1D()
{
  _dim = DIM;
  _info = "C + Q operator in (x) 1D ";

}



// x[0] = v/vth = x

  
void RF1D::evalIC(double* values, const double* x) const
{

   values[0] = 10.0*1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-( (x[0]-0.0)*(x[0]-0.0) )/0.2 ) ;

}


//populate drift vector  
void RF1D::evalDrift(double* values, const double* x) const
  {

  long double Gp =( erf(x[0]) - 0.2e1*x[0] 
		     *exp(  -std::pow(x[0],2) )
		     *std::pow(M_PI,(-0.1e1/0.2e1)) )
     /std::pow(x[0],2)/0.2e1;
   
   long double Ge =(erf(lelp*x[0]) - 0.2e1*lelp*x[0]
		    *exp(-std::pow(lelp*x[0],2))
		    *std::pow(M_PI,(-0.1e1/0.2e1))  )
     /std::pow(lelp,2)/std::pow(x[0],2)/0.2e1; 
   
   long double gamma =  (   (Cp/vth)/x[0]    )*( erf(x[0]) - Gp   ) 
     + ( (Ce/vth)/x[0]  )*( erf(lelp*x[0]) - Ge   );  
   
   long double beta = ((Cp/vth)/x[0])*Gp + ((Ce/vth)/x[0])*Ge;
   long double alpha  = -Cp*std::pow(lp,2)*( 1.0 + m/mp )*Gp
     -Ce*std::pow(le,2)*(1.0  + m/mee )*Ge 
     + Cp/(2.0*std::pow(vth*x[0],2))*( erf(x[0]) -  Gp  ) 
     + Ce/(2.0*std::pow(vth*x[0],2))*( erf(lelp*x[0])  - Ge  );
   
   
   values[0]= alpha + 2.0*K/x[0];


    
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


// The diffsion Give rise to Huge velocity diffusion
void RF1D::evalDiffusion(double* values, const double* x) const
  {
    // Ordering of the variables
    // r=x[0], x = x[0], xi = x[1]
    
    long double Gp =( erf(x[0]) - 0.2e1*x[0] 
		      *exp(  -std::pow(x[0],2) )
		      *std::pow(M_PI,(-0.1e1/0.2e1)) )
      /std::pow(x[0],2)/0.2e1;
    
    long double Ge =(erf(lelp*x[0]) - 0.2e1*lelp*x[0]
		     *exp(-std::pow(lelp*x[0],2))
		     *std::pow(M_PI,(-0.1e1/0.2e1))  )
      /std::pow(lelp,2)/std::pow(x[0],2)/0.2e1; 
    
   
   long double gamma =  (   (Cp/vth)/x[0]    )*( erf(x[0]) - Gp   ) 
     + ( (Ce/vth)/x[0]  )*( erf(lelp*x[0]) - Ge   );  
   
   long double beta = ((Cp/vth)/x[0])*Gp + ((Ce/vth)/x[0])*Ge;
        

   values[0] = std::sqrt(beta  + 2.0*K);


  }


void RF1D::evalSource(double* val, const double* x) const {

  val[0] = std::abs((1.0/4.0)*K*nh*sqrt(2)*std::pow(m,2) 
	       *sqrt(m/(kb*T))*exp(-(1.0/2.0)*m*std::pow(x[0],2)/(kb*T))
	       *(-3*kb*T+m*std::pow(x[0],2) )
	       /(std::pow(M_PI, 3.0/2.0)*std::pow(kb*T, 3))
	       *std::pow(x[0],2)); // Jacobian
  ///std::pow(vth, 2)
	       

  
};

// Code snippet needed for object initalization
// Function always needed!
extern "C" {
  Operator *load()  { return new RF1D;  }
}



