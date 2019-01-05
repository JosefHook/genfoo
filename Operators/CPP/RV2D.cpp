//
// Implementation of 1D RF diffusion on radial fluxsurfaces
// with a T(r) = T_0( 1 - (r/r_0)^2 ) this gives a source term 
// will vary over r.
 
// $Id$ 
// Drift and diffusion coefficient file for 
// the Ornstein-Uhlenbeck process in 3D 
// By Josef Höök Copyright 2010 (C) all rights reserved
//
#include <cmath>
#include "Operator.hpp"
#include "RV2DParams.hpp"

#define DIM 2

#include <iostream>


// The operator class 
class RV2D : public Operator
{

public:
  
  RV2D();
  void evalDrift(double* val, const double* x) const;
  void evalDiffusion(double* val, const double* x)const;
  void evalIC(double* val, const double* x) const ;
  void evalSource(double* val, const double* x) const;
  double T(const double* r) const;

};

// Default constructor
RV2D::RV2D()
{
  _dim = DIM;
  _info = "C + Q operator in (x) 1D ";

}

// Temperature function
double RV2D::T(const double* r) const {
  double r0 = 1.0;
  return T0*(1.0 - std::pow(r[0]/r0,2.0));
}



// x[0] = v/vth = x

  
void RV2D::evalIC(double* values, const double* x) const
{

  values[0] = 10.0*1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-( (x[0]-0.0)*(x[0]-0.0) )/0.2 )*std::pow(x[0], 2.0);

}


//populate drift vector  
void RV2D::evalDrift(double* values, const double* x) const
  {
    long double Tp, Te;
    Tp = T(&x[1]);
    Te = Tp;
    long double lelp = sqrt( mee*Tp/(mp*Te) );
    long double loglambda = std::log( CONV*3/(2*Z*Zp*std::pow(q,3))*sqrt( std::pow(kb*T(&x[1]),3)/(M_PI*np)) );
    long double lp = sqrt(mp/(2*kb*Tp));
    long double le = sqrt(mee/(2*kb*Te));
    long double vth = 1/lp;
    long double p1 = q/m;
    long double p2 = std::pow(Z,2)*loglambda;
    long double p3 = q/eps0;
    long double p4 = np*q/eps0;
    long double p5 = ne*q/eps0;
    long double Cp = (8.0*M_PI*p2*std::pow(Zp,2))*(p1*p4)*(p1*p3)/( std::pow(4.0*M_PI,2));
    long double Ce = (8.0*M_PI*p2*std::pow(Ze,2))*(p1*p5)*(p1*p3)/( std::pow(4.0*M_PI,2));

    
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
    
    
    values[0]= alpha/std::pow(x[0],2.0) + 2.0*K/std::pow(x[0],3.0);
    values[1] = values[1];

    
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
void RV2D::evalDiffusion(double* values, const double* x) const
  {
    long double Tp, Te;
    Tp = T(&x[1]);
    Te = Tp;
    long double lelp = sqrt( mee*Tp/(mp*Te) );
    long double loglambda = std::log( CONV*3/(2*Z*Zp*std::pow(q,3))*sqrt( std::pow(kb*T(&x[1]),3)/(M_PI*np)) );
    long double lp = sqrt(mp/(2*kb*Tp));
    long double le = sqrt(mee/(2*kb*Te));
    long double vth = 1/lp;
    long double p1 = q/m;
    long double p2 = std::pow(Z,2)*loglambda;
    long double p3 = q/eps0;
    long double p4 = np*q/eps0;
    long double p5 = ne*q/eps0;
    long double Cp = (8.0*M_PI*p2*std::pow(Zp,2))*(p1*p4)*(p1*p3)/( std::pow(4.0*M_PI,2));
    long double Ce = (8.0*M_PI*p2*std::pow(Ze,2))*(p1*p5)*(p1*p3)/( std::pow(4.0*M_PI,2));


    
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
        

   values[0] = std::sqrt(beta  + 2.0*K)/x[0];
   values[1] = values[1];

  }

//
// x[0] = v, x[1] = r
//
void RV2D::evalSource(double* val, const double* x) const {


  long double Tp, Te;
  Tp = T(&x[1]);
  Te = Tp;
  long double lelp = sqrt( mee*Tp/(mp*Te) );
  long double loglambda = std::log( CONV*3/(2*Z*Zp*std::pow(q,3))*sqrt( std::pow(kb*T(&x[1]),3)/(M_PI*np)) );
  long double lp = sqrt(mp/(2*kb*Tp));
  long double le = sqrt(mee/(2*kb*Te));
  long double vth = 1/lp;
  long double p1 = q/m;
  long double p2 = std::pow(Z,2)*loglambda;
  long double p3 = q/eps0;
  long double p4 = np*q/eps0;
  long double p5 = ne*q/eps0;
  long double Cp = (8.0*M_PI*p2*std::pow(Zp,2))*(p1*p4)*(p1*p3)/( std::pow(4.0*M_PI,2));
  long double Ce = (8.0*M_PI*p2*std::pow(Ze,2))*(p1*p5)*(p1*p3)/( std::pow(4.0*M_PI,2));
  
  val[0] = std::abs((1.0/4.0)*K*nh*sqrt(2)*std::pow(m,2) 
		    *sqrt(m/(kb*T(&x[1])))*exp(-(1.0/2.0)*m*std::pow(x[0],2)/(kb*T(&x[1])))
		    *(-3*kb*T(&x[1])+m*std::pow(x[0],2) )
		    /(std::pow(M_PI, 3.0/2.0)*std::pow(kb*T(&x[1]), 3))
		    *std::pow(x[0],2)); // Jacobian
  
  val[1] = val[1];
  
}


// Code snippet needed for object initalization
// Function always needed!
extern "C" {
  Operator *load()  { return new RV2D;  }
}



