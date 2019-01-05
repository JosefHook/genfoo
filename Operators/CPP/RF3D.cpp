// $Id$ 
// Drift and diffusion coefficient file for 
// the Ornstein-Uhlenbeck process in 3D 
// By Josef Höök Copyright 2010 (C) all rights reserved
//
#include <cmath>
#include "Operator.hpp"
#include "RF3DParams.hpp"

#define DIM 3



// The operator class 
class RF3D : public Operator
{

public:
  
  RF3D();
  void evalDrift(double* val, const double* x) const;
  void evalDiffusion(double* val, const double* x)const;
  void evalIC(double* val, const double* x) const ;
  void evalSource(double* val, const double* x) const;


};



// Default constructor
RF3D::RF3D()
{
  _dim = DIM;
  _info = "C + Q operator in (r,v, xi) ";

}


// x[0] = r
// x[1] = v/vth = x
// x[2] = xi

  
void RF3D::evalIC(double* values, const double* x) const
{
  
  //  values[0] = 20.0*1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-( (x[0]-0.0)*(x[0]-0.0) + (x[1]-0.0)*(x[1]-0.0) + (x[2]-0.0)*(x[2]-0.0))/0.2 ) ; 
  
  
   values[0] = 10.0*1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-(  (x[0]-0.0)*(x[0]-0.0) + (x[1]-0.0)*(x[1]-0.0) )/0.2 ) ;


}


// Drift to (x,y,z)=(5,5,5)
//populate drift vector  
void RF3D::evalDrift(double* values, const double* x) const
  {

  long double Gp =( erf(x[0]) - 0.2e1*x[0] 
		     *exp(  -std::pow(x[0],2) )
		     *std::pow(M_PI,(-0.1e1/0.2e1)) )
     /std::pow(x[0],2)/0.2e1;
   
   long double Ge =(erf(lelp*x[0]) - 0.2e1*lelp*x[0]
		    *exp(-std::pow(lelp*x[0],2))
		    *std::pow(M_PI,(-0.1e1/0.2e1))  )
     /std::pow(lelp,2)/std::pow(x[0],2)/0.2e1; 
    
   long double gamma = ((Cp/vth)/x[1])*( erf(x[1]) - Gp   ) 
     + ((Ce/vth)/x[1])*( erf(x[1]) - Ge   );
   
   long double beta = ((Cp/vth)*x[1])*Gp + ((Ce/vth)*x[1])*Ge;
   
   long double alpha = -Cp*std::pow(lp,2.0)*( 1 + m/mp )*Gp 
     -Ce*std::pow(le,2.0)*( 1 + m/mee )*Ge 
     + 1/( 2.0*vth*x[1] )*gamma;
   
       
   long double q01 = K;
   
     
   
   // Pitch-Angle scattering
   //long double gamma = ((Cp*lp*lp*lp)/(x[0]))*( erf(x[0]) - Gp   )
   //  + ((Ce*lp*lp*lp)/(x[0]))*( erf(x[0]) - Ge   );
   

   values[0] =  1/x[0]; //A1 
   values[1]= 1/(std::pow(vth, 2.0))*(   alpha + q01/(vth*x[1])*( 1 + std::pow(x[2],2.0))  ); //A2
   values[2]= 1/(std::pow(vth, 2.0))
     *(  q01/(std::pow(x[1],2.0))*x[2]*(1 - 3.0*std::pow(x[2],2.0))
	 - 1/(2.0*std::pow(x[1],2.0))*x[2]*gamma );//A3


    
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

void RF3D::evalDiffusion(double* values, const double* x) const
  {
    // Ordering of the variables
    // r=x[0], x = x[1], xi = x[2]
    
    long double Gp =( erf(x[1]) - 0.2e1*x[1] 
		      *exp(  -std::pow(x[1],2) )
		      *std::pow(M_PI,(-0.1e1/0.2e1)) )
      /std::pow(x[1],2)/0.2e1;
    
    long double Ge =(erf(lelp*x[1]) - 0.2e1*lelp*x[1]
		     *exp(-std::pow(lelp*x[1],2))
		     *std::pow(M_PI,(-0.1e1/0.2e1))  )
      /std::pow(lelp,2)/std::pow(x[1],2)/0.2e1; 
    
    
    long double beta = ((Cp/vth)*x[1])*Gp + ((Ce/vth)*x[1])*Ge;
    
    long double gamma = ((Cp/vth)/x[1])*( erf(x[1]) - Gp   ) 
      + ((Ce/vth)/x[1])*( erf(x[1]) - Ge   );
    
    
    long double q01 = K;


    long double g = 1/(std::pow(vth,2.0))*(  beta  + 4.0*q01 
					     - 2.0*q01*(  1.0 + std::pow(x[2],2.0)  ) ); 
    long double k = 2.0/(std::pow(vth, 2.0))*q01/x[1];
    long double o = ( std::pow(x[2],2)  -1 )*x[2];
    long double m = 1/(std::pow(vth,2.0))*( gamma/( 2.0*std::pow(x[1],2.0)  ) 
					    *( 1.0 - std::pow(x[2], 2.0) )  
					    + 2.0*q01/( std::pow(x[1],2.0) )
					    *std::pow(x[2],2)*( 1.0 - std::pow(x[2], 2.0) ) ); 

    values[0] = std::sqrt(2.0);       // B11 
    values[1] = 0.0;                  // B12
    values[2] = 0.0;                  // B13

    values[3] = 0.0;                  // B21
    values[4] = sqrt(g);          // B22 
    values[5] = 0.0;                  // B23
    
    values[6] = 0.0;                  // B31 
    values[7] = k*o/std::sqrt(g);       // B32 
    values[8] = std::sqrt(( m*g -k*k*o*o ) /g );  // B33     





  }


void RF3D::evalSource(double* val, const double* x) const {

  val[0] = 0.0;
  
};

// Code snippet needed for object initalization
// Function always needed!
extern "C" {
  Operator *load()  { return new RF3D;  }
}



