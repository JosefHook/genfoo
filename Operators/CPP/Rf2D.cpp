// $Id$ 
// Drift and diffusion coefficient file for 
// the Ornstein-Uhlenbeck process in 3D 
// By Josef Höök Copyright 2010 (C) all rights reserved
//
#include <cmath>
#include "Operator.hpp"
#include "RF2DParams.hpp"

#define DIM 2

#include <iostream>

// The operator class 
class RF2D : public Operator
{

public:
  
  RF2D();
  void evalDrift(double* val, const double* x) const;
  void evalDiffusion(double* val, const double* x)const;
  void evalIC(double* val, const double* x) const ;
  void evalSource(double* val, const double* x) const;


};



// Default constructor
RF2D::RF2D()
{
  _dim = DIM;
  _info = "C + Q operator in (r,v, xi) ";

}



// x[0] = v/vth = x
// x[1] = xi

  
void RF2D::evalIC(double* values, const double* x) const
{


  
  //    values[0] = 10.0*1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-( (x[0]-0.0)*(x[0]-0.0) + (x[1]-0.0)*(x[1]-0.0))/0.2 ) ;
 
   values[0] = 10.0*1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-( (x[0]-0.0)*(x[0]-0.0) )/0.2 ) ;

 //  values[0] = std::pow(sin(x[0]),2)*std::pow(sin(x[1]),2);
  //values[0] = nh*std::pow(m/(2.0*M_PI*kb*T),3.0/2.0)*exp(- m*std::pow(x[0],2)/(2*kb*T) )
  //  *1/(2.0*M_PI)*sqrt(1/(2.0*M_PI))*exp(-std::pow(x[1],2)/0.2 );
 //  values[0] = nh*std::pow(m/(2.0*M_PI*kb*T),3.0/2.0)*exp(- m*std::pow(x[0],2)/(2*kb*T) );
  
}


// Drift to (x,y,z)=(5,5,5)
//populate drift vector  
void RF2D::evalDrift(double* values, const double* x) const
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
   
	//   long double alpha = -Cp*std::pow(lp,2.0)*( 1.0 + m/mp )*Gp 
	//  -Ce*std::pow(le,2.0)*( 1.0 + m/mee )*Ge 
	// + 1/( 2.0*vth*x[0] )*gamma;
	
        long double alpha  = -Cp*std::pow(lp,2)*( 1.0 + m/mp )*Gp
	  -Ce*std::pow(le,2)*(1.0  + m/mee )*Ge 
	+ Cp/(2.0*std::pow(vth*x[0],2))*( erf(x[0]) -  Gp  ) 
	+ Ce/(2.0*std::pow(vth*x[0],2))*( erf(lelp*x[0])  - Ge  );
	
   long double q01 = K;
   
     
   
   // Pitch-Angle scattering
   //long double gamma = ((Cp*lp*lp*lp)/(x[0]))*( erf(x[0]) - Gp   )
   //  + ((Ce*lp*lp*lp)/(x[0]))*( erf(x[0]) - Ge   );
   

   values[0]= 1/(std::pow(vth, 2))*(   alpha + q01/(vth*x[0])*( 1 + std::pow(x[1],2))  ); //A2
   values[1]= 1/(std::pow(vth, 2))*(  q01/(std::pow(x[0],2))*x[1]*(1 - 3.0*std::pow(x[1],2))
					- 1/(2.0*std::pow(x[0],2))*x[1]*gamma );//A3



   /*
   std::cout	      << "DRIFT:  alpha: " << alpha
		      << " Beta: " << beta 
		      << " Gamma: " << gamma 
		      << " A1: " << values[0] 
		      << " A2: " << values[1] 
		      << std::endl;

   */
    
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
void RF2D::evalDiffusion(double* values, const double* x) const
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
        
    long double q01 = K;

   //   long double gamma = ((Cp*lp*lp*lp)/(x[0]))*( erf(x[0]) - Gp   )
   //+ ((Ce*lp*lp*lp)/(x[0]))*( erf(x[0]) - Ge   );
   
   // B = Matrix(2, 2, {
   // (1, 1) = g(x, xi),
   // (1, 2) = k(x)*o(xi),
   // (2, 1) = k(x)*o(xi),
   //  (2, 2) = m(x, xi)
   // })
   // g(x, xi ) = \frac{1}{v_{th}^2}\left(\langle (\Delta v_{||})^2\rangle  + 4 q_{0,1} - 2q_{0,1}(1+\xi^2)\right)
   // k(x) = \frac{2}{v_{th}^2}   \frac{q_{0,1}}{ x}
   // o(xi) = (2\xi^2 -1)\xi  
   // m(x,xi) = \frac{1}{v_{th}^2} \left( \frac{\langle (\Delta v_{\perp})^2\rangle}{2  x^2 }( 1 - \xi^2) +   \frac{2 q_{0,1}}{ x^2}(1 + \xi^2) \right)



    // 
    // Shape of Sigma matrix : 
    // Matrix(2, 2, {
    // (1, 1) = sqrt(g(x, xi))
    // (1, 2) = k(x)*o(xi)/sqrt(g(x, xi))
    // (2, 1) = 0
    // (2, 2) = sqrt(  (m(x, xi)*g(x, xi)-k(x)^2*o(xi)^2)  /g(x, xi) )
    //  })



    /*
    long double g = 1/(std::pow(vth,2.0))*(  beta  + 4.0*q01 
					     - 2.0*q01*(  1.0 + std::pow(x[1],2.0)  ) ); 
    long double k = 2.0/(std::pow(vth, 2.0))*q01/x[0];
    long double o = ( std::pow(x[1],2)  -1 )*x[1];
    long double m = 1/(std::pow(vth,2.0))*( gamma/( 2.0*std::pow(x[0],2.0)  ) 
					    *( 1.0 - std::pow(x[1], 2.0) )  
					    + 2.0*q01/( std::pow(x[0],2.0) )
					    *std::pow(x[1],2)*( 1.0 - std::pow(x[1], 2.0) ) ); 

					    */


    long double g = 1/(std::pow(vth,2.0))*(  beta  + 4.0*q01 
					     - 2.0*q01*(  1.0 + std::pow(x[1],2.0)  ) ); 
    long double k = 2.0/(std::pow(vth, 2.0))*q01/x[0];
    long double o = ( std::pow(x[1],2)  -1 )*x[1];
    long double m = 1/(std::pow(vth,2.0))*( gamma/( 2.0*std::pow(x[0],2.0)  ) 
					    *( 1.0 - std::pow(x[1], 2.0) )  
					    + 2.0*q01/( std::pow(x[0],2.0) )
					    *std::pow(x[1],2)*( 1.0 - std::pow(x[1], 2.0) ) ); 
    
    long double m1 = gamma/( 2.0*std::pow(x[0],2.0) )*( 1.0 - std::pow(x[1], 2.0) );
    long double m2 = 2.0*q01/( std::pow(x[0],2.0) )*std::pow(x[1],2)*( 1.0 - std::pow(x[1], 2.0) );


    long double g1 = 1/(std::pow(vth,2.0))*(  beta  );

    long double g2 = 1/(std::pow(vth,2.0))*(  4.0*q01 
						 - 2.0*q01*(  1.0 + std::pow(x[1],2.0)  ) ); 



    // Factorization B = L^TL 
   
      values[0] = sqrt(g);     // B11
      values[1] = k*o/sqrt(g); // B12
      values[2] = 0.0; //B21 
      values[3] = sqrt(( m*g -k*k*o*o ) /g ); //B22
   



    // Cholesky factorisation B = LL^T
      /*
    values[0] = sqrt(g);     // B11
    values[1] = 0;           // B12
    values[2] = k*o/sqrt(g);              // B21
    values[3] = sqrt(( m*g -k*k*o*o ) /g ); // B22
      */

    /*
    std::cout	      << "DIFF: Beta: " << beta
		      << " Gamma: " << gamma 
		      << " g: " << g 
		      << " g1: " << g1 
		      << " g2: " << g2
		      << " k: " << k
		      << " o: " << o 
		      << " m: " << m 
		      << " m1: " << m1/(std::pow(vth,2.0)) 
		      << " m2: " << m2/(std::pow(vth,2.0)) 
		      << " v11: " << values[0] 
		      << " v12: " << values[1] 
		      << " v21: " << values[2] 
      		      << " v22: " << values[3] 
		      << std::endl;

    */


  }


void RF2D::evalSource(double* val, const double* x) const {


  //K = Power/(3.0*m*nh) =>
  // k*nh = Power/3m , where m = 2*mp
  // vth^2 = 2*kb*Tp/mp

  val[0] = 2.0*K*nh*exp(-std::pow(x[1],2))*(-3.0 + 2.0*std::pow(x[1],2) )
    /(std::pow(M_PI, 3.0/2.0)*std::pow(vth, 5));
  
  
  //  Normalized with /sqrt(nh)
  // val[0] = std::abs(2.0*K*std::pow(nh,0.5)*exp(-std::pow(x[1],2))*(-3.0 + 2.0*std::pow(x[1],2) )
  // 		    /(std::pow(M_PI, 3.0/2.0)*std::pow(vth, 5)));
  

  /*
  val[0] = 2.0*exp(-std::pow(x[0],2))*(-3.0 + 2.0*std::pow(x[0],2) )
    /(std::pow(M_PI, 3.0/2.0));
  */


  //Unormalized
  /*  val[0] = std::abs((1.0/4.0)*K*nh*sqrt(2)*std::pow(m,2) 
		      *sqrt(m/(kb*T))*exp(-(1.0/2.0)*m*std::pow(x[0],2)/(kb*T))
		    *(-3*kb*T+m*std::pow(x[0],2) )
		    /(std::pow(M_PI, 3.0/2.0)*std::pow(kb*T, 3))
		    *std::pow(x[0],2)); // Jacobian
		    */       
  ///std::pow(vth, 2)
	       


  
};

// Code snippet needed for object initalization
// Function always needed!
extern "C" {
  Operator *load()  { return new RF2D;  }
}



