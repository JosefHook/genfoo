//
//
// Monte Carlo Base class
//
//


#include "Particle.hpp"
#include "mArray.hpp"
#include "Operator.hpp"
#include "RandomNumberGenerator.hpp"
#include "ReflectedBoundaryCondition.hpp"
#include "MonteCarlo.hpp"

#include "H5PartIO.hpp"
#include <boost/math/special_functions/sign.hpp>

// Assemble MonteCarlo pieces
int MonteCarlo::assemble( Operator *ops, ReflectedBoundaryCondition *bc,
			  RandomNumberGenerator *rng ) {
  
    _drift = new double[ops->getDim()];
    _diff = new double[(ops->getDim()*ops->getDim())];

  _ops = ops; 
  _bc = bc;
  _dW = rng;
  return 0;
  
  }


//
// Adaptive time-step
//
// 
int MonteCarlo::adaptivestep( Particle *p, mArray *dt ) {
  

  int dim = _ops->getDim();
  
  double *ptr = &(p->data())[0]; 

  long long p_row = p->size1();

  for(int i=0; i<dim; i++) {
    _drift[i] = 0.0;
    for(int j=0; j<dim; j++)
      _diff[i*dim + j] = 0.0;
  }
  

  for(int i=0; i<p_row; i++) {
    
   
    // The Particle type is ordered in row major
    // m x n  dimensional matrix M_ij  => (i x n + j)
    // n = columns
    //  int col = p.size2();
    // We want pointers to the first row for each column
    
    _ops->evalDrift(_drift, &ptr[i*dim]);  
    _ops->evalDiffusion(_diff, &ptr[i*dim]);  
    
    
    for(int l=0; l<dim; l++) {
      (*p) (i,l) +=  _drift[l]*(*dt)(i,0);
      for(int j=0; j<dim; j++)
	(*p) (i,l) += _diff[l*dim+j]*sqrt((*dt)(i,0))*( (*_dW)()); 
      
    }
    
    
    
    
    
 
  }
     // Apply boundary conditions
    _bc->apply(p);
  
  
  return 0;
  
}



//
// Regular time-step
//
// 
int MonteCarlo::step( Particle *p, double dt ) {
  

  int dim = _ops->getDim();
  
  double *ptr = &(p->data())[0]; 

  long long p_pos = p->pos();

  for(int i=0; i<dim; i++) {
    _drift[i] = 0.0;
    for(int j=0; j<dim; j++)
      _diff[i*dim + j] = 0.0;
 }
  
  for(int i=0; i<p_pos; i++) {
    
    
    // The Particle type is ordered in row major
    // m x n  dimensional matrix M_ij  => (i x n + j)
    // n = columns
    //  int col = p.size2();
    // We want pointers to the first row for each column
    
    _ops->evalDrift(_drift, &ptr[i*dim]);  
    _ops->evalDiffusion(_diff, &ptr[i*dim]);  
    
    for(int l=0; l<dim; l++) {
      (*p) (i,l) +=  _drift[l]*dt;
      for(int j=0; j<dim; j++)
	(*p) (i,l) += _diff[l*dim+j]*sqrt(dt)*( (*_dW)()); 
      
    }
    
    
    
    
  }
  // Apply boundary conditions
    _bc->apply(p);
  

     
  return 0;
  
}



// Sample Initial condition using Metropolis-Hastings
int MonteCarlo::init(Particle *p) {

  int dim = _ops->getDim();  
  // Zero out Particle matri
  double *ptr = &(p->data())[0]; 


 for(int i = 0; i < p->size1 (); i++ ) 
   for(int j=0; j<  p->size2() ; j++)
     (*p)(i,j) = 0.0;
 
 // boost::uniform_01<> ud; 
 boost::mt19937 u_seed;
 boost::uniform_01<boost::mt19937> U(u_seed);
 

//boost::variate_generator<boost::mt19937&, boost::uniform_01<> >U(rng, ud);
  double burnin = 100;
  double *xp,*xt, Pxp, Pxt; 
  double a,r; 
  double dl = 1e-2;
  
  xp = new double[p->size2()];
  xt = new double[p->size2()];

  H5PartIO file("GenFoo-MetropolisHastings.h5part", p->size1());
  file << *p; 

 for(int b=0; b< burnin; b++) 
   { 
     
     

     for(int i = 0; i < p->size1(); i++ ) 
       {
	 for(int j=0; j<  p->size2(); j++) 
	   {
	     
	   // Create a proposal step x' 
	   // Q(x', x^t)
	     xp[j]  = (*p)(i,j) +    sqrt(dl)*((*_dW)());
	     xt[j] = (*p) (i,j);
	  
	   }


	     // Evaluate P(x')/P(x^t) 
	     _ops->evalIC(&Pxp, xp);
	     _ops->evalIC(&Pxt, xt);
	    
	     
	     // since Q(x', x^t) is the normal distribution 
	     // Q(x', x^t)/Q(x^t, x') = 1
	     
	     a = std::min(Pxp/Pxt, 1.0);
	     r = U();
	     //  std::cout << "xp: " << xp[0] << xp[1] 
	     //	       << " xt: " << xt[0] << xt[0]
	     //	       << "a: " << a 
	     //	       << " Rand: " << r << std::endl;

	 for(int j=0; j<  p->size2(); j++) 
	   {
 
	     if( r<= a ) {
	       // Keep new x'
	       (*p) (i,j) = xp[j];
	     } else {
	       // Keep old pos;
	     }


	   }
	 
       }
	

     // Apply boundary conditions
     _bc->apply(p);
	     

     file << *p; 
 
   }


 file.close();
 delete xp, xt;
 return 0; 
}




// Sample Initial condition using Metropolis-Hastings
int MonteCarlo::initsource( Particle *p,  double dl) {

  int dim = _ops->getDim();  
  // Zero out Particle matri
  double *ptr = &(p->data())[0]; 
  boost::mt19937 u_seed;
  boost::uniform_01<boost::mt19937> U(u_seed);
  
  
//boost::variate_generator<boost::mt19937&, boost::uniform_01<> >U(rng, ud);
  double burnin = 100;
  double *xp,*xt, Pxp, Pxt; 
  double a,r; 
  //  double dl = 1e-2;
  
  xp = new double[p->size2()];
  xt = new double[p->size2()];

  H5PartIO file("GenFoo-InitSource.h5part", p->size1());
  file << *p; 

 for(int b=0; b< burnin; b++) 
   { 
     


     for(int i = 0; i < p->size1(); i++ ) 
       {
	 for(int j=0; j<  p->size2(); j++) 
	   {
	     
	   // Create a proposal step x' 
	   // Q(x', x^t)
	     xp[j]  = (*p)(i,j) +   sqrt(dl)*((*_dW)());
	     xt[j] = (*p) (i,j);
	   }
	     // Evaluate P(x')/P(x^t) 
	     _ops->evalSource(&Pxp, xp);
	     _ops->evalSource(&Pxt, xt);
	     
	     // since Q(x', x^t) is the normal distribution 
	     // Q(x', x^t)/Q(x^t, x') = 1
	     
	     a = std::min(Pxp/Pxt, 1.0);
	     r = U();


	 for(int j=0; j<  p->size2(); j++) 
	   {
 
	     if( r<= a ) {
	       // Keep new x'
	       (*p) (i,j) = xp[j];
	     
	     } else {
	       // Keep old pos;
	     }
	     
	   }
	 
       }
     // Apply boundary conditions
     _bc->apply(p);
;
     
     std::cout << "Init source sample iteration : " << b
	       << ". Burnin expected at : " << burnin << std::endl;
     file << *p; 
 
   }


 file.close();
 delete xp, xt;
 return 0; 
}






// Sample Initial condition using Metropolis-Hastings
int MonteCarlo::initsource( Particle *p,  double dl, int burnin) {

  int dim = _ops->getDim();  
  // Zero out Particle matri
  double *ptr = &(p->data())[0]; 
  boost::mt19937 u_seed;
  boost::uniform_01<boost::mt19937> U(u_seed);
  
  

  double *xp,*xt, Pxp, Pxt; 
  double a,r; 
  //  double dl = 1e-2;
  
  xp = new double[p->size2()];
  xt = new double[p->size2()];

  H5PartIO file("GenFoo-InitSource.h5part", p->size1());
  file << *p; 

 for(int b=0; b< burnin; b++) 
   { 
     


     for(int i = 0; i < p->size1(); i++ ) 
       {
	 for(int j=0; j<  p->size2(); j++) 
	   {
	     
	   // Create a proposal step x' 
	   // Q(x', x^t)
	     xp[j]  = (*p)(i,j) +   sqrt(dl)*((*_dW)());
	     xt[j] = (*p) (i,j);
	   }
	     // Evaluate P(x')/P(x^t) 
	     _ops->evalSource(&Pxp, xp);
	     _ops->evalSource(&Pxt, xt);
	     
	     // since Q(x', x^t) is the normal distribution 
	     // Q(x', x^t)/Q(x^t, x') = 1
	     
	     a = std::min(Pxp/Pxt, 1.0);
	     r = U();


	 for(int j=0; j<  p->size2(); j++) 
	   {
 
	     if( r<= a ) {
	       // Keep new x'
	       (*p) (i,j) = xp[j];
	     
	     } else {
	       // Keep old pos;
	     }
	     
	   }
	 
       }
     // Apply boundary conditions
     _bc->apply(p);
;
     
     std::cout << "Init source sample iteration : " << b
	       << ". Burnin expected at : " << burnin << std::endl;
     file << *p; 
 
   }


 file.close();
 delete xp, xt;
 return 0; 
}





//
// USED IN MAIN
//
// Sample Initial condition using Metropolis-Hastings
int MonteCarlo::initsource( WeightedParticle *wp,  double dl, int burnin) {

  Particle *p; 
  mArray *w; 
  int dim = _ops->getDim();  
  int sxt, sxp;
  double a,r; 
  double *ptr;  
  double *xp,*xt, Pxp, Pxt; 
  boost::mt19937 u_seed;
  boost::uniform_01<boost::mt19937> U(u_seed);
  
  
  p = wp->getParticle();
  w = wp->getWeight();
  ptr = &(p->data())[0]; 

  
  xp = new double[p->size2()];
  xt = new double[p->size2()];

  H5PartIO file("GenFoo-InitSource.h5part", p->size1());
  file << *wp; 

 for(int b=0; b< burnin; b++) 
   { 
     


     for(int i = 0; i < p->size1(); i++ ) 
       {
	 for(int j=0; j<  p->size2(); j++) 
	   {
	     
	   // Create a proposal step x' 
	   // Q(x', x^t)
	     xp[j]  = (*p)(i,j) +   sqrt(dl)*((*_dW)());
	     xt[j] = (*p) (i,j);

	   }
	     // Evaluate P(x')/P(x^t) 
	     _ops->evalSource(&Pxp, xp);
	     _ops->evalSource(&Pxt, xt);

	     // Store signs
	     sxp = boost::math::sign(Pxp);
	     sxt = boost::math::sign(Pxt);
	     Pxp = std::abs(Pxp);
	     Pxt = std::abs(Pxt);

	    
	     
	     // since Q(x', x^t) is the normal distribution 
	     // Q(x', x^t)/Q(x^t, x') = 1
	     
	     a = std::min(Pxp/Pxt, 1.0);
	     r = U();


	 for(int j=0; j<  p->size2(); j++) 
	   {
 
	     if( r<= a ) {
	       // Keep new x'
	       (*p) (i,j) = xp[j];
	       (*w) (i,0) = sxp;
	       	       
	     } else {
	       // Keep old pos;
	       (*w) (i,0) = sxt; // uggly 
	     }
	     
	   }
	 
       }
     // Apply boundary conditions
     _bc->apply(p);

     
     std::cout << "Init source sample iteration : " << b
	       << ". Burnin expected at : " << burnin << std::endl;
     file << *wp; 
 
   }


 file.close();
 delete xp, xt;
 return 0; 
}






//
// Sample source term from Offline-generated point cloud
//
int MonteCarlo::sample(Particle *p, Particle *s, double np) {

  
  long s_row = s->size1();
  long s_col = s->size2();
  long p_row = p->size1();
  long p_col = p->size2();
  
  long  p_start_pos = p->pos();
  long  s_start_pos = s->pos();


  // If array full return
  if( p_start_pos >= p_row) {
    return 0;
  }

  // If for some reason we try to sample 
  // more than allowed return 1
  if( p_start_pos > p_row) {
    std::cout << "Cannot sample more. Array full" << std::endl;
    return 1;
  }



  // TODO FIX SO THAT FRACTIONAL PARTI IS TREATED


  // If np < 1 create one particle with P(X<np) by
  // sampling from uniform distribution
  if(np<1.0) {
    boost::mt19937 u_seed;
    boost::uniform_01<boost::mt19937> U(u_seed);
  
    if(np<U()) 
      np=1.0;
    else {
      std::cout << "Returned without particle realization" << std::endl;
      return 0;
    }

  }


  

  

  // Pick np last particles from samples 
  int I = std::min(p_start_pos+(int)np, p_row);

  for (int i=p_start_pos; i<I; i++) 
    for (int k=0; k<p_col; k++) {
      /*      	      std::cout << "p("<< i <<","<< k << ") : " <<(*p)(i, k) 
		      << "  s("<< i <<","<< k << ") : " <<(*s)(i, k) 
		      << std::endl;*/
      (*p)(i, k) = (*s)(i, k);      
    }
  p->set_pos(p_start_pos+np);
  s->set_pos(s_start_pos+np);
  /*
  std::cout << "s_start_pos: " << s_start_pos << std::endl;
  std::cout << "np " << np << std::endl;
  std::cout << "p->pos() " << p->pos() << std::endl;
  std::cout << "s->pos() " << s->pos() << std::endl;
  */
}






//
// Sample source term from Offline-generated point cloud
//
int MonteCarlo::sample(WeightedParticle *wp, WeightedParticle *ws, double np) {

  Particle *p = wp->getParticle();
  Particle *s = ws->getParticle();
  mArray *p_weight = wp->getWeight();
  mArray *s_weight = ws->getWeight();


  long s_row = s->size1();
  long s_col = s->size2();
  long p_row = p->size1();
  long p_col = p->size2();
  
  long  p_start_pos = p->pos();
  long  s_start_pos = s->pos();


  // If array full return
  if( p_start_pos >= p_row) {
    return 0;
  }

  // If for some reason we try to sample 
  // more than allowed return 1
  if( p_start_pos > p_row) {
    std::cout << "Cannot sample more. Array full" << std::endl;
    return 1;
  }



  // TODO FIX SO THAT FRACTIONAL PARTI IS TREATED


  // If np < 1 create one particle with P(X<np) by
  // sampling from uniform distribution
  if(np<1.0) {
    boost::mt19937 u_seed;
    boost::uniform_01<boost::mt19937> U(u_seed);
  
    if(np<U()) 
      np=1.0;
    else {
      std::cout << "Returned without particle realization" << std::endl;
      return 0;
    }

  }


  

  

  // Pick np last particles from samples 
  int I = std::min(p_start_pos+(int)np, p_row);

  for (int i=p_start_pos; i<I; i++) 
    for (int k=0; k<p_col; k++) {
      /*      	      std::cout << "p("<< i <<","<< k << ") : " <<(*p)(i, k) 
		      << "  s("<< i <<","<< k << ") : " <<(*s)(i, k) 
		      << std::endl;*/
      (*p)(i, k) = (*s)(i, k);     
      (*p_weight)(i,0) = (*s_weight)(i,0);

    }
  p->set_pos(p_start_pos+np);
  s->set_pos(s_start_pos+np);
  /*
  std::cout << "s_start_pos: " << s_start_pos << std::endl;
  std::cout << "np " << np << std::endl;
  std::cout << "p->pos() " << p->pos() << std::endl;
  std::cout << "s->pos() " << s->pos() << std::endl;
  */
}


