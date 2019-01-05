//
// Implementation file for ReflectedBoundaryCondition class for MC particles
// by Josef Höök, joh@kth.se, Copyright ... 2010
//
#include "Particle.hpp"
#include "ReflectedBoundaryCondition.hpp"



int ReflectedBoundaryCondition::apply(Particle *p) {
  int row = p->size1(); 
  int col = p->size2();
  // Number of Cols decide the dimensionality of the problem.
  
  // SpaceBoundaries are ordered as
  // [x_min, y_min, ..., x_max, y_max, ...]
  // _border[0...n]
  // Reflect space

  // Iterate over particles
  for(int i=0; i<row; i++) {

    for (int n=0; n< col; n++)  {


    // Reflect at lower boundaries
      //    std::cout << "i: " << i 
      //	      << " n: " << n 
      //	      << " _border[n]: " << _border[n] 
      //	      << " _border[col+n]: " << _border[col+n] 
      //	      << "  col: " << col 
      //	      << "  row: " << row
      //       << std::endl;

      int iterator = 0;
      while(true) {
	// Reflect lower bound
	if(  (*p) (i,n) < _border[n] ) {
	  (*p) (i,n) = 2.0*_border[n] -(*p) (i,n);
	  // Upper bounds
	} else if( (*p) (i,n) > _border[col+n] ) {  
	  
	  (*p) (i,n) = 2.0*_border[col+n] -(*p) (i,n);
	  
	} else
	  break;

	if(iterator>1e3) {
	  break;
	} 
	iterator++;
      }     
      
      // Reflect at upper boundaries
      //while(true) {
	// Reflect
      //	if(  ) {
      //	} else
      //	  break;
      // }     
    }
    
  }





}



