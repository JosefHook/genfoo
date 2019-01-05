//
// Test Ball resampling 
// on two Gaussians


#include <iostream>
#include "BallResample.hpp"
#include "mArray.hpp"
#include "Particle.hpp"

#include <boost/random.hpp>



using namespace std;

int main(void)
{

  int N = 50;
  int D = 2;
  RandomNumberGenerator *rng = new RandomNumberGenerator();
  boost::mt19937 u_seed;
  boost::uniform_01<boost::mt19937> U(u_seed);
  
  BallResample *br = new BallResample(3, 0.50, rng);
  Particle *p = new Particle(N,D);
  mArray *weight = new mArray(N,1);
  weight->clear();

  // populate particles and weights
 for(int i=0; i< N; i++) {
   for(int j=0; j<D; j++) 
     (*p)(i,j) =U();
   (*weight)(i,0) = -2.0*(double)(i%2) +1.0;
 }
  



 p->set_pos(p->size1());
  WeightedParticle *wp = new WeightedParticle(p,weight);     

  // Print out ball-resmple
 for(int i=0; i< N; i++) {
   cout << "( "<<  (*p)(i,0); 
   for(int j=1; j<D; j++) 
     cout << ", " << (*p)(i,j); 
   
   cout << " )\t\t(" << (*weight)(i,0) 	<< ") " << endl; 
 }


 br->resample(wp);
 
 cout << " -----------------AFTER RESAMPLING -------------------------" << endl;
  // Print out ball-resmple
 for(int i=0; i< N; i++) {
   cout << "( "<<  (*p)(i,0); 
   for(int j=1; j<D; j++) 
     cout << ", " << (*p)(i,j); 
   
   cout << " )\t\t(" << (*weight)(i,0) 	<< ") " << endl; 
 }



  return 0;

}
