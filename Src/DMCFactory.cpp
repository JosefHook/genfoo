//
// GenFoo-dMC 
// A delta-f generalized Fokker-Planck solver
// By Josef Höök, joh@kth.se copyright (C) all rights reserved
// 
#include <iostream>
#include "Params.hpp"
#include "MonteCarlo.hpp"
#include "Particle.hpp"
#include "mArray.hpp"
#include "WeightedParticle.hpp"
#include "H5PartIO.hpp"
#include "DatIO.hpp"
#include "Operator.hpp"
#include "DynamicLoader.hpp"
#include "BallResample.hpp"
#include "ReflectedBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"
#ifdef _OPENMP_
#include <omp.h>
#endif 

using namespace std;

//
// Main 
//
int
DMCFactory( Params pm )
{

  double dt = pm.getDT();
  double END = pm.getTimeEnd();
  double t = pm.getTimeBegin();
    int N = pm.getNP();
  int D = pm.getDimension();
  double dN = ((double) N)/ceil(END/dt);

  
  cout << "dN: " << dN << endl;
  cout << "NParticles: " << N << endl;

  //
  // Temporary security check, exit if number of particle are > 1e12
  //
  if(N>1e12) {
    cout << "To many particles!, max:1e12, exiting" << endl;
    return(1);
      }

  
  //
  // Load operators dynamically from library
  //
  DynamicLoader dl;
  if(dl.Load(   ((std::string) pm.getLibraryFile()).c_str()   ))
    exit(1);

  Operator *ops;
  Operator *(*load)() = (Operator*(*)())(dl.bootstrap());
  // Create new instance of Operator object
  ops = (*load)();



  // Pass XML doc to operator
  ops->registerXMLDoc(pm.getXMLDoc());



  //
  // Check if dimensions match between parameter file and library 
  //
  if(ops->getDim()!=D) {
    std::cout << "Dimension missmatch between library and parameters!"
	      << std::endl;
    throw "Dimension missmatch between library and parameters!";
  }


  //
  // Define particle arrays
  //
  Particle *p = new Particle(N,D);
  Particle *samples = new Particle(N,D);
  p->clear();
  samples->clear();


  //
  // Define Randomnumbergenerator
  //
 RandomNumberGenerator *rng = new RandomNumberGenerator();


 //
 // Define Monte Carlo solver
 //
  MonteCarlo *mc = new MonteCarlo();


 
 //
 // Setup reflected boundary conditions 
 //


 //
 // Setup reflected boundary conditions 
 //
  std::vector<double> boarder(pm.getBoundaryMinValues());
  std::vector<double> tmp(pm.getBoundaryMaxValues());
  
  boarder.insert( boarder.end(),
		  tmp.begin(),
		  tmp.end() );  

 ReflectedBoundaryCondition *rbc = 
   new ReflectedBoundaryCondition(boarder);
 
 //
 // Assemble the MC factory
 //
 mc->assemble(ops, rbc, rng);
 

 // Setup some time-variables
  int t_i=0;
 int collected = 0;
  
 //
 // Create weight arrays
 //
 mArray *sample_weight = new mArray(N,1);
 sample_weight->clear();
 mArray *weight = new mArray(N,1);
 weight->clear();

 //
 // Define Weighted particles
 //
 WeightedParticle *ws = new WeightedParticle(samples, sample_weight); 
 WeightedParticle *wp = new WeightedParticle(p,weight);     
 
 //
 // Wanted behavior of BallResample:
 // BallResample *br  = new BallResample(nneighbours,reduce, rng);
 // BallResample *br  = new BallResample(rng);
 // br->resample(wp);
 // 
 BallResample *br = new BallResample(3, 0.50, rng);

 //
 // Pre-sample the source term
 //
 mc->initsource(ws, 1e-3, 1000);

   
   
 //
 // Open output file, (H5Part file format)
 //
 H5PartIO file("dMCFactoryOut.h5part", N);
 DatIO datfile("dMCFactoryOut");
  
 /////////////////////////////////////////////////////
 //  time-step
 ////////////////////////////////////////////////////

  while(t<END)
    {

   
      //
      // Sample dN new particles from ws to wp.
      // 
      mc->sample(wp,ws, dN);

      //
      // Perform MC-step and apply reflected boundary condition
      //
      mc->step(p,dt);

      //
      // Flush particles
      //
      file.setNumParticles(p->pos());
      file << (*wp);   

      datfile<< (*wp);

      //
      // Write something
      //
      cout << "iter: " << t_i << " t: " << t 
	   << " T: " << END << " dt: " << dt 
	   << " NP: " << p->pos()
	   << " N: " << N <<endl;
            

      //
      // Resample particles
      //

      //if((t_i % 100) ==0) {

      /*
      if ( (( t_i % 500)==0 ) && t_i >1)  {	
	int pos_pre = p->pos();
	cout << "----------------------------------------" << endl;
	cout << "Pre pos: " << pos_pre << endl;	
	//
	//	
	int maxn = 30;    // Max neigbors
	double red = 0.3;
	double radius = 0.9;
	br->reduce(wp,radius, maxn, red );
	cout << " Pos : " << p->pos() << endl;
	// Pickup new pointers old ones have been destroyed in reduce.
	p  = wp->getParticle(); 
	weight = wp->getWeight();
	cout << " Pos new pointer : " << p->pos() << endl;
	

      }
      */
      t +=dt;
      t_i++;
    }
  
  file.close();
  
  delete mc;
  
  return 0;

}

