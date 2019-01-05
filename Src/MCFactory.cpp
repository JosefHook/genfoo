//
// GenFoo-MC 
// A GENeralized FOkker-Planck sOlver -
// A daptive Monte Carlo 
//  By Josef Höök, joh@kth.se copyright (C) all rights reserved
// 
#include <iostream>
#include "Params.hpp"
#include "MonteCarlo.hpp"
#include "Particle.hpp"
#include "H5PartIO.hpp"
#include "DatIO.hpp"
#include "Operator.hpp"
#include "DynamicLoader.hpp"
#include "ReflectedBoundaryCondition.hpp"
#include "RandomNumberGenerator.hpp"

// MPI related 
#include <mpi.h>


using namespace std;

int
MCFactory(Params pm) 
{

  double dt = pm.getDT();
  double END = pm.getTimeEnd();
  double t = pm.getTimeBegin();
  
  int N = pm.getNP();
  int D = pm.getDimension();


  //
  // Operators
  // load operators dynamically from library

  DynamicLoader dl;
  if(dl.Load(  ((std::string) pm.getLibraryFile()).c_str() ))
    exit(1);

  Operator *ops;
  Operator *(*load)() = (Operator*(*)())(dl.bootstrap());
  // Create new instance of Operator object
  ops = (*load)();

  // Pass XML doc to operator
  ops->registerXMLDoc(pm.getXMLDoc());


  if(ops->getDim()!=D) {

    std::cout << "Dimension missmatch between library and parameters!"
	      << std::endl;
    throw "Dimension missmatch between library and parameters!";
  }
  Particle *p = new Particle(N,D);
  p->clear();
  p->set_pos(p->size1());


  //
  // Initialize MPI 
  //
  MPI::Init();

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
  std::vector<double> boarder(pm.getBoundaryMinValues());
  std::vector<double> tmp(pm.getBoundaryMaxValues());
  
  boarder.insert( boarder.end(),
		  tmp.begin(),
		  tmp.end() );  
  

  ReflectedBoundaryCondition *rbc = 
    new ReflectedBoundaryCondition(boarder);
  
    //cf.SpaceBoundaries);
 
 //
 // Assemble the MC factory
 //
 mc->assemble(ops, rbc, rng);
 



    mc->init(p); // initialize distribution provided in ops
    cout << "Initialization.....done " << endl;
    
    

  H5PartIO file("MCFactoryOut.h5part", N);
  H5PartIO filef("MCFactoryOutFinalT.h5part", N);
  DatIO datfile("MCFactoryOut");
  
  int t_i=0;
  file << (*p);
  
  while(t<END)
    {
      
      mc->step(p,dt);
      file << (*p);
      datfile << (*p); 
  
      if( (t_i %100 )==0) {
	cout << "iter: " << t_i << " t: " << t 
	     << " T: " << END << " dt: " << dt <<endl;
      }
      t +=dt;
      t_i++;
      
    }





  filef << (*p);
  // Fast hack
  file << (*p);
  
  filef.close();

  file.close();
  

  delete mc;
  delete p;
  delete rng;

  return 0;

}

