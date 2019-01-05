//
// By Josef Höök, 2010-10-07
//
// Interface class for boundary conditions
//
#pragma once
#ifndef BOUNDARYCONDITION_H
#define BOUNDARYCONDITION_H



class BoundaryCondition 
{

public:
  virtual int apply(Particle *p) const { 
    std::cout << " apply() must be implemented " <<std::endl; 
    return(1);
  } 
  virtual ~BoundaryCondition() {}
  
};
  


#endif
