// Implements reflected Boundary condition
// for MC particles
#pragma once
#ifndef REFLECTEDBOUNDARYCONDITION_H
#define REFLECTEDBOUNDARYCONDITION_H



//Todo add support for general : public BoundaryCondition
class ReflectedBoundaryCondition 
{

public: 
  ReflectedBoundaryCondition();
  ReflectedBoundaryCondition(std::vector<double> br) : _border(br) {}

   int apply(Particle *p);
  private:
  std::vector<double> _border;

};


#endif
