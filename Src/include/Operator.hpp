// By Josef Höök, 2010-05-31
// $Id$
//
// 
// Operator interface template class
//
//
// Example usage:
// class Foo : public Operator<float, 3> { ... };
// Defines and implements a three dimensional problem 
// in float precision
// 
#pragma once
#ifndef OPERATOR_H
#define OPERATOR_H


#include <string>
#include <libxml/parser.h>

#define UNUSED(x) ((void)x)


class Operator
{

 public:
  // here we should have function pointers
  virtual void evalDrift(double* val, const double* x) const = 0;
  virtual void evalDiffusion(double* val, const double* x) const = 0;
  virtual void evalIC(double* val, const double* x) const = 0;
  virtual void evalSource(double* val, const double* x) const  = 0;
  virtual void registerXMLDoc(xmlDocPtr doc) { UNUSED(doc); }
  virtual std::string getInfo() const { return _info; };
  int getDim() const { return _dim; }
  int getRandomPerParticle() const { return _rrp; }

protected:

  int _dim, _rrp;
  std::string _info;


};



#endif
