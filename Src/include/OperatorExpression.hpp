//
// Small template class for casting Operator types to Expressions
//

#pragma once
#ifndef OPERATOREXPRESSION_H
#define OPERATOREXPRESSION_H



using namespace dolfin;
//#include "Operator.hpp"



// Drift Expression Operator wrapper
class DriftExpressionScalar : public Expression
{
public:
  DriftExpressionScalar(Operator *ops) : _ops(ops) {
  }
  void eval(Array<double>& values, const Array<double>& x) const {
    _ops->evalDrift(&values[0], &x[0]);
  }
private:
  Operator *_ops; 
};


// Diffusion Expression Operator wrapper
class DiffusionExpressionScalar : public Expression
{
 public:
  DiffusionExpressionScalar(Operator *ops) : _ops(ops) {
  }
  void eval(Array<double>&  values, const Array<double>& x) const {
    
    _ops->evalDiffusion(&values[0], &x[0]);
    
    
  }
private:
  Operator *_ops; 
};





// Initial Condition Expression Operator wrapper
template <int T> class ICExpression : public Expression
{
 public:
  ICExpression(Operator *ops) : _ops(ops) {

  }
  void eval(Array<double>& values, const Array<double>& x) const {
    _ops->evalIC(&values[0], &x[0]);

  }
 private:
  Operator *_ops; 
};




// Drift Expression Operator wrapper
template <int T> class DriftExpression : public Expression
{
 public:
    DriftExpression(Operator *ops) : Expression(T), _ops(ops) {
    }
  void eval(Array<double>& values, const Array<double>& x) const {
           _ops->evalDrift(&values[0], &x[0]);
  }
 private:
  Operator *_ops; 
};


// Diffusion Expression Operator wrapper
template <int T, int D> class DiffusionExpression : public Expression
{
 public:
 
     DiffusionExpression(Operator *ops) : Expression(T,D), _ops(ops) {
  }
  void eval(Array<double>& values, const Array<double>& x) const {

    _ops->evalDiffusion(&values[0], &x[0]);
    }
 private:
  Operator *_ops; 
};







class OperatorExpression 
{
  
public:
  OperatorExpression(Operator *ops) {


    std::cout << " Dim in Operator Expression is " << ops->getDim() << std::endl;
    switch(ops->getDim())
      {
      case 1 : 
	ic = new ICExpression<1>(ops) ;
	dr = new DriftExpressionScalar(ops);
	di = new DiffusionExpressionScalar(ops);
	break;
      case 2 :
	ic = new ICExpression<2>(ops);
	dr = new DriftExpression<2>(ops);
	di = new DiffusionExpression<2,2>(ops);
	break;
      case 3 :
	ic = new ICExpression<3>(ops);
	dr = new DriftExpression<3>(ops);
	di = new DiffusionExpression<3,3>(ops);
	break;
      default:
	ic = new ICExpression<1>(ops);
	dr = new DriftExpressionScalar(ops);
	di = new DiffusionExpressionScalar(ops);
      }
  }      
  
  Expression* getIC() { return ic;}
  Expression* getDrift() { return dr; }
  Expression* getDiffusion() { return di;}   
  

  ~OperatorExpression() {
	delete ic;
	delete dr; 
	delete di;
  }      
  
  
  
private:
  Expression *ic;
  Expression *dr;
  Expression *di;
    
};






#endif
