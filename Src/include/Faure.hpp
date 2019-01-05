#ifndef FAURE_H
#define FAURE_H
#pragma once

#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;


//
// Implements Faure sequence, with Linear scrambling according to 
// Matousek, Tezuka
//
class Faure
{
public: 
  int init(int dim, int N); // Initialize sequence and scrambling
  int init(int dim, int N, int maxdig); // Initialize sequence and scrambling
  int scramble(void);                   // Scramble sequence
  int scramble(int maxdig);             // Scramble sequence

  matrix* rand(M); // Return N quasi random numbers   
  
private:
  int s;                // Dimension
  int atmost;           // Largest value of N
  int maxdig;           // Maximum number of scrambled digits 
  int qs;               // Smallest prime >= s
  int nextn;            // The indx number of the next quasi-vector 
  int testn;
  int hisum;
  double rqs;           // 1/qs;
  double pgtemp; 
  double ztemp;
  matrix *quasi;        // Private hold of the last quasi sequence
  matrix *scoef;        // Scrambling generating matrices
  matrix *coef;         // Binomial coefficient  



};




#endif

