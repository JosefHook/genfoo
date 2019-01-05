//
// Data Input Output 
// Saves data in readable text format: time weight x_1  x_2 x_3 ... x_n  
//
#include "DatIO.hpp"
#include <fstream>
#include <sstream>
DatIO::DatIO(std::string filename) {
  _filename = filename;
  _time_step = 0;
}

void DatIO::setNumParticles(long long n, int D) {
  if(n<=0)
    _np= 1;
  else 
    _np = n;
}


void DatIO::operator<<(const WeightedParticle &wp) {

  
  Particle *p;
  mArray *w;
  std::ostringstream file_str;
  // Build complete filename
  file_str << _filename << "."<< _time_step << ".dat";
  std::ofstream file_h;
  file_h.open((file_str.str()).c_str());
   
  
  p = wp.getParticle();
  w = wp.getWeight();
  long long row = p->size1();
  int col = p->size2();


  // Use p.pos if it has been set
  if(p->pos() > 0 ) 
    row = p->pos();

  // print header
  file_h << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  file_h << "% Particle data file for iteration: " << _time_step << "\n";
  file_h << "% Data format: x_1 x_2 ... x_d weight \n";
  file_h << "% Number of particles: " << row << "\n";
  file_h << "% Particle dimension: " << col << "\n";
  file_h << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";



  for(long long k=0; k<row; k++) {    
    // Store position  
    for(int l=0; l<col; l++) {
      file_h << (*p)(k,l) << " ";
    }
    // Store weight
    file_h << (*w)(k,0) << "\n"; 
  }

  file_h.close();
  _time_step++;
  
}







void DatIO::operator<<(const Particle &p) {

  
  std::ostringstream file_str;
  // Build complete filename
  file_str << _filename << "."<< _time_step << ".dat";
  std::ofstream file_h;
  file_h.open((file_str.str()).c_str());
   
  
  long long row = p.size1();
  int col = p.size2();


  // Use p.pos if it has been set
  if(p.pos() > 0 ) 
    row = p.pos();

  // print header
  file_h << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  file_h << "% Particle data file for iteration: " << _time_step << "\n";
  file_h << "% Data format: x_1 x_2 ... x_d  \n";
  file_h << "% Number of particles: " << row << "\n";
  file_h << "% Particle dimension: " << col << "\n";
  file_h << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";

  for(long long k=0; k<row; k++) {    
    // Store position  
    for(int l=0; l<col; l++) {
      file_h << p(k,l) << " ";
    }
    file_h << "\n"; 
  }

  file_h.close();
  _time_step++;
  
}













void DatIO::setIterationStep(int t) {
  _time_step = t;
}







