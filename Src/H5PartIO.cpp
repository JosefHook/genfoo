#include <iostream>
#include "H5PartIO.hpp"
#include "mArray.hpp"
//#include <H5Part.h>
#include <boost/lexical_cast.hpp>



H5PartIO::H5PartIO(std::string filename) {


  /* Open an HDF5 file for writing */
  _filename = filename;
#ifdef SERIAL
  _file_h = H5PartOpenFile(_filename.c_str(),H5PART_WRITE);
#else
  _file_h = H5PartOpenFileParallel(_filename.c_str(),H5PART_WRITE, 
				   MPI_COMM_WORLD);
#endif
  //  if(H5PartFileIsValid(_file_h)) {
  //  std::cout << "Could not open file for writing" << std::endl;
  // }
  _time_steps=0;
  _np=0;


}


H5PartIO::H5PartIO(std::string filename, long long n) {


  /* Open an HDF5 file for writing */
  _filename = filename;
#ifdef SERIAL
  _file_h = H5PartOpenFile(_filename.c_str(),H5PART_WRITE);
#else
  _file_h = H5PartOpenFileParallel(_filename.c_str(),H5PART_WRITE, 
				   MPI_COMM_WORLD);
#endif
  _time_steps=0;
  if(n<=0) 
    _np = 1;
  else
    _np = n;
  H5PartSetNumParticles(_file_h,_np);
  
}



H5PartIO::H5PartIO(std::string filename, long long n, int D) {


  /* Open an HDF5 file for writing */
  _filename = filename;
  _file_h = H5PartOpenFile(_filename.c_str(),H5PART_WRITE);
  _time_steps=0;
  if(n<=0) 
    _np = 1;
  else
    _np = n;
  H5PartSetNumParticlesStrided(_file_h,_np,D);
  
}






int H5PartIO::file(std::string filename) {


  /* Open an HDF5 file for writing */
  _filename = filename;
  _file_h = H5PartOpenFile(_filename.c_str(),H5PART_WRITE);
  // if(H5PartFileIsValid(_file_h)) {
  //  std::cout << "Could not open file for writing" << std::endl;
  //  return 1;
  //}
  // Set default values
  _time_steps=0;
  _np=0;
  
  return 0;

}


void H5PartIO::setTimeSteps(int n) {
  if(n<=0) 
    _time_steps = 0;
  else
    _time_steps = n;
  H5PartSetStep(_file_h,_time_steps);
}
  

void H5PartIO::setNumParticles(long long n, int D) {

 if(n<=0) 
    _np = 1;
  else
    _np = n;
 H5PartSetNumParticlesStrided(_file_h,_np, D);

}  



void H5PartIO::setNumParticles(long long n) {

 if(n<=0) 
    _np = 1;
  else
    _np = n;
 H5PartSetNumParticles(_file_h,_np);

}  



  
// Assume for now that p(3)
void H5PartIO::operator<<(const Particle &p) {




  const double *ptr;
  ptr = &(p.data()[0]); // acess the data
  long long row = p.size1();
  int col = p.size2();
  //  std::string space[3] = {"x", "y", "z"};
  const char *space[] = {"x", "y", "z"};

  // Use p.pos if it has been set
  if(p.pos() > 0 ) 
    row = p.pos();



  // The Particle type is ordered in row major
  // m x n  dimensional matrix M_ij  => (i x n + j)
  // i=0, j=0 => x_0 , i=0, j=1 => y_0 ...
  // n = columns
  //  int col = p.size2();
  // We want pointers to the first row for each column


  // Write out data independent of dimensions
  // Note H5Part IO library needs keywords "x", "y" and "z"
  // which is very stupid indeed!

  // Ok so this code is TO general for H5Part (leave it for now)!!
  //  std::string base;
  //  for(int l=0; l<p.size2(); l++) {
  //  if(l>=3)  
  //    base = space[l%3] + boost::lexical_cast<std::string>(l);
  //    else
  //    base = space[l];
  //    H5PartWriteDataFloat64(_file_h,base.c_str(),&ptr[l*col]);
  //}
  // H5Part striding does not seem to work.
  // We therefore need to split each particle to individual arrays
  //

  double *data = new double[(int)row];

  H5PartSetStep(_file_h,_time_steps);  
  // Paraview requires a 3D array no matter what!
  for(int l=0; l<3; l++) {
    
    if(l<col) {

      for(int k=0; k<row; k++)
	data[k] = p(k,l);

          //     &ptr[l*col]

      //      H5PartWriteDataFloat64(_file_h,space[l].c_str(),data);
      H5PartWriteDataFloat64(_file_h,space[l],data);
      
      // Fill out the rest of the space with zeros !
    }
    else {
      double *zero = new double[p.size1()];
      for(int m=0; m< p.size1(); m++)
  	zero[m] = 0.0;
      
      //H5PartWriteDataFloat64(_file_h,space[l].c_str(),zero);
      H5PartWriteDataFloat64(_file_h,space[l],zero);
      delete zero;
    }
    
   }
  delete data;
  

  

  // Automatic increment time step
  _time_steps++;

}







// Assume for now that p(3)
void H5PartIO::operator<<(const WeightedParticle &wp) {

  H5PartSetStep(_file_h,_time_steps);



  const char *pspace[] = {"x", "y", "z"};
  const char *wspace[] = {"w"};
  Particle *p;
  mArray *w;
  
  p = wp.getParticle();
  //w = wp.getWeight();
  w = wp.getWeight();
  long long row = p->size1();
  int col = p->size2();


  // Use p.pos if it has been set
  if(p->pos() > 0 ) 
    row = p->pos();

  
  // The Particle type is ordered in row major
  // m x n  dimensional matrix M_ij  => (i x n + j)
  // i=0, j=0 => x_0 , i=0, j=1 => y_0 ...
  // n = columns
  //  int col = p.size2();
  // We want pointers to the first row for each column


  // Write out data independent of dimensions
  // Note H5Part IO library needs keywords "x", "y" and "z"
  // which is very stupid indeed!

  // Ok so this code is TO general for H5Part (leave it for now)!!
  //  std::string base;
  //  for(int l=0; l<p.size2(); l++) {
  //  if(l>=3)  
  //    base = space[l%3] + boost::lexical_cast<std::string>(l);
  //    else
  //    base = space[l];
  //    H5PartWriteDataFloat64(_file_h,base.c_str(),&ptr[l*col]);
  //}
  // H5Part striding does not seem to work.
  // We therefore need to split each particle to individual arrays
  //

  double *pdata = new double[row];
  double *wdata = new double[row];
 
  
  // Paraview requires a 3D array no matter what!
  // Loop over space and weights assuming 
  // WeightedParticle having the ordering: x,y,z, wx, wy,wz

  // Loop over particles
  for(int l=0; l<3; l++) {
    
    if(l<col) {

      for(int k=0; k<row; k++) {
	pdata[k] = (*p)(k,l);
	if(l==0)
	  wdata[k] = (*w)(k,l);
      }

      if(l==0)
	H5PartWriteDataFloat64(_file_h,wspace[l],wdata);
      
      H5PartWriteDataFloat64(_file_h,pspace[l],pdata);
    }      // Fill out the rest of the space with zeros !
    else {
      double *zero = new double[p->size1()];
      for(int m=0; m< p->size1(); m++)
  	zero[m] = 0.0;
      H5PartWriteDataFloat64(_file_h,pspace[l],zero);
      delete zero;
    }
    
  }





  delete pdata, wdata;
  

  

  // Automatic increment time step
  _time_steps++;

}


int H5PartIO::close(void) {
  
  H5PartCloseFile(_file_h);
  return 0;
}

