//
// Testing of the ANN library
// 
#include <iostream>
#include "BallResample.hpp"
#include "Particle.hpp"

#include <cmath>

#include <boost/numeric/bindings/traits/std_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/lapack/lapack.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/ublas/io.hpp> 
#include <boost/random.hpp>

boost::mt19937 u_seed;
boost::uniform_01<boost::mt19937> U(u_seed);


using namespace std;

//
// Output newp , New points
// Output neww , New weights for above points
// Input oldp,   Old points
// Input oldw,   Old weights
// Input nnIdx,  ANN KD-tree index array. 
// Input radius, Ball-radius
// Input qp, Query-point
// Input qw, Query-point weight
// Input n_ball, N- partilces in ball = length of nnIdx
// Input n_reduced, reduced N
//
// Observe allocations must be done outside function
//
// Only works for 2D cases for now.
//
int BallResample::calcParticleWeights(ANNpointArray newp,
				      mArray *neww,
				      ANNpointArray oldp,
				      mArray *oldw, 
				      ANNidxArray nnIdx,
				      double radius,
				      ANNpoint qp,
				      double qw,
				      int n_ball, 
				      int n_reduced,
				      long gIdx,
				      int dim) {
  
  namespace lapack = boost::numeric::bindings::lapack;
  namespace ublas =  boost::numeric::ublas;
  
  ublas::matrix<double, ublas::column_major>  *G, *K,*GK, *ID;
  ublas::vector<double> *lweight;
  // Calculate new weights from 
  // G*a = K*b where G is a Gram matrix 
  // for new particles with gaussian distributed positions
  // and K is the old Gram matrix.
  // a re the new weights
  // b is the old weight
  
  if(n_ball < n_reduced)  {
    
    cout << " Not enough  particles to reduce" << endl;
    return 1;
  }


  ANNpoint pt, ptn;
  int ncol = std::max(n_ball, n_reduced);
  int nrow;

  cout << " In educe, dim: " << dim << endl;
  // Allocation of memory follows the formula below, 
  //
  // 1D case is special
  // 1D: 0  + 2*1 + 1 = 3
  // 2D: 1  + 2*2 + 1 = 6
  // 3D: 3*1 + 2*3 + 1 =  10 
  // 4D: 3*2 + 2*4 + 1 = 15  
  // 5D: 10 + 2*5 +1 = 21
  // General case: 
  // Sd = (dim-1)/2 ( 2*(dim-1) -(dim-2) )  => (dim-1)*(dim-1) - 1/2 (dim-1)*(dim-2)  
  //
  // Total size of matrix is : Sd + 2*dim +1 
  // 2D: 1*1 -1/2*1*0 = 1 + 2*dim +1 = 6
  // 3D: 2*2 -1/2*2*1 = 3  + 2*dim +1 = 10
  // 4D: 3*3 -3*2/2 =   6  + 2*dim +1 = 15   
  // 5D: 4*4 -1/2 * 4*3 = 10 + 2*dim +1 = 21 
  //
  if(dim==1) 
    nrow = 3;
  else 
    nrow = (dim-1)*(dim-1) - 1/2*(dim-1)*(dim-2) + 2*dim +1; 

  cout << "In reduce, nrow: " << nrow << endl; 
  G = new ublas::matrix<double, ublas::column_major> (nrow,n_reduced);
  K = new ublas::matrix<double, ublas::column_major> (nrow,n_ball);
  GK = new ublas::matrix<double, ublas::column_major> (n_reduced,n_ball);
  lweight = new ublas::vector<double>(n_ball);
    

  for(int j=0; j<n_ball; j++) {
    if(j<n_ball) {

	// Pointer to particle j
	pt = oldp[nnIdx[j]];

	(*K)(0,j) = 1.0; 	
	
	// Position
	for(int i=0; i<dim; i++) {
	  (*K)(i+1,j) = pt[i] - qp[i];      // x
#ifdef DBUG_INDX	      
	  cout << "K(" << i+1 << ", " << j << ") = pt[" << i <<"]"<< endl; 
#endif
	  // Break exit function if not a number
	  if(std::isnan(pt[i]) || std::isnan(qp[i]))
	    return 1;
	}
	
	// Second moment
	for(int i=dim; i<2*dim; i++) {
	  (*K)(i+1,j) = ((*K)(i-dim +1 ,j))*((*K)(i-dim+1,j));  // x*x
#ifdef DBUG_INDX	      
	  cout << "K(" << i+1 << ", " << j << ") = K(" 
	       << i-dim+1 << ", " << j << ")*K("
	       << i-dim+1 << ", " << j << ")" << endl; 
#endif
	}

	if(dim>1) {	
	  // Covariance
	  for(int i=1; i<dim; i++) {
	    for(int k=i; k<dim; k++) {
	      (*K)(2*dim+k+i-1,j) = ((*K)(i,j))*((*K)(k+1,j));  

#ifdef DBUG_INDX	      
	      cout << "K(" << 2*dim +k + i-1 << ", " << j << ") =K(" 
		   << i << ", " << j << ")*K("
		   << k+1 << ", " << j << ")" << endl; 
#endif


	    }
	  }
	}




	// 1D
	//
	//  (*K)(1,j) = pt[0] - qp[0];      // x
	//  0 
	//  (*K)(2,j) = ((*K)(1,j))*((*K)(1,j));  // x*x
	//
	// 2D 
	//  (*K)(1,j) = pt[0] - qp[0];      // x
	//  (*K)(2,j) = pt[1] - qp[1];      // y
	//
	//
	//  (*K)(3,j) = ((*K)(1,j))*((*K)(1,j));  // x*x
	//  (*K)(4,j) = ((*K)(2,j))*((*K)(2,j));  // y*y
	//	
	//  (*K)(5,j) = ((*K)(1,j))*((*K)(2,j));  // x*y 
	//
	//
	// 3D
	//
	//  (*K)(1,j) = pt[0] - qp[0];      // x
	//  (*K)(2,j) = pt[1] - qp[1];      // y
	//  (*K)(3,j) = pt[1] - qp[1];      // z
	//
	//  (*K)(4,j) = ((*K)(1,j))*((*K)(1,j));  // x*x
	//  (*K)(5,j) = ((*K)(2,j))*((*K)(2,j));  // y*y
	//  (*K)(6,j) = ((*K)(3,j))*((*K)(3,j));  // z*z
	//	
	//  (*K)(7,j) = ((*K)(1,j))*((*K)(2,j));  // x*y 
	//  (*K)(8,j) = ((*K)(1,j))*((*K)(3,j));  // x*z 
	//
	//  (*K)(9,j) = ((*K)(2,j))*((*K)(3,j));  // y*z 
	//
	//
	// 4D
	//
	//  (*K)(1,j) = pt[0] - qp[0];      // x
	//  (*K)(2,j) = pt[1] - qp[1];      // y
	//  (*K)(3,j) = pt[1] - qp[1];      // z
	//  (*K)(4,j) = pt[1] - qp[1];      // q
	//
	//  (*K)(5,j) = ((*K)(1,j))*((*K)(1,j));  // x*x
	//  (*K)(6,j) = ((*K)(2,j))*((*K)(2,j));  // y*y
	//  (*K)(7,j) = ((*K)(3,j))*((*K)(3,j));  // z*z
	//  (*K)(8,j) = ((*K)(4,j))*((*K)(4,j));  // q*q
	//
	//  (*K)(9,j) = ((*K)(1,j))*((*K)(2,j));  // x*y 
	//  (*K)(10,j) = ((*K)(1,j))*((*K)(3,j));  // x*z 
	//  (*K)(11,j) = ((*K)(1,j))*((*K)(4,j));  // x*q 
	//
	//  (*K)(12,j) = ((*K)(2,j))*((*K)(3,j));  // y*z 
	//  (*K)(13,j) = ((*K)(2,j))*((*K)(4,j));  // y*q 
	//
	//  (*K)(14,j) = ((*K)(3,j))*((*K)(4,j)); // z*q 
	//	
	// 5D
	//
	//  (*K)(1,j) = pt[0] - qp[0];      // x
	//  (*K)(2,j) = pt[1] - qp[1];      // y
	//  (*K)(3,j) = pt[1] - qp[1];      // z
	//  (*K)(4,j) = pt[1] - qp[1];      // q
	//  (*K)(5,j) = pt[1] - qp[1];      // p
	//
	//  (*K)(6,j) = ((*K)(1,j))*((*K)(1,j));  // x*x
	//  (*K)(7,j) = ((*K)(2,j))*((*K)(2,j));  // y*y
	//  (*K)(8,j) = ((*K)(3,j))*((*K)(3,j));  // z*z
	//  (*K)(9,j) = ((*K)(4,j))*((*K)(4,j));  // q*q
	//  (*K)(10,j) = ((*K)(4,j))*((*K)(4,j));  // p*p
	//
	//  (*K)(11,j) = ((*K)(1,j))*((*K)(2,j));  // x*y 
	//  (*K)(12,j) = ((*K)(1,j))*((*K)(3,j));  // x*z 
	//  (*K)(13,j) = ((*K)(1,j))*((*K)(4,j));  // x*q 
	//  (*K)(14,j) = ((*K)(1,j))*((*K)(5,j));  // x*p 
	//
	//  (*K)(15,j) = ((*K)(2,j))*((*K)(3,j));  // y*z 
	//  (*K)(16,j) = ((*K)(2,j))*((*K)(4,j));  // y*q 
	//  (*K)(17,j) = ((*K)(2,j))*((*K)(5,j));  // y*p 
	//
	//  (*K)(18,j) = ((*K)(3,j))*((*K)(4,j)); // z*q 
	//  (*K)(19,j) = ((*K)(3,j))*((*K)(4,j)); // z*p 
	// 
	//  (*K)(20,j) = ((*K)(4,j))*((*K)(5,j)); // q*p 
	//

	
	// Populate weight vector 
	(*lweight)(j) = (*oldw)(j,0);

    }
    
    
    // Populate the LHS
    if(j<n_reduced) {



      // Store new points in newp, 
      ptn = newp[gIdx + j];
      (*G)(0,j) = 1.0; 	

#ifdef DBUG_INDX	      
      cout << " ----------- GINDEX-------------------------" << endl; 
#endif            
      // Position
      for(int i=0; i<dim; i++) {
	(*G)(i+1,j) = (*K)(i+1,j)  + 1e-1*U();      // pos
	ptn[i] = (*G)(i+1,j)+ qp[i];          // Store new pos
#ifdef DBUG_INDX	      
	cout << "G(" << i+1 << ", " << j << ") =K(" << i+1 << ", " << j << ")"<< endl; 
#endif
      }
      
      // Second moment
      for(int i=dim; i<2*dim; i++) {
	(*G)(i+1,j) = ((*G)(i-dim +1 ,j))*((*G)(i-dim+1,j));  // pos*pos
#ifdef DBUG_INDX	      
	cout << "G(" << i+1 << ", " << j << ") =G(" 
	     << i-dim+1 << ", " << j << ")* G("
	     << i-dim+1 << ", " << j << ")" << endl; 
#endif
      }
      
      if(dim>1) {	
	// Covariance
	for(int i=1; i<dim; i++) {
	  for(int k=i; k<dim; k++) {
	    (*G)(2*dim+k+i-1,j) = ((*G)(i,j))*((*G)(k+1,j));  
#ifdef DBUG_INDX	      
	    cout << "G(" << 2*dim +k + i-1 << ", " << j << ") =G(" 
		 << i << ", " << j << ")* G("
		 << k+1 << ", " << j << ")" << endl; 
#endif


	  }
	}
      }

#ifdef DBUG_INDX	      
      cout <<"------------------------ NEW G ---------------- " << endl;
      cout << "G: " << *G << endl;
#endif
      //// Randomly sample a new position ( OBS Uniform distribution ) 
      
      //      (*G)(0,j)  = 1.0;     
      //// draw new particles around the query point
      ////double r, phi;
      ////      r = U();
      ////phi = U();
      ////      (*G)(1,j) = radius*r*cos( 2.0*M_PI*phi ) + qp[0] ;// //1.0/std::sqrt(2.0)*radius*U(); //x
      ////  (*G)(2,j) = radius*r*sin( 2.0*M_PI*phi )  + qp[1]; // y
      // (*G)(1,j) = (*K)(1,j)  + 1e-1*U();
      // (*G)(2,j) = (*K)(2,j)  + 1e-1*U();
      
      //      (*G)(3,j) = ((*G)(1,j))*((*G)(1,j));  // x*x
      // (*G)(4,j) = ((*G)(2,j))*((*G)(2,j));  // y*y
      // (*G)(5,j) = ((*G)(1,j))*((*G)(2,j));  // x*y
 
      //// Store new points in newp, 
      //      ptn = newp[gIdx + j];
      // ptn[0] = (*G)(1,j)+ qp[0];
      // ptn[1] = (*G)(2,j) + qp[1];

#ifdef DBUG_INDX	      
      cout <<"------------------------ OLD G ---------------- " << endl;
      cout << "G: " << *G << endl;
#endif
 
    }
    
  }

#ifdef DBUG_INDX	      
    cout << "n_ball: "<< n_ball << endl; 
    cout << " K: " << *K << endl;

    cout << "lweight: " << *lweight << endl;
    cout << "n_recuded: "<< n_reduced << endl; 
    cout << " G: " << *G << endl;
#endif
  
  lapack::optimal_workspace workspace;

  int v1; 
  int max_s;
  v1 = std::max(nrow, n_reduced);
  max_s = std::max(1, v1); 
  ID = new ublas::matrix<double, ublas::column_major> (ublas::identity_matrix<double>(max_s, max_s));
  ublas::vector<double> *alpha = new ublas::vector<double>(max_s);
  lapack::gelss(*G,*ID,workspace);


  // Result is written in K and  have dimension K(M x P) given
  // G (N x M) and K(N x P).
  // Pseudo inverse
  ID->resize(n_reduced,nrow);
  (*GK) = prod(*ID,*K);      	
  // Multiply with the weight
  (*alpha) = prod(*GK,*lweight);      

  double maxa=0;
  double mina=0;
  // Todo: fix for speed
  for(int o=0; o<n_reduced; o++) {
    std::cout << " alpha(o) " << (*alpha)(o) << std::endl;
    maxa=std::max(maxa, (*alpha)(o));
    mina=std::min(mina, (*alpha)(o));

    (*neww)(o+gIdx,0) = (*alpha)(o);
  }
  //
  // Safeguard against large weights
  // return error if to large or small
  // 
  if(maxa>100 || mina<-100) 
    return 1;


  
  delete lweight, alpha;
  delete G, K, GK, ID;
  return 0;
}
	




//
// Output: New weightedParticle
// Input: Weighted Particles array 
//
//

int BallResample::reduce(WeightedParticle *owp, 
			 double radius, int maxneighbors, double reducefactor )
{
  
  
  // General variables
  mArray *newweights; // New weight array 
  mArray *oldweights = owp->getWeight(); // Old weight array 
  mArray *tmpweights;
  mArray *wp;
  Particle *op = owp->getParticle();
  Particle *np;  // New particle weights
  
  int retval = 0;

  // GENERAL PARAMETERS
  double qw;                     // Query weight
  int dim;                // NR- dimensions
  long TotSize;
  long Size;             // Max particles in inital array.
  long ReducedSize;        // Reduced particle size  
  long gIdx  = 0;             // Global index

    
  // ANN PARAMETERS
  double err = 0;              // KD-tree error
  
  // ANN  variables
  ANNpointArray dataPts, tmpdataPts;
  ANNpointArray newdataPts;
  ANNpoint queryPt, pt, tpt;
  ANNidxArray nnIdx;
  ANNdistArray dists;
  ANNkd_tree* kdTree;
  ANNdist eps =radius;
  int nRet;  
  int newnRet;
  long QueryPos;
  double remainder = 0; // Keep  track of remainder  

  bool noluck = false;

  dim =  op->size2(); 
  TotSize = op->size1();
  Size = op->pos();  //op->size1();   
  ReducedSize = round(Size*reducefactor);
  
  nnIdx = new ANNidx[maxneighbors];
  dists = new ANNdist[maxneighbors];

  // Allocate memory for local particle array 
  queryPt = annAllocPt(dim);
  dataPts = annAllocPts(Size, dim);
  
  
  // Copy data points 
  for(int i=0; i<Size; i++) {
    pt = dataPts[i]; 
    for(int j=0; j<dim; j++) 
      pt[j] = (*op)(i,j);
  }


  // Allocate memory for weight arrays
  newweights = new mArray(ReducedSize,1);
  newweights->clear();

  // Allocate memory for new data point array 
  newdataPts = annAllocPts( ReducedSize , dim);

  
  //
  // Loop over particles 
  //
  //
  while( gIdx< ReducedSize) { 
	 
    
    // Randomly pick a  query - point 
    QueryPos = floor((Size-1)*U());
    queryPt = dataPts[QueryPos]; //(*op)( QueryPos, j); 
    // Build kd-tree and query.
    kdTree = new ANNkd_tree(dataPts, Size, dim);
    nRet= kdTree->annkFRSearch(queryPt, std::pow(eps,2.0), maxneighbors, nnIdx, dists, err);
 

    nRet = std::min(nRet, maxneighbors);
    
    // Define new number of particles
    newnRet =  round( (double) nRet*(double)reducefactor  + (double)remainder);
    remainder = abs( newnRet  - (double) nRet* (double)reducefactor - (double) remainder  );
    cout << "New number of particles: " << newnRet << endl; 
    cout << " Remainder : " << remainder << endl;
    
      

    if(nRet>1 && (newnRet + gIdx < ReducedSize)) {
      
      // Print diagnostics
      cout << "Query-point: ( " << queryPt[0];
      for(int j=1; j<dim; j++)
	cout << ", " << queryPt[j];
      cout << " ) "<< endl;
      cout << " Max distance: " << eps << endl;
      cout << "N points found: " << nRet << endl;
      cout << "NN: i\tIdx[i]\tdistance  \t point" << endl;
      
      //
      // Loop over result
      // 
      for(int i=0; i<nRet; i++) {
	dists[i] = sqrt(dists[i]);
	pt = dataPts[nnIdx[i]];
	cout << i << "\t" << nnIdx[i] << "\t" << dists[i];
	cout << "  \t ( "<< pt[0] ;
	for(int j=1; j<dim; j++) {
	  cout << ", " << pt[j] ;
	}
	cout << " ) " << endl; 

	
      }
      
      
      
      if(newnRet == 0) {
	cout << "To few particles found. " << endl;
	break;
      }
 
      /////////////////////////////////// DIAGNOSTICS //////////////////////////////////
      //
      // Loop over result
      //
      cout << "--------- BEFORE CALL TO REDUCE --------------------" << endl;
      for(int i=0; i<nRet; i++) {
	pt = dataPts[nnIdx[i]];
	cout << "i: " <<i << " ( "<< pt[0] ;
	for(int j=1; j<dim; j++) {
	  cout << ", " << pt[j] ;
	}
	cout << " ) old weights: " <<  (*oldweights)(i,0) << endl; 
	
      }
      cout << "--------- BEFORE CALL TO REDUCE --------------------" << endl;
      /////////////////////////////////// DIAGNOSTICS //////////////////////////////////

      // Calc weights
      retval = calcParticleWeights(newdataPts, newweights, 
				   dataPts, oldweights, nnIdx, 
				   radius,  queryPt , qw,
				   nRet , newnRet, gIdx, dim  );
      



      /////////////////////////////////// DIAGNOSTICS //////////////////////////////////
      //
      // Loop over result
      //
      cout << "--------- AFTER CALL TO REDUCE --------------------" << endl;
      for(int i=0; i<newnRet; i++) {
	pt = newdataPts[gIdx + i];
	cout << "i: " <<i << " ( "<< pt[0] ;
	for(int j=1; j<dim; j++) {
	  cout << ", " << pt[j] ;
	}
	cout << " ) new weights: " <<  (*newweights)(gIdx + i,0) << endl; 
	
      }
      cout << "--------- AFTER CALL TO REDUCE --------------------" << endl;
      /////////////////////////////////// DIAGNOSTICS //////////////////////////////////

      
      if(retval == 0 ) { 
      	//
	// Mark and  indexed particles in dataPts
	// 
	for(int i=0; i<nRet; i++) 
	  dataPts[nnIdx[i]] = NULL;	
	// Mark  Query position
	dataPts[QueryPos] = NULL;

	tmpdataPts = annAllocPts( Size-nRet, dim);	
	tmpweights = new mArray(Size-nRet, 1);
	// Reduce dataPts and oldweights
	int ii = 0;
	for(int i=0; i<Size; i++) {
	  pt = dataPts[i];
	  // If unmarked particle, copy values
	  if(pt != NULL) {
	    tpt = tmpdataPts[ii];
	    for(int j = 0; j<dim ; j++) {
	      tpt[j] = pt[j];
	    }
	    (*tmpweights)(ii,0) = (*oldweights)(i,0);
	    ii+=1;
	  } 
	}
	
	
	// Replace dataPts with reduced, tmpdatapts 
	annDeallocPts(dataPts);
	dataPts = tmpdataPts;
	
	// Replace oldweights with reduced, tmpweights
	delete oldweights;
	oldweights = tmpweights;
#ifdef DBUG3
	cout << " DataPts reduced size: " << ii << "  from Size: " << Size << endl;
	cout << " newDataPts increased size: " << gIdx+newnRet << "  from Size: " << gIdx << endl;
	cout << " newnRet : " << newnRet << endl;
	cout << " Size-newnRet: "<< Size - nRet  << endl;
	cout << "nRet " << nRet << endl;
#endif	
    
    
	// Update global index
	gIdx += newnRet; 
      Size -= nRet;
      } else {
	noluck = true;
      }
      //
      // If no particles found with KD-tree , copy over Query-point and continue!
      // 
#ifdef DBUG2     
      cout << "----------------- END IF CASE -------------------" << endl;
      cout << "-- nRet: " << nRet << " (newnRet + gIdx < ReducedSize): " << 
	(newnRet + gIdx < ReducedSize)<< "-----"<<endl;
#endif
      
    } else { 
      noluck = true;
    }
    
    
    
    if(noluck) {

#ifdef DBUG1     
      cout << "-------------------------------------------------------------" << endl;
      cout << "-------------------------------------------------------------" << endl;

      cout << "Did not find any neighbors. ";
      cout << "Trying with a new point " << endl;
      cout << "Global index in reduced data: " << gIdx << endl;
      cout << "reduced data size: " << ReducedSize << endl;
      cout << "Query Pos: " << QueryPos<< endl;
      pt = dataPts[QueryPos];
      cout << " QoeryPoint (" << pt[0] << ", " << pt[1] << " )"<< endl;
      cout << "data size: " << Size << endl;
      cout << " dim : " << dim << endl;
      cout << " -------- DATA -------------" << endl;
      for(int i=0; i<Size; i++) {
	 pt = dataPts[i];
	cout << "  \t ( "<< pt[0] ;
	for(int j=1; j<dim; j++) {
	  cout << ", " << pt[j] ;
	}
	cout << " ) " << endl;
      }
#endif
      







      // Store Query Position and weight
      pt =  dataPts[QueryPos];
      tpt = newdataPts[gIdx];
      for(int j=0; j<dim; j++){
	tpt[j] = pt[j];
      }
      (*newweights)(gIdx,0) = (*oldweights)(QueryPos,0);
   
      gIdx++;


#ifdef DBUG1     
      //////////////////////////////// DIAGNOSTICS /////////////////////////////

    cout << " -------- NEW DATA -------------" << endl;
      for(int i=0; i<ReducedSize; i++) {
	 pt = newdataPts[i];
	cout << "  \t ( "<< pt[0] ;
	for(int j=1; j<dim; j++) {
	  cout << ", " << pt[j] ;
	}
	cout << " ) " << endl;
      }



      cout << " -------- QUERY RESULT  -------------" << endl;

      // Print diagnostics
      cout << "Query-point: ( " << queryPt[0];
      for(int j=1; j<dim; j++)
	cout << ", " << queryPt[j];
      cout << " ) "<< endl;
      cout << " Max distance: " << eps << endl;
      cout << "N points found: " << nRet << endl;
      cout << "NN: i\tIdx[i]\tdistance  \t point" << endl;
      
      //
      // Loop over result
      // 
      for(int i=0; i<nRet; i++) {
	dists[i] = sqrt(dists[i]);
	pt = dataPts[nnIdx[i]];
	cout << i << "\t" << nnIdx[i] << "\t" << dists[i];
	cout << "  \t ( "<< pt[0] ;
	for(int j=1; j<dim; j++) {
	  cout << ", " << pt[j] ;
	}
	cout << " ) " << endl; 
	
      }
      //////////////////////////////// DIAGNOSTICS /////////////////////////////

#endif


      // Mark  Query position
      dataPts[QueryPos] = NULL;
      tmpdataPts = annAllocPts( Size-1, dim);
      tmpweights = new mArray(Size-1, 1);
	
      // Reduce dataPts
      int ii = 0;
      for(int i=0; i<Size; i++) {
	pt = dataPts[i];
	// If unmarked particle, copy values
	if(pt != NULL) {
	  tpt = tmpdataPts[ii];
	  for(int j = 0; j<dim ; j++) {
	    tpt[j] = pt[j];
	  }
	  (*tmpweights)(ii,0) = (*oldweights)(i,0);
	  ii+=1;
	} 
      }
      
      // Replace dataPts with reduced, tmpdatapts 
      annDeallocPts(dataPts);
      dataPts = tmpdataPts;
      
      // Replace oldweights with reduced, tmpweights
      delete oldweights;
      oldweights = tmpweights;
      
      Size--;


#ifdef DBUG2     
      cout << "----------------- END ELSE CASE -------------------" << endl;
      cout << "-- nRet: " << nRet << " (newnRet + gIdx < ReducedSize): " << 
	(newnRet + gIdx < ReducedSize)<< "-----"<<endl;
#endif
      noluck = false;
    }
    

   	// Remove Tree and start over
	delete kdTree;
	annClose();
 
  }


  // There might still be some particles left which needs to be added.
  long finalSize = ReducedSize + Size;

  // Copy newdatPts to new WeightedParticle array.
  //  Allocate memory for new particles
  np = new Particle(TotSize,dim);
  np->set_pos(finalSize);
  wp = new mArray(TotSize, 1);
  wp->clear();

  for(long i = 0; i<ReducedSize; i++) {
    pt = newdataPts[i];
    for(int j = 0; j<dim ; j++)
      (*np)(i,j) = pt[j];
    (*wp)(i,0) = (*newweights)(i,0);
  }


  for(int i=0; i<Size; i++) {
    pt = dataPts[i];
    for(int j = 0; j<dim ; j++)
      (*np)(ReducedSize + i,j) = pt[j];
  (*wp)(ReducedSize + i,0) = (*oldweights)(i,0);
  }
  
  delete newweights;
  delete oldweights;
  delete op; 

  owp->setParticle(np);
  owp->setWeight(wp);

  cout << "gIdx : " << gIdx << endl;
  cout << "finalSize : " << finalSize << endl;

  // Free up some memory
  delete newdataPts;  
  delete nnIdx, dists;
  
  return 0;
}












