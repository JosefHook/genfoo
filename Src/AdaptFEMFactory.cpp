//
// AdaptFEMFactory by Josef Höök, 2012 
// Copyright (c) all rights reserved
//
#include <dolfin.h>
//#include <glib.h>  
#include <boost/algorithm/string/trim.hpp>
#include "FokkerPlanck1D.h"
#include "FokkerPlanck2D.h"
#include "FokkerPlanck3D.h"
#include "DualFokkerPlanck2D.h"
#include "Kernel2D.h"



#include <dolfin/adaptivity/ErrorControl.h>

//#include "include/Coefficients.hpp"
#include "DynamicLoader.hpp"
#include "Params.hpp"
#include "Operator.hpp"
#include "OperatorExpression.hpp"
using namespace dolfin;


void ShowAFactoryInfo()
{


  cout << "     Using Adaptive Finite Element backend version 1.1            " << endl;


}

////////////////////////////////////////////////
// MAIN
////////////////////////////////////////////////

int AdaptFEMFactory(Params pm)
{

  
  // General info
  ShowAFactoryInfo();

  DynamicLoader dl;
  
  
  if(  dl.Load( ((std::string) pm.getLibraryFile()).c_str() )  )
    exit(1);
  
  Operator *ops; 
  Operator *(*load)() = (Operator*(*)())(dl.bootstrap());
  // Create new instance of coefficient object


  
// Call operator constructor with XML file 
  ops = (*load)();  

 // Pass XML doc to operator
  ops->registerXMLDoc(pm.getXMLDoc());



  cout << "Info: " << ops->getInfo() << endl;
  cout << "Problem defined in R^" << ops->getDim() << endl;

  Expression* ui;
  Expression* drift;
  Expression* diffusion;

  OperatorExpression oe(ops);
  
  ui        = oe.getIC();
  drift     = oe.getDrift();
  diffusion = oe.getDiffusion();


  cout << " Value rank of ui: " << ui->value_rank() << endl;
  cout << " Value rank of drift: " << drift->value_rank() << endl;
  cout << " Value rank of diffusion: " << diffusion->value_rank() << endl;
 


  cout << "Drift rank: " << drift->value_rank() << endl;
  cout << "Diffusion rank: " << diffusion->value_rank() << endl;



  if((ops->getDim()<1) || (ops->getDim()>3) ) {
    cout << "Wrong dimension in coefficient library. Received: " 
	 << ops->getDim()  << " should be [1,2,3]." << endl; 
    exit(1);
  }
    
  
  // Test if library file has same dimension as config file
  if(ops->getDim() != pm.getDimension() ) {
    cout << "Space dimension differs between library- and configfile." << endl; 
    cout << "In "<< pm.getLibraryFile() << " space dimension is:" 
	 << ops->getDim() << endl;
    cout << "In " << pm.getParameterFile() <<" space dimension is:" 
	 << pm.getDimension() << endl;
    exit(1);
  }
  unsigned int dim;
  void *a, *L, *R;
  void *dual_a, *dual_L;
 

  FunctionSpace *V, *dual_V;  
  Mesh *mesh;

  // Domain information
  std::vector<int> grid = pm.getGridDimension();
  std::vector<double> bmax = pm.getBoundaryMaxValues();
  std::vector<double> bmin = pm.getBoundaryMinValues();
      

 
  switch(ops->getDim())
    {
      
    case 1: {
      cout << "One dimensional problem detected" << endl;

      
      mesh = new Interval(grid[0], bmin[0], bmax[0]);
      //mesh = new Interval(1000, 1e-9,1e5);

      dim = mesh->topology().dim();
      // Create Function space
      V = new FokkerPlanck1D::FunctionSpace(*mesh);
      a = new FokkerPlanck1D::BilinearForm(*V,*V);
      // Define variational forms
      L = new FokkerPlanck1D::LinearForm(*V);
      break;
    }
    case 2: {
      cout << "Two dimensional problem detected" << endl;
     
      // Define solution space 
      //mesh = new Rectangle(-3.0,-3.0,3.0,3.0,20,20);
      mesh = new Rectangle(bmin[0],
			   bmin[1],
			   bmax[0],
			   bmax[1],
			   grid[0],
			   grid[1]);
      dim = mesh->topology().dim();
      // Create Function space
      V = new FokkerPlanck2D::FunctionSpace(*mesh);
      dual_V = new DualFokkerPlanck2D::FunctionSpace(*mesh);
      
      // Define variational forms
      a = new FokkerPlanck2D::BilinearForm(*V,*V);
      L = new FokkerPlanck2D::LinearForm(*V);


      dual_a = new DualFokkerPlanck2D::BilinearForm(*dual_V,*dual_V);
      dual_L = new DualFokkerPlanck2D::LinearForm(*dual_V);

      break;
    }
    case 3: {
      cout << "Three dimensional problem detected" << endl;

      // Define solution space 
      //mesh = new Box(-3.0,-3.0,-1.0,3.0,3.0,1.0,10,10,20);
      mesh = new Box(bmin[0],
		     bmin[1],
		     bmin[2],
		     bmax[0],
		     bmax[1],
		     bmax[2],
		     grid[0],
		     grid[1],
		     grid[2]);
      dim = mesh->topology().dim();
      // Create Function space
      V = new FokkerPlanck3D::FunctionSpace(*mesh);
      
      // Define variational forms
      a = new FokkerPlanck3D::BilinearForm(*V,*V);
      L = new FokkerPlanck3D::LinearForm(*V);
      break;
    }
    default:
      cout << "Invalid dimensionality: " << ops->getDim();
      cout << " valid values [1,2,3]" << endl;

      delete ops;
      exit(1);
      
    }
    


  
  //   DirichletBoundary db;
  // Constant zero(*mesh,0.0);
  // DirichletBC bc(*V,zero,db);


  Function u1(*V);
  Function u0(*V); 
  Function dual_u1(*dual_V);
  Function dual_u0(*dual_V);
 Constant dual_ui(1.);


  double END = pm.getTimeEnd();
  double dt  = pm.getDT(); 
  double t =  pm.getTimeBegin(); ; 
  Constant DT( dt);


 // Linear system
  Matrix A;
  Vector b;
  Vector r;
  
  switch(ops->getDim())
    {
      
      
    case 1: {
      ((FokkerPlanck1D::LinearForm*)L)->u0 = *ui;
      ((FokkerPlanck1D::LinearForm*)L)->dt = DT;
      ((FokkerPlanck1D::LinearForm*)L)->A = *drift;
      ((FokkerPlanck1D::LinearForm*)L)->B = *diffusion;      
      ((FokkerPlanck1D::BilinearForm*)a)->A = *drift;
      ((FokkerPlanck1D::BilinearForm*)a)->B = *diffusion;
      ((FokkerPlanck1D::BilinearForm*)a)->dt = DT;
      
      assemble(A,*((FokkerPlanck1D::BilinearForm*)a));
      break;
      
    } 
    case 2: {
      
      ((FokkerPlanck2D::LinearForm*)L)->u0 = *ui;
      ((FokkerPlanck2D::LinearForm*)L)->dt = DT;
      ((FokkerPlanck2D::LinearForm*)L)->A = *drift;
      ((FokkerPlanck2D::LinearForm*)L)->B = *diffusion;      
      ((FokkerPlanck2D::BilinearForm*)a)->A = *drift;
      ((FokkerPlanck2D::BilinearForm*)a)->B = *diffusion;
      ((FokkerPlanck2D::BilinearForm*)a)->dt = DT;
      
      assemble(A,*((FokkerPlanck2D::BilinearForm*)a));
      break; 
      
    }
    case 3: {
      
      ((FokkerPlanck3D::LinearForm*)L)->u0 = *ui;
      ((FokkerPlanck3D::LinearForm*)L)->dt = DT;
      ((FokkerPlanck3D::LinearForm*)L)->A = *drift;
      ((FokkerPlanck3D::LinearForm*)L)->B = *diffusion;      
      ((FokkerPlanck3D::BilinearForm*)a)->A = *drift;
      ((FokkerPlanck3D::BilinearForm*)a)->B = *diffusion;
      ((FokkerPlanck3D::BilinearForm*)a)->dt = DT;
      
      assemble(A,*((FokkerPlanck3D::BilinearForm*)a));
      break;
    }
      
}

  
  // Apply Boundary Condition
  // bc.apply(A);

  //std::string *str = new std::string(myConfig.Solver);
  //boost::algorithm::trim(*str);


  GenericLinearSolver *kry;
  
  if(pm.getSolver() == "LU") 
    {
      kry = new LUSolver();
    }
  else 
    {
      kry = new KrylovSolver(pm.getSolver(), "sor");
    }
  // 
  cout << "Using solver: " << pm.getSolver() << endl;



  File file("AdaptFEMFactoryOut.pvd");
File dual_file("DualAdaptFEMFactoryOut.pvd");	 

  u0 = *ui;

   //   u0 /= u0.sum();
   file << u0;
 

  // Define L2 error
  //  L2Error1D::LinearForm *l2 = new L2Error1D::LinearForm(*V);
            
  //   l2->A = *drift; 
  //   l2->B = *diffusion;
  //   l2->U = u0;
  //    l2->dt = DT;  

  //----------------------------------------
  // Timestep
  //----------------------------------------
  Progress p("Time-stepping");
   set_log_level(2);
  //  plot(u0);


   //
   // Idea of adaptive method
   // Solve primals solution on initial mesh to T
   // for each time-step calculate residuals (if possible)
   // 
   // At the final time T solve the dual equation backward in time
   // with the Goal function as source
   // Calculate cell and facet residuals with this information
   //



//
// Forward Solution
//

  while ( t<=END)
    {
      // Assemble RHS and solve
      switch(ops->getDim())
	{
	case 1: {
	  assemble(b,*((FokkerPlanck1D::LinearForm*)L));		
	  break;
	}
	case 2: {
	  assemble(b,*((FokkerPlanck2D::LinearForm*)L));
		
	  break;
	}
	case 3: {
	  assemble(b,*((FokkerPlanck3D::LinearForm*)L));
      
	  break;
	}
	}
	
#ifdef FENICS_1_0
      kry->solve(A, *(u1.vector()), b);	
#else
      kry->solve(A, u1.vector(), b);
#endif
      //Vector u1v = u1.vector();
      
      

      //u1 /= u1v.sum();
      // Save result
      file << u1;
      


      // Calculate error by inserting the calculated solution in the equations 
      // a(u,v) = l(u0,v)
      // Error || u - u_h || <= C || (0.5B^TB)^{-1/2}*[   ] ||
      //  l2->u = u1;
      //l2->U = u0;
#ifdef FENICS_1_0
        cout << " Residual FENICS:  " << residual(A, *(u1.vector()), b) << endl;	
#else
	cout << " Residual FENICS:  " << residual(A, u1.vector(), b) << endl;
#endif 
      // Transfer data
      switch(ops->getDim())
	{
	case 1: {
	  ((FokkerPlanck1D::LinearForm*)L)->u0 = u1;
	  break;
	}
	case 2: {
	  ((FokkerPlanck2D::LinearForm*)L)->u0 = u1;
	  break;
	}
	case 3: {
	  ((FokkerPlanck3D::LinearForm*)L)->u0 = u1;	  
	  break;
	}
	}

      t += dt;
      p = t/END;

    }


cout << "Forward solution done. Contiuing with backward solution" << endl;
//
// Backward solution
//
// Currently only supports 2D
//


//
// Calculate moment of the distribution function
//

Kernel2D::Functional *kernel = new Kernel2D::Functional(*mesh);

kernel->f = u1;
Constant   one(1.0);
kernel->K = one;



double kres =assemble(*kernel);

// Print kres
cout << "Kernel value is: " << kres << endl;

dual_u0 = Constant(kres);
// Initial condition the final distribution of f
 ((DualFokkerPlanck2D::LinearForm*)dual_L)->u0 = dual_u0;
      ((DualFokkerPlanck2D::LinearForm*)dual_L)->dt = DT;
      ((DualFokkerPlanck2D::LinearForm*)dual_L)->A = *drift;
      ((DualFokkerPlanck2D::LinearForm*)dual_L)->B = *diffusion;
      ((DualFokkerPlanck2D::BilinearForm*)dual_a)->A = *drift;
      ((DualFokkerPlanck2D::BilinearForm*)dual_a)->B = *diffusion;
      ((DualFokkerPlanck2D::BilinearForm*)dual_a)->dt = DT;

      assemble(A,*((DualFokkerPlanck2D::BilinearForm*)dual_a));

dual_file << dual_u0;

double tau = END;
while(tau>=0) {


     assemble(b,*((DualFokkerPlanck2D::LinearForm*)dual_L));
 
#ifdef FENICS_1_0
kry->solve(A, *(dual_u1.vector()), b);
#else
kry->solve(A, dual_u1.vector(), b);
#endif

dual_file << dual_u1;

((DualFokkerPlanck2D::LinearForm*)dual_L)->u0 = dual_u1;




tau -= dt;

}





      // destruct 
      switch(ops->getDim())
	{
	case 1: {
	  delete ((FokkerPlanck1D::BilinearForm*)a);
	  delete ((FokkerPlanck1D::LinearForm*)L);

	  break;
	}
	case 2: {
	  delete ((FokkerPlanck2D::BilinearForm*)a);
	  delete ((FokkerPlanck2D::LinearForm*)L);
	  break;
	}
	case 3: {
	  ((FokkerPlanck3D::BilinearForm*)a);
	  ((FokkerPlanck3D::LinearForm*)L);	  
	  break;
	}
	}
      

      delete V;

      delete mesh;
      //delete str;
      delete kry;

      
}

