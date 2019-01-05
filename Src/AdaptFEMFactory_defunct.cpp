//
// Adaptive FEM Factory by Josef Höök, 2012 
// Copyright (c) all rights reserved
//
#include <dolfin.h>

#include <boost/algorithm/string/trim.hpp>
//#include "include/L2Error1D.h"
#include "AdaptFokkerPlanck1D.h"
#include "AdaptFokkerPlanck2D.h"
#include "AdaptFokkerPlanck3D.h"
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
  if(dl.Load(   ((std::string) pm.getLibraryFile()).c_str() ))
    exit(1);
  
  Operator *ops; 
  Operator *(*load)() = (Operator*(*)())(dl.bootstrap());
  // Create new instance of coefficient object
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
  void *a, *L;
  FunctionSpace *V;  
  Mesh *mesh;

  void *M;



  // Domain information
  std::vector<int> grid = pm.getGridDimension();
  std::vector<double> bmax = pm.getBoundaryMaxValues();
  std::vector<double> bmin = pm.getBoundaryMinValues();
      

 
  switch(ops->getDim())
    {
      
    case 1: {
      cout << "One dimensional problem detected" << endl;

      // Define solution space 
      // TODO Generalize:
      // For Ornstein
      //mesh = new Interval(500, -3.0,3.0);
      // For RF-heating
      
      
      mesh = new Interval(grid[0], bmin[0], bmax[0]);
      //mesh = new Interval(1000, 1e-9,1e5);

      dim = mesh->topology().dim();
      // Create Function space
      V = new AdaptFokkerPlanck1D::BilinearForm::TrialSpace(*mesh);
      M = new AdaptFokkerPlanck1D::GoalFunctional(*mesh);
      // Define variational forms
      a = new AdaptFokkerPlanck1D::BilinearForm(*V,*V);
      L = new AdaptFokkerPlanck1D::LinearForm( *V);

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
      V = new AdaptFokkerPlanck2D::BilinearForm::TrialSpace(*mesh);
      M = new AdaptFokkerPlanck2D::GoalFunctional(*mesh);
      
      
      // Define variational forms
      a = new AdaptFokkerPlanck2D::BilinearForm(*V,*V);
      L = new AdaptFokkerPlanck2D::LinearForm(*V);
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
      V = new AdaptFokkerPlanck3D::BilinearForm::TrialSpace(*mesh);
      M = new AdaptFokkerPlanck3D::GoalFunctional(*mesh);
      
      // Define variational forms
      a = new AdaptFokkerPlanck3D::BilinearForm(*V,*V);
      L = new AdaptFokkerPlanck3D::LinearForm(*V);
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
  //  Function u1(*(AdaptFokkerPlanck2D::BilinearForm::TrialSpace*)V);
  //Function u0(*(AdaptFokkerPlanck2D::BilinearForm::TrialSpace*)V);



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
      ((AdaptFokkerPlanck1D::LinearForm*)L)->u0 = *ui;
      ((AdaptFokkerPlanck1D::LinearForm*)L)->dt = DT;
      ((AdaptFokkerPlanck1D::LinearForm*)L)->A = *drift;
      ((AdaptFokkerPlanck1D::LinearForm*)L)->B = *diffusion;      
      ((AdaptFokkerPlanck1D::BilinearForm*)a)->A = *drift;
      ((AdaptFokkerPlanck1D::BilinearForm*)a)->B = *diffusion;
      ((AdaptFokkerPlanck1D::BilinearForm*)a)->dt = DT;
      
      //      assemble(A,*((AdaptFokkerPlanck1D::BilinearForm*)a));
      break;
      
    } 
    case 2: {
      
      ((AdaptFokkerPlanck2D::LinearForm*)L)->u0 = *ui;
      ((AdaptFokkerPlanck2D::LinearForm*)L)->dt = DT;
      ((AdaptFokkerPlanck2D::LinearForm*)L)->A = *drift;
      ((AdaptFokkerPlanck2D::LinearForm*)L)->B = *diffusion;      
      ((AdaptFokkerPlanck2D::BilinearForm*)a)->A = *drift;
      ((AdaptFokkerPlanck2D::BilinearForm*)a)->B = *diffusion;
      ((AdaptFokkerPlanck2D::BilinearForm*)a)->dt = DT;
      
      //      assemble(A,*((AdaptFokkerPlanck2D::BilinearForm*)a));
      break; 
      
    }
    case 3: {
      
      ((AdaptFokkerPlanck3D::LinearForm*)L)->u0 = *ui;
      ((AdaptFokkerPlanck3D::LinearForm*)L)->dt = DT;
      ((AdaptFokkerPlanck3D::LinearForm*)L)->A = *drift;
      ((AdaptFokkerPlanck3D::LinearForm*)L)->B = *diffusion;      
      ((AdaptFokkerPlanck3D::BilinearForm*)a)->A = *drift;
      ((AdaptFokkerPlanck3D::BilinearForm*)a)->B = *diffusion;
      ((AdaptFokkerPlanck3D::BilinearForm*)a)->dt = DT;
      
      //      assemble(A,*((AdaptFokkerPlanck3D::BilinearForm*)a));
      break;
    }
      
}

  

  GenericLinearSolver *kry;
  
  if(pm.getSolver() == "LU") 
    {
      kry = new LUSolver();
    }
  else 
    {
      kry = new KrylovSolver(pm.getSolver(), "sor");

      kry->parameters["monitor_convergence"] = true;
      kry->parameters["error_on_nonconvergence"] = true; //false;
    }
  // 
  cout << "Using solver: " << pm.getSolver() << endl;



  File file("AdaptiveFEMFactoryOut.pvd");
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
  while ( t<=END)
    {

      // Adaptive TOL 
      double tol = 1e-5;
      


      // Assemble RHS and solve
      switch(ops->getDim())
	{
	case 1: {
	  //	  assemble(b,*((AdaptFokkerPlanck1D::LinearForm*)L));		

	  LinearVariationalProblem problem(  *((AdaptFokkerPlanck1D::BilinearForm*)   a),*((AdaptFokkerPlanck1D::LinearForm*)L), u1);
	  AdaptiveLinearVariationalSolver solver(problem);
	  //solver.parameters("error_control")("dual_variational_solver")["linear_solver"] = "gmres";
	  solver.solve(tol,*((AdaptFokkerPlanck1D::GoalFunctional*)M ) );


	  break;
	}
	case 2: {
	    //	  assemble(b,*((AdaptFokkerPlanck2D::LinearForm*)L));

	    LinearVariationalProblem problem(*((AdaptFokkerPlanck2D::BilinearForm*)a),*((AdaptFokkerPlanck2D::LinearForm*)L), u1);
	  AdaptiveLinearVariationalSolver solver(problem);

	  //solver.parameters("error_control")("dual_variational_solver")["linear_solver"] = "gmres";
	  solver.solve(tol,*((AdaptFokkerPlanck2D::GoalFunctional*)M)   );
	  
		
	  break;
	}
	case 3: {
	    //	  assemble(b,*((AdaptFokkerPlanck3D::LinearForm*)L));

	  LinearVariationalProblem problem(*((AdaptFokkerPlanck3D::BilinearForm*)a),*((AdaptFokkerPlanck3D::LinearForm*)L), u1);
	  AdaptiveLinearVariationalSolver solver(problem);
	  //	  solver.parameters("error_control")("dual_variational_solver")["linear_solver"] = "gmres";
	  solver.solve(tol,*((AdaptFokkerPlanck3D::GoalFunctional*)M)  );
      
      
	  break;
	}
	}


      //           kry->solve(A, u1.vector(), b);
      
    


      //      solve(A == L, u1, bc, tol, M); 

      //u1 /= u1v.sum();
      // Save result
       file << u1.fine();
      // Use the things below in 1.0-dev and onward
      //          file << u1.leaf_node();
      

      // Calculate error by inserting the calculated solution in the equations 
      // a(u,v) = l(u0,v)
      // Error || u - u_h || <= C || (0.5B^TB)^{-1/2}*[   ] ||
      //  l2->u = u1;
      //l2->U = u0;
      //assemble(r, *l2);
      //real residual = sqrt(fabs(r.sum()));
      //cout << " Residual error:  " << residual << endl;
 
      // Transfer data
      switch(ops->getDim())
	{
	case 1: {
	  ((AdaptFokkerPlanck1D::LinearForm*)L)->u0 = u1;
	  break;
	}
	case 2: {
	  ((AdaptFokkerPlanck2D::LinearForm*)L)->u0 = u1;
	  break;
	}
	case 3: {
	  ((AdaptFokkerPlanck3D::LinearForm*)L)->u0 = u1;	  
	  break;
	}
	}

      t += dt;
            p = t/END;


    }

      // Transfer data
      switch(ops->getDim())
	{
	case 1: {
	  delete ((AdaptFokkerPlanck1D::BilinearForm*)a);
	  delete ((AdaptFokkerPlanck1D::LinearForm*)L);

	  break;
	}
	case 2: {
	  delete ((AdaptFokkerPlanck2D::BilinearForm*)a);
	  delete ((AdaptFokkerPlanck2D::LinearForm*)L);
	  break;
	}
	case 3: {
	  ((AdaptFokkerPlanck3D::BilinearForm*)a);
	  ((AdaptFokkerPlanck3D::LinearForm*)L);	  
	  break;
	}
	}
      

      delete V;

      delete mesh;
      delete kry;

      
}

