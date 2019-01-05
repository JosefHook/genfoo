//
// GenFoo Main file
// Generalized Fokker-Planck Solver 
// By Josef Höök, joh@kth.se, copyright (C), all rights reserved
//

#include <iostream>
#include "JIT.hpp"
#include "Params.hpp"


using namespace std;


void ShowInfo()
{

  cout << "                       GenFoo                            " << endl;
  cout << "           Generalized Fokker-Planck Solver.              " << endl;
  cout << "Version 1.2, 2012, by Josef Höök Copyright (C), joh@kth.se "<< endl;

}


// Define factories
#ifdef BUILD_MC
int MCFactory( Params pm );
#endif
#ifdef BUILD_DMC
int DMCFactory( Params pm );
#endif
#ifdef BUILD_FEM
int FEMFactory( Params pm );
int AdaptFEMFactory( Params pm );
#endif

#include <boost/program_options.hpp>
namespace po = boost::program_options;



int 
main(int argc, char *argv[]) 
{




  
  // General info
  ShowInfo();


  // Parse command line
  int opt;
  po::options_description desc("Allowed options");
  desc.add_options()
    ("help", "produce help message")
    ("parameter-file,p", po::value<std::string>(), "parameter file");


  if(argc==1) {
    std::cout << desc << std::endl;
    return 1;
  }

  
  po::positional_options_description p;
  p.add("parameter-file", -1);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).
	    options(desc).positional(p).run(), vm);
  po::notify(vm);
  if (vm.count("parameter-file"))
    {
      std::cout << vm["parameter-file"].as<std::string>() << std::endl;
      
    } else {
    std::cout << "No arguments specified on the command line" << std::endl;
  }
  
  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }
  


  std::string pm_file; 
  pm_file = vm["parameter-file"].as<std::string>() ;
  
  Params pm(pm_file, argc, argv); 

  std::string method  = pm.getFactory();

  // Test parameter functions


   //int foo = pm.validate();
  //cout << "Validate main: " << foo << endl;
   cout << " In Main: " << method << endl;
   cout << " getDimension(): " << pm.getDimension() << endl;

   cout << "getLibraryFile(): " << pm.getLibraryFile() << endl;
   cout << "getParameterFile(): " << pm.getParameterFile() << endl;
   cout << "getXMLDoc(): " << pm.getXMLDoc() << endl;
   cout << "getBoundaryMinValues(): " << endl; 
   pm.getBoundaryMinValues();
   cout << "getBoundaryMaxValues(): " << endl;
   pm.getBoundaryMaxValues();
   cout << " getGridDimension() " << endl;
   pm.getGridDimension();

   cout << "getTimeBegin(): " << pm.getTimeBegin() << endl;
   cout << "getTimeEnd(): " << pm.getTimeEnd() << endl;
   cout << "getDT(): " << pm.getDT() << endl;
   

   cout << "getSolver(): " << pm.getSolver() << endl;

   // Investigate if the Factory is an .so file or and .xml file
   string op_file = pm.getLibraryFile();
   
   // find ending .so or .xml 


   string file_ending = op_file.substr(op_file.find_last_of(".")+1);

   cout << "file ending :" << file_ending << ":" <<endl;
   if(!file_ending.compare("xml")) {
     cout << "Operator defined in XML. Calling JIT" << endl;  
     // Invoke the JIT compiler if .xml file
     JIT jit(op_file);
     if(jit.compileOp()) {
       std::cerr << "Failed to invoke Just-In-Time compilation of :"
		 << pm.getLibraryFile() << std::endl; 
       // End program
       return -1; 
     } else {
       std::cout << "JIT compilation successful" << endl;
	string path = "./" + jit.getLibraryFile();
       pm.setLibraryFile(path);
       
     }
     

   }
   

   




  if(method=="MonteCarlo") {
#ifdef BUILD_MC
    std::cout << "Using Monte Carlo solver" << std::endl;
    MCFactory(pm);

#else
    std::cout << " Monte Carlo solver not included in this build" << std::endl;
#endif

  }
  else if(method=="QMC") {
    std::cout << "Using quasi-Monte Carlo solver" << std::endl;
    std::cout << " NOTE IMPLEMENTED YET " << std::endl;
  }

  //    InitMC(myArguments, myConfig);

  else if(method == "DeltaFMonteCarlo") {

#ifdef BUILD_DMC
    std::cout << "Using delta-f Monte Carlo solver" << std::endl;
    DMCFactory( pm );
#else
    std::cout << "delta-f Monte Carlo solver not included in this build" << std::endl;
#endif





  }
    //    InitdMC(myArguments, myConfig);
  else if(method == "FEM") {

#ifdef BUILD_FEM
    std::cout << "Using FEM solver" << std::endl;
    FEMFactory( pm );
    
#else
    std::cout << "FEM solver not included in this build" << std::endl;
#endif


  }

  else if(method == "AFEM") {

#ifdef BUILD_ADAPTIVE_FEM
    std::cout << "Using adaptive FEM solver" << std::endl;
    AdaptFEMFactory( pm );
    
#else
    std::cout << "Adaptive FEM solver not included in this build" << std::endl;
#endif




  }
  else {
    std::cout << " Could not find the requested factory backend " << pm.getFactory() << std::endl;

  }




return 0;
} // End main
