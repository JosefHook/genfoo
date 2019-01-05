//
// A  Dynamic library loader
//
//
#include "DynamicLoader.hpp"
#include <iostream>
#include <dlfcn.h>
#include <sys/stat.h>






int DynamicLoader::Load(char *file)
{


        // Check if library exists
      struct stat st;
      if(stat(file, &st) ==0) 
	{
	  std::cout << "Loading coefficients from: "  << file << std::endl;
	  
	  _object = NULL;
	  _object = dlopen(file, RTLD_NOW);
	  _error = dlerror();
	  
	  if(_object == NULL)
	    {
	      std::cout << "Error loading library " << _error << std::endl;
	      return 1;
	      
	    }
	  
	  
	  
	  
	}
      else 
	{
	  std::cout << "Library does not exists" << std::endl;
	  return 1;
	}      
      
      
      
      
      
      
      return 0;
}







int DynamicLoader::Load(const char *file)
{


        // Check if library exists
      struct stat st;
      if(stat(file, &st) ==0) 
	{
	  std::cout << "Loading coefficients from: "  << file << std::endl;
	  
	  _object = NULL;
	  _object = dlopen(file, RTLD_NOW);
	  _error = dlerror();
	  
	  if(_object == NULL)
	    {
	      std::cout << "Error loading library " << _error << std::endl;
	      return 1;
	      
	    }
	  
	  
	  
	  
	}
      else 
	{
	  std::cout << "Library does not exists" << std::endl;
	  return 1;
	}      
      
      
      
      
      
      
      return 0;
}




int DynamicLoader::Load(int argc, char *argv[])
{


  // Arguments specified.
  // Go ahead and check library
  if(argc>=2) 
    {
      
      // Check if library exists
      struct stat st;
      if(stat(argv[1], &st) ==0) {
	std::cout << "Loading operators from: "  << argv[1] << std::endl;

	_object = NULL;
	_object = dlopen(argv[1], RTLD_NOW);
	_error = dlerror();
	
	if(_object == NULL)
	  {
	    std::cout << "Error loading library " << _error << std::endl;
	    return 1;

	  }
	

	

      }
      else {
	std::cout << "Library does not exists" << std::endl;
	return 1;
      }      
     

    } 
  else 
    {
      std::cout << "Could not parse arguments" << std::endl;
      return 1;
    }
  


  return 0;
}



void* DynamicLoader::getObject()
{

  return _object;

}



void* DynamicLoader::bootstrap()
{
  return dlsym(_object, "load");
}
