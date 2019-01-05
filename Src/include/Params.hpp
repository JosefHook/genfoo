/**
 * Header file for Params.cpp
 *
 * Copyright (c) 2012 Josef Höök
 * joh@kth.se
 */
#ifndef _PARAMS_H
#define _PARAMS_H


#include <libxml/tree.h>

#include <libxml/parser.h> // 4 XSD validation 
#include <libxml/xmlschemas.h> // 4 XSD validation 

#include <string>
#include <vector>
class Params 
{

public:  
  
   Params();


 // ITM compatible constructor
  //Params(xmlDocPtr *params, void *itm) : _xml_params(params), _operator_p(itm) {} 


  // Main constructor
  Params(std::string pm_file, int argc, char *argv[]);
  
  std::string getFactory(void);
  std::string getParameterFile(void) { return _parameter_file; }
  std::string getLibraryFile(void); 
  void setLibraryFile( std::string jit_file ); 
  xmlDocPtr getXMLDoc(void) { return _xml_params; }
  
  int getDimension(void);
  int getNP(void); // Number of particles
  std::vector<double> getBoundaryMinValues(void);
  std::vector<double> getBoundaryMaxValues(void);
  std::vector<int> getGridDimension(void);
  
  
  double getTimeEnd(void);
  double getTimeBegin(void);
  double getDT(void);
  std::string getSolver(void);

  int getargc() { return _argc; }
  char** getargv() { return _argv; }
  // Validate schema file
  int validate(void); 
  
 
private: 


  std::string _solver;
  std::string _parameter_file;
  std::string _library_file;
  std::vector<int> _grid_dim;
  std::vector<double> _sb_minvalues;
  std::vector<double> _sb_maxvalues;
  double _timeBegin;
  double _timeEnd;
  double _dt; 
  xmlDocPtr _xml_params;
  
  int _argc; 
  char **_argv;
  void *_operator_p;  // a general object parameter for operator
  
  const char static _schema[];
};


#endif 
