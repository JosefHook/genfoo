//
// A Just-In-Time compiler class for compilation 
// op operators defined in XML
// 
// Input: XML file 
// Output: A shared library  (.so file) with stubs
//         for runtime loading.
//
//
#pragma once
#ifndef _JIT_H
#define _JIT_H


// XSLT
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>

#include <string>


class JIT 
{

public:
  JIT();
  JIT( std::string operator_file );
  void registerOp( std::string operator_file );
  std::string getLibraryFile(void);
  int compileOp( void );
  //int compileOp( xmlDocPtr xsl );
  
private:

  const static char _xsl[];   // XSLT document 
  const static char  _schema[]; // Schema file
  //  xmlDocPtr _op;
  std::string _op_xml_file; 
  

};


#endif
