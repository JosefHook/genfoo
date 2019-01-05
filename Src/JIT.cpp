//
// Functions for JIT class
//
// By Josef Höök and Thomas Johnson
// Copyright (C) all rights reserved 
//

#include "JIT.hpp"

// Pipe and fork
#include <unistd.h>
#include <sys/wait.h>
#include <string.h>

// XSLT
#include <libxslt/transform.h>
#include <libxslt/xsltutils.h>


// Embedded XLST file
#include "xml/PRE_OPERATOR_XSL.hpp"

// Embedded schema file 
#include "xml/PRE_OPERATOR_XSD.hpp"


#include <iostream>


/*--------------------------------------------*/
// Define XSLT document 
const  char JIT::_xsl[] = PRE_OPERATOR_XSL_STR;

// Define XSD schema
const  char JIT::_schema[] = PRE_OPERATOR_XSD_STR;

// load and parse operator file
 JIT::JIT(std::string operator_file) : _op_xml_file(operator_file) {}



// Register and parse operator file
void JIT::registerOp( std::string operator_file) {

  _op_xml_file = operator_file;
}

//
// Transforms operator xml file to .cpp code 
// and compile it to a shared library using g++
//
int JIT::compileOp(void) {
  int child_status;
  int pipefd[2];
  pid_t cpid;



  
    // XSLT part open file and parse
    xsltStylesheetPtr cur = NULL;
    xmlDocPtr doc, res;
    xmlSubstituteEntitiesDefault(1);
    xmlLoadExtDtdDefaultValue = 1;
    // cur = xsltParseStylesheetFile((const xmlChar *)"operator.xsl");
    xmlDocPtr xslDoc = xmlReadMemory(_xsl, strlen(_xsl), "null.xml", NULL,0);
    cur = xsltParseStylesheetDoc(xslDoc);
    
    if(cur==NULL) {
      std::cerr << "Failed to  parse stylesheet" << std::endl;
      xsltCleanupGlobals();
      xmlCleanupParser();
      return 1;
    }

    doc = xmlParseFile(_op_xml_file.c_str());

    if(doc==NULL) {
      std::cerr << "Failed to parse xml file" << std::endl;
      xsltFreeStylesheet(cur);
      xsltCleanupGlobals();
      xmlCleanupParser();
      return 1;
    }

    res = xsltApplyStylesheet(cur, doc, NULL);

    if(res==NULL) {
      std::cerr << "Failed to apply stylesheet" << std::endl;
      xsltFreeStylesheet(cur);
      xmlFreeDoc(doc);
      xsltCleanupGlobals();
      xmlCleanupParser();
      return 1;
    }



  // Create a pipe
  if (pipe(pipefd) == -1) {
    std::cerr <<"Failed to create pipe"<< std::endl;
    return 1;
  }
  // Fork a child process
  cpid = fork();
  if (cpid == -1) {
    std::cerr << "Failed to fork" << std::endl;
    return 1;
  }
  
  if (cpid == 0) {    // Child process 
    close(pipefd[1]); // close write end of pipe
    // duplicate file descriptor and connect it to stdin
    dup2(pipefd[0],0);

    // Write a copy of Operator.hpp to a temporary directory and link 
    // it to g++
    
    // compile the code to a shared library file
    execlp("g++", "g++", "-xc++","-fPIC", 
	   "-shared", 
	   "-I/usr/include/libxml2", "-lxml2",
	   "-I../Src/include", 
	   "-o", "libXMLDefinedOperator.so", "-", NULL);
    // If execlp failed this code will be executed
    std::cerr << "Call to g++ failed" << std::endl;
    close(pipefd[0]);
    return 1;
  } else {            // Parent process 
    close(pipefd[0]); //Close read end of pipe
    

    // Open a write stream to pipe
    FILE *stream = fdopen(pipefd[1], "w");
    xsltSaveResultToFile(stream, res, cur);	   // Write to pipe
    close(pipefd[1]);          // send EOF
    // Wait for child to finish
    wait(&child_status);
   
    // Debug send code to STDOUT
      xsltSaveResultToFile(stdout, res, cur);	   // Write to stdout

    // If compilation failed print out the code
    if( WEXITSTATUS(child_status) == 1) {
          std::cout << "Compilation failed! " << std::endl;
      std::cout << "-----------------------------------------------"
		<< std::endl;
      std::cout << " Code sent to g++ -----------------------------" 
		<< std::endl;  
      std::cout << "-----------------------------------------------"
		<< std::endl;
      
      xsltSaveResultToFile(stdout, res, cur);	   // Write to stdout
      
      std::cout << "-----------------------------------------------"
		<< std::endl;
      
      std::cout << " Code sent to g++ -----------------------------" 
		<< std::endl; 
      std::cout << "-----------------------------------------------"
		<< std::endl;


    xsltFreeStylesheet(cur);
    xmlFreeDoc(res);
    xmlFreeDoc(doc);
    xsltCleanupGlobals();
    xmlCleanupParser();
    return 1;
      
    }

    xsltFreeStylesheet(cur);
    xmlFreeDoc(res);
    xmlFreeDoc(doc);
    xsltCleanupGlobals();
    xmlCleanupParser();



  } // End of parent process
  
 return 0;
}


std::string JIT::getLibraryFile(void) {
  
  std::string str("libXMLDefinedOperator.so");
  return str;


}



