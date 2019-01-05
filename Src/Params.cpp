/**
 * Params.cpp 
 * Copyright (c) Josef Höök all rights reserved
 * joh@kth.se
 */


//#include "ITMObject.hpp"
#include "Params.hpp"
#include <fstream>
#include <iostream>
#include <boost/regex.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>

// Embedded schema file
#include "xml/PRE_PARAMETER_SCHEMA_XSD.hpp"


#define DBUG2

//------------------------------
// Parameter schema rules
//------------------------------


const char Params::_schema[] = PRE_PARAMETER_SCHEMA_XSD_STR;


// Argc, argv constructor 
Params::Params(std::string pm_file,int argc, char *argv[])
{


  // Now read the xml parameter file and local store it 
  _xml_params = xmlReadFile(pm_file.c_str(), NULL, 0); 
  if(_xml_params == NULL) {
    std::cerr << " Failed to parse parameter file: " << pm_file << std::endl;
    exit(1);
  }

  // Validate XML file
  if(validate() ==0) 
    exit(1);


  // local store pointer to arguments
  _argc = argc;
  _argv = argv; 
  _parameter_file = pm_file;
  _library_file.clear();


} // Params(argc,...)





// ---------------------------------------------------
// Validate XML file against the _schema XSD string
// ---------------------------------------------------
int Params::validate(void)
{

  xmlDocPtr schema_ptr = xmlReadMemory(_schema, strlen(_schema), "schema.xml", NULL,0);
  if (schema_ptr == NULL) {
    std::cerr << " Could not read schema. Exiting! " << std::endl;
    exit(1);
    
  } else {
#ifdef DBUG1
    std::cout << " Read schema file into xmlDocPtr " << std::endl;
#endif
  }
  xmlSchemaParserCtxtPtr parser_ctxt = xmlSchemaNewDocParserCtxt(schema_ptr);
  if (parser_ctxt == NULL) {
    std::cerr << " Could not create a XML parser. Exiting!" << std::endl;
    xmlFreeDoc(schema_ptr);
    exit(1);
  }
  xmlSchemaPtr schema = xmlSchemaParse(parser_ctxt);
  if (schema == NULL) {
    std::cerr << " Schema is not valid. Exiting! " << std::endl;
    xmlSchemaFreeParserCtxt(parser_ctxt);
    xmlFreeDoc(schema_ptr);
    exit(1);
  }
  xmlSchemaValidCtxtPtr valid_ctxt = xmlSchemaNewValidCtxt(schema);
  if (valid_ctxt == NULL) {
    std::cerr << " Validation context failed. Exiting!" << std::endl;
    
    xmlSchemaFree(schema);
    xmlSchemaFreeParserCtxt(parser_ctxt);
    xmlFreeDoc(schema_ptr);
    exit(1);
  }
  int retval = (xmlSchemaValidateDoc(valid_ctxt, _xml_params) == 0);
  xmlSchemaFreeValidCtxt(valid_ctxt);
  xmlSchemaFree(schema);
  xmlSchemaFreeParserCtxt(parser_ctxt);
  xmlFreeDoc(schema_ptr);
  
  //
  // If error in validation write the schema file to parameter_schema.xsd 
  //
  if(retval == 0) {
    std::ofstream fout;
    fout.open("parameter_schema.xsd");    
    std::string sc(_schema);
    fout << sc << std::endl;
    fout.close();
    std::cerr << "Compare input XML with Schema: parameter_schema.xsd" 
	      << std::endl;
  }
  // Return 1 if success
  return retval ? 1 : 0; 
  
}




//
// Retrieve library file
//
std::string Params::getLibraryFile(void) {

  if(_library_file.empty()) {
    
    xmlNodePtr cur;
    std::string lib;
    lib.clear();
    
    cur = xmlDocGetRootElement(_xml_params);
    
    // Traverse XML tree and search for <library>VALUE</library>
    cur = cur->xmlChildrenNode;
    xmlChar *key;
    while (cur != NULL) {
      // Step down 
      if((!xmlStrcmp(cur->name, (const xmlChar *)"library"))){
	cur = cur->xmlChildrenNode;
      }
      
      if((!xmlStrcmp(cur->name, (const xmlChar *)"filename"))){
	key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
	lib = std::string((const char*)key);
	boost::algorithm::trim(lib);
	xmlFree(key);
	// local store lib file
	_library_file = lib;
	return lib;
      }
      cur = cur->next;
    }
 
    
  } else {
    return _library_file;
  }
  
}


//
// set library file
//
void Params::setLibraryFile(std::string jit_file) {  
  _library_file = jit_file;
}





//
// Retrieve factory from parameter file or the command line
//
std::string Params::getFactory(void) 
{


  std::string factory;
  // Parse the XML file and search for factory tag
  xmlNodePtr cur;
  cur = xmlDocGetRootElement(_xml_params);
  
  if(cur==NULL) {
    std::cerr << " Empty XML input file. Exiting " << std::endl;
    xmlFreeDoc(_xml_params);
    exit(1);
  }
  
// Traverse XML tree and search for <factory>VALUE</factory>
  cur = cur->xmlChildrenNode;
  xmlChar *key;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"factory"))){
      cur = cur->xmlChildrenNode; // move down and search for <select>
    }
    if((!xmlStrcmp(cur->name, (const xmlChar *)"select"))){
      key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
      //      std::cout << "Factory found: " << key << std::endl;
      factory = std::string((const char*)key);
      boost::algorithm::trim(factory);
      xmlFree(key);
      return factory;
    }

   cur = cur->next;
  }
  
  //
  // Factory not found in XML file search on the command line
  //
  //  if(factory.empty()) {
  // }
}



//---------------------------------------
// Retrieve dimension from XML param file
//---------------------------------------
int Params::getDimension(void) {
  int dim;
  // Parse the XML file and search for factory tag
  xmlNodePtr cur;
  cur = xmlDocGetRootElement(_xml_params);
  if(cur==NULL) {
    std::cerr << " Empty XML input file. Exiting " << std::endl;
    xmlFreeDoc(_xml_params);
    exit(1);
  }
  
// Traverse XML tree and search for <factory>VALUE</factory>
  cur = cur->xmlChildrenNode;
  xmlChar *key;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"domain"))){
      cur = cur->xmlChildrenNode;
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"dimensions"))){
      key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
      std::string local_str((const char*) key);
      boost::algorithm::trim(local_str);
      dim = boost::lexical_cast<int>(local_str);
      return dim;
    }

   cur = cur->next;
  }
  

  std::cerr << " Could not find dimension element. Exiting!" << std::endl;


}








double Params::getTimeBegin(void ) {
  xmlNodePtr cur;
  cur = xmlDocGetRootElement(_xml_params);
  // Traverse XML tree and search for <domain>...</domain>
  cur = cur->xmlChildrenNode;
  xmlChar *key;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"domain"))){
      cur = cur->xmlChildrenNode;
    }
    // Search for <timeBegin>
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"timeBegin"))){
      key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
      std::string local_str((const char*) key);
      boost::algorithm::trim(local_str);
      _timeBegin = boost::lexical_cast<double>(local_str);
      break;
   }
    cur = cur->next;
 }
  return _timeBegin;
} 


double Params::getTimeEnd(void ) {
  xmlNodePtr cur;
  cur = xmlDocGetRootElement(_xml_params);
  // Traverse XML tree and search for <domain>...</domain>
  cur = cur->xmlChildrenNode;
  xmlChar *key;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"domain"))){
      cur = cur->xmlChildrenNode;
    }
    // Search for <timeBegin>
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"timeEnd"))){
      key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
      std::string local_str((const char*) key);
      boost::algorithm::trim(local_str);
      _timeEnd = boost::lexical_cast<double>(local_str);
      break;
   }
    cur = cur->next;
 }
  return _timeEnd;
} 

double Params::getDT(void ) {
  xmlNodePtr cur;
  cur = xmlDocGetRootElement(_xml_params);
  // Traverse XML tree and search for <domain>...</domain>
  cur = cur->xmlChildrenNode;
  xmlChar *key;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"domain"))){
      cur = cur->xmlChildrenNode;
    }
    // Search for <timeBegin>
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"dt"))){
      key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
      std::string local_str((const char*) key);
      boost::algorithm::trim(local_str);
      _dt = boost::lexical_cast<double>(local_str);
      break;
   }
    cur = cur->next;
 }
  return _dt;
} 


//---------------------------
// getSolver...
//--------------------------
std::string Params::getSolver(void ) {
  xmlNodePtr cur;
  cur = xmlDocGetRootElement(_xml_params);
  // Traverse XML tree and search for <factory>...</factory>
  cur = cur->xmlChildrenNode;
  xmlChar *key;
  bool found_select = false;
  std::string local_str;
  local_str.clear();

  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"factory"))){
      cur = cur->xmlChildrenNode;
    }
    // Search for <select>
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"select"))){
      key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
      local_str.assign((const char*) key, strlen((const char*)key) );
    }

    // Factory in select found now search for solver in subtree
    if(!local_str.empty()) {
      if((!xmlStrcmp(cur->name, (const xmlChar*)local_str.c_str()))) {
	cur = cur->xmlChildrenNode;
	found_select = true;
      }
    }

    /*
    if(local_str == "FEM") {
	// FEM Solver selected step down in FEM tree
	if ((!xmlStrcmp(cur->name, (const xmlChar *)"FEM"))){
	  cur = cur->xmlChildrenNode;
	}
	found_select = true;
    } 


    if(local_str == "AFEM") {
	// Adaptive FEM Solver selected step down in FEM tree
	if ((!xmlStrcmp(cur->name, (const xmlChar *)"AFEM"))){
	  cur = cur->xmlChildrenNode;
	}
	found_select = true;
    } 

    if(local_str == "MonteCarlo") {
	// MonteCarlo Solver selected step down in  tree
	if ((!xmlStrcmp(cur->name, (const xmlChar *)"MonteCarlo"))){
	  cur = cur->xmlChildrenNode;
	}
	found_select = true;
    } 


    if(local_str == "DeltaFMonteCarlo") {
	// MonteCarlo Solver selected step down in tree
	if ((!xmlStrcmp(cur->name, (const xmlChar *)"DeltaFMonteCarlo"))){
	  cur = cur->xmlChildrenNode;
	}
	found_select = true;
	} 
*/

    // If select factory found continue extracting solver
    if(found_select) {
      // Search for <select>
      if ((!xmlStrcmp(cur->name, (const xmlChar *)"solver"))){
	key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
	std::string local_str((const char*) key);
	boost::algorithm::trim(local_str);
	_solver = local_str; 
	break;
      }
      
    }
    
    cur = cur->next;
  }
  return _solver;
} 



//---------------------------------------
// Retrieve dimension from XML param file
//---------------------------------------
int Params::getNP(void) {
  int np;
  // Parse the XML file and search for factory tag
  xmlNodePtr cur;
  cur = xmlDocGetRootElement(_xml_params);
  
  // Traverse XML tree and search for <factory>VALUE</factory>
  cur = cur->xmlChildrenNode;
  xmlChar *key;
  while (cur != NULL) {
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"factory"))){
      cur = cur->xmlChildrenNode;
    }
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"MonteCarlo"))){
      cur = cur->xmlChildrenNode;
    }


    if ((!xmlStrcmp(cur->name, (const xmlChar *)"nParticles"))){
      key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
      std::string local_str((const char*) key);
      boost::algorithm::trim(local_str);
      np  =(int) boost::lexical_cast<double>(local_str);
      return np;
    }

   cur = cur->next;
  }
  

  std::cerr << "Could not find nParticles tag in XML!" << std::endl;
  return -1; 

}




//-----------------------------------------
// Retrieve boundary min values from XML file
//-----------------------------------------
std::vector<double> Params::getBoundaryMinValues(void) {

  //
  // If first time called
  // populate local variable _sb 
  // 

  xmlNodePtr cur = xmlDocGetRootElement(_xml_params);
 
  // Traverse XML tree and search for <domain></domain>
  cur = cur->xmlChildrenNode;
  xmlChar *key;
  while (cur != NULL) {
    // Walk down if true
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"domain"))){
      cur = cur->xmlChildrenNode;
    }
    // Walk down if true
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"boundingbox"))){
      cur = cur->xmlChildrenNode; 
    }
    // Extract minValues
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"minValues"))) {

      key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
      // Expect string on the form, " double double double ..."
      // Split the string and lexical cast each double
      std::string mykey((const char*)key);
#ifdef DBUG1
      std::cout << "myKey: " << mykey << std::endl << std::flush; 
#endif
      boost::regex rexp("\\s+");
      boost::sregex_token_iterator ti(mykey.begin(), mykey.end(), rexp, -1);
      boost::sregex_token_iterator  tj;

      // Loop over tokens " 0.2 2 12 ..." and drop " "
      while(ti !=tj) 
	{
	  std::string local_str(*ti++);
#ifdef DBUG1
	  std::cout << "local_str : " << local_str << std::endl << std::flush;
#endif
	  if(!local_str.empty()) { 	  
	    boost::algorithm::trim(local_str);
	    _sb_minvalues.push_back(boost::lexical_cast<double>(local_str));
	  }
	}
#ifdef DBUG2 
      for( std::vector<double>::iterator it = _sb_minvalues.begin(); it != _sb_minvalues.end(); ++it)
	std::cout << "MinValues are: " << *it  << std::endl; 
#endif
      break;
    }
    cur = cur->next;
  }
  
  return _sb_minvalues;
}



//-----------------------------------------
// Retrieve boundary max values from XML file
//-----------------------------------------
std::vector<double> Params::getBoundaryMaxValues(void) {

  //
  // If first time called
  // populate local variable _sb 
  // 

  xmlNodePtr cur = xmlDocGetRootElement(_xml_params);
 
  // Traverse XML tree and search for <domain></domain>
  cur = cur->xmlChildrenNode;
  xmlChar *key;
  while (cur != NULL) {
    // Walk down if true
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"domain"))){
      cur = cur->xmlChildrenNode;
    }
    // Walk down if true
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"boundingbox"))){
      cur = cur->xmlChildrenNode; 
    }
    // Extract minValues
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"maxValues"))) {

      key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
      // Expect string on the form, " double double double ..."
      // Split the string and lexical cast each double
      std::string mykey((const char*)key);
#ifdef DBUG1
      std::cout << "myKey: " << mykey << std::endl << std::flush; 
#endif
      boost::regex rexp("\\s+");
      boost::sregex_token_iterator ti(mykey.begin(), mykey.end(), rexp, -1);
      boost::sregex_token_iterator  tj;

      // Loop over tokens " 0.2 2 12 ..." and drop " "
      while(ti !=tj) 
	{
	  std::string local_str(*ti++);
#ifdef DBUG1
	  std::cout << "local_str : " << local_str << std::endl << std::flush;
#endif
	  if(!local_str.empty()) { 	 
	    boost::algorithm::trim(local_str); 
	    _sb_maxvalues.push_back(boost::lexical_cast<double>(local_str));
	  }
	}
#ifdef DBUG2
      for( std::vector<double>::iterator it = _sb_maxvalues.begin(); it != _sb_maxvalues.end(); ++it) 
	std::cout << "MaxValues are: " << *it  << std::endl; 
#endif

      break;
    }
    cur = cur->next;
  }
  
  return _sb_maxvalues;
}



//-----------------------------------------
// Retrieve grid dimension from XML file
//-----------------------------------------
std::vector<int> Params::getGridDimension(void) {


  xmlNodePtr cur = xmlDocGetRootElement(_xml_params);
  // Traverse XML tree and search for <factory></factory>
    cur = cur->xmlChildrenNode;
  xmlChar *key;

  while (cur != NULL) {
    // Walk down if true
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"factory"))){
      cur = cur->xmlChildrenNode;
    }
    // Walk down if true


    
    if ((!xmlStrcmp(cur->name, (const xmlChar *)(getFactory()).c_str()  ))){
      cur = cur->xmlChildrenNode; 
    }
    // Extract gridSize values
    if ((!xmlStrcmp(cur->name, (const xmlChar *)"gridSize"))) {
      key = xmlNodeListGetString(_xml_params, cur->xmlChildrenNode, 1);
      // Expect string on the form, " int int int ..."
      // Split the string and lexical cast each int
      std::string mykey((const char*)key);
#ifdef DBUG1
      std::cout << "myKey: " << mykey << std::endl << std::flush; 
#endif
      boost::regex rexp("\\s+");
      boost::sregex_token_iterator ti(mykey.begin(), mykey.end(), rexp, -1);
      boost::sregex_token_iterator  tj;
      // Loop over tokens " 2 2 12 ..." and drop " "
      while(ti !=tj) 
	{
	  std::string local_str(*ti++);
#ifdef DBUG1
	  std::cout << "local_str : " << local_str << std::endl << std::flush;
#endif
	  if(!local_str.empty()) {	 
	    boost::algorithm::trim(local_str);
	    _grid_dim.push_back(boost::lexical_cast<int>(local_str));
	  }
	}
#ifdef DBUG2
      for( std::vector<int>::iterator it = _grid_dim.begin(); it != _grid_dim.end(); ++it) 
	std::cout << "Grid dim : " << *it  << std::endl; 
#endif

      break;
    }
    cur = cur->next;
  }
  return _grid_dim;
}


