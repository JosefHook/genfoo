//
// 
// This code converts XML file to a hpp file 
// with the text embedded in preprocessor MACRO 
//
//
// Input:  EmbedXML file.xml
//
// Output: file.hpp 
// 
// The MACRO name is: M_FILE_STR
// Ex: usage
// const static char[] Params::_file_str = "M_FILE_STR";
//
//

#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;

// 
// Rewrite \ to \\
// and 
// " to \"  
//  
string fix_string(string inp) 
{
  string ret;
  ret.clear();
  int size = inp.size();
  // Loop over characters
  for(int i=0; i<size; i++) 
    {
      
      switch(inp[i])
	{
	  // add extra backslash
	case '\\':
	  ret += "\\\\";
	  break;

	case '"':
	  ret += "\\\"";
	  break;

	default:
	  ret += inp[i];
	  break;
	}
      
    }
  ret += "\\n";	
  //ret += " &#xA; ";
  return ret;
}

int main(int argc, char *argv[]) 
{
  
  // argv 
  if(argc!=2) {
    cout << " No file specified on the command line " << endl;
    return 1;
  } 
  
  string filename(argv[1]);
  cout << "Reading file:" << argv[1] << endl;
  
  ifstream fin( argv[1], ios::in); 

  // search for ending
  size_t pos;
  pos = filename.find('.');
  
  // Extract substring
  string name;
  //  name = filename.substr(0,pos);



  // Replace dot with underscore
  filename.replace(pos, 1, "_");
  name = filename; 

  transform(name.begin(), name.end(),name.begin(), ::toupper);
  cout << "Preprocessor variable: PRE_" << name << "_STR" << endl; 


  string ofilename;
  ofilename = "PRE_" + name + ".hpp";

  /*struct stat buf;
  if (stat(filename.c_str(), &buf) != -1) 
    {
      cout << "File: " << filename << " exists. quiting!"; 
      
    }
  */


  cout << "Writing to: " << ofilename << endl;;
  ofstream fout( ofilename.c_str(), ios::out);
   
  
  
  

  // f.open(fname);
  
  string b;

  fout << "#define PRE_" << name <<"_STR \"";
  while(fin.good()) {
    if ( getline( fin, b ) )
      {
	fout <<  fix_string(b) << " ";
      }
    
  }
  fout <<"\"";
  fout << endl;
  
  fin.close();
  fout.close();
  cout << "done." << endl;
  
  
  

  return 0;
}

