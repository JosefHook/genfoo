//
// A Dynamic loader class
//

#ifndef _DYN_LOAD
#define _DYN_LOAD

class DynamicLoader 
{


public:

  int Load(char *file);
  int Load(const char *file);

  int Load(int argc, char *argv[]);
  void* getObject();
  void* bootstrap();


private:

  void* _object;
  char* _error;

};


#endif
