#
# CMAKE build file for GenFoo 
# 
#


cmake_minimum_required(VERSION 2.6)


if(NOT GenFoo_SOURCE_DIR)
message(FATAL_ERROR "CMake must be run from genfoo root directory" )
endif()

#get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#set(CMAKE_CXX_FLAGS "-ggdb")

#########################
# Build executable
#########################
#add_executable(EmbedDoc  EmbedDoc.cpp )

# Try an explicit build 
message("Building EmbedDoc with compiler: ${CMAKE_CXX_COMPILER}") 

execute_process(COMMAND ${CMAKE_CXX_COMPILER} EmbedDoc.cpp -o EmbedDoc
    WORKING_DIRECTORY ${GenFoo_SOURCE_DIR}/Misc 
    RESULT_VARIABLE retval )

  if(retval)
    message(FATAL_ERROR " Failed to compile Misc/EmbedDoc.cpp. Try to manually compile it!")
  else()
    message(" EmbedDoc compiled successfully")
  endif()