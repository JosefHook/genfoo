#
# Top level CMAKE build file for GenFoo 
# 
# GenFoo binary is built in Src/
#
cmake_minimum_required(VERSION 2.6)


project(GenFoo)
option(BUILD_MC_BACKEND "Build Monte Carlo backend"  off)
option(BUILD_DMC_BACKEND "Build delta-f Monte Carlo backend"  off)
option(BUILD_FEM_BACKEND "Build Finite Element backend"  on)
option(BUILD_ADAPTIVE_FEM_BACKEND "Build adaptive Finite Element backend"  off)


option(BUILD_FOR_ITM "Build as ITM (Integrated Tokamak Modelling) compatible library " off)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")


#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${GenFoo_SOURCE_DIR}/Build)
#message("SOURCEDIR", ${GenFoo_SOURCE_DIR})

enable_language( Fortran )

# Build EmbedDoc binary in Misc/
# which is used to embed XSL docs
# into the source in Src
add_subdirectory(Misc)

add_subdirectory(Src)
add_subdirectory(Operators)





