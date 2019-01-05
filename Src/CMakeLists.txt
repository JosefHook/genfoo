#
# CMAKE build file for GenFoo 
# 
#

cmake_minimum_required(VERSION 2.6)

if(NOT GenFoo_SOURCE_DIR)
message(FATAL_ERROR "CMake must be run from genfoo root directory" )
endif()

include_directories(${GenFoo_SOURCE_DIR}/Src/include)
include_directories(${GenFoo_SOURCE_DIR}/Src/xml)
include_directories(${GenFoo_SOURCE_DIR}/Deps/include)
include_directories(${GenFoo_SOURCE_DIR}/Deps/include/boost-numeric-bindings)


add_definitions(-DFENICS_1_0)

get_property(inc_dirs DIRECTORY PROPERTY INCLUDE_DIRECTORIES)
#message("inc_dirs = ${inc_dirs}")

#
# Before doing anything compile the XSL and XML files in xml/ 
# to .hpp code
#
add_subdirectory(xml)

message("XML to hpp transform done")
set(CMAKE_CXX_FLAGS "-ggdb")

##########################
# Common source files
# Build a list "GSource" 
# with source files
##########################
list(APPEND GSource DynamicLoader Main Params JIT) 

# Build list of source files based on ccmake configuration options


##########################
# MC Factory 
##########################
if(BUILD_MC_BACKEND)
message("Including Monte Carlo factory")
add_definitions(-DBUILD_MC)
list(APPEND GSource DatIO  H5PartIO  MCFactory MonteCarlo ReflectedBoundaryCondition)
endif()

##########################
# DMC Factory 
##########################
if(BUILD_DMC_BACKEND)
message("Including delta-f Monte Carlo factory")
add_definitions(-DBUILD_DMC)

  if(BUILD_MC_BACKEND)
    list(APPEND GSource BallResample DMCFactory)
  else()
    list(APPEND GSource  BallResample DatIO  H5PartIO DMCFactory MonteCarlo ReflectedBoundaryCondition)
  endif()
endif()
##########################
# FEM Factory 
##########################
if(BUILD_FEM_BACKEND)
message("Including FEM factory")
#list(APPEND GSource FEMFactory AdaptFEMFactory)
list(APPEND GSource FEMFactory )
add_definitions(-DBUILD_FEM)

endif()


##########################
# FEM Factory 
##########################
if(BUILD_ADAPTIVE_FEM_BACKEND)
message("Including adaptive FEM factory")
#list(APPEND GSource FEMFactory AdaptFEMFactory)
list(APPEND GSource AdaptFEMFactory )
add_definitions(-DBUILD_ADAPTIVE_FEM)

endif()


# Print out the source list
#message(${GSource})




####################### 
# Find libraries
####################### 

 find_package(dolfin REQUIRED)
if(DOLFIN_FOUND)
  include_directories(${DOLFIN_INCLUDE_DIRS})
endif(DOLFIN_FOUND)

find_package(LibXml2 REQUIRED)
if(LIBXML2_FOUND)
  include_directories(${LIBXML2_INCLUDE_DIR})
endif(LIBXML2_FOUND)


find_package(LibXslt REQUIRED)
if(LIBXSLT_FOUND)
  include_directories(${LIBXSLT_INCLUDE_DIR})
endif(LIBXSLT_FOUND)


find_package(dl REQUIRED)
if(DL_FOUND)
  include_directories(${DL_INCLUDE_DIRS})
endif(DL_FOUND)

find_package(MPI REQUIRED)
if(MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif(MPI_FOUND)


#find_package(Boost 1.36 COMPONENTS filesystem program_options system thread iostreams  regex REQUIRED)
find_package(Boost 1.36 COMPONENTS program_options regex REQUIRED)
if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif(Boost_FOUND)


if(BUILD_MC_BACKEND OR BUILD_DMC_BACKEND)

find_package(LAPACK)
if(LAPACK_FOUND)
  include_directories(${LAPACK_INCLUDE_DIRS})
endif(LAPACK_FOUND)


endif()

find_package(ZLIB)
if(ZLIB_FOUND)
  include_directories(${ZLIB_INCLUDE_DIRS})
endif(ZLIB_FOUND)

if(BUILD_MC_BACKEND OR BUILD_DMC_BACKEND)

find_package(HDF5)
if(HDF5_FOUND)
  include_directories(${HDF5_INCLUDE_DIRS})
endif(HDF5_FOUND)


endif()

#########################
# Build executable
#########################
add_executable(GenFoo  ${GSource} )





####################### 
# Link libraries
####################### 

if(DOLFIN_FOUND)
  target_link_libraries(GenFoo ${DOLFIN_LIBRARIES})
endif(DOLFIN_FOUND)


if(LIBXML2_FOUND)
  target_link_libraries(GenFoo ${LIBXML2_LIBRARIES})
endif(LIBXML2_FOUND)


if(LIBXSLT_FOUND)
  target_link_libraries(GenFoo ${LIBXSLT_LIBRARIES})
endif(LIBXSLT_FOUND)


if(DL_FOUND)
  target_link_libraries(GenFoo ${DL_LIBRARIES})
endif(DL_FOUND)


if(Boost_FOUND)
  target_link_libraries(GenFoo ${Boost_LIBRARIES})
endif(Boost_FOUND)

if(BUILD_MC_BACKEND OR BUILD_DMC_BACKEND)


if(LAPACK_FOUND)
  target_link_libraries(GenFoo ${LAPACK_LIBRARIES})
endif(LAPACK_FOUND)

if(MPI_FOUND)
  target_link_libraries(GenFoo ${MPI_LIBRARIES})
endif(MPI_FOUND)


endif()

if(ZLIB_FOUND)
  target_link_libraries(GenFoo ${ZLIB_LIBRARIES})
endif(ZLIB_FOUND)



if(BUILD_MC_BACKEND OR BUILD_DMC_BACKEND)

if(HDF5_FOUND)
  target_link_libraries(GenFoo ${HDF5_LIBRARIES})
endif(HDF5_FOUND)

endif()


if(MPI_FOUND)
  target_link_libraries(GenFoo ${MPI_LIBRARIES})
endif()

if(MPI_COMPILE_FLAGS)
   set_target_properties(GenFoo PROPERTIES 
     COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(GenFoo PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()




###############################################
# Add glib library (needed for cmd line parsing
###############################################
#   add_library(imp_glib STATIC IMPORTED)
#   # point the imported target at the real file
#   set_property(TARGET imp_glib PROPERTY
#                IMPORTED_LOCATION  ${GenFoo_SOURCE_DIR}/Deps/lib/libglib-2.0.a)
#   target_link_libraries(GenFoo imp_glib)


########################################## 
# ANN library needed for delta-f MC
##########################################
if(BUILD_DMC_BACKEND)
 # Add ANN library  
   add_library(imp_ann STATIC IMPORTED)
   # point the imported target at the real file
   set_property(TARGET imp_ann PROPERTY
                IMPORTED_LOCATION  ${GenFoo_SOURCE_DIR}/Deps/lib/libANN.a)

   target_link_libraries(GenFoo imp_ann)

endif()


########################################## 
# h5part needed for Monte Carlo solvers
##########################################
if(BUILD_MC_BACKEND OR BUILD_DMC_BACKEND)
 # Add h5part library 
   add_library(imp_h5part STATIC IMPORTED)
   # point the imported target at the real file
   set_property(TARGET imp_h5part PROPERTY
                IMPORTED_LOCATION  ${GenFoo_SOURCE_DIR}/Deps/lib/libH5Part.a)

   target_link_libraries(GenFoo imp_h5part)

endif()



# Install binary in Build dir
install(TARGETS GenFoo DESTINATION  ${GenFoo_SOURCE_DIR}/Build)
