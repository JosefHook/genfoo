#
# CMAKE build file for GenFoo 
# 
#


include_directories("${GenFoo_SOURCE_DIR}/Src/include")
include_directories("${GenFoo_SOURCE_DIR}/Operators/CPP/PlasmaPhysics")

set(CMAKE_CXX_FLAGS "-ggdb")

#
# Operator class depend on XML library
#
find_package(LibXml2 REQUIRED)
if(LIBXML2_FOUND)
  include_directories(${LIBXML2_INCLUDE_DIR})
endif(LIBXML2_FOUND)



# Build share library 

add_library( 2DX SHARED 2DX.cpp)
add_library( 3DX SHARED 3DX.cpp)
add_library( DuffingVanDerPol SHARED DuffingVanDerPol.cpp)
add_library( Mutation SHARED Mutation.cpp)
add_library( NorthSouthFlowOnCircle SHARED NorthSouthFlowOnCircle.cpp)
add_library( NorthSouthFlowOnSphere SHARED NorthSouthFlowOnSphere.cpp)
add_library( OrnUhl1D SHARED OrnUhl1D.cpp)
add_library( OrnUhl2D SHARED OrnUhl2D.cpp)
add_library( OrnUhl3D SHARED OrnUhl3D.cpp)
add_library( Rf1D SHARED Rf1D.cpp)
add_library( RF1D SHARED RF1D.cpp)
add_library( Rf2D SHARED Rf2D.cpp)
add_library( RF2D SHARED RF2D.cpp)
add_library( RV2D SHARED RV2D.cpp RV2DParams.hpp)
add_library( RF3D SHARED RF3D.cpp)


if(LIBXML2_FOUND)
  target_link_libraries(2DX ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries(3DX ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( DuffingVanDerPol ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( Mutation ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( NorthSouthFlowOnCircle ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( NorthSouthFlowOnSphere ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( OrnUhl1D ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( OrnUhl2D ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( OrnUhl3D ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( Rf1D ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( RF1D ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( Rf2D ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( RF2D ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( RV2D ${LIBXML2_INCLUDE_LIBRARIES})
  target_link_libraries( RF3D ${LIBXML2_INCLUDE_LIBRARIES})
endif(LIBXML2_FOUND)

