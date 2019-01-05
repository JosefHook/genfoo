# - Try to find DL
# Once done, this will define
#
# Once done this will define
#  DL_FOUND - System has dl
#  DL_INCLUDE_DIRS - The dl include directories
#  DL_LIBRARIES - The libraries needed to use dl
#  DL_DEFINITIONS - Compiler switches required for using dl

find_path(DL_INCLUDE_DIR dlfcn.h
  HINTS /usr/include)


find_library(DL_LIBRARY NAMES dl libdl
             HINTS /usr/lib)

set(DL_LIBRARIES ${DL_LIBRARY} )
set(DL_INCLUDE_DIRS ${DL_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set DL_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(dl  DEFAULT_MSG
                                  DL_LIBRARY DL_INCLUDE_DIR)

mark_as_advanced(DL_INCLUDE_DIR DL_LIBRARY )