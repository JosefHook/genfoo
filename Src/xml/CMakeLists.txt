#
# Before doing anything compile the XSL and XML files in xml/ 
# to .hpp code
#
#execute_process(COMMAND ls 
#  WORKING_DIRECTORY ${GenFoo_SOURCE_DIR}/Src/xml 
#  RESULT_VARIABLE e_res 
#  OUTPUT_VARIABLE e_out )
#message(STATUS "EmbedDoc: '${e_res}' : '${e_out}' ")


if(NOT GenFoo_SOURCE_DIR)
message(FATAL_ERROR "CMake must be run from genfoo root directory" )
endif()

file(GLOB files 
  RELATIVE ${GenFoo_SOURCE_DIR}/Src/xml
   "*.x*" )
foreach(filename ${files})
  message("Processing ${filename}")
  execute_process(COMMAND ${GenFoo_SOURCE_DIR}/Misc/EmbedDoc ${filename}
    WORKING_DIRECTORY ${GenFoo_SOURCE_DIR}/Src/xml )

endforeach()



