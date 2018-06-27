cmake_minimum_required( VERSION 3.0 )

project( ONCOLIB VERSION 1.0 )

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )

set( python_mod_src
  src/utils.cpp
  src/clonetree.cpp
  src/mutclonetree.cpp
  src/basetree.cpp
  src/migrationtree.cpp
  src/migrationgraph.cpp
  src/cell.cpp
  src/simulation.cpp
  src/basematrix.cpp
  src/readmatrix.cpp
  src/generatemixture.cpp
  src/frequencymatrix.cpp
  src/enumeratemutationtrees.cpp
  src/spruce/perfectphylograph.cpp
  src/spruce/perfectphylomatrix.cpp
  src/spruce/perfectphylotree.cpp
  src/spruce/realmatrix.cpp
  src/spruce/matrix.cpp
  src/spruce/realtensor.cpp
  src/spruce/rootedcladisticancestrygraph.cpp
  src/spruce/rootedcladisticnoisyancestrygraph.cpp
  src/spruce/rootedcladisticenumeration.cpp
  src/spruce/rootedcladisticnoisyenumeration.cpp
  src/spruce/solution.cpp
  src/spruce/solutionset.cpp
  src/spruce/statetree.cpp
  src/spruce/tensor.cpp
)

set( python_mod_hdr
  src/cell.h
  src/simulation.h
  src/beta_distribution.hpp
  src/clonetree.h
  src/mutclonetree.h
  src/basetree.h
  src/migrationtree.h
  src/migrationgraph.h
  src/basematrix.h
  src/readmatrix.h
  src/frequencymatrix.h
  src/generatemixture.h
  src/utils.h
  src/enumeratemutationtrees.h
  src/spruce/matrix.h
  src/spruce/realmatrix.h
  src/spruce/realtensor.h
  src/spruce/rootedcladisticancestrygraph.h
  src/spruce/rootedcladisticnoisyancestrygraph.h
  src/spruce/rootedcladisticenumeration.h
  src/spruce/rootedcladisticnoisyenumeration.h
  src/spruce/solution.h
  src/spruce/solutionset.h
  src/spruce/statetree.h
  src/spruce/tensor.h
)

set( sequence_hdr
  src/beta_distribution.hpp
  src/clonetree.h
  src/mutclonetree.h
  src/basetree.h
  src/utils.h
  src/generatemixture.h
  src/basematrix.h
  src/frequencymatrix.h
  src/readmatrix.h
)

set( sequence_src
  src/utils.cpp
  src/clonetree.cpp
  src/mutclonetree.cpp
  src/basetree.cpp
  src/sequencemain.cpp
  src/basematrix.cpp
  src/frequencymatrix.cpp
  src/readmatrix.cpp
)

set( mix_hdr
  src/beta_distribution.hpp
  src/clonetree.h
  src/mutclonetree.h
  src/basetree.h
  src/utils.h
  src/generatemixture.h
  src/basematrix.h
  src/frequencymatrix.h
  src/readmatrix.h
)

set( mix_src
  src/utils.cpp
  src/clonetree.cpp
  src/mutclonetree.cpp
  src/basetree.cpp
  src/mixmain.cpp
  src/generatemixture.cpp
  src/basematrix.cpp
  src/frequencymatrix.cpp
  src/readmatrix.cpp
)

set( graph2dot_hdr
  src/migrationgraph.h
  src/basetree.h
  src/clonetree.h
  src/utils.h
)

set( graph2dot_src
  src/utils.cpp
  src/clonetree.cpp
  src/basetree.cpp
  src/migrationgraph.cpp
  src/graph2dotmain.cpp
)

set( tree2dot_hdr
  src/beta_distribution.hpp
  src/clonetree.h
  src/mutclonetree.h
  src/basetree.h
  src/utils.h
  src/basematrix.h
  src/frequencymatrix.h
  src/readmatrix.h
)

set( tree2dot_src
  src/utils.cpp
  src/clonetree.cpp
  src/mutclonetree.cpp
  src/basetree.cpp
  src/tree2dotmain.cpp
  src/basematrix.cpp
  src/readmatrix.cpp
  src/frequencymatrix.cpp
)

set( tree2freqs_hdr
  src/beta_distribution.hpp
  src/clonetree.h
  src/mutclonetree.h
  src/basetree.h
  src/utils.h
  src/basematrix.h
  src/frequencymatrix.h
  src/readmatrix.h
)

set( tree2freqs_src
  src/utils.cpp
  src/clonetree.cpp
  src/mutclonetree.cpp
  src/basetree.cpp
  src/tree2freqsmain.cpp
  src/basematrix.cpp
  src/frequencymatrix.cpp
  src/readmatrix.cpp
)

set( downsample_hdr
  src/basematrix.h
  src/readmatrix.h
  src/frequencymatrix.h
  src/utils.h
)

set( downsample_src
  src/utils.cpp
  src/downsamplemain.cpp
  src/basematrix.cpp
  src/readmatrix.cpp
  src/frequencymatrix.cpp
)

set( simulate_hdr
  src/cell.h
  src/simulation.h
  src/beta_distribution.hpp
  src/clonetree.h
  src/mutclonetree.h
  src/basetree.h
  src/migrationtree.h
  src/migrationgraph.h
  src/utils.h
  src/basematrix.h
  src/readmatrix.h
  src/frequencymatrix.h
)

set( simulate_src
  src/utils.cpp
  src/clonetree.cpp
  src/mutclonetree.cpp
  src/basetree.cpp
  src/migrationtree.cpp
  src/migrationgraph.cpp
  src/cell.cpp
  src/simulatemain.cpp
  src/simulation.cpp
  src/basematrix.cpp
  src/readmatrix.cpp
  src/frequencymatrix.cpp
)

set( reads2freqs_src
  src/utils.cpp
  src/reads2freqsmain.cpp
  src/basematrix.cpp
  src/readmatrix.cpp
  src/frequencymatrix.cpp
)

set( reads2freqs_hdr
  src/basematrix.h
  src/readmatrix.h
  src/frequencymatrix.h
  src/utils.h
)

set( cluster_hdr
  src/bronkerbosch.h
  src/utils.h
  src/basematrix.h
  src/readmatrix.h
  src/frequencymatrix.h
  src/cluster/cluster.h
)

set( cluster_src
  src/bronkerbosch.cpp
  src/utils.cpp
  src/basematrix.cpp
  src/frequencymatrix.cpp
  src/readmatrix.cpp
  src/cluster/clustermain.cpp
  src/cluster/cluster.cpp
)

set( precluster_hdr
  src/utils.h
  src/basematrix.h
  src/frequencymatrix.h
)

set( precluster_src
  src/cluster/preclustermain.cpp
  src/utils.cpp
  src/basematrix.cpp
  src/frequencymatrix.cpp
)

set( enumerate_hdr
  src/utils.h
  src/basematrix.h
  src/frequencymatrix.h
  src/readmatrix.h
  src/clonetree.h
  src/mutclonetree.h
  src/enumeratemutationtrees.h
  src/spruce/matrix.h
  src/spruce/realmatrix.h
  src/spruce/realtensor.h
  src/spruce/rootedcladisticancestrygraph.h
  src/spruce/rootedcladisticnoisyancestrygraph.h
  src/spruce/rootedcladisticenumeration.h
  src/spruce/rootedcladisticnoisyenumeration.h
  src/spruce/solution.h
  src/spruce/solutionset.h
  src/spruce/statetree.h
  src/spruce/tensor.h
)

set( enumerate_src
  src/enumeratemain.cpp
  src/utils.cpp
  src/basematrix.cpp
  src/frequencymatrix.cpp
  src/readmatrix.cpp
  src/enumeratemutationtrees.cpp
  src/basetree.cpp
  src/clonetree.cpp
  src/mutclonetree.cpp
  src/spruce/perfectphylograph.cpp
  src/spruce/perfectphylomatrix.cpp
  src/spruce/perfectphylotree.cpp
  src/spruce/realmatrix.cpp
  src/spruce/matrix.cpp
  src/spruce/realtensor.cpp
  src/spruce/rootedcladisticancestrygraph.cpp
  src/spruce/rootedcladisticnoisyancestrygraph.cpp
  src/spruce/rootedcladisticenumeration.cpp
  src/spruce/rootedcladisticnoisyenumeration.cpp
  src/spruce/solution.cpp
  src/spruce/solutionset.cpp
  src/spruce/statetree.cpp
  src/spruce/tensor.cpp
)

# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

if (CMAKE_BUILD_TYPE MATCHES "^[Rr]elease")
  # check if Doxygen is installed
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
      # set input and output files
      set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
      set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

      # request to configure the file
      configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
      message("Doxygen build started")

      # note the option ALL which allows to build the docs together with the application
      add_custom_target( doc_doxygen ALL
          COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
          WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
          COMMENT "Generating API documentation with Doxygen"
          VERBATIM )
  else (DOXYGEN_FOUND) 
    message("Doxygen need to be installed to generate the doxygen documentation")
  endif (DOXYGEN_FOUND)
endif()

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )

if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unused-local-typedefs" )
endif()

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: Debug Release Profile" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

set( GitCommand1 "symbolic-ref HEAD 2> /dev/null | cut -b 12-" )
execute_process( COMMAND git symbolic-ref HEAD  WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
COMMAND tail -c -7 COMMAND tr -d "\\n" OUTPUT_VARIABLE GitBranch )
set( GitCommand2 "log --pretty=format:\\\"%h\\\" -1" )
execute_process( COMMAND git log --pretty=format:%h WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
COMMAND head -n 1  COMMAND tr -d "\\n" OUTPUT_VARIABLE GitHashtag )

#set( Boost_USE_STATIC_LIBS ON )
find_package(PythonLibs)
if( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  find_package( Boost REQUIRED COMPONENTS thread system filesystem python2.7 )
else( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  find_package( Boost REQUIRED COMPONENTS thread system filesystem python )
endif()
MESSAGE( STATUS "found boost:"
"\nINCLUDE: ${Boost_INCLUDE_DIRS}"
"\nLIB: ${Boost_LIBRARIES}"
)

FIND_PACKAGE(PythonInterp REQUIRED)
FIND_PACKAGE(PythonLibsNew REQUIRED)
FIND_PACKAGE(Armadillo REQUIRED)

include_directories( "${LIBLEMON_ROOT}/include" "src" ${Boost_INCLUDE_DIRS} )
link_directories( "${LIBLEMON_ROOT}/lib" )

# custom macro with most of the redundant code for making a python example module
MACRO( addPythonExe _name _srccpp )
ADD_EXECUTABLE(${_name} ${_srccpp} ${python_mod_src} ${python_mod_hdr} )
# make the pyd library link against boost_numpy python and boost
TARGET_INCLUDE_DIRECTORIES(${_name} PUBLIC "${PROJECT_SOURCE_DIR}/src" ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} ${ARMADILLO_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${_name} ${CommonLibs} boost_python2.7)
ENDMACRO()

MACRO( addPythonMod _name _srccpp )
PYTHON_ADD_MODULE(${_name} ${_srccpp} ${python_mod_src} ${python_mod_hdr} )
# make the pyd library link against boost_numpy python and boost
TARGET_INCLUDE_DIRECTORIES(${_name} PUBLIC "${PROJECT_SOURCE_DIR}/src" ${PYTHON_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS} ${GUROBI_INCLUDE_DIR} ${ARMADILLO_INCLUDE_DIRS})
TARGET_LINK_LIBRARIES(${_name} ${CommonLibs})
ENDMACRO()

set( CommonLibs
  m
  pthread
  emon
  ${Boost_LIBRARIES}
  ${ARMADILLO_LIBRARIES}
)

set( CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG" )

add_executable( simulate ${simulate_src} ${simulate_hdr} )
target_link_libraries( simulate ${CommonLibs} )

add_executable( mix ${mix_src} ${mix_hdr} )
target_link_libraries( mix ${CommonLibs} )

add_executable( sequence ${sequence_src} ${sequence_hdr} )
target_link_libraries( sequence ${CommonLibs} )

add_executable( downsample ${downsample_src} ${downsample_hdr} )
target_link_libraries( downsample ${CommonLibs} )

add_executable( tree2freqs ${tree2freqs_src} ${tree2freqs_hdr} )
target_link_libraries( tree2freqs ${CommonLibs} )

add_executable( reads2freqs ${reads2freqs_src} ${reads2freqs_hdr} )
target_link_libraries( reads2freqs ${CommonLibs} )

add_executable( tree2dot ${tree2dot_src} ${tree2dot_hdr} )
target_link_libraries( tree2dot ${CommonLibs} )

add_executable( graph2dot ${graph2dot_src} ${graph2dot_hdr} )
target_link_libraries( graph2dot ${CommonLibs} )

add_executable( cluster ${cluster_src} ${cluster_hdr} )
target_link_libraries( cluster ${CommonLibs} )

add_executable( enumerate ${enumerate_src} ${enumerate_hdr} )
target_link_libraries( enumerate ${CommonLibs} )

add_executable( precluster ${precluster_src} ${precluster_hdr} )
target_link_libraries( precluster ${CommonLibs} )

addPythonMod( oncolib src/python.cpp )
