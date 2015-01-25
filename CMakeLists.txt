################################################################################
project( "TOPCODER-ABCSpeedup" CXX )

cmake_minimum_required( VERSION 2.8.3 )

set( CMAKE_ECLIPSE_VERSION 3.8 )
set( _ECLIPSE_VERSION 3.8 )

set( CMAKE_CXX_COMPILER_ARG1 -std=c++11 ) ### for Eclipse's discovery extra arguments

#add_definitions( -O0 -ggdb -std=c++11 -Wall -pedantic )
#add_definitions( -O2 -std=c++11 -Wall -g -pedantic )
add_definitions( -O2 -std=c++11 -Wall -g -pedantic -fopenmp )
#add_definitions( -O2 -std=c++0x -ggdb -D__GXX_EXPERIMENTAL_CXX0X__ )

################################################################################

include_directories(
    src
)

################################################################################

add_executable( clonify0 src/clonify0.cc )
add_executable( clonify src/clonify.cc )
add_executable( main src/main.cpp )

target_link_libraries( main gomp)
target_link_libraries( clonify gomp)
target_link_libraries( clonify0 gomp)

################################################################################
