cmake_minimum_required(VERSION 3.15)
project(p2pDDSketch)

set(CMAKE_CXX_STANDARD 14)

add_executable(p2pDDSketch main.cpp ddsketch.cc ddsketch.h graph.cpp graph.h error.cpp error.h)

FIND_PACKAGE( Boost COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )
target_link_libraries(p2pDDSketch LINK_PUBLIC igraph ${Boost_LIBRARIES} )
