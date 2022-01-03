# only for cmake --version >= 3.5.1
cmake_minimum_required(VERSION 3.5.1)

# project name
project(MDRun)

# creates the variable EXEC and sets it to hello
set(EXEC MDRun)

# set the C++11 standard
set(CMAKE_CXX_STANDARD 11)
add_definitions(-std=c++17 -O3 -DUSE_DOUBLES -DUSE_OPENMP -fopenmp -I/home/jonathan/programs/cpp/headers/eigen-eigen-5a0156e40feb)
# I../includes
include_directories(includes)
#include_directories(/home/jonathan/programs/cpp/headers/eigen-eigen-5a0156e40feb)

# puts all .cpp files inside src to the SOURCES variable
file(GLOB SOURCES src/*.cpp)

# compiles the files defined by SOURCES to generante the executable defined by EXEC
add_executable(${EXEC} ${SOURCES})


SET_TARGET_PROPERTIES(${EXEC} PROPERTIES COMPILE_FLAGS "-fopenmp -O3 "
                      LINK_FLAGS "-fopenmp")
