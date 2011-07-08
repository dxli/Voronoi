cmake_minimum_required(VERSION 2.6)
project(mps)

#There are lots of scripts with cmake
#for finding external libraries. 
#see /usr/local/share/cmake-2.6/Modules/Find*.cmake for more examples
#find_package(GLUT)
#find_package(OpenGL)
set(CMAKE_MODULE_PATH ".")
find_package(Grace)
set(CORELIBS ${GRACE_LIBRARY})
add_library(myGrace myGrace.cpp)
add_library(mkRandom mkRandom.cpp)
add_library(atom atom.cpp)
add_library(layer layer.cpp)
set(CORELIBS ${CORELIBS} myGrace mkRandom atom layer)


set(CMAKE_CXX_FLAGS "-g -O0 -march=native -fomit-frame-pointer -Wall")
set(PROGRAMS mps0)
foreach(program ${PROGRAMS})
  add_executable(${program} ${program}.cpp)
  target_link_libraries(${program} ${CORELIBS} ${GRACE_LIBRARY} myGrace mkRandom atom layer)
endforeach(program)

