CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(stator_gen)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++11")

ADD_EXECUTABLE(stator_gen stator_gen.cpp)
