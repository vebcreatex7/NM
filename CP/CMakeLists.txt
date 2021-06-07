cmake_minimum_required(VERSION 3.16.3)

project(cp)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall -ggdb")

add_executable(cp src/main.cpp src/calculation.cpp src/matrix.cpp)