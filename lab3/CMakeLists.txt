cmake_minimum_required(VERSION 3.16.3)

project(lab3)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall -ggdb")

add_executable(lab3 src/main.cpp src/matrix.cpp src/calculation.cpp)
