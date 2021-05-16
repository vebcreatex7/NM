cmake_minimum_required(VERSION 3.16.3)

project(lab2)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall")
add_executable(lab2 src/main.cpp src/calculation.cpp src/matrix.cpp)