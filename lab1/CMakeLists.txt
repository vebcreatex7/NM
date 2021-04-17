cmake_minimum_required(VERSION 3.16.3)

project(lab1)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "-Wall")
add_executable(lab1 src/matrix.cpp src/main.cpp src/calculation.cpp)

