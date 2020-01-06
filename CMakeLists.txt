cmake_minimum_required(VERSION 3.15)
project(mpi)
find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})

set(CMAKE_CXX_STANDARD 14)

add_executable(mpi main.cpp)