cmake_minimum_required(VERSION 3.15)
project(test LANGUAGES C CXX)

find_package(MPI REQUIRED) 
find_package(OpenMP)

add_executable(test test.cpp)
# if the program is c, use MPI_C
target_link_libraries(test MPI::MPI_CXX OpenMP::OpenMP_CXX)