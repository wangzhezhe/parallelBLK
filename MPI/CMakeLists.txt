cmake_minimum_required(VERSION 3.15)
project(ipgather LANGUAGES C CXX)

find_package(MPI REQUIRED) 

add_executable(ipgather ipgather.cpp)
target_link_libraries(ipgather MPI::MPI_CXX)

add_executable(mpi_isend mpi_isend.cpp)
target_link_libraries(mpi_isend MPI::MPI_CXX)