cmake_minimum_required(VERSION 3.14)
project(mpiclientserver LANGUAGES C CXX)
set(CMAKE_CXX_STANDARD 11)

find_package(MPI REQUIRED) 


add_executable(server server.c)
target_link_libraries(server MPI::MPI_C)

add_executable(client client.c)
target_link_libraries(client MPI::MPI_C)

add_executable(serverdor serverdor.c)
target_link_libraries(serverdor MPI::MPI_C)

add_executable(cliente cliente.c)
target_link_libraries(cliente MPI::MPI_C)

add_executable(manager manager.c)
target_link_libraries(manager MPI::MPI_C)

add_executable(worker worker.c)
target_link_libraries(worker MPI::MPI_C)