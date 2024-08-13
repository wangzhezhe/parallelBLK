cmake_minimum_required(VERSION 3.15)
project(SimpleAdd LANGUAGES C CXX)

find_package(OpenCL REQUIRED)

add_executable(SimpleAdd SimpleAdd.cpp)
target_link_libraries(SimpleAdd OpenCL::OpenCL)
set(CL_HPP_TARGET_OPENCL_VERSION 220)