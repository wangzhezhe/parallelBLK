cmake_minimum_required(VERSION 3.15)
project(OCLHelloWorld LANGUAGES C CXX)

find_package(OpenCL REQUIRED)

add_executable(OCLHelloWorld OCLHelloWorld.cpp)
target_link_libraries(OCLHelloWorld OpenCL::OpenCL)
set(CL_HPP_TARGET_OPENCL_VERSION 220)