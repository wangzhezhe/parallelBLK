
#define CL_HPP_TARGET_OPENCL_VERSION 220

#ifdef __APPLE__
#include <OpenCL/cl.hpp>
#else
#include <CL/opencl.hpp>
#endif

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

int main()
{
    // get all platforms (drivers), e.g. NVIDIA
    std::vector<cl::Platform> all_platforms;
    cl::Platform::get(&all_platforms);

    if (all_platforms.size() == 0)
    {
        std::cout << " No platforms found. Check OpenCL installation!\n";
        exit(1);
    }
    cl::Platform default_platform = all_platforms[0];
    std::cout << "Using platform: " << default_platform.getInfo<CL_PLATFORM_NAME>() << "\n";

    // get default device (CPUs, GPUs) of the default platform
    // TODO differences between plaform and device?
    std::vector<cl::Device> all_devices;
    default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
    if (all_devices.size() == 0)
    {
        std::cout << " No devices found. Check OpenCL installation!\n";
        exit(1);
    }
    std::cout << "Num of detected devices is " << all_devices.size() << std::endl;

    // // use device[1] because that's a GPU; device[0] is the CPU
    cl::Device default_device = all_devices[0];
    std::cout << "Using device: " << default_device.getInfo<CL_DEVICE_NAME>() << "\n";

    // a context is like a "runtime link" to the device and platform;
    // i.e. communication is possible
    cl::Context context({default_device});

    // create the program that we want to execute on the device
    cl::Program::Sources sources;

    // calculates for each element; C = A + B
    // compile the source code firstly, this is similar to the opengl
    std::string kernel_code =
        "   void kernel simple_add(global const int* A, global const int* B, global int* C, "
        "                          global const int* N) {"
        "       int ID, Nthreads, n, ratio, start, stop;"
        "       ID = get_global_id(0);"
        "       Nthreads = get_global_size(0);"
        "       n=N[0];" 
        "       for(int index=ID; index<n; index=index+Nthreads){"
        "           C[index]=A[index]+B[index];}"
        "   }";
    sources.push_back({kernel_code.c_str(), kernel_code.length()});

    cl::Program program(context, sources);
    if (program.build({default_device}) != CL_SUCCESS)
    {
        std::cout << "Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << std::endl;
        exit(1);
    }

    // apparently OpenCL only likes arrays ...
    // N holds the number of elements in the vectors we want to add
    cl_int N[1] = {100};
    const cl_int n = 100;

    // create buffers on device (allocate space on GPU)
    // this is just like the memory buffer
    cl::Buffer buffer_A(context, CL_MEM_READ_WRITE, sizeof(cl_int) * n);
    cl::Buffer buffer_B(context, CL_MEM_READ_WRITE, sizeof(cl_int) * n);
    cl::Buffer buffer_C(context, CL_MEM_READ_WRITE, sizeof(cl_int) * n);
    cl::Buffer buffer_N(context, CL_MEM_READ_ONLY, sizeof(cl_int));

    // create things on here (CPU), the data is located on CPU here
    cl_int A[n], B[n];
    for (cl_int i = 0; i < n; i++)
    {
        A[i] = i;
        B[i] = n - i - 1;
    }
    // create a queue (a queue of commands that the GPU will execute)
    // this commandQueue might be similar to the opengl pipeline
    // the difference is that in opengl, the pipeline is fixed
    // here, we have more flexibility to customize that pipeline
    cl::CommandQueue queue(context, default_device);

    // push write commands to queue
    // we first write these thress
    queue.enqueueWriteBuffer(buffer_A, CL_TRUE, 0, sizeof(int) * n, A);
    queue.enqueueWriteBuffer(buffer_B, CL_TRUE, 0, sizeof(int) * n, B);
    queue.enqueueWriteBuffer(buffer_N, CL_TRUE, 0, sizeof(int), N);

    // specifying the kernel, it seems that execution the kernel is the last step of the execution queue
    // adding kernel program, adding queue, adding NUll, what does te NDRange(10) specifies the number of threads we want to run
    cl::Kernel simple_add(program, "simple_add");


    //Attention! 
    //There is issue to get resut
    //if it fails to setting input parameter properly
    //for example, setting 3 parameter but the kernel uses 4 parameters will cause error
    simple_add.setArg(0, buffer_A);
    simple_add.setArg(1, buffer_B);
    simple_add.setArg(2, buffer_C);
    simple_add.setArg(3, buffer_N);

    // executing the functor
    queue.enqueueNDRangeKernel(simple_add, cl::NullRange, cl::NDRange(10), cl::NullRange);

    // another way is to use
    // cl::make_kernel<cl::Buffer&, cl::Buffer&, cl::Buffer&> simple_add(cl::Kernel(program, "simple_add"));
    // cl::EnqueueArgs eargs(queue, cl::NullRange, cl::NDRange(10), cl::NullRange);
    // simple_add(eargs, buffer_A, buffer_B, buffer_C).wait();
    queue.finish();

    cl_int C[n];
    // read result from GPU to here
    cl_int err = queue.enqueueReadBuffer(buffer_C, CL_TRUE, 0, sizeof(cl_int) * n, C);
    if (err!=CL_SUCCESS){
        std::cout << "Err in enqueueReadBuffer" << std::endl;
    }


    std::cout << "result: {";
    for (int i = 0; i < n; i++)
    {
        std::cout << C[i] << " ";
    }
    std::cout << "}" << std::endl;

    return 0;
}