// matrix vector multiplecation
#include <chrono>
#include <cstdlib>
#include <iostream>
#include <vector>
using namespace std::chrono;

#define NUM_THREADS_PERBLOCK 128

// the macro to check the cudaAPI return code
#define cudaCheck(error)                                                      \
  if (error != cudaSuccess) {                                                 \
    printf("Fatal error: %s at %s:%d\n", cudaGetErrorString(error), __FILE__, \
           __LINE__);                                                         \
    exit(1);                                                                  \
  }

// the kernel operation executed by each thread
// the m and n is the size of input
__global__ void MVfunc(float* dMatrix, float* dVector, float* dAns, size_t m,
                       size_t n) {
  // for each thread
  // there is no memory colescing issue here
  size_t gindex = threadIdx.x + (blockIdx.x * blockDim.x);

  // it is necessary to know the location of the element this thread access
  // in the matrix
  size_t rowIndex = gindex / n;
  size_t columnIndx = gindex % n;

  // if (gindex == 0) {
  //  printf("device debug: %lf\n", dMatrix[1 * n + 1]);
  //}

  // printf(
  //    "threadIdx.x %d threadIdx.y %d threadIdx.z %d blockIdx.x %d blockIdx.y "
  //    "%d blockIdx.z %d gridDim.x %d "
  //    "gridDim.y %d gridDim.z %d\n",
  //    threadIdx.x, threadIdx.y, threadIdx.z, blockIdx.x, blockIdx.y,
  //    blockIdx.z, gridDim.x, gridDim.y, gridDim.z);

  // do nothing if the id is large than allocated one
  if (gindex >= m * n) {
    return;
  }

  // get particular data, do the multiplication
  // the index is not equals to the gindex
  // error in this way:
  // float a1 = dMatrix[gindex];
  // we need to map the thread id into the correct data position
  // one thread id is matched to one data position
  float a1 = dMatrix[gindex];
  float a2 = dVector[columnIndx];
  float v = a1 * a2;

  //printf("gindex %ld rowIndex %ld columnIndx %ld a1 %lf a2 %lf\n", gindex,
  //       rowIndex, columnIndx, a1, a2);
  // put it into the dAns vector
  // using atomicAdd to avoid the race condition
  atomicAdd(&dAns[rowIndex], v);
}

// do the memory copy and the data partition
void gpuMV(std::vector<float>& matrix, std::vector<float>& vector,
           std::vector<float>& ans, size_t m, size_t n) {
  float* dMatrix;
  float* dVector;
  float* dAns;

  // set the device id
  cudaSetDevice(0);

  // allocate memory on gpu
  cudaCheck(cudaMalloc((void**)&dMatrix, m * n * sizeof(float)));
  cudaCheck(cudaMalloc((void**)&dVector, n * sizeof(float)));
  cudaCheck(cudaMalloc((void**)&dAns, m * sizeof(float)));
  // set the memory to zero
  cudaCheck(cudaMemset(dAns, 0, m * sizeof(float)));

  // std::cout << "debug input  " << matrix[1 * n + 1] << std::endl;

  // init the data on the device
  // copy the data from the cpu into the device

  auto memcpy1start = high_resolution_clock::now();

  cudaCheck(cudaMemcpy(dMatrix, matrix.data(), m * n * sizeof(float),
                       cudaMemcpyHostToDevice));
  cudaCheck(cudaMemcpy(dVector, vector.data(), n * sizeof(float),
                       cudaMemcpyHostToDevice));
  auto memcpy1end = high_resolution_clock::now();
  auto memcpy1DurationMicro =
      duration_cast<microseconds>(memcpy1end - memcpy1start);

  std::cout << "memcpy to device time " << memcpy1DurationMicro.count()
            << std::endl;

  // do computation
  size_t elementsNum = m * n;

  auto kernelstart = high_resolution_clock::now();

  size_t blockSize =
      (elementsNum + NUM_THREADS_PERBLOCK - 1) / NUM_THREADS_PERBLOCK;

  // check blockSize, to make sure it is less then the limitation

  MVfunc<<<blockSize, NUM_THREADS_PERBLOCK>>>(dMatrix, dVector, dAns, m, n);

  auto kernelstop = high_resolution_clock::now();
  auto kernelDurationMicro =
      duration_cast<microseconds>(kernelstop - kernelstart);
  std::cout << "kernel time " << kernelDurationMicro.count() << std::endl;

  auto memcpy2start = high_resolution_clock::now();

  // copy results back, set it into the vector direactly
  cudaCheck(cudaMemcpy((float*)ans.data(), dAns, m * sizeof(float),
                       cudaMemcpyDeviceToHost));

  auto memcpy2end = high_resolution_clock::now();
  auto memcpy2Macro = duration_cast<microseconds>(memcpy2end - memcpy2start);

  std::cout << "memcpy to host time " << memcpy2Macro.count() << std::endl;
}

std::vector<float> cpuMV(std::vector<float>& matrix, std::vector<float>& vector,
                         size_t m, size_t n) {
  std::vector<float> ans;
  // the position of (0,0) is at the top left corner
  float tempv = 0;

  for (size_t i = 0; i < m; i++) {
    tempv = 0;
    for (size_t j = 0; j < n; j++) {
      tempv = tempv+(matrix[i * n + j] * vector[j]);
    }
    ans.push_back(tempv);
  }
  return ans;
}

int main(int argc, char** argv) {
  // int matrix size m*n
  size_t m = 100;
  size_t n = 200;

  // matrix (it is unnecessary to use the insert vector here)
  std::vector<float> matrix(m * n);

  // vector is n*1
  std::vector<float> vector(n);

  // init matrix and vector
  srand(static_cast<unsigned>(time(0)));

  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      // matrix[i * m + j] =
      //    static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
      matrix[i * n + j] = (i * 10 + j) * 0.1;
    }
  }

  // init vector
  for (size_t j = 0; j < n; j++) {
    // vector[j] = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
    vector[j] = j * 0.1;
  }

  // std::cout << "debug input 1  " << matrix[1 * n + 1] << std::endl;

  // ans is m*1
  std::vector<float> ansCPU;
  std::vector<float> ansGPU(m);

  auto cpustart = high_resolution_clock::now();
  ansCPU = cpuMV(matrix, vector, m, n);
  auto cpuend = high_resolution_clock::now();
  auto cputime = duration_cast<microseconds>(cpuend - cpustart);
  std::cout << "cputime time " << cputime.count() << std::endl;

  std::cout << "gputime time: " << std::endl;
  gpuMV(matrix, vector, ansGPU, m, n);

  // compare the difference
  // the max diff is around 0.06 in this case
  float epsilon = 0.1;
  for (int j = 0; j < m; j++) {
    float diff = fabs(ansCPU[j] - ansGPU[j]);
    if (diff > epsilon) {
      std::cout << "error " << j << " " << ansCPU[j] << ", " << ansGPU[j]
                << " diff " << diff << std::endl;
    }
  }

  return 0;
}