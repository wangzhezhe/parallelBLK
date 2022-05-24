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
  size_t globalthreadIndex = threadIdx.x + (blockIdx.x * blockDim.x);
  // elem num is more than thread number in all blocks

  size_t elemid = 0;
  size_t loopCount = 0;

  // it becomes to the naive case if there is enough numbers of block
  // and every thread access one data point
  // if (globalthreadIndex < 10) {
  //  printf("loopCount1 %ld", loopCount);
  //}

  extern __shared__ float s[];

  // each thread copy a chunk of mem from global into shared
  // this is equals to ceil operation
  int chunkSize = ((n + blockDim.x - 1) / blockDim.x);
  //if(globalthreadIndex==0){
  //  printf("debug chunkSize %d blockDim.x %d\n",chunkSize, blockDim.x);
  //}
  // each thread load at most two elements
  for (size_t partition = 0; partition < chunkSize; partition++) {
    size_t baseIndex = partition * blockDim.x;
    size_t elemIndex = baseIndex + threadIdx.x ;
    if (elemIndex < n) {
      s[elemIndex] = dVector[elemIndex];
      //printf("globalthreadIndex %ld elemIndex %ld s[elemIndex] %lf dVector[elemIndex] %lf\n",globalthreadIndex, elemIndex,s[elemIndex],dVector[elemIndex]);
    }
  }

  __syncthreads();  // wait for each thread to copy its elemenet

  size_t elemAllBlocks = gridDim.x * blockDim.x;

  while (true) {
    // if (globalthreadIndex < 10) {
    //  printf("loopCount2 %ld", loopCount);
    //}

    // add parathesis explicitly
    elemid = (loopCount * elemAllBlocks) + globalthreadIndex;
    // if (globalthreadIndex < 10) {
    //  printf(
    //      "globalthreadIndex %ld elemid %ld loopCount %ld gridDim.x %d blockDim.x "
    //      "%d\n",
    //      globalthreadIndex, elemid, loopCount, gridDim.x, blockDim.x);
    //}

    if (elemid >= m * n) {
      break;
    }

    size_t rowIndex = elemid / n;
    size_t columnIndx = elemid % n;

    // float a1 = dMatrix[rowIndex * n + columnIndx];
    // the dMatrix is accessed once for each thread
    // we do not need to put it into the shared memory here
    float a1 = dMatrix[elemid];
    //if(globalthreadIndex<100){
    //  printf("blockid %d globalthreadIndex %ld s[columnIndx] %lf dVector[columnIndx] %lf\n",blockIdx.x, globalthreadIndex,s[columnIndx],dVector[columnIndx]);
    //}
    float a2 = s[columnIndx];
    float v = a1 * a2;

    //TODO, the dAns can also be computed in the shared memory

    // printf("n gindex %ld rowindex %ld columnIndx %ld a1 %lf a2 %lf\n",
    // gindex,
    //       rowindex, columnIndx, a1, a2);
    // put it into the dAns vector
    // using atomicAdd to avoid the race condition
    atomicAdd(&dAns[rowIndex], v);

    loopCount++;
  }
  return;
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
  auto kernelstart = high_resolution_clock::now();

  // size_t blockSize =
  //     (elementsNum + NUM_THREADS_PERBLOCK - 1) / NUM_THREADS_PERBLOCK;

  size_t blockSize = 64;
  // the shared memory size is n*sizeof(float) that stores vector
  MVfunc<<<blockSize, NUM_THREADS_PERBLOCK, n * sizeof(float)>>>(
      dMatrix, dVector, dAns, m, n);

  // MVfunc<<<(elementsNum + NUM_THREADS_PERBLOCK - 1) / NUM_THREADS_PERBLOCK,
  //         NUM_THREADS_PERBLOCK>>>(dMatrix, dVector, dAns, m, n);

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
      tempv = tempv + (matrix[i * n + j] * vector[j]);
      // printf("i %ld j %ld m %f v %f tempv %f\n", i, j, matrix[i * n + j],
      //       vector[j], tempv);
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
  // the max diff is around 0.25 in this case
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