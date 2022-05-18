// test how to map the id into the thread for different cases
// try different block id and thred id in this case

#include <stdio.h>

__global__ void threadid1d() {
  // index
  printf("Hello from block %d %d %d, thread %d %d %d\n", blockIdx.x, blockIdx.y,
         blockIdx.z, threadIdx.x, threadIdx.y, threadIdx.z);
  // use the blockid plus the thread id
  // we can calculate the 1d thread
}

__global__ void threadid2d() {
  // index
  printf("Hello from block %d %d %d, thread %d %d %d\n", blockIdx.x, blockIdx.y,
         blockIdx.z, threadIdx.x, threadIdx.y, threadIdx.z);
}

// show the 1d 2d and 3d layout?
int main() {
  const int n = 80;
  int blocksize = 8;  // value usually chosen by tuning and hardware constraints
  int nblocks = n / blocksize;  // val
  // this is interpreted as 1d
  //threadid1d<<<nblocks, blocksize>>>();
  //cudaDeviceSynchronize();
  //the x*y*z in the dim3 structure should equals to the number of blocks or the blocksize
  threadid2d<<<dim3(8,10,1), dim3(2,2,2)>>>();
  cudaDeviceSynchronize();

  return 0;
}