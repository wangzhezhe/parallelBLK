#include <stdio.h>

//variable that resides on the device
__device__ const char *STR = "HELLO WORLD!";
const char STR_LENGTH = 12;

__global__ void hello()
{
    //the char could be outputed by different threads 
	printf("%c\n", STR[threadIdx.x % STR_LENGTH]);
}

int main(void)
{
	int num_threads = STR_LENGTH;
	int num_blocks = 1;
	hello<<<num_blocks,num_threads>>>();
	cudaDeviceSynchronize();

	return 0;
}
