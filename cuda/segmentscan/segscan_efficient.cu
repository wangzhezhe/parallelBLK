#include <stdio.h>
#include "utils.h"

__global__ void segscan_kernel(int spanSize, int workPerwarp, int listlen, int *gpuInput, int *gpuInputBuff, int *gpuFlag, int *gpuType, int stride)
{

    int blockSize = blockDim.x;
    int threadidKernal = blockIdx.x * blockSize + threadIdx.x;

    //printf("warp num %d span size %d\n", warpNum,spanSize);

    int warpId = threadidKernal / spanSize;
    int laneId = threadidKernal % spanSize;

    //int beg = gbase * spanSize + threadidKernal;
    int beg = workPerwarp * warpId + laneId;
    int end = min(listlen, beg + workPerwarp);

    int spanIndexInWarp = 0;
    int i;

    for (i = beg, spanIndexInWarp = 0; i < end; i += spanSize, spanIndexInWarp++)
    {
        //printf("global id %d warpId %d spanIndexInWarp %d taskValue %d flag %d type %d\n",
        //       i, warpId, spanIndexInWarp, gpuInput[i], gpuFlag[i], gpuType[i]);

        //get value from buffer

        if (i >= stride && gpuType[i] == gpuType[i - stride])
        {
            gpuInput[i] = gpuInputBuff[i] + gpuInputBuff[i - stride];
        }
    }

    memcpy(gpuInputBuff, gpuInput, sizeof(int) * listlen);
}

int main(void)
{

    printf("seg scan test\n");
    //set up the configuration
    int blockNum = 2;
    int blockSize = 64;
    int eleNum = 1300;

    int spanSize = 32;
    int warpNum = blockNum * blockSize / spanSize;
    int workPerwarp;
    int reminder;
    int extra;
    int listlen = eleNum;
    if (listlen % warpNum == 0)
    {
        workPerwarp = listlen / warpNum;
    }
    else
    {
        reminder = listlen % warpNum;
        if (reminder % warpNum == 0)
        {
            extra = reminder / warpNum;
        }
        else
        {
            extra = (reminder / warpNum) + 1;
        }

        workPerwarp = extra + (listlen / warpNum);
    }

    //input vector

    int i, j;

    int *input = (int *)malloc(sizeof(int) * eleNum);
    int *flag = (int *)malloc(sizeof(int) * eleNum);
    int *type = (int *)malloc(sizeof(int) * eleNum);

    int segLen=11;

    for (i = 0, j = 1; i < eleNum; i++, j++)
    {
        input[i] = 1;
        if (j % segLen == 0)
        {
            //the last element for every segment is zero
            flag[i] = 1;
        }
        type[i] = i / segLen;
    }

    //gpu allocation
    int *gpuInput;
    int *gpuInputBuff;
    int *gpuFlag;
    int *gpuType;

    cudaMallocManaged((void **)&gpuInput, sizeof(int) * eleNum);
    cudaMemcpy(gpuInput, input, sizeof(int) * eleNum, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&gpuInputBuff, sizeof(int) * eleNum);
    cudaMemcpy(gpuInputBuff, input, sizeof(int) * eleNum, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&gpuFlag, sizeof(int) * eleNum);
    cudaMemcpy(gpuFlag, flag, sizeof(int) * eleNum, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&gpuType, sizeof(int) * eleNum);
    cudaMemcpy(gpuType, type, sizeof(int) * eleNum, cudaMemcpyHostToDevice);

    printf("init configuration ok\n");

    int stride = 1;
    setTime();
    for (stride = 1; stride <= eleNum; stride *= 2)
    {
        segscan_kernel<<<blockNum, blockSize>>>(spanSize, workPerwarp, eleNum, gpuInput, gpuInputBuff, gpuFlag, gpuType, stride);
    }


    cudaMemcpy(input, gpuInput, sizeof(int) * eleNum, cudaMemcpyDeviceToHost);

    for (i = 0; i < eleNum; i++)
    {
        printf("checkValue index %d value %d\n", i, input[i]);
        //check the flag and get the last element for the segmentation, this value is useful for some cases
    }

    //free operations
    int finishTime=getTime();
    printf("Took %d ms.\n",finishTime);

    cudaDeviceSynchronize();

    return 0;
}
