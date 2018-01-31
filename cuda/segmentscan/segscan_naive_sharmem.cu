#include <stdio.h>
#include "utils.h"

__global__ void segscan_kernel(int spanSize, int workPerwarp, int listlen, int *gpuInput, int *gpuInputBuff, int *gpuType)
{

    int blockSize = blockDim.x;
    int threadidKernal = blockIdx.x * blockSize + threadIdx.x;
    int threadidBlock = threadIdx.x;
    //printf("warp num %d span size %d\n", warpNum,spanSize);

    int warpId = threadidKernal / spanSize;
    int laneId = threadidKernal % spanSize;

    //int WarpIdinBlock = threadidBlock / spanSize;

    //int beg = gbase * spanSize + threadidKernal;
    int beg = workPerwarp * warpId + laneId;
    int end = min(listlen, (warpId+1)*workPerwarp);

    int spanIndexInWarp = 0;
    int i;

    //copy 32 elements to share mem

    int stride = 1;

    __shared__ int sharTask[1024];
    __shared__ int sharTaskBuff[1024];

    int sharBegGlobal;
    int sharMemIndex;
    int sharEnd;

    for (i = beg, spanIndexInWarp = 0; i < end; i += spanSize, spanIndexInWarp++)
    {

        sharBegGlobal = warpId * workPerwarp + spanIndexInWarp * spanSize;
        //sharMemIndex=WarpIdinBlock*spanSize+laneId;
        //sharEnd=sharBeg+spanSize;
        sharMemIndex = warpId * spanSize + laneId;

        sharTask[sharMemIndex] = gpuInput[i];
        sharTaskBuff[sharMemIndex] = gpuInput[i];

        printf("tid %d BegInGlobal %d sharMemIndex %d WarpId %d threadidKernal %d workPerwarp %d\n", i, sharBegGlobal, sharMemIndex, warpId, threadidKernal, workPerwarp);

        for (stride = 1; stride <= spanSize; stride *= 2)
        {
            //printf("global id %d warpId %d spanIndexInWarp %d taskValue %d flag %d type %d\n",
            //       i, warpId, spanIndexInWarp, gpuInput[i], gpuFlag[i], gpuType[i]);

            if (laneId >= stride && gpuType[i] == gpuType[i - stride] )
            {

                sharTask[sharMemIndex] = sharTaskBuff[sharMemIndex] + sharTaskBuff[sharMemIndex - stride];
            }



            __syncthreads();
            atomicExch(&sharTaskBuff[sharMemIndex], sharTask[sharMemIndex]);
        }

        //put the element into the global mem
        atomicExch(&gpuInputBuff[i], sharTask[sharMemIndex]);
    }

    //memcpy(gpuInput,gpuInputBuff, sizeof(int) * listlen);
}

int main(void)
{

    printf("seg scan test\n");
    //set up the configuration
    int blockNum = 1;
    int blockSize = 64;
    int eleNum = 131;

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

    printf("warpNum %d workPerwarp %d\n", warpNum, workPerwarp);
    //input vector

    int i, j;

    int *input = (int *)malloc(sizeof(int) * eleNum);
    int *flag = (int *)malloc(sizeof(int) * eleNum);
    int *type = (int *)malloc(sizeof(int) * eleNum);

    int segLen = 12;

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
    //last element should be 1
    flag[i-1] = 1;

    //gpu allocation
    int *gpuInput;
    int *gpuInputBuff;
    int *gpuFlag;
    int *gpuType;

    cudaMallocManaged((void **)&gpuInput, sizeof(int) * eleNum);
    cudaMemcpy(gpuInput, input, sizeof(int) * eleNum, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&gpuInputBuff, sizeof(int) * eleNum);
    cudaMemcpy(gpuInputBuff, input, sizeof(int) * eleNum, cudaMemcpyHostToDevice);

    //cudaMallocManaged((void **)&gpuFlag, sizeof(int) * eleNum);
    //cudaMemcpy(gpuFlag, flag, sizeof(int) * eleNum, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&gpuType, sizeof(int) * eleNum);
    cudaMemcpy(gpuType, type, sizeof(int) * eleNum, cudaMemcpyHostToDevice);

    printf("init configuration ok\n");

    int stride = 1;
    setTime();

    segscan_kernel<<<blockNum, blockSize>>>(spanSize, workPerwarp, eleNum, gpuInput, gpuInputBuff, gpuType);

    cudaMemcpy(input, gpuInputBuff, sizeof(int) * eleNum, cudaMemcpyDeviceToHost);

    //add to the last position of the segment

    //check every spanSize, check the end of every workWarp
    for (i = -1; i < eleNum; i++)
    {
        if (i == -1)
        {
            continue;
        }
        if ((i + 1) % workPerwarp % spanSize == 0 || (i + 1) % workPerwarp == 0 )
        {

            printf("index i %d flag %d\n", i, flag[i]);

            //if the flag is not 1, add to the next position with flag equal to 1
            if (flag[i] != 1 && i > 0)
            {
                for (j = i; j < eleNum; j++)
                {
                    if (flag[j] == 1)
                    {
                        printf("index j %d flag %d\n", j, flag[j]);
                        input[j] = input[j] + input[i];
                        break;
                    }
                }
            }
        }
    }

    for (i = 0; i < eleNum; i++)
    {
        printf("checkValue index %d value %d\n", i, input[i]);
        //check the flag and get the last element for the segmentation, this value is useful for some cases
    }

    //free operations
    int finishTime = getTime();
    printf("Took %d ms.\n", finishTime);

    cudaDeviceSynchronize();

    return 0;
}
