#include <vector>
#include <iostream>
#include <sys/time.h>
#include "stdlib.h"

#include "utils.h"
#include "stdlib.h"
#include "cuda_error_check.cuh"
#include "initial_graph.cuh"
#include "parse_graph.cuh"

#define SSSP_INF 1073741824

int compareDest(const void *e1, const void *e2)
{
    const SdwEdge *edge1 = (SdwEdge *)e1;
    const SdwEdge *edge2 = (SdwEdge *)e2;

    if (edge1->dst <= edge2->dst)
    {
        return -1;
    }
    return 1;
}

__global__ void segscan_kernel_min(int spanSize, int vertexNum, int workPerwarp, int listlen, SdwEdge *gpuElist, unsigned int *gpuDisCur, unsigned int *gpuDisCurBuff, int *gpuFlag, int *gpuChange, int iterindex)
{

    int blockSize = blockDim.x;
    int threadidKernal = blockIdx.x * blockSize + threadIdx.x;
    //printf("warp num %d span size %d\n", warpNum,spanSize);

    int warpId = threadidKernal / spanSize;
    int laneId = threadidKernal % spanSize;

    //int WarpIdinBlock = threadidBlock / spanSize;

    //int beg = gbase * spanSize + threadidKernal;
    int beg = workPerwarp * warpId + laneId;
    int end = min(listlen, (warpId + 1) * workPerwarp);

    //int spanIndexInWarp = 0;
    int i;

    //copy 32 elements to share mem

    //int stride = 1;

    // shareMem is used for element of distance
    //__shared__ unsigned int sharTask[1024];
    //__shared__ unsigned int sharTaskBuff[1024];

    //int sharBegGlobal;
    //int sharMemIndex;

    int src, dst, weight, tempDist, tmpOld;

    for (i = beg; i < end; i += spanSize)
    {
        src = gpuElist[i].src;
        dst = gpuElist[i].dst;
        weight = gpuElist[i].weight;

        tempDist = gpuDisCur[src] + weight;

        //if(dst==1 && i<32){
        //    printf("index i %d src %d dst %d weight %d gpuDisCur[dst] %d\n",i,src,dst,weight,gpuDisCur[dst]);
        //}
        if (tempDist < gpuDisCur[dst])
        {
            tmpOld = gpuDisCur[dst];
            atomicMin(&gpuDisCur[dst], tempDist);

            if (tmpOld != gpuDisCur[dst])
            {
                atomicExch(gpuChange, 1);
                //printf("dst %d old %d new %d\n",dst,tmpOld,gpuDisCur[dst]);
            }
        }
    }

    //memcpy(gpuInput,gpuInputBuff, sizeof(int) * listlen);
}

void pullerSortByDst(std::vector<initial_vertex> *peeps, int blockSize, int blockNum)
{
    if (blockSize % 32 != 0)
    {
        printf("blockSize should be the multiple of 32\n");
        exit(1);
    }
    printf("start puller, sorted by dest\n");
    setTime();
    //Do all the things here!
    int i, j;
    int nbLen;
    //input parameter is a inverse adjacent list, transfer it into csv file
    int vertexNum = peeps->size();
    int edgeNum = 0;
    for (i = 0; i < vertexNum; i++)
    {
        nbLen = (*peeps)[i].nbrs.size();
        edgeNum = edgeNum + nbLen;
    }
    //printf("vertex num %d edge number %d\n", vertexNum, edgeNum);

    //std::vector<SdwEdge*> edgeList;
    SdwEdge *edgeList = (SdwEdge *)malloc(sizeof(SdwEdge) * edgeNum);

    if (edgeList == NULL)
    {
        printf("malloc fail");
        exit(1);
    }

    int edgeIndex = 0;
    for (i = 0; i < vertexNum; i++)
    {
        nbLen = (*peeps)[i].nbrs.size();
        for (j = 0; j < nbLen; j++)
        {
            edgeList[edgeIndex].dst = i;
            edgeList[edgeIndex].src = (*peeps)[i].nbrs[j].srcIndex;
            edgeList[edgeIndex].weight = (*peeps)[i].nbrs[j].edgeValue.weight;
            edgeIndex++;
        }
    }

    //sort
    qsort(edgeList, edgeNum, sizeof(SdwEdge), compareDest);

    //for (i = 0; i < 128; i++)
    //{

    //    printf("edge index %d src %d dst %d weight %d\n", i,edgeList[i].src, edgeList[i].dst, edgeList[i].weight);

    //printf("src (%d) dst (%d) wieght (%d) flag %d\n", edgeList[i].src, edgeList[i].dst, edgeList[i].weight, flagEdge[i]);
    //}

    int *flagEdge = (int *)malloc(sizeof(int) * edgeNum);

    //check after sorting
    for (i = 0; i < edgeNum - 1; i++)
    {
        if (edgeList[i].dst != edgeList[i + 1].dst)
        {

            flagEdge[i] = 1;
        }

        // printf("src (%d) dst (%d) wieght (%d) flag %d\n", edgeList[i].src, edgeList[i].dst, edgeList[i].weight, flagEdge[i]);
    }
    //last one
    if (edgeList[i].dst != edgeList[i - 1].dst)
    {
        flagEdge[i] = 1;
    }

    unsigned int *DisCur = (unsigned int *)malloc(sizeof(unsigned int) * edgeNum);
    unsigned int *DisCurBuff = (unsigned int *)malloc(sizeof(unsigned int) * edgeNum);

    //attention to the diff of element number here
    unsigned int *finalDis = (unsigned int *)malloc(sizeof(unsigned int) * vertexNum);

    DisCur[0] = 0;
    DisCurBuff[0] = 0;

    for (i = 1; i < edgeNum; i++)
    {
        DisCur[i] = SSSP_INF;
        DisCurBuff[i] = SSSP_INF;
    }

    DisCur[edgeList[0].dst] = 0;
    DisCurBuff[edgeList[0].dst] = 0;

    //check init dist

    //init the parameters on GPU

    SdwEdge *gpuElist;
    unsigned int *gpuDisCur;
    unsigned int *gpuDisCurBuff;
    int *gpuFlag;

    cudaMallocManaged((void **)&gpuElist, sizeof(SdwEdge) * edgeNum);
    cudaMemcpy(gpuElist, edgeList, sizeof(SdwEdge) * edgeNum, cudaMemcpyHostToDevice);
    //memcpy(gpuElist, edgeList, sizeof(SdwEdge) * edgeNum);

    cudaMallocManaged((void **)&gpuDisCur, sizeof(unsigned int) * vertexNum);
    //memcpy(gpuDisCur, DisCur, sizeof(unsigned int) * vertexNum);
    cudaMemcpy(gpuDisCur, DisCur, sizeof(unsigned int) * vertexNum, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&gpuDisCurBuff, sizeof(unsigned int) * vertexNum);
    //memcpy(gpuDisCur, DisCur, sizeof(unsigned int) * vertexNum);
    cudaMemcpy(gpuDisCurBuff, DisCur, sizeof(unsigned int) * vertexNum, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&gpuFlag, sizeof(int) * edgeNum);
    cudaMemcpy(gpuFlag, flagEdge, sizeof(int) * edgeNum, cudaMemcpyHostToDevice);

    //gpu task
    int *gpuChange;
    cudaMallocManaged((void **)&(gpuChange), sizeof(int));

    int iteIndex = 0;
    int change = 0;

    int spanSize = 32;
    int warpNum = blockNum * blockSize / spanSize;
    int workPerwarp;
    int reminder;
    int extra;
    int listlen = edgeNum;
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

    struct timespec gpuStart, gpuEnd;
    double gpuTotal = 0;

    printf("bn %d bs %d workPerwarp %d\n", blockNum, blockSize, workPerwarp);

    for (iteIndex = 0; iteIndex < edgeNum; iteIndex++)
    {
        //printf("iteration %d\n",iteIndex);

        change = 0;
        cudaMemcpy(gpuChange, &change, sizeof(int), cudaMemcpyHostToDevice);

        //TODO shareMem <check if shareMem is larger than limitation>

        //int stride = 1;
        //for (stride = 1; stride <= edgeNum; stride *= 2)
        //{
        //type could be represent by the index of node
        //    pulling_kernel_seg<<<blockNum, blockSize>>>(gpuElist, spanSize, vertexNum, workPerwarp, edgeNum, gpuDisCur, gpuDisCurBuff, gpuFlag, stride, gpuChange);
        //}

        //__global__ void segscan_kernel_min(int spanSize, int vertexNum, int workPerwarp, int listlen, SdwEdge *gpuElist, unsigned int *gpuDisCur, unsigned int *gpuDisCurBuff, int *gpuFlag, int *gpuChange)

        clock_gettime(CLOCK_MONOTONIC_RAW, &gpuStart);

        segscan_kernel_min<<<blockNum, blockSize>>>(spanSize, vertexNum, workPerwarp, edgeNum, gpuElist, gpuDisCur, gpuDisCurBuff, gpuFlag, gpuChange, iteIndex);

        clock_gettime(CLOCK_MONOTONIC_RAW, &gpuEnd);

        gpuTotal = gpuTotal + (double)1000 * (gpuEnd.tv_sec - gpuStart.tv_sec) + (double)(gpuEnd.tv_nsec - gpuStart.tv_nsec) / 1000000;

        //copy buffer to host

        //cudaMemcpy(DisCur, gpuDisCur, sizeof(unsigned int) * edgeNum, cudaMemcpyDeviceToHost);
        //filter the last element in segment

        //check change

        cudaMemcpy(&change, gpuChange, sizeof(int), cudaMemcpyDeviceToHost);

        if (change == 0 && iteIndex > 0)
        {

            printf("iteration no change, ineration num %d\n", iteIndex);
            break;
        }
    }

    cudaMemcpy(DisCur, gpuDisCur, sizeof(unsigned int) * vertexNum, cudaMemcpyDeviceToHost);

    for (i = 0; i < vertexNum; i++)
    {
        //printf("test value %d\n",newDisCur[i]);
        (*peeps)[i].vertexValue.distance = DisCur[i];
    }

    cudaFree(gpuElist);
    cudaFree(gpuDisCur);

    printf("The computation kernel time on GPU is %f milli-seconds\n", gpuTotal);
    std::cout << "The total computation time is " << getTime() << " milli-seconds.\n";
}