#include <vector>
#include <iostream>

#include "unistd.h"
#include "utils.h"
#include "cuda_error_check.cuh"
#include "initial_graph.cuh"
#include "parse_graph.cuh"

#define SSSP_INF 1073741824

int CompareSrc(const void *e1, const void *e2)
{
    const SdwEdge *edge1 = (SdwEdge *)e1;
    const SdwEdge *edge2 = (SdwEdge *)e2;

    if (edge1->src <= edge2->src)
    //if (edge1->dst <= edge2->dst)
    {
        return -1;
    }
    return 1;
}

/*
(gpuTPElist, gpuTplFlag, TPElen, gpuDisCur, gpuDisPrev, gpuChange, spanSize, iterIndex, blockNum, warpNum, workPerwarp);
*/

__global__ void neighborHandling_kernel(
    SdwEdge *gpuElist, int *gpuTplFlag, int TPElen, unsigned int *gpuDisCur, unsigned int *gpuDisPrev, int *gpuChange, int spanSize, int iterIndex, int blockNum, int warpNum, int workPerwarp)
{

    int blockSize = blockDim.x;
    int threadidKernal = blockIdx.x * blockSize + threadIdx.x;
    //printf("warp num %d span size %d\n", warpNum,spanSize);

    int warpId = threadidKernal / spanSize;
    int laneId = threadidKernal % spanSize;

    //int WarpIdinBlock = threadidBlock / spanSize;

    //int beg = gbase * spanSize + threadidKernal;
    int beg = workPerwarp * warpId + laneId;
    int end = min(TPElen, (warpId + 1) * workPerwarp);

    int i;
    int src, dst, weight;

    //printf("tid %d workload %d warpid %d beg %d end %d\n",threadidKernal,workPerwarp,warpId,beg,end);
    for (i = beg; i < end; i += spanSize)
    {
        src = gpuElist[i].src;
        dst = gpuElist[i].dst;
        weight = gpuElist[i].weight;

        atomicMin(&gpuDisCur[dst], gpuDisCur[src] + weight);

        if (gpuDisCur[dst] != gpuDisPrev[dst])
        {

            //printf("iteration %d index %d prev %d cur %d\n",iterIndex,i,gpuDisPrev[dst],gpuDisCur[dst]);
            atomicExch(&gpuDisPrev[dst], gpuDisCur[dst]);
            atomicExch(gpuChange, 1);
            gpuTplFlag[dst] = 1;
        }

    }
}

void pullerSortBySrcTPE(std::vector<initial_vertex> *peeps, int blockSize, int blockNum)
{

    printf("start tpe, pullerSortBySrcTPE, sort by src\n");

    if (blockSize % 32 != 0)
    {
        printf("blockSize should be the multiple of 32\n");
        exit(1);
    }

    struct timespec sortStart, sortEnd;
    clock_gettime(CLOCK_MONOTONIC_RAW, &sortStart);
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
    SdwEdge *TPEdge = (SdwEdge *)malloc(sizeof(SdwEdge) * edgeNum);

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
    qsort(edgeList, edgeNum, sizeof(SdwEdge), CompareSrc);

    //check after sorting
    /*
    for (i = 0; i < edgeNum; i++)
    {
        if (i < 64)
        {
            printf("index %d src (%d) dst (%d) wieght (%d)\n", i, edgeList[i].src, edgeList[i].dst, edgeList[i].weight);
        }
    }
    */
    unsigned int *DisCur = (unsigned int *)malloc(sizeof(unsigned int) * vertexNum);
    DisCur[0] = 0;
    for (i = 1; i < vertexNum; i++)
    {
        DisCur[i] = SSSP_INF;
    }

    DisCur[edgeList[0].src]=0;

    //init the parameters on GPU

    SdwEdge *gpuTPElist;
    unsigned int *gpuDisCur;
    unsigned int *gpuDisPrev;

    int *gpuTplFlag;
    int *tplFlag;
    tplFlag = (int *)malloc(sizeof(int) * vertexNum);

    cudaMallocManaged((void **)&gpuTPElist, sizeof(SdwEdge) * edgeNum);
    cudaMemcpy(gpuTPElist, edgeList, sizeof(SdwEdge) * edgeNum, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&gpuDisCur, sizeof(unsigned int) * vertexNum);
    cudaMemcpy(gpuDisCur, DisCur, sizeof(unsigned int) * vertexNum, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&gpuDisPrev, sizeof(unsigned int) * vertexNum);
    cudaMemcpy(gpuDisPrev, DisCur, sizeof(unsigned int) * vertexNum, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&gpuTplFlag, sizeof(int) * vertexNum);
    cudaMemcpy(gpuTplFlag, tplFlag, sizeof(int) * vertexNum, cudaMemcpyHostToDevice);

    //gpu task
    int *gpuChange;
    cudaMallocManaged((void **)&(gpuChange), sizeof(int));

    int iterIndex = 0;
    int change = 0;

    int spanSize = 32;
    int warpNum = blockNum * blockSize / spanSize;
    int workPerwarp;
    int reminder;
    int extra;
    int TPElen = edgeNum;
    int newEdgeLen;
    int tempSrc;
    int tempDst;

    struct timespec gpuStart, gpuEnd;
    struct timespec filterStart, filterEnd;

    double gpuTotal = 0;
    double filterTotal=0;

    clock_gettime(CLOCK_MONOTONIC_RAW, &sortEnd);
    double sortTime;

    sortTime = (double)1000 * (sortEnd.tv_sec - sortStart.tv_sec) + (double)(sortEnd.tv_nsec - sortStart.tv_nsec) / 1000000;

    printf("initial time before iteration is %f milli-seconds\n",sortTime);

    for (iterIndex = 0; iterIndex < edgeNum; iterIndex++)
    {

        change = 0;

        //printf("iteration index %d\n", iterIndex);
        cudaMemcpy(gpuChange, &change, sizeof(int), cudaMemcpyHostToDevice);

        for (i = 0; i < vertexNum; i++)
        {
            tplFlag[i] = 0;
        }

        cudaMemcpy(gpuTplFlag, tplFlag, sizeof(int) * vertexNum, cudaMemcpyHostToDevice);

        //trick gsharemem space array

        if (TPElen % warpNum == 0)
        {
            workPerwarp = TPElen / warpNum;
        }
        else
        {
            reminder = TPElen % warpNum;
            if (reminder % warpNum == 0)
            {
                extra = reminder / warpNum;
            }
            else
            {
                extra = (reminder / warpNum) + 1;
            }

            workPerwarp = extra + (TPElen / warpNum);
        }

        //printf("iteration %d new edge length %d workPerWarpNew %d\n", iterIndex, TPElen, workPerwarp);

        clock_gettime(CLOCK_MONOTONIC_RAW, &gpuStart);

        neighborHandling_kernel<<<blockNum, blockSize>>>(
            gpuTPElist, gpuTplFlag, TPElen, gpuDisCur, gpuDisPrev, gpuChange, spanSize, iterIndex, blockNum, warpNum, workPerwarp);

        clock_gettime(CLOCK_MONOTONIC_RAW, &gpuEnd);

        gpuTotal = gpuTotal + (double)1000 * (gpuEnd.tv_sec - gpuStart.tv_sec) + (double)(gpuEnd.tv_nsec - gpuStart.tv_nsec) / 1000000;

        //copy TPEflag to cpu

        cudaMemcpy(tplFlag, gpuTplFlag, sizeof(int) * vertexNum, cudaMemcpyDeviceToHost);

        //filter out the updated edge , TPE and the length of TPEdge

        /*
        for (i = 0; i < vertexNum; i++)
        {
            if (tplFlag[i] == 1)
          {
                printf("iteration %d change vertex %d\n", iterIndex, i);
            }
        }
        */

        clock_gettime(CLOCK_MONOTONIC_RAW, &filterStart);

        //update the TPSE

        newEdgeLen = 0;
        for (i = 0; i < edgeNum; i++)
        {
            tempSrc = edgeList[i].src;
            if (tplFlag[tempSrc] == 1)
            {
                memcpy(TPEdge + newEdgeLen, edgeList + i, sizeof(SdwEdge));
                newEdgeLen++;

                tempDst=edgeList[i].dst;

                if(tempDst > tempSrc){
                    tplFlag[tempDst] = 1;
                }
            }
        }

        //check newEdgeLen

        /*
        for (i = 0; i < newEdgeLen; i++)
        {
            printf("TPE index %d src %d dst %d weight %d\n", i, TPEdge[i].src, TPEdge[i].dst, TPEdge[i].weight);
        }
        */

        cudaMemcpy(gpuTPElist, TPEdge, sizeof(SdwEdge) * newEdgeLen, cudaMemcpyHostToDevice);

        TPElen = newEdgeLen;


        clock_gettime(CLOCK_MONOTONIC_RAW, &filterEnd);

        filterTotal = filterTotal + (double)1000 * (filterEnd.tv_sec - filterStart.tv_sec) + (double)(filterEnd.tv_nsec - filterStart.tv_nsec) / 1000000;



        cudaMemcpy(&change, gpuChange, sizeof(int), cudaMemcpyDeviceToHost);

        if (change == 0)
        {
            printf("iteration no change, ineration num %d\n",iterIndex);
            break;
        }
    }

    cudaMemcpy(DisCur, gpuDisCur, sizeof(unsigned int) * vertexNum, cudaMemcpyDeviceToHost);

    //update the vector
    for (i = 0; i < vertexNum; i++)
    {
        (*peeps)[i].vertexValue.distance = DisCur[i];
    }

    cudaFree(gpuTPElist);
    cudaFree(gpuDisCur);

    printf("The computation kernel time on GPU is %f milli-seconds\n",gpuTotal);
    printf("The filter Operation on CPU is %f milli-seconds\n",filterTotal);

    std::cout << "The total computation time is " << getTime() << " milli-seconds.\n";
}
