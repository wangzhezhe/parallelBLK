#include <vector>
#include <iostream>
#include "stdlib.h"

#include "utils.h"
#include "stdlib.h"
#include "cuda_error_check.cuh"
#include "initial_graph.cuh"
#include "parse_graph.cuh"

#define SSSP_INF 1073741824

int compare(const void *e1, const void *e2)
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

bool disChange(unsigned int *dp, unsigned int *dc, int len)
{
    int i;
    for (i = 0; i < len; i++)
    {
        if (dp[i] != dc[i])
        {
            return true;
        }
    }
    return false;
}

__global__ void pulling_kernel(SdwEdge *gpuElist, int listlen, unsigned int *gpuDisCur, int *gpuChange, int spanSize, int blockNum, int workPerwarp,int warpNum){

    int blockSize = blockDim.x;
    int threadidKernal = blockIdx.x * blockSize + threadIdx.x;
    
    //printf("warp num %d span size %d\n", warpNum,spanSize);

    int warpId= threadidKernal / spanSize;
    int laneId = threadidKernal % spanSize;

    //int beg = gbase * spanSize + threadidKernal;
    int beg=workPerwarp*warpId+laneId;
    int end = min(listlen,beg+workPerwarp);

    int src,dst,weight,tempDist,tmpOld;
    int i;
    
    //printf("tid %d workload %d warpid %d beg %d end %d\n",threadidKernal,workPerwarp,warpId,beg,end);
    for (i = beg; i < end; i += spanSize)
    {
         src = gpuElist[i].src;
         dst = gpuElist[i].dst;
         weight = gpuElist[i].weight;

         tempDist = gpuDisCur[src] + weight;

        //if(src==0 && i<32){
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
}

void pullerSortBySrc(std::vector<initial_vertex> *peeps, int blockSize, int blockNum)
{
    if (blockSize % 32 != 0)
    {
        printf("blockSize should be the multiple of 32\n");
        exit(1);
    }
    printf("start puller, sorted by src\n");
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
    qsort(edgeList, edgeNum, sizeof(SdwEdge), compare);

    //check after sorting
    //for (i = 0; i < edgeNum; i++)
    //{
    //    printf("src (%d) dst (%d) wieght (%d)\n", edgeList[i].src, edgeList[i].dst, edgeList[i].weight);
    //}

    unsigned int *DisCur = (unsigned int *)malloc(sizeof(unsigned int) * vertexNum);
    unsigned int *newDisCur = (unsigned int *)malloc(sizeof(unsigned int) * vertexNum);
    DisCur[0] = 0;
    for (i = 1; i < vertexNum; i++)
    {
        DisCur[i] = SSSP_INF;
    }

    DisCur[edgeList[0].src]=0;

    //check init dist

    //for (i = 0; i < vertexNum; i++)
    //{
    //    printf("index %d dist %d\n", i, finalDist[i]);
    //}

    //init the parameters on GPU

    SdwEdge *gpuElist;
    unsigned int *gpuDisCur;

    cudaMallocManaged((void **)&gpuElist, sizeof(SdwEdge) * edgeNum);
    cudaMemcpy(gpuElist, edgeList, sizeof(SdwEdge) * edgeNum, cudaMemcpyHostToDevice);
    //memcpy(gpuElist, edgeList, sizeof(SdwEdge) * edgeNum);

    cudaMallocManaged((void **)&gpuDisCur, sizeof(unsigned int) * vertexNum);
    //memcpy(gpuDisCur, DisCur, sizeof(unsigned int) * vertexNum);
    cudaMemcpy(gpuDisCur, DisCur, sizeof(unsigned int) * vertexNum, cudaMemcpyHostToDevice);

    //copy in every iteration

    //gpu task
    int *gpuBase;
    int *gpuChange;
    cudaMallocManaged((void **)&(gpuBase), sizeof(int));
    cudaMallocManaged((void **)&(gpuChange), sizeof(int));
    int iteIndex = 1;
    int change = 0;
    
    int spanSize = 32;
    int warpNum=blockNum*blockSize/spanSize;
    int workPerwarp;
    int reminder;
    int extra;
    int listlen=edgeNum;
    if (listlen % warpNum == 0)
    {
        workPerwarp = listlen / warpNum;
    }
    else
    {
        reminder=listlen % warpNum;
        if(reminder % warpNum==0){
            extra=reminder/warpNum;
        }else{
            extra=(reminder/warpNum)+1;
        }
        
        workPerwarp = extra+(listlen / warpNum);
    }

    struct timespec gpuStart, gpuEnd;
    double gpuTotal = 0;
    
    for (iteIndex = 0; iteIndex < edgeNum; iteIndex++)
    //for (iteIndex = 0; iteIndex < 1000; iteIndex++)
    {
        //for (iteIndex=0; iteIndex<2; iteIndex++){
        //printf("debug point 0\n");
        //gpu task
        //memset(distFlag,0,sizeof(int)*vertexNum);
        change = 0;

        //int *change=(int*)malloc(sizeof(int));
        //*change=0;
        //printf("iteration index %d\n", iteIndex);
        cudaMemcpy(gpuChange, &change, sizeof(int), cudaMemcpyHostToDevice);

        //trick gsharemem space array

        //1 size for edgeList  1*blockSize
        //2 size for distPrev  1*blockSize
        //3 size for distCur   1*blockSize
        //printf("debug point 1\n");
        //int sharSize = blockSize * (sizeof(SdwEdge) + sizeof(unsigned int) * 2);
        //int sharSize = blockSize * (sizeof(SdwEdge));
        
        //int sharSize = workPerwarp * (sizeof(SdwEdge));
        
        //TODO shareMem <check if shareMem is larger than limitation>
        int sharSize=0;

        if (blockNum * blockSize < spanSize)
        {
            spanSize = blockNum * blockSize;
        }

        //printf("iteration %d curr base %d\n",i,base);
        
        //pulling_kernel<<<blockNum, blockSize, sharSize>>>(gpuElist, edgeNum, gpuDisCur, base, gpuChange, spanSize,iteIndex);
        
        clock_gettime(CLOCK_MONOTONIC_RAW, &gpuStart);

        pulling_kernel<<<blockNum, blockSize, sharSize>>>(gpuElist, edgeNum, gpuDisCur, gpuChange, spanSize, blockNum, workPerwarp, warpNum);

        clock_gettime(CLOCK_MONOTONIC_RAW, &gpuEnd);

        gpuTotal = gpuTotal + (double)1000 * (gpuEnd.tv_sec - gpuStart.tv_sec) + (double)(gpuEnd.tv_nsec - gpuStart.tv_nsec) / 1000000;

        
        //pulling_kernel_testshar<<<blockNum, blockSize, sharSize>>>(gpuElist, edgeNum, gpuDisCur, base, gpuChange, spanSize);
        //printf("debug point 3\n");
        //get return value, check if change
        cudaMemcpy(&change, gpuChange, sizeof(int), cudaMemcpyDeviceToHost);

        //cudaDeviceSynchronize();

        //printf("iteration number %d gpuchange %d\n",iteIndex,change);

        //for(i=0;i<vertexNum;i++){
        //    printf("iterIndex %d change %d, final index %d dis %d\n",iteIndex,change,i,DisCur[i]);
        //}
        //bool ifchange=disChange(DisPre,DisCur,vertexNum);
        if (change == 0)
        {
           printf("iteration no change, ineration num %d\n",iteIndex);
           //cudaMemcpy(newDisCur, gpuDisCur, sizeof(unsigned int) * vertexNum, cudaMemcpyDeviceToHost);
           break;
        }
    }

    //check result

    //for(i=0;i<vertexNum;i++){
    //       printf("index %d dis %d\n",i,DisCur[i]);
    //}

    cudaMemcpy(newDisCur, gpuDisCur, sizeof(unsigned int) * vertexNum, cudaMemcpyDeviceToHost);
    //update the vector
    for (i = 0; i < vertexNum; i++)
    {
        //printf("test value %d\n",newDisCur[i]);
        (*peeps)[i].vertexValue.distance = newDisCur[i];
    }

    cudaFree(gpuElist);
    cudaFree(gpuDisCur);

    printf("The computation kernel time on GPU is %f milli-seconds\n",gpuTotal);
    std::cout << "The total computation time is " << getTime() << " milli-seconds.\n";
}