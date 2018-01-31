#include "genresult.cuh"
#include <sys/time.h>
#include "spmv_graph.cu"

/* Put your own kernel(s) here*/
/*input parameters:
mat: sparse matrix M in coo format
vec: input vector V in coo frmat
res: the multiplication of M*V
blocksize: means the number of threads per block, num_threads
blockNum: means the number of block per grid, num_blocks
*/

#define MAXTHREADNUMX 1024
#define MAXTHREADNUMY 1024
#define MAXGRIDNUMX 2147483647

__global__ void scan_kernal_opt_graph(MatrixInfo *gmat, MatrixInfo *gvec, float gresult[], float gresultbuf[], int spanSize, int workPerwarp)
{

    int blockSize = blockDim.x;
    int threadidKernal = blockIdx.x * blockSize + threadIdx.x;
    //int threadidBlock = threadIdx.x;

    int warpId = threadidKernal / spanSize;
    int laneId = threadidKernal % spanSize;

    //int WarpIdinBlock = threadidBlock / spanSize;

    int listlen = gmat->nz;
    //int beg = gbase * spanSize + threadidKernal;
    int beg = workPerwarp * warpId + laneId;
    int end = min(listlen, (warpId + 1) * workPerwarp);

    int spanIndexInWarp = 0;
    int i;

    //printf("threaidInKernel %d laneid %d\n", threadidKernal, laneId);

    //copy 32 elements to share mem

    int stride = 1;

    //__shared__ float sharTask[1024];
    //__shared__ float sharTaskBuff[1024];

    //int sharBegGlobal;
    //int sharMemIndex = laneId;
    int vecIndex;
    //int sharEnd;

    //printf("testmat mat index %d rnum %d cnum %d val %f\n", threadidKernal, gmat->rIndex[threadidKernal], gmat->cIndex[threadidKernal], gmat->val[threadidKernal]);

    for (i = beg, spanIndexInWarp = 0; i < end; i += spanSize, spanIndexInWarp++)
    {

        //sharBegGlobal = warpId * workPerwarp + spanIndexInWarp * spanSize;
        //sharMemIndex=WarpIdinBlock*spanSize+laneId;
        //sharEnd=sharBeg+spanSize;
        //sharMemIndex = warpId * spanSize + laneId;

        vecIndex = gmat->cIndex[i];
        //sharTask[sharMemIndex] = 1.0 * (gmat->val[i]) * (gvec->val[vecIndex]);
        gresult[i] = 1.0 * (gmat->val[i]) * (gvec->val[vecIndex]);
        gresultbuf[i] = gresult[i];
        //sharTaskBuff[sharMemIndex] = sharTask[sharMemIndex];

        for (stride = 1; stride <= spanSize; stride *= 2)
        {

            if (laneId >= stride && gmat->rIndex[i] == gmat->rIndex[i - stride])
            {

                //sharTask[sharMemIndex] = sharTaskBuff[sharMemIndex] + sharTaskBuff[sharMemIndex - stride];
                gresult[i] = gresultbuf[i] + gresultbuf[i - stride];
            }

            //atomicExch(&sharTaskBuff[sharMemIndex], sharTask[sharMemIndex]);
            //atomicExch(&gresultbuf[i], gresult[i]);
            gresultbuf[i] = gresult[i];
            //__syncthreads();
        }

        //put the element into the global mem
        //atomicExch(&gresult[i], sharTask[sharMemIndex]);
    }
}

void matrixSortQuickOptGraph(MatrixInfo *mat, int left, int right)
{
    int i = left, j = right;
    int pivot = mat->rIndex[(left + right) / 2];

    /* partition */
    while (i <= j)
    {
        while (mat->rIndex[i] < pivot)
            i++;
        while (mat->rIndex[j] > pivot)
            j--;
        if (i <= j)
        {
            matSwap(mat, i, j);
            i++;
            j--;
        }
    };

    /* recursion */
    if (left < j)
        matrixSortQuickOpt(mat, left, j);
    if (i < right)
        matrixSortQuickOpt(mat, i, right);
}

int fillFlagGraph(MatrixInfo *mat, int *flag)
{
    printf("row num %d\n", mat->M);
    int len = mat->nz;
    int i;
    int fCount = 0;
    for (i = 0; i < len - 1; i++)
    {
        if (mat->rIndex[i] != mat->rIndex[i + 1])
        {
            //printf("cur fcount %d\n",fCount);
            flag[i] = 1;
            fCount++;
        }
    }
    //flag of last element should always be 1
    flag[len - 1] = 1;

    return 0;
}

//assume the input data is sorted by row number
void transformToNewGraph(
    MatrixInfo *oldM, MatrixInfo *newM, MatrixInfo *oldvec, MatrixInfo *newvec, int *posvc, int *posmr, int *posorir, int *compactlen)
{
    int i;
    int nzlen = oldM->nz;
    int lenvec = newvec->M;
    int oldmrow, oldmcol, newmrow;
    int oldvrow;
    int *vflag = (int *)calloc(lenvec, sizeof(int));
    //init flag to label if the vector element changed the position
    //the initial element in newM and newV should same with the old one
    int vcount = 0;
    for (i = 0; i < nzlen; i++)
    {
        oldvrow = oldM->cIndex[i];

        //if current elem is changed previously, continue, else, change position
        if (vflag[oldvrow] != 1)
        {
            newvec->val[vcount] = oldvec->val[oldvrow];
            posvc[oldvrow] = vcount;
            posorir[vcount] = oldvrow;
            vflag[oldvrow] = 1;
            vcount++;
        }

        //update matrix column

        oldmcol = oldM->cIndex[i];
        oldmrow = oldM->rIndex[i];
        newM->cIndex[i] = posvc[oldmcol];

        //update matrix row
        if (i == 0)
        {
            newmrow = 0;
        }
        else
        {

            if (oldM->rIndex[i] == oldM->rIndex[i - 1])
            {
                newmrow = newM->rIndex[i - 1];
            }
            else
            {
                newmrow = newM->rIndex[i - 1] + 1;
            }
        }
        newM->rIndex[i] = newmrow;
        posmr[newmrow] = oldmrow;
    }

    *compactlen = (vcount - 1);
    free(vflag);
}

void transformBackGraph(float *reorderresult, float *origresult, int *posorir, int compactlen)
{
    int i;
    int origpos;
    for (i = 0; i < compactlen; i++)
    {

        origpos = posorir[i];
        if (reorderresult[i] != 0)
        {
            origresult[origpos] = reorderresult[i];
        }
    }
}

int getMulDesignByGraph(MatrixInfo *mat, MatrixInfo *vec, MatrixInfo *finalres, int blockSize, int blockNum)
{
    /*Allocate*/

    struct timespec start, end;
    struct timespec prestart, preend;
    struct timespec kerstart, kerend;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    clock_gettime(CLOCK_MONOTONIC_RAW, &prestart);
    /*Your own magic here!*/

    printf("**getMulDesign optimization for graph based algorithm**\n");

    int i, j;

    matrixSortQuickOpt(mat, 0, mat->nz - 1);

    printf("sort ok\n");

    //check sort
    //for (i = 0; i < mat->nz; i++)
    //{
    //    printf("index (%d) row (%d) colum (%d) value (%f)\n", i, mat->rIndex[i], mat->cIndex[i], mat->val[i]);
    //}

    //check flag

    //mat transform
    //copy the mat

    MatrixInfo *reordermat;
    MatrixInfo *reordervec;
    reordermat = (MatrixInfo *)malloc(sizeof(MatrixInfo));

    reordermat->rIndex = (int *)malloc(sizeof(int) * mat->nz);
    reordermat->cIndex = (int *)malloc(sizeof(int) * mat->nz);
    reordermat->val = (float *)malloc(sizeof(float) * mat->nz);

    printf("step0 \n");

    memcpy(reordermat->rIndex, mat->rIndex, sizeof(int) * mat->nz);
    memcpy(reordermat->cIndex, mat->cIndex, sizeof(int) * mat->nz);
    memcpy(reordermat->val, mat->val, sizeof(float) * mat->nz);

    reordermat->M = mat->M;
    reordermat->N = mat->N;
    reordermat->nz = mat->nz;

    reordervec = (MatrixInfo *)malloc(sizeof(MatrixInfo));
    //reordervec->rIndex = (int *)malloc(sizeof(int) * vec->nz);
    //reordervec->cIndex = (int *)malloc(sizeof(int) * vec->nz);
    reordervec->val = (float *)malloc(sizeof(float) * vec->nz);

    memcpy(reordervec->val, vec->val, sizeof(float) * vec->nz);

    printf("step1 \n");

    reordervec->rIndex = NULL;

    reordervec->cIndex = NULL;

    reordervec->M = vec->M;
    reordervec->N = vec->N;
    reordervec->nz = vec->nz;
    printf("step2 mat->nz %d\n", mat->nz);

    int *rowchange = (int *)malloc(sizeof(int) * mat->nz);
    //int *vaccess = (int *)malloc(sizeof(int) * mat->nz);
    int *reordervIndex = (int *)malloc(sizeof(int) * mat->M);
    int *oldvtonew = (int *)malloc(sizeof(int) * mat->M);
    //int *posmr = (int *)malloc(sizeof(int) * mat->nz);
    //int *posorir = (int *)calloc(sizeof(int), mat->nz);

    //struct DataItem **hashCindex = (DataItem **)malloc(mat->nz * sizeof(struct DataItem *));
    //(struct DataItem *)hashCindex[mat->nz];
    //return 1;
    //for (i = 0; i < mat->nz; i++)
    //{
    //    hashCindex[i] = NULL;
    //}

    //struct DataItem *graph[mat->nz];
    //struct DataItem **graph = (DataItem **)malloc(mat->nz * sizeof(struct DataItem *));
    //for (i = 0; i < mat->nz; i++)
    //{
    //    graph[i] = NULL;
    //}

    //init the graph
    vector<struct HeadNode *> graph;
    for (i = 0; i < mat->nz; i++)
    {

        vector<struct DataItem *> ev;
        HeadNode *node = new (HeadNode);
        node->edgeV = ev;
        node->maxweight = 0;
        node->flag = 0;
        graph.push_back(node);
    }

    int *lastFlag = (int *)calloc(((mat->nz) + 1), sizeof(int));
    int r = fillFlag(reordermat, lastFlag);

    if (r == -1)
    {
        printf("failed to fill flag\n");
    }

    printf("mattraverse begin\n");
    int newrownum;
    matTraverse(mat, reordermat, &newrownum, rowchange);

    printf("mattraverse ok\n");

    //printf("row chage\n");
    //for (i = 0; i < mat->nz; i++)
    //{
    //    printf("change new %d old %d\n", i, rowchange[i]);
    // }
    printf("newrownum %d\n", newrownum);
    //printf("vaccess\n");
    //for (i = 0; i < mat->nz; i++)
    //{
    //    printf("valindex %d\n", vaccess[i]);
    //}

    //check the reordering result
    int eleNum = mat->nz;

    //partition

    int spanSize = 32;
    //int spanSize = 4;
    int warpNum = (blockNum * blockSize) / spanSize;

    printf("warp num %d\n", warpNum);
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

    //construct map
    struct DataItem *maxEdge = (DataItem *)malloc(sizeof(DataItem));
    maxEdge->weight = 0;
    printf("addVbasedWarp begin\n");
    addVbasedWarp(mat->cIndex, graph, spanSize, workPerwarp, warpNum, maxEdge);

    //vdisplay(graph);
    printf("map construct ok\n");

    sortgraph(graph);
    //vdisplay(graph);
    printf("sort ok\n");

    //get the reordered column, record the transorm position (new to old)
    //remember the cache line number, find new after each cache line full

    vgetReorderedVindex(graph, reordervIndex, spanSize, mat->N);

    printf("check new v index\n");
    for (i = 0; i < mat->N; i++)
    {
        //printf("index %d vale %d\n",i,reordervIndex[i]);

        oldvtonew[reordervIndex[i]] = i;
        reordervec->val[i] = vec->val[reordervIndex[i]];
        //printf("oldvindex %d newvindex %d\n", reordervIndex[i], i);
    }

    //use transform position array to update the column index

    //for (i = 0; i < mat->nz; i++)
    //{
    //    printf("index %d vindex %d\n", i, oldvtonew);
    //}

    vmodifymatColumn(reordermat->cIndex, oldvtonew, mat->nz);


    //for (i = 0; i < reordermat->nz; i++)
    //{
    //    printf("index (%d) row (%d) colum (%d) val (%f) flag %d\n", i, reordermat->rIndex[i], reordermat->cIndex[i], reordermat->val[i], lastFlag[i]);
    //}

    //use reordered v to caculate

    //use transform for row change the result back

    //printf("map sort ok");

    //printf("graph\n");
    //display(graph, mat->nz);

    //printf("reversegraph\n");
    //display(reversegraph, mat->nz);

    //printf("maxedge\n");
    //printf("s %d e %d weight %d\n", maxEdge->start, maxEdge->end, maxEdge->weight);

    //printf("getReorderedVindex begin\n");

    //getReorderedVindex(graph, reversegraph, maxEdge, reordervIndex, mat->N);

    //printf("check reorder vector\n");

    //modify mat column
    //modifymatColumn(reordermat, hashCindex, oldvtonew);

    //check graph and reverse graph

    //for (i = 0; i < reordermat->nz; i++)
    //{
    //    printf("index (%d) row (%d) colum (%d) val (%f) flag %d\n", i, reordermat->rIndex[i], reordermat->cIndex[i], reordermat->val[i], lastFlag[i]);
    //}

    //void transformToNew(MatrixInfo *oldM, MatrixInfo *newM, MatrixInfo *oldvec, MatrixInfo *newvec,int *posvc, int *posmr)

    //int compactlen;

    //transformToNew(mat, reordermat, vec, reordervec, posvc, posmr, posorir, &compactlen);

    //init variable on gpu
    MatrixInfo *gpu_mat;
    MatrixInfo *gpu_vec;
    float *gpu_result;
    float *gpu_result_buf;
    //int *gpu_flag;
    cudaMallocManaged((void **)&gpu_mat, sizeof(MatrixInfo));

    cudaMallocManaged((void **)&gpu_vec, sizeof(MatrixInfo));

    cudaMallocManaged((void **)&(gpu_mat->rIndex), sizeof(int) * (mat->nz));
    //cudaMemcpy(gpu_mat->rIndex, mat->rIndex, sizeof(int) * mat->nz, cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_mat->rIndex, reordermat->rIndex, sizeof(int) * mat->nz, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&(gpu_mat->cIndex), sizeof(int) * (mat->nz));
    //cudaMemcpy(gpu_mat->cIndex, mat->cIndex, sizeof(int) * mat->nz, cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_mat->cIndex, reordermat->cIndex, sizeof(int) * mat->nz, cudaMemcpyHostToDevice);

    cudaMallocManaged((void **)&(gpu_mat->val), sizeof(float) * (mat->nz));
    //cudaMemcpy(gpu_mat->val, mat->val, sizeof(float) * mat->nz, cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_mat->val, reordermat->val, sizeof(float) * mat->nz, cudaMemcpyHostToDevice);

    //gpu_mat->M = mat->M;
    //gpu_mat->N = mat->N;
    //gpu_mat->nz = mat->nz;

    gpu_mat->M = reordermat->M;
    gpu_mat->N = reordermat->N;
    gpu_mat->nz = reordermat->nz;

    //vector assignment
    gpu_vec->rIndex = NULL;
    gpu_vec->cIndex = NULL;

    //only val is useful here, other value is null for vector
    cudaMallocManaged((void **)&(gpu_vec->val), sizeof(float) * (vec->nz));
    //cudaMemcpy(gpu_vec->val, vec->val, sizeof(float) * vec->nz, cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_vec->val, reordervec->val, sizeof(float) * vec->nz, cudaMemcpyHostToDevice);

    //result
    cudaMallocManaged((void **)&gpu_result, sizeof(float) * (mat->nz));
    cudaMallocManaged((void **)&gpu_result_buf, sizeof(float) * (mat->nz));
    //cudaMallocManaged((void **)&gpu_flag, sizeof(int) * (mat->nz));

    float *result = (float *)calloc(mat->nz, sizeof(float));

    //still need to initialise the gpu memory even if the memory is just allocated
    //cudaMemcpy(gpu_flag, lastFlag, mat->nz * sizeof(int), cudaMemcpyHostToDevice);
    //cudaMemcpy(gpu_result, result, mat->nz * sizeof(float), cudaMemcpyHostToDevice);
    //cudaMemcpy(gpu_result_buf, result, mat->nz * sizeof(float), cudaMemcpyHostToDevice);

    cudaMemcpy(gpu_result, result, reordermat->nz * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(gpu_result_buf, result, reordermat->nz * sizeof(float), cudaMemcpyHostToDevice);

    //kernel function
    printf("scan begin\n");

    clock_gettime(CLOCK_MONOTONIC_RAW, &preend);
    clock_gettime(CLOCK_MONOTONIC_RAW, &kerstart);

    scan_kernal_opt<<<blockNum, blockSize>>>(gpu_mat, gpu_vec, gpu_result, gpu_result_buf, spanSize, workPerwarp);

    printf("scan finishing\n");
    cudaMemcpy(result, gpu_result, mat->nz * sizeof(float), cudaMemcpyDeviceToHost);

    clock_gettime(CLOCK_MONOTONIC_RAW, &kerend);

    printf("start sweeping operation...\n");

    //printf("raw result\n");
    //for (i = 0; i < mat->nz + 1; i++)
    //{
    //    printf("check result index %d val %f\n", i, result[i]);
    //}

    //sweep operation (time consuming)
    for (i = -1; i < eleNum; i++)
    {
        //printf("copy value i (%d)\n",i);
        if (i == -1)
        {
            continue;
        }
        if ((i + 1) % workPerwarp % spanSize == 0 || (i + 1) % workPerwarp == 0)
        {

            //printf("index i %d flag %d\n", i, flag[i]);

            //if the flag is not 1, add to the next position with flag equal to 1
            if (lastFlag[i] != 1 && i > 0)
            {
                for (j = i; j < eleNum; j++)
                {
                    if (lastFlag[j] == 1)
                    {
                        //printf("index j %d flag %d\n", j, flag[j]);
                        result[j] = result[j] + result[i];
                        break;
                    }
                }
            }
        }
    }

    float *sweepresult = (float *)calloc(mat->M + 1, sizeof(float));

    float *origresult = (float *)calloc(mat->M + 1, sizeof(float));

    //printf("check elem in reordermat again\n");
    //for (i = 0; i < reordermat->nz; i++)
    //{
    //    printf("index (%d) row (%d) colum (%d) val (%f) flag %d\n", i, reordermat->rIndex[i], reordermat->cIndex[i], reordermat->val[i], lastFlag[i]);
    //}

    int finalCount;

    for (i = 0; i < reordermat->nz + 1; i++)
    {
        if (lastFlag[i] == 1)
        {
            //attention some row is all zero, skip those rows
            finalCount = reordermat->rIndex[i];
            sweepresult[finalCount] = result[i];
        }
    }

    printf("**transfer back operation**\n");
    // tranfer to originlal row and the compace row number
    transformBackGraph(sweepresult, origresult, rowchange, mat->N);

    for (i = 0; i < reordermat->M + 1; i++)
    {

        finalres->val[i] = origresult[i];
    }

    //cudaDeviceSynchronize(); // this line is needed for the completeness of the program
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    printf("Your total Time: %lu milli-seconds\n", 1000 * (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000);
    printf("Your preprocess Time: %lu milli-seconds\n", 1000 * (preend.tv_sec - prestart.tv_sec) + (preend.tv_nsec - prestart.tv_nsec) / 1000000);
    printf("Your kernel Time: %lu milli-seconds\n", 1000 * (kerend.tv_sec - kerstart.tv_sec) + (kerend.tv_nsec - kerstart.tv_nsec) / 1000000);

    free(result);
    free(lastFlag);

    //be careful for bus error if do the following operation first
    cudaFree(gpu_mat->rIndex);
    cudaFree(gpu_mat->cIndex);
    cudaFree(gpu_mat->val);

    cudaFree(gpu_mat);
    cudaFree(gpu_vec);

    //return 1 if everything is ok

    return 1;
}
