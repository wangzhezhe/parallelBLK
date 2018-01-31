#include "genresult.cuh"
#include <sys/time.h>

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

__global__ void scan_kernal_opt(MatrixInfo *gmat, MatrixInfo *gvec, float gresult[], float gresultbuf[], int spanSize, int workPerwarp)
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

void matrixSortQuickOpt(MatrixInfo *mat, int left, int right)
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

int fillFlag(MatrixInfo *mat, int *flag)
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
void transformToNew(
    MatrixInfo *oldM, MatrixInfo *newM, MatrixInfo *oldvec, MatrixInfo *newvec, int *posvc, int *posmr, int *posorir, int *newvcount,int *compactlen)
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
        if (oldmrow == 449588)
        {
            printf("new row %d \n", newmrow);
        }
        newM->rIndex[i] = newmrow;
        posmr[newmrow] = oldmrow;
    }

    *newvcount = vcount;
    *compactlen = newmrow+1;
    free(vflag);
}

void transformBack(float *reorderresult, float *origresult, int *posorir, int compactlen)
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

int getMulDesign(MatrixInfo *mat, MatrixInfo *vec, MatrixInfo *finalres, int blockSize, int blockNum, int reorder)
{
    /*Allocate*/

    struct timespec start, end;
    struct timespec prestart, preend;
    struct timespec kerstart, kerend;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    clock_gettime(CLOCK_MONOTONIC_RAW, &prestart);
    /*Your own magic here!*/
    if (reorder == 1)
    {
        printf("**getMulDesign optimization based on first touch packing**\n");
    }

    int i, j;

    matrixSortQuickOpt(mat, 0, mat->nz - 1);

    printf("sort ok\n");

    //check sort
    //for (i = 0; i < 100; i++)
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

    reordervec->rIndex = NULL;

    reordervec->cIndex = NULL;

    reordervec->M = vec->M;
    reordervec->N = vec->N;
    reordervec->nz = vec->nz;

    int *posvc = (int *)malloc(sizeof(int) * mat->nz);
    int *posmr = (int *)malloc(sizeof(int) * mat->nz);
    int *posorir = (int *)calloc(sizeof(int), mat->nz);

    //add flag (last element of every line) line number
    int *lastFlag = (int *)calloc(((mat->nz) + 1), sizeof(int));
    int r = fillFlag(reordermat, lastFlag);

    if (r == -1)
    {
        printf("failed to fill flag\n");
    }

    //for (i = 0; i < reordermat->nz; i++)
    //{
    //   printf("index (%d) row (%d) colum (%d) val (%f) flag %d\n", i, reordermat->rIndex[i], reordermat->cIndex[i], reordermat->val[i], lastFlag[i]);
    //}

    //void transformToNew(MatrixInfo *oldM, MatrixInfo *newM, MatrixInfo *oldvec, MatrixInfo *newvec,int *posvc, int *posmr)

    int compactlen;
    int newvcount;
    if (reorder == 1)
    {
        transformToNew(mat, reordermat, vec, reordervec, posvc, posmr, posorir, &newvcount,&compactlen);
    }

    //check the reordering result
    /*
    printf("reordered mat:\n");
    for (i = 0; i < mat->nz; i++)
    {

        printf("r %d c %d v %f\n", reordermat->rIndex[i], reordermat->cIndex[i], reordermat->val[i]);
    }

    printf("vector:\n");
    for (i = 0; i < vec->nz; i++)
    {

        printf("r %d val %f\n", i, reordervec->val[i]);
    }

    printf("vector mapping\n");

    for (i = 0; i < compactlen; i++)
    {
        printf("newindex %d originalindex %d\n", i, posorir[i]);
    }

    float *origresult = (float *)malloc(sizeof(float) * 10);

    transformBack(reordervec->val, origresult, posorir, compactlen);

    printf("after transformBack\n");

    for (i = 0; i < compactlen; i++)
    {
        printf("index %d val %f\n", i, origresult[i]);
    }
*/

    int eleNum = mat->nz;

    //partition

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
    //void transformBack(float *orderresult, float *origresult, int *posorir, int compactlen)

    if (reorder == 1)
    {
        printf("**transfer back operation**\n");
        transformBack(sweepresult, origresult, posmr, compactlen);
        //transformBack(sweepresult, origresult, posmr, mat->nz);
    }
    for (i = 0; i < reordermat->M + 1; i++)
    {
        if (reorder == 1)
        {
            finalres->val[i] = origresult[i];
        }
        else
        {
            finalres->val[i] = sweepresult[i];
        }
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
