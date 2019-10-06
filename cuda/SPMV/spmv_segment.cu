#include "genresult.cuh"
#include <sys/time.h>

#define MAXTHREADNUMX 1024
#define MAXTHREADNUMY 1024
//the max number of integer block number is 2147483647 (upperlimitation of integer)
#define MAXGRIDNUMX 2147483647

#define SCANNUMBER 32
#define LOADROWNUM 1024

typedef struct
{
    int row;
    int index;
} Guard;

//scan_kernal<<<blockNum, blockSize>>>(gpu_mat, gpu_vec, gpu_result, gpu_flag, gpu_base);
//shared memory -> only one pointer could access
__global__ void scan_kernal(MatrixInfo *gmat, MatrixInfo *gvec, float gresult[], int gflag[], int const gpubase)
{
    //try to modify the gresult and gflag to see the correcness of return
    int gbase = gpubase;
    //printf("test scan ... scan kernal gbase %d\n",gbase);
    int blockSize = blockDim.x;
    //int threadGlobal = gbase + blockIdx.x * blockSize + threadIdx.x;
    //int threadidKernal = blockIdx.x * blockSize + threadIdx.x;
    int threadidBlock = threadIdx.x;
    //printf("curr gbase %d, curr blockId %d, curr threadidBlock %d, curr threadidKernal %d\n",gbase,blockIdx.x,threadidBlock,threadidKernal);
    //gresult[gbase]=1.0;
    //gflag[gbase]=1;

    //init the shared memory
    //from 0-blockSize-1 store original data from the global area
    //from blockSize-2*blockSize-1 the data after scan add
    /*
    //1 size for shareMem 2*blockSize
    //2 size for mat row  1*blockSize
    //3 size for mat val  1*blockSize
    //4 size for vec val  1*blockSize
    */
    extern __shared__ float shareMem[];
    //initiallise half of the sharMem
    int blockBase = blockIdx.x * blockSize;
    float *shareAdd = (float *)shareMem;
    memcpy(shareAdd, gmat->val + gbase + blockBase, sizeof(float) * blockSize);

    int *matrow = (int *)&shareAdd[2 * blockSize];
    memcpy(matrow, gmat->rIndex + gbase + blockBase, sizeof(int) * blockSize);

    int *matcol = (int *)&matrow[blockSize];
    memcpy(matcol, gmat->cIndex + gbase + blockBase, sizeof(int) * blockSize);

    float *matval = (float *)&matcol[blockSize];
    memcpy(matval, gmat->val + gbase + blockBase, sizeof(float) * blockSize);

    //value of vector is not continuous
    float *vecval = (float *)&matval[blockSize];
    //init in parallel way
    //printf("threadidBlock %d, vector %f\n",threadidBlock,gvec->val[matcol[threadidBlock]]);
    vecval[threadidBlock] = gvec->val[matcol[threadidBlock]];

    //if (threadGlobal == debugtid1 || threadGlobal == debugtid2 || matrow[threadidBlock]==debugrnum)
    //{
    //printf("gbase %d blockid %d tidblock %d, tidkernal %d , tidglobal %d,  sharadd %f, matrow %d, matcol %d matval %f, vecval %f\n",
    //         gbase, blockIdx.x, threadidBlock, threadidKernal, threadGlobal, shareAdd[threadidBlock], matrow[threadidBlock], matcol[threadidBlock], matval[threadidBlock], vecval[threadidBlock]);
    //}
    int span = 1;
    float product = 0;
    float addValue = 0;
    //compute
    product = 1.0 * shareAdd[threadidBlock] * vecval[threadidBlock];
    //shareAdd[threadidBlock] = product;
    atomicExch(&shareAdd[threadidBlock], product);
    __syncthreads();

    //if (threadGlobal == debugtid1 || threadGlobal == debugtid2 || matrow[threadidBlock]==debugrnum)
    //{
    //    printf("debug step1 id %d threadidBlock %f \n", threadGlobal, shareAdd[threadidBlock]);
    //}

    for (span = 1; span < blockSize; span *= 2)
    {
        //add operation
        //second part of the if condifiton is used to make ture the row are in same row
        if (threadidBlock >= span && matrow[threadidBlock] == matrow[threadidBlock - span])
        {
            addValue = shareAdd[threadidBlock] + shareAdd[threadidBlock - span];
            atomicExch(&shareAdd[blockSize + threadidBlock], addValue);

            //if (threadGlobal == debugtid1 || threadGlobal == debugtid2 || matrow[threadidBlock]==debugrnum)
            //{
            //    printf("debug step2 id %d threadidBlock %f \n", threadGlobal, shareAdd[threadidBlock]);
            //}
        }
        else
        {
            //same with original one, do not add
            atomicExch(&shareAdd[blockSize + threadidBlock], shareAdd[threadidBlock]);
        }

        //if (threadGlobal == debugtid1 || threadGlobal == debugtid2 || matrow[threadidBlock]==debugrnum)
        //{
        //    printf("debug step3 id %d threadidBlock %f \n", threadGlobal, shareAdd[threadidBlock]);
        //}

        __syncthreads();

        shareAdd[threadidBlock] = shareAdd[blockSize + threadidBlock];

        //if (threadGlobal == debugtid1 || threadGlobal == debugtid2 || matrow[threadidBlock]==debugrnum)
        //{
        //    printf("debug step4 id %d threadidBlock %f \n", threadGlobal, shareAdd[threadidBlock]);
        //}

        __syncthreads();
    }

    //computed the final result, label at the end position of the block
    //check curr shareADD
    //if (threadGlobal == debugtid1 || threadGlobal == debugtid2 || matrow[threadidBlock]==debugrnum)
    //{
    //    printf("debugbeforegresult %f globalid %d\n", shareAdd[threadidBlock], threadGlobal);
    //}
    atomicExch(&gresult[gbase + blockBase + threadidBlock], shareAdd[threadidBlock]);
    atomicExch(&gflag[gbase + blockBase + blockSize - 1], 1);
    //if (threadGlobal == debugtid1 || threadGlobal == debugtid2 || matrow[threadidBlock]==debugrnum)
    //{
    //    printf("debug globalid %d gbase %d gflag label %d gflag value %d\n", threadGlobal, gbase, gbase + blockBase + blockSize - 1, gflag[gbase + blockBase + blockSize - 1]);
    //}
}

// refer to the csr version in paper
//Efficient Sparse Matrix-Vector Multiplication on CUDA
__global__ void
/*
const int num_rows, //number of rows
const int *ptr, //ptr,indices,data represent the matrix stored by csr format
const int *indices, 
const float *data,
const float *x, //x represent the vetor to take the multiplication operation
float *y //y represent the final value after multiplication operation
*/
spmv_csr_vector_kernel(const int num_rows,
                       const int *ptr,
                       const int *indices,
                       const float *data,
                       const float *x,
                       float *y)
{
    __shared__ float vals[];
    // global thread index
    // both blocknumber and blocksize determine the total threads that run in paralllel
    int thread_id = blockDim.x * blockIdx.x + threadIdx.x; 
    // global warp index
    int warp_id = thread_id / 32;                          
    // thread index within the warp
    int lane = thread_id & (32 - 1);                      
    
    // one warp per row
    int row = warp_id;
    if (row < num_rows)
    {
        //start and end could make sure the start and the end data for this array
        int row_start = ptr[row];
        int row_end = ptr[row + 1];
        // compute running sum per thread

        vals[threadIdx.x] = 0;

        //attention, the row_end could make sure 
        //do multiplication for every element in this row
        for (int jj = row_start + lane; jj < row_end; jj += 32){
            //why this is not lane?   shared memroy 0-31 belong to warp1 row1 32-64 belong to warp2 row2
            //different rows partition shared memory into different rows
            vals[threadIdx.x] += data[jj] * x[indices[jj]];
        }


        // parallel reduction in shared memory
        // one warp has 32 thread
        // only consider first part of thread in warp
        //every thread in a warp will do the same operation every time +1 +2 +4...
        //every time of calculation
        //one specific row will be calculated out
        //it's a little bit tricky for this part of code
        //because of the if divergency, the thread with id > 16 will not excute the following code
        //we assuem that row number is large
        //every iteration , the value in current warp will be added into the y[row]
        //thread with a small id will excute multiple steps
        //for lane16, the thread with lane > 16 will wait here, same idea for the following parts
        //for lane = 0, every if condition will be excuted
        if (lane < 16)
            vals[threadIdx.x] += vals[threadIdx.x + 16];
        if (lane < 8)
            vals[threadIdx.x] += vals[threadIdx.x + 8];
        if (lane < 4)
            vals[threadIdx.x] += vals[threadIdx.x + 4];
        if (lane < 2)
            vals[threadIdx.x] += vals[threadIdx.x + 2];
        if (lane < 1)
            vals[threadIdx.x] += vals[threadIdx.x + 1];
        // first thread writes the result
        if (lane == 0)
            y[row] += vals[threadIdx.x];
    }
}

void matSwap(MatrixInfo *mat, int i, int j)
{
    int tempr;
    int tempc;
    float tempv;

    tempr = mat->rIndex[j];
    mat->rIndex[j] = mat->rIndex[i];
    mat->rIndex[i] = tempr;
    //change cIndex
    tempc = mat->cIndex[j];
    mat->cIndex[j] = mat->cIndex[i];
    mat->cIndex[i] = tempc;
    //change val
    tempv = mat->val[j];
    mat->val[j] = mat->val[i];
    mat->val[i] = tempv;
}

int partition(MatrixInfo *mat, int low, int high)
{
    int pivot = mat->rIndex[high]; // pivot
    int i = (low - 1);             // Index of smaller element

    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (mat->rIndex[j] <= pivot)
        {
            i++; // increment index of smaller element
            matSwap(mat, i, j);
        }
    }
    matSwap(mat, i + 1, high);
    return (i + 1);
}

void matrixSortQuick(MatrixInfo *mat, int low, int high)
{

    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
           at right place */
        int pi = partition(mat, low, high);
        printf("qsort %d\n", pi);
        // Separately sort elements before
        // partition and after partition
        matrixSortQuick(mat, low, pi - 1);
        matrixSortQuick(mat, pi + 1, high);
    }
}
void matrixSortQuick2(MatrixInfo *mat, int left, int right)
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
        matrixSortQuick2(mat, left, j);
    if (i < right)
        matrixSortQuick2(mat, i, right);
}

void matrixSort(MatrixInfo *mat)
{
    int elemNum = mat->nz;

    //sort according the mat->M
    int i, j;
    int tempr;
    int tempc;
    float tempv;
    for (i = 0; i < elemNum; i++)
    {
        for (j = i + 1; j < elemNum; j++)
        {
            if (mat->rIndex[j] < mat->rIndex[i])
            {
                tempr = mat->rIndex[j];
                mat->rIndex[j] = mat->rIndex[i];
                mat->rIndex[i] = tempr;
                //change cIndex
                tempc = mat->cIndex[j];
                mat->cIndex[j] = mat->cIndex[i];
                mat->cIndex[i] = tempc;
                //change val
                tempv = mat->val[j];
                mat->val[j] = mat->val[i];
                mat->val[i] = tempv;
            }
        }
    }
    return;
}

int getGuardRecord(MatrixInfo *mat, Guard **gRecord)
{
    int elem_num = mat->nz;
    int i;
    int gindex = 0;
    for (i = 1; i < elem_num; i++)
    {
        if (mat->rIndex[i] != mat->rIndex[i - 1])
        {

            (*gRecord)[gindex].index = i - 1;
            (*gRecord)[gindex].row = mat->rIndex[i - 1];
            //printf("ri %d ri-1 %d guard %d i-1 %d gindex %d\n",mat->rIndex[i],mat->rIndex[i-1],(*gRecord)[gindex].index,i-1,gindex);
            gindex++;
        }
    }
    //last one
    (*gRecord)[gindex].index = elem_num - 1;
    (*gRecord)[gindex].row = mat->rIndex[elem_num - 1];
    return (gindex + 1);
}

//glen is the length of the guard
//return cloest guard position
int getCloseGuardIndex(int currIndex, Guard *guard, int glen)
{

    int i = 0;
    int r = -2;
    for (i = 0; i < glen - 1; i++)
    {
        if (currIndex == guard[i].index || currIndex == guard[i + 1].index)
        {
            return -1;
        }
        else
        {
            if (currIndex > guard[i].index && currIndex < guard[i + 1].index)
            {
                return guard[i + 1].index;
            }
        }
    }

    printf("debugclose currIndex %d\n", currIndex);
    return r;
}

//make sure if the index equal to the element in guard
//range the element in guard to make sure if the index is in the array
int checkGuard(int index, int *guard, int len)
{

    int i = 0;
    int *temp = guard;

    for (i = 0; i < len; i++)
    {
        if (index == *temp)
        {
            return 1;
        }
    }
    return 0;
}

/*input parameters:
mat: sparse matrix M in coo format
vec: input vector V in coo frmat
res: the multiplication of M*V

*/
int getMulScan(MatrixInfo *mat, MatrixInfo *vec, MatrixInfo *finalres, int blockSize, int blockNum)
{

    struct timespec prestart, preend;

    clock_gettime(CLOCK_MONOTONIC_RAW, &prestart);

    //check parameters
    printf("mulScan matrix m n nz %d %d %d\n", mat->M, mat->N, mat->nz);
    printf("mulScan vector m n nz %d %d %d\n", vec->M, vec->N, vec->nz);

    int errCode = 1;
    if (mat == NULL || vec == NULL || finalres == NULL)
    {
        errCode = 2;
        printf("getMulScan err, mat , vec, res could not be null\n");
        return errCode;
    }

    if (blockSize > (MAXTHREADNUMX * MAXTHREADNUMX))
    {
        errCode = 2;
        printf("parameters err, blockSize (%d) is larger than limitation\n", blockSize);
        return errCode;
    }

    //if (mat->nz > LOADROWNUM)
    //{
    //TODO divide logic into different parts for large size matrix
    //    errCode = 2;
    //    printf("matrix size (%d) is larger than upper limitation\n", mat->nz);
    //    return errCode;
    // }

    printf("non-zero number %d\n", mat->nz);

    //original mat value
    int i;
    for (i = 0; i < mat->nz; i++)
    {
        //printf("indexoriginal %d r %d c %d val %f\n",i, mat->rIndex[i], mat->cIndex[i], mat->val[i]);
        if (mat->rIndex[i] < 0 || mat->cIndex[i] < 0)
        {
            printf("matrix load fail: \nindexoriginal %d r %d c %d val %f\n", i, mat->rIndex[i], mat->cIndex[i], mat->val[i]);
            exit(1);
        }
    }

    /*
    sorting by column
    */

    //matrixSort(mat);
    //matrixSortQuick(mat,0,mat->nz-1);
    matrixSortQuick2(mat, 0, mat->nz - 1);

    /*
    summary the matrix, get guardRecord, length= mat->M*sizeof(int)
    */

    Guard *guard_record = (Guard *)malloc((mat->M + 1) * sizeof(Guard));

    int glength = getGuardRecord(mat, &guard_record);

    //check guard

    //for (i = 0; i < glength; i++)
    //{
    //    printf("index %d guardPosition %d guard rowid %d\n", i, guard_record[i].index,guard_record[i].row);
    //}

    //check sort

    //for (i = 0; i < mat->nz; i++)
    //{
    //    printf("index %d r %d c %d val %f\n", i, mat->rIndex[i], mat->cIndex[i], mat->val[i]);
    //}

    printf("Sorting by row number finished\n");

    /*Allocate things...*/

    //malloc innter pointer value
    //refer https://github.com/parallel-forall/code-samples/blob/master/posts/unified-memory/dataElem_um.cu
    // void** is really important here to prevent the error like all cuda device is busy

    MatrixInfo *gpu_mat;
    MatrixInfo *gpu_vec;
    float *gpu_result;
    int *gpu_flag;
    int *gpu_base;

    cudaMallocManaged((void **)&gpu_mat, sizeof(MatrixInfo));
    cudaMallocManaged((void **)&gpu_vec, sizeof(MatrixInfo));

    cudaMallocManaged((void **)&(gpu_mat->rIndex), sizeof(int) * (mat->nz));
    memcpy(gpu_mat->rIndex, mat->rIndex, sizeof(int) * mat->nz);

    cudaMallocManaged((void **)&(gpu_mat->cIndex), sizeof(int) * (mat->nz));
    memcpy(gpu_mat->cIndex, mat->cIndex, sizeof(int) * mat->nz);

    cudaMallocManaged((void **)&(gpu_mat->val), sizeof(float) * (mat->nz));
    memcpy(gpu_mat->val, mat->val, sizeof(float) * mat->nz);

    gpu_mat->M = mat->M;
    gpu_mat->N = mat->N;
    gpu_mat->nz = mat->nz;

    //vector assignment
    gpu_vec->rIndex = NULL;
    gpu_vec->cIndex = NULL;

    //only val is useful here, other value is null for vector
    cudaMallocManaged((void **)&(gpu_vec->val), sizeof(float) * (vec->nz));
    memcpy(gpu_vec->val, vec->val, sizeof(float) * vec->nz);

    //result

    cudaMallocManaged((void **)&gpu_result, sizeof(float) * (mat->nz));
    cudaMallocManaged((void **)&gpu_flag, sizeof(int) * (mat->nz));

    int *flag = (int *)calloc(mat->nz, sizeof(int));
    float *result = (float *)calloc(mat->nz, sizeof(float));

    //still need to initialise the gpu memory even if the memory is just allocated
    memcpy(gpu_flag, flag, mat->nz * sizeof(int));
    memcpy(gpu_result, result, mat->nz * sizeof(float));

    //base
    cudaMallocManaged((void **)&(gpu_base), sizeof(int));

    /* Start Execution */

    printf("allocation things on gpu ok\n");
    printf("block size (%d),block number (%d)\n", blockSize, blockNum);

    clock_gettime(CLOCK_MONOTONIC_RAW, &preend);

    struct timespec start, end;

    clock_gettime(CLOCK_MONOTONIC_RAW, &start);

    int all_elem_num = mat->nz;
    int workload = blockNum * blockSize;
    int iteration = all_elem_num / workload;
    int reminder_load = all_elem_num % workload;

    printf("workload %d iteration %d reminder %d\n", workload, iteration, reminder_load);

    //int *base = (int *)malloc(sizeof(int));
    //base[0] = 0;
    int base = 0;
    //trick gsharemem space array

    //1 size for shareMem 2*blockSize
    //2 size for mat row  1*blockSize
    //3 size for mat val  1*blockSize
    //4 size for vec val  1*mat->M

    int sharSize = blockSize * (3 * sizeof(float) + 2 * sizeof(int)) + blockSize * (sizeof(float));
    //int rowNum=mat->M;
    i = 0;
    for (i = 0; i < iteration; i++)
    {
        base = workload * i;
        //cudaMemcpy(gpu_base, base, sizeof(int), cudaMemcpyHostToDevice);
        //printf("input base iteration %d\n",base);
        scan_kernal<<<blockNum, blockSize, sharSize>>>(gpu_mat, gpu_vec, gpu_result, gpu_flag, base);
        //cudaDeviceSynchronize();
    }
    if (reminder_load > 0)
    {
        int newblocknum = reminder_load / blockSize;
        int newreminder = reminder_load % blockSize;
        printf("newblocknum %d newreminder %d\n", newblocknum, newreminder);

        base = workload * i;

        if (newblocknum > 0)
        {
            //printf("input base iteration %d\n",base);
            scan_kernal<<<newblocknum, blockSize, sharSize>>>(gpu_mat, gpu_vec, gpu_result, gpu_flag, base);
        }
        if (newreminder > 0)
        {
            //less than one block
            int newbase = base + newblocknum * blockSize;
            int reminderSize = newreminder * (3 * sizeof(float) + 2 * sizeof(int)) + newreminder * (sizeof(float));
            printf("input reminder iteration %d\n", newbase);
            scan_kernal<<<1, newreminder, reminderSize>>>(gpu_mat, gpu_vec, gpu_result, gpu_flag, newbase);
        }

        //gpu_base = base;
        //cudaMemcpy(gpu_base, base, sizeof(int), cudaMemcpyHostToDevice);
    }

    //sum up from guard, flag, result
    //copy result and flag from gpu
    cudaMemcpy(flag, gpu_flag, mat->nz * sizeof(int), cudaMemcpyDeviceToHost);
    cudaMemcpy(result, gpu_result, mat->nz * sizeof(float), cudaMemcpyDeviceToHost);

    //check flag, result, ok here, the result could be write back into the cpu
    //for (i = 0; i < mat->nz; i++)
    //{
    //    printf("index %d flag %d  raw result %f\n", i, flag[i], result[i]);
    //}

    //TODO improve the parallelism ability here

    int closeIndex = 0;
    for (i = 0; i < all_elem_num; i++)
    {
        if (flag[i] == 1)
        {

            //closeIndex = getCloseGuardIndex(i, guard_record, mat->M);
            closeIndex = getCloseGuardIndex(i, guard_record, glength);
            //printf("debug curr index %d closeIndex %d\n",i,closeIndex);
            if (closeIndex >= 0)
            {
                result[closeIndex] = result[closeIndex] + result[i];
            }
            if (closeIndex == -2)
            {
                printf("bug in get cloestguard\n");
                exit(1);
            }
        }
    }

    Guard *temp = guard_record;
    int nzindex = 0;
    int finalindex = 0;
    int padnum = 0;
    for (i = 0; i < all_elem_num; i++)
    {
        if (i == temp[nzindex].index)
        {
            if (nzindex == 0)
            {
                if (temp[nzindex].row > 0)
                {
                    padnum = temp[nzindex].row;
                    while (padnum > 0)
                    {
                        finalres->val[finalindex] = 0;
                        finalindex++;
                        padnum--;
                    }
                }
            }
            else if (nzindex >= 1)
            {
                padnum = temp[nzindex].row - temp[nzindex - 1].row;
                while ((padnum - 1) > 0)
                {
                    finalres->val[finalindex] = 0;
                    finalindex++;
                    padnum--;
                }
            }
            finalres->val[finalindex] = result[i];
            finalindex++;

            //add zero here
            nzindex++;
        }
    }

    //print the final result res
    //for (i = 0; i < mat->M; i++)
    //{
    //    printf("final lnum %d result %f\n", i, finalres->val[i]);
    //}

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);

    printf("Preprocess Time: %lu milli-seconds\n", 1000 * (preend.tv_sec - prestart.tv_sec) + (preend.tv_nsec - prestart.tv_nsec) / 1000000);
    printf("Kernel Time: %lu milli-seconds\n", 1000 * (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1000000);

    /*Deallocate, please*/
    //cudaDeviceReset();
    free(guard_record);
    free(result);
    free(flag);

    //be careful for bus error if do the following operation first
    //cudaFree(gpu_mat->rIndex);
    //cudaFree(gpu_mat->cIndex);
    //cudaFree(gpu_mat->val);

    cudaFree(gpu_mat);
    cudaFree(gpu_vec);

    //return 1 if everything is ok
    //printf("curr errCode %d\n", errCode);
    return errCode;
}
