#include <stdio.h>
#include <string.h>
#include "spmv.cuh"
#include "genresult.cuh"
#include "mmio.h"

void logError(const char *errArg, const char *eMsg)
{
	if (eMsg != NULL)
		printf("Error: %s\n", eMsg);
	if (errArg != NULL)
		printf("Error found near: '%s'\n", errArg);
	puts("USAGE: spmv -mat [matrix file] -ivec [vector file] -alg [segment|design] -blksize [blocksize] -blknum [blocknum]");
	puts("Where the order of the parameters and string case do not matter");
	puts("And the algorithms are:");
	puts("     - segment = simple segment based scan approach");
	puts("     - design = your own design implementation");
}

typedef enum { CMDLN_ARG_NULL,
			   CMDLN_ARG_RTYPE = 0,
			   CMDLN_ARG_MAT = 1,
			   CMDLN_ARG_VEC = 2,
			   CMDLN_ARG_ALG = 4,
			   CMDLN_ARG_BLOCK = 8,
			   CMDLN_ARG_BLOCKNUM = 16,
			   CMDLN_ARG_ERR = 32 } CmdLnArg;

CmdLnArg getArgType(const char *argv)
{
	//path of the matrix file
	if (strcasecmp(argv, "-mat") == 0)
		return CMDLN_ARG_MAT;
	//path of the vector file
	else if (strcasecmp(argv, "-ivec") == 0)
		return CMDLN_ARG_VEC;
	//algorithm
	//segment means simple segment based scan approach
	//design means own designed implementation
	else if (strcasecmp(argv, "-alg") == 0)
		return CMDLN_ARG_ALG;
	//You can have at most 65535 blocks in one dimension
	else if (strcasecmp(argv, "-blksize") == 0)
		return CMDLN_ARG_BLOCK;

	else if (strcasecmp(argv, "-blknum") == 0)
		return CMDLN_ARG_BLOCKNUM;

	else if (strcasecmp(argv, "-rtype") == 0)
	{
		printf("get rtype info\n");
		return CMDLN_ARG_RTYPE;
	}

	else
		return CMDLN_ARG_ERR;
}

typedef enum { ALG_SEGMENT,
			   ALG_DESIGN,
			   ALG_OPTFIRST,
			   ALG_OPTGRAPH
			    } AlgType;

int populateAlgType(const char *argv, AlgType *toPop)
{
	if (strcasecmp(argv, "segment") == 0)
	{
		*toPop = ALG_SEGMENT;
		return 1;
	}
	else if (strcasecmp(argv, "design") == 0)
	{
		*toPop = ALG_DESIGN;
		return 1;
	}
		else if (strcasecmp(argv, "optfirst") == 0)
	{
		*toPop = ALG_OPTFIRST;
		return 1;
	}		
	else if (strcasecmp(argv, "optgraph") == 0)
	{
		*toPop = ALG_OPTGRAPH;
		return 1;
	}
	else
	{
		return 0;
	}
}

int doSpmv(MatrixInfo *mat, MatrixInfo *vec, MatrixInfo *res, AlgType how, int blockSize, int blockNum, int rtype)
{
	int errCode = 0;
	switch (how)
	{
	case ALG_SEGMENT:
		errCode = getMulScan(mat, vec, res, blockSize, blockNum);
		return errCode;
	case ALG_DESIGN:
		errCode = getMulDesign(mat, vec, res, blockSize, blockNum, 0);
		//errCode =1 if everything is ok here
		return errCode;

	case ALG_OPTFIRST:
		errCode = getMulDesign(mat, vec, res, blockSize, blockNum, 1);
		//errCode =1 if everything is ok here
		return errCode;
	case ALG_OPTGRAPH:
		errCode = getMulDesignByGraph(mat, vec, res, blockSize, blockNum);
		//errCode =1 if everything is ok here
		return errCode;

	default:
		return errCode;
	}
}

/* a function that verifies the output with the provided sample solution
 * the function will print the total number of incorrect rows */
int verify(const int nz, const int M, const int *rIndex, const int *cIndex, const float *val, const float *vec, const float *res)
{

	float *correct = (float *)malloc(sizeof(float) * M);
	memset(correct, 0, sizeof(float) * M);

	/* get the correct output vector */

	for (int i = 0; i < nz; ++i)
	{
		correct[rIndex[i]] += val[i] * vec[cIndex[i]];
	}

	int o = 0; // the total number of incorrect rows, initialized to 0

	for (int i = 0; i < M; ++i)
	//for (int i = 0; i < 20; ++i)
	{
		float l = correct[i] > 0 ? correct[i] : -1 * correct[i];
		float m = res[i] > 0 ? res[i] : -1 * res[i];
		float k = l - m > 0 ? l - m : m - l;
		float rel = k / l;
		if (rel > .01)
		{
			o++;
			printf("Yours - %f, correct - %f, Relative error - %f the rownumber %d\n", res[i], correct[i], rel, i);
		}
	}

	return o;
}

int main(int argc, char **argv)
{
	//printf("debug argc (%d)\n", argc);
	//the command line check is original
	//first element is progress name
	//second element -mat third element is the path of matrix
	//11 commands requires 5 input parameters

	if (argc != 11)
	{
		logError(NULL, NULL);
		return 1;
	}


	//This is so that the arguments can be presented in any order with the blocksize defaulting to 1024
	int cumArgs = CMDLN_ARG_NULL;
	CmdLnArg argOrder[10];
	int i;
	for (i = 1; i < argc; i += 2)
	{
		CmdLnArg currArg = getArgType(argv[i]);

		if (currArg == CMDLN_ARG_ERR)
		{
			logError(argv[i], "Invalid or duplicate argument.");
			return 1;
		}
		else
		{
			argOrder[i / 2] = currArg; //May the truncation be ever in our favor.
			cumArgs |= currArg;
		}
	}

	if (!(31 & cumArgs))
	{
		logError(NULL, "Missing arguments!");
		return 1;
	}

	char *mFile, *vFile;
	int rtype;
	AlgType algo; //Si, debe ser algo!
	int blockSize;
	int blockNum;

	for (i = 0; i < (argc - 1) / 2; i++)
	{
		switch (argOrder[i])
		{
		case CMDLN_ARG_RTYPE:
			if (sscanf(argv[i * 2 + 2], "%d", &rtype) != 1 || rtype < 0)
			{
				logError(argv[i * 2 + 2], "rtype must be 0 (no reorder) or 1 (first access) or 2 (graph)");
				return 1;
			}
			break;
		case CMDLN_ARG_ALG:
			if (!populateAlgType(argv[i * 2 + 2], &algo))
			{
				logError(argv[i * 2 + 2], "Unsupported algorithm");
				return 1;
			}
			break;
		case CMDLN_ARG_MAT:
			mFile = argv[i * 2 + 2];
			break;
		case CMDLN_ARG_VEC:
			vFile = argv[i * 2 + 2];
			break;
		case CMDLN_ARG_BLOCK:
			if (sscanf(argv[i * 2 + 2], "%d", &blockSize) != 1 || blockSize <= 0)
			{
				logError(argv[i * 2 + 2], "Block size must be a positive integer (greater than 0)");
				return 1;
			}
			break;
		case CMDLN_ARG_BLOCKNUM:
			if (sscanf(argv[i * 2 + 2], "%d", &blockNum) != 1 || blockNum <= 0)
			{
				logError(argv[i * 2 + 2], "Block num must be a positive integer (greater than 0)");
				return 1;
			}
			break;

		default:
			puts("Logic is literally broken. This should never be seen!");
		}
	}

	printf("Reading matrix from %s\n", mFile);
	MatrixInfo *matrix = read_file(mFile);
	if (matrix == NULL)
	{
		logError(mFile, "Error regarding matrix file.");
		return 1;
	}

	printf("Reading vector from %s\n", vFile);
	printf("the value of N: (%d)\n", matrix->N);
	MatrixInfo *vector = read_vector_file(vFile, matrix->N);
	if (vector == NULL)
	{
		logError(mFile, "Error regarding vector file.");
		return 1;
	}
	//give product the initialised space and let the pointer be null
	MatrixInfo *product = initMatrixResult(matrix->M, blockSize);
	cudaError_t err;
	if (doSpmv(matrix, vector, product, algo, blockSize, blockNum, rtype) != 1 || (err = cudaDeviceSynchronize()) != cudaSuccess || !writeVect(product, "output.txt"))
	{

		printf("\x1b[31m%s\x1b[0m\n", cudaGetErrorString(err));
		logError(NULL, "Failed to produce output");
	}
	else
	{
		printf("start Verifying...\n");

		int o = verify(matrix->nz, matrix->M, matrix->rIndex, matrix->cIndex, matrix->val, vector->val, product->val);
		//int o=1;
		if (o > 0)
		{
			printf("\t %d Error rows found \n", o);
			puts("Test failed!");
			//freeMatrixInfo(matrix);
			//freeMatrixInfo(vector);
			//freeMatrixInfo(product);
			return 1;
		}
	}

	//freeMatrixInfo(matrix);
	//freeMatrixInfo(vector);
	//freeMatrixInfo(product);

	puts("Test passed!");
	return 0;
}
