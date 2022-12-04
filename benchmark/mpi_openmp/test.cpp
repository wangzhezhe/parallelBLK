#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

int main(int argc, char *argv[])
{

    MPI_Init(&argc, &argv);

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char name[MPI_MAX_PROCESSOR_NAME];
    int resultlength;
    MPI_Get_processor_name(name, &resultlength);
    int num_threads = 0;
#pragma omp parallel default(shared)
    {
        num_threads = omp_get_num_threads();
    }

    printf("\n----------Processor name: %s, MPI Ranks: %d, Currnet rank %d, OpenMP Threads: %d----------\n", name, size, rank, num_threads);
    MPI_Finalize();
    return 0;
}