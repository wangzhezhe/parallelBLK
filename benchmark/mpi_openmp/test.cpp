#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <omp.h>

#include <unistd.h>
#include <sys/syscall.h>
#include <sys/types.h>

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
   
    unsigned cpu, node;
 
    // Get current CPU core and NUMA node via system call
    // Note this has no glibc wrapper so we must call it directly
    syscall(SYS_getcpu, &cpu, &node, NULL);
 
    // Display information
    printf("Rank %d is running on CPU core %u and NUMA node %u.\n\n", rank, cpu, node);
 
    
    MPI_Finalize();
    return 0;
}