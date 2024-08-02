#include <stdio.h>
#include "mpi.h"
#include <iostream>

////////////
// MPI_Isend
////////////
//
// int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,
//              MPI_Comm comm, MPI_Request *request)
//
// This example uses MPI_Isend to do a non-blocking send of information from the root process to a destination process.
// The destination process is set as a variable in the code and must be less than the number of processes started.
//
// example usage:
//		compile: mpicc -o mpi_isend mpi_isend.c
//		run: mpirun -n 4 mpi_isend
//

int main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);
    int rank, size;
    int tag, destination, count;
    int buffer = 345678; // value to send

    tag = 1234;
    destination = 1; // destination process
    count = 1;       // number of elements in buffer

    MPI_Status status;
    // this request is a kind of the handler
    // for processing the request
    MPI_Request request = MPI_REQUEST_NULL;

    MPI_Comm_size(MPI_COMM_WORLD, &size); // number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // rank of current process

    if (size != 2)
    {
        std::cout << "only works for the case where there are 2 processes" << std::endl;
        exit(0);
    }

    if (rank == 0)
    {
        printf("rank 0 send value to processor %d:\n", destination);
        MPI_Isend(&buffer, count, MPI_INT, destination, tag, MPI_COMM_WORLD, &request); // non blocking send to destination process
    }

    if (rank == destination)
    {
        MPI_Irecv(&buffer, count, MPI_INT, 0, tag, MPI_COMM_WORLD, &request); // destination process receives
    }

    //the code blocks and wait here
    MPI_Wait(&request, &status); // bloks and waits for destination process to receive data

    printf("MPI process %d received value %d from rank %d, with tag %d and error code %d.\n",
           rank,
           buffer,
           status.MPI_SOURCE,
           status.MPI_TAG,
           status.MPI_ERROR);

    if (rank == 0)
    {
        printf("processor %d sent %d\n", rank, buffer);
    }

    if (rank == destination)
    {
        printf("processor %d got %d\n", rank, buffer);
    }

    MPI_Finalize();

    return 0;
}