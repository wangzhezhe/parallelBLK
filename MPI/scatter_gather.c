
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

#define N 60

int main(int argc, char *argv[])
{
    int i, my_id, num_procs, num_elem;
    int array[N], array_final[N];
    int *array_recv;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    // Write the code to complete the exercise below.
    // if thread number is 11 (not divisible by 60)
    int num_task_for_everythread = 10;
    array_recv = (int *)malloc(sizeof(int) * num_task_for_everythread);
    //init the original array
    if (my_id == 0)
    {
        for (i = 0; i < N; i++)
        {
            array[i] = i;
            printf("init index %d value %d\n", i, array[i]);
        }
    }

    MPI_Scatter(
        array,
        num_task_for_everythread,
        MPI_INT,
        array_recv,
        num_task_for_everythread,
        MPI_INT,
        0,
        MPI_COMM_WORLD);

    for (i = 0; i < num_task_for_everythread; i++)
    {
        printf("thread id %d index %d recieve value %d\n", my_id, i, array_recv[i]);

        //update operation
        array_recv[i] = array_recv[i] + 1;
    }



    MPI_Gather(
        array_recv,
        num_task_for_everythread,
        MPI_INT,
        array_final,
        num_task_for_everythread,
        MPI_INT,
        0,
        MPI_COMM_WORLD);

    if (my_id == 0)
    {
        for (i = 0; i < N; i++)
        {

            printf("init index %d final value %d\n", i, array_final[i]);
        }
    }
    MPI_Finalize();
    return 0;
}
