#include <stdio.h>
#include <string.h>
#include <mpi.h>

main(int argc, char *argv[])
{
    int p_rank;
    int p_num;
    int source;
    int dest;
    int tag = 0;
    char message[100];
    char process_name[100];
    int name_len;
    MPI_Status status;

    //MPI_Init(&argc, &argv);
    MPI_Init(NULL, NULL);

    MPI_Comm_rank(MPI_COMM_WORLD, &p_rank);

    MPI_Comm_size(MPI_COMM_WORLD, &p_num);

    if (p_rank != 0)
    {
        //create message
        sprintf(message, "Greetings from the process %d", p_rank);
        dest = 0;
        tag = p_rank;
        MPI_Send(message, strlen(message) + 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
        MPI_Get_processor_name(process_name, &name_len);

        printf("%s send message\n", process_name);
    }
    else
    { //if the process with rank 0
        for (source = 1; source < p_num; source++)
        {
            MPI_Get_processor_name(process_name, &name_len);
            //tag=source-1;
            //MPI_Recv(message, 100, MPI_CHAR, source, tag , MPI_COMM_WORLD, &status);
            MPI_Recv(message, 100, MPI_CHAR, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
            printf("process id %d in processor %s get message : %s, comm src rank %d comm tag %d\n",
                   p_rank, process_name, message, status.MPI_SOURCE, status.MPI_SOURCE);
        }
    }

    //shut down
    MPI_Finalize();
}