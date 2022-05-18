#include "mpi.h"
#include <stdio.h>

int main(int argc, char *argv[])
{  
   MPI_Init(&argc, &argv);

   MPI_Comm intercomm; 
   MPI_Comm_get_parent(&intercomm);

   int sendbuf[2]; // redundant for worker.
   int recvbuf;

   MPI_Scatter(sendbuf, 1, MPI_INT, &recvbuf, 1, MPI_INT, 0, intercomm);
   printf("recvbuf = %d\n", recvbuf);

   MPI_Finalize();
   return 0;
}

/*
Rank 0 [Sat May 22 07:06:41 2021] [c3-0c2s4n2] Fatal error in MPI_Comm_spawn: Other MPI error, error stack:
MPI_Comm_spawn(144).........: MPI_Comm_spawn(cmd="worker_program", argv=(nil), maxprocs=2, MPI_INFO_NULL, root=0, MPI_COMM_SELF, intercomm=0x7fffffff4498, errors=(nil)) failed
MPID_Comm_spawn_multiple(70): Function MPID_Comm_spawn_multiple not implemented
srun: error: nid00722: task 0: Aborted
srun: Terminating job step 42763955.0
*/