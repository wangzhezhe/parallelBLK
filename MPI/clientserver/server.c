// https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node248.htm
// http://silveiraneto.net/estudos/mpi-client-server/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "mpi.h"

int main(int argc, char *argv[]) {
  int size, again;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  printf("size %d\n", size);

  MPI_Comm client;
  MPI_Status status;
  char port_name[MPI_MAX_PORT_NAME];
  int MAX_DATA = 10;

  double buf[MAX_DATA];

  if (size != 1) {
    printf("error, server too big");
    exit(-1);
  }
  MPI_Open_port(MPI_INFO_NULL, port_name);
  printf("server available at %s\n", port_name);
  while (1) {
    MPI_Comm_accept(port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &client);
    again = 1;
    while (again) {
      MPI_Recv(buf, MAX_DATA, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, client,
               &status);
      switch (status.MPI_TAG) {
        case 0:
          MPI_Comm_free(&client);
          MPI_Close_port(port_name);
          MPI_Finalize();
          return 0;
        case 1:
          MPI_Comm_disconnect(&client);
          again = 0;
          break;
        case 2: /* do something */
          printf("do sth/n");
          sleep(1);
        default:
          /* Unexpected message type */
          MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
  }
}