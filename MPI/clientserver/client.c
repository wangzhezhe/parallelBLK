#include "mpi.h"
#include <string.h>
int main(int argc, char **argv) {
  MPI_Comm server;
  int MAX_DATA=10;
  double buf[MAX_DATA];
  char port_name[MPI_MAX_PORT_NAME];

  MPI_Init(&argc, &argv);
  strcpy(port_name, argv[1]); /* assume server's name is cmd-line arg */

  MPI_Comm_connect(port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &server);
  
  int ifloop=1;
  while (ifloop) {
    int tag = 2; /* Action to perform */
    MPI_Send(buf, MAX_DATA, MPI_DOUBLE, 0, tag, server);
    ifloop = 0;
  }
  MPI_Send(buf, 0, MPI_DOUBLE, 0, 1, server);
  MPI_Comm_disconnect(&server);
  MPI_Finalize();
  return 0;
}