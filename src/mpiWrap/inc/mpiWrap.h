#ifndef MPI_WRAP_H_
#define MPI_WRAP_H_

#include <mpi.h>

void init_MPI(int &mpi_np, int &mpi_r);
void finalize_MPI();

#endif // !MPI_WRAP_H_