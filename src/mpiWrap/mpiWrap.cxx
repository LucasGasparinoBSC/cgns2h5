#include "mpiWrap.h"

void init_MPI(int &mpi_np, int &mpi_r)
{
    // Initialize MPI env
    MPI_Init(NULL, NULL);

    // Get mpi_nprocs and mpi_rank
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_np);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_r);
}

void finalize_MPI()
{
    // Finalize MPI env
    MPI_Finalize();
}