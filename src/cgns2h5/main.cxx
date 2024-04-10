// C/C++ headers
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cstdint>

// Lib headers
#include "mpiWrap.h"

int main( int argc, char *argv[] )
{
    // Initialize MPI env
    int mpi_nprocs, mpi_rank;
    init_MPI(mpi_nprocs, mpi_rank);
    return 0;
}