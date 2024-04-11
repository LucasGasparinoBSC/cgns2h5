// C/C++ headers
#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cstdint>

// MPI header
#include <mpi.h>

// Parallel CGNS header
#include <pcgnslib.h>

// HDF5 headers
#include <hdf5.h>
#include <hdf5_hl.h>

// GPU headers
#include <cuda.h>
#include <cuda_runtime.h>
#include <nvToolsExt.h>
#include <openacc.h>

// Element Type processor
#include "cgnsElemInfo.h"

int main( int argc, char *argv[] )
{
    // Initialize MPI
    MPI_Init( NULL, NULL );

    // Set mpi_nprocs and mpi_rank
    int mpi_nprocs, mpi_rank;
    MPI_Comm_size( MPI_COMM_WORLD, &mpi_nprocs );
    MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );

    // for now, support runs with only 1 process
    if ( mpi_nprocs != 1 )
    {
        if ( mpi_rank == 0 )
        {
            std::cerr << "Error: This program only supports runs with 1 process" << std::endl;
        }
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }

    // Check if the number of arguments is correct
    if ( argc != 3 )
    {
        if ( mpi_rank == 0 )
        {
            std::cerr << "Too few arguments!" << std::endl;
            std::cerr << "Usage: " << argv[0] << " <input.cgns> <output.h5>" << std::endl;
        }
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }

    // Input should be a CGNS file
    if ( std::string(argv[1]).find(".cgns") == std::string::npos )
    {
        if ( mpi_rank == 0 )
        {
            std::cerr << "Error: Input file should be a CGNS file" << std::endl;
        }
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }

    // Output should be an HDF5 file
    if ( std::string(argv[2]).find(".h5") == std::string::npos )
    {
        if ( mpi_rank == 0 )
        {
            std::cerr << "Error: Output file should be an HDF5 file" << std::endl;
        }
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }

    // Reader:

    // Set the MPI communicator for parallel CGNS
    cgp_mpi_comm(MPI_COMM_WORLD);

    // Set input and output file names from args
    std::string input_cgns = argv[1];
    std::string output_h5 = argv[2];

    // Use parallel CGNS to open the input file
    int cgns_file;
    if ( cgp_open( input_cgns.c_str(), CG_MODE_READ, &cgns_file ) != CG_OK )
    {
        if ( mpi_rank == 0 )
        {
            std::cerr << "Error: Cannot open CGNS file " << input_cgns << std::endl;
        }
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }

    // Simple version: assume 1 base with 1 zone
    const int n_bases = 1;
    const int idx_Base = 1;
    const int n_zones = 1;
    const int idx_Zone = 1;

    // Get zone information
    cgsize_t zone_size[1][3];
    char zone_name[33];
    if ( cg_zone_read( cgns_file, idx_Base, idx_Zone, zone_name, (cgsize_t*)zone_size ) ) {
        if ( mpi_rank == 0 )
        {
            std::cerr << "Error: Cannot read zone information" << std::endl;
        }
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }
    int npoin = zone_size[0][0]; // Number of grid nodes
    int nelem = zone_size[0][1]; // Number of elements
    int nface = zone_size[0][2]; // TODO: Need better name for this

    // Allocate data for the coordinates
    float *x = new float[npoin];
    float *y = new float[npoin];
    float *z = new float[npoin];
    cgsize_t irmin, irmax, istart, iend;

    // Lower and upper range indexes
    irmin = 1;
    irmax = npoin;

    // Read the coordinates using serial CGNS
    cg_coord_read( cgns_file, idx_Base, idx_Zone, "CoordinateX",
                   CGNS_ENUMV(RealSingle), &irmin, &irmax, x );
    cg_coord_read( cgns_file, idx_Base, idx_Zone, "CoordinateY",
                   CGNS_ENUMV(RealSingle), &irmin, &irmax, y );
    cg_coord_read( cgns_file, idx_Base, idx_Zone, "CoordinateZ",
                   CGNS_ENUMV(RealSingle), &irmin, &irmax, z );

    // Get the number of sections
    int nsections;
    cg_nsections( cgns_file, idx_Base, idx_Zone, &nsections );
    if ( mpi_rank == 0 ) printf("Number of sections: %d\n", nsections);

    // Read the element connectivity
    int nbound;
    int iparent_flag;
    char section_name[33];
    cgsize_t connec[nelem][27], iparent_data;
    CGNS_ENUMT(ElementType_t) itype;
    if ( mpi_rank == 0 ) printf("Reading section data...\n");
    for (int idx_sec = 1; idx_sec <= nsections; idx_sec++)
    {
        cg_section_read( cgns_file, idx_Base, idx_Zone, idx_sec, section_name,
                         &itype, &istart, &iend, &nbound, &iparent_flag );
        if ( mpi_rank == 0 )
        {
            printf("Section %d: %s\n", idx_sec, section_name);
            printf("  Element type: %s\n", ElementTypeName[itype]);
            printf("  Start: %d\n", istart);
            printf("  End: %d\n", iend);
            printf("  Number of bounds: %d\n", nbound);
            printf("  Parent flag: %d\n", iparent_flag);
            int nnode, porder, eldim;
            getElemInfo(itype, nnode, porder, eldim);
        }
        if ( itype == CGNS_ENUMV(HEXA_27) )
        {
            cg_elements_read( cgns_file, idx_Base, idx_Zone, idx_sec, connec[0], &iparent_data );
            if ( mpi_rank == 0 ) printf("  Parent data: %d\n", iparent_data);
        }
    }

    // Print the first 5 elements connectivity
    if ( mpi_rank == 0 )
    {
        printf("First 5 elements connectivity:\n");
        for (int i = 0; i < 5; i++)
        {
            printf("Element %d: ", i);
            for (int j = 0; j < 27; j++)
            {
                printf("%d ", connec[i][j]);
            }
            printf("\n");
        }
    }

    // Writer:

    // Create a HDF5 file with parallel I/O
    hid_t fileOut = H5Fcreate( output_h5.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    // Close the HDF5 file
    H5Fclose( fileOut );

    // Finalize MPI
    MPI_Finalize();
    return 0;
}