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
#include "Conversor.h"

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
    const int idx_Base = 1;
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
    uint64_t npoin = zone_size[0][0]; // Number of grid nodes
    uint64_t nelem = zone_size[0][1]; // Number of elements
    uint64_t nface = zone_size[0][2]; // TODO: Need better name for this
    if ( mpi_rank == 0 )
    {
        printf("Zone name: %s\n", zone_name);
        printf("Number of points: %ld\n", npoin);
        printf("Number of elements: %ld\n", nelem);
        printf("Number of faces: %ld\n", nface);
    }

    // Allocate data for the coordinates (assuming 3D)
    double *x = new double[npoin];
    double *y = new double[npoin];
    double *z = new double[npoin];
    double *xyz = new double[npoin*3];
    cgsize_t irmin, irmax, istart, iend;

    // Lower and upper range indexes
    irmin = 1;
    irmax = npoin;

    // Read the coordinates using serial CGNS
    cg_coord_read( cgns_file, idx_Base, idx_Zone, "CoordinateX",
                   CGNS_ENUMV(RealDouble), &irmin, &irmax, xyz );
    cg_coord_read( cgns_file, idx_Base, idx_Zone, "CoordinateY",
                   CGNS_ENUMV(RealDouble), &irmin, &irmax, xyz+npoin );
    cg_coord_read( cgns_file, idx_Base, idx_Zone, "CoordinateZ",
                   CGNS_ENUMV(RealDouble), &irmin, &irmax, xyz+npoin*2 );

    // Get the number of sections
    int nsections;
    cg_nsections( cgns_file, idx_Base, idx_Zone, &nsections );
    if ( mpi_rank == 0 ) printf("Number of sections: %d\n", nsections);

    // Read the element connectivity
    int nnode, nbound, porder, eldim;
    int sec_dims_nnode[nsections];
    int iparent_flag;
    char section_name[33];
    cgsize_t* connec = new cgsize_t[nelem];
    cgsize_t iparent_data;
    CGNS_ENUMT(ElementType_t) itype;
    cgsize_t* connecSOD2D = (cgsize_t*)malloc(nelem*sizeof(cgsize_t));
    if ( mpi_rank == 0 ) printf("Reading section data...\n");
    for (int idx_sec = 1; idx_sec <= nsections; idx_sec++)
    {
        // Read the rection
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
        }
        sec_dims_nnode[idx_sec] = nnode;

        // Reallocate memory based on the number of nodes
        getElemInfo(itype, nnode, porder, eldim);
        connec = (cgsize_t*)realloc(connec, nelem*nnode*sizeof(cgsize_t));
        connecSOD2D = (cgsize_t*)realloc(connecSOD2D, nelem*nnode*sizeof(cgsize_t));

        // Extract the element base type, without number of nodes, to a string
        std::string elemTypeStr = ElementTypeName[itype];
        elemTypeStr = elemTypeStr.substr(0, elemTypeStr.find("_"));

        // If the base is a HEXA, read the connectivity and convert to SOD2D
        if ( elemTypeStr == "HEXA" )
        {
            cg_elements_read( cgns_file, idx_Base, idx_Zone, idx_sec, connec, &iparent_data );
            if ( mpi_rank == 0 ) printf("  Parent data: %d\n", iparent_data);
            // Conversion to SO2D format:
            Conversor conv;
            conv.convert2sod_HEXA(porder, nelem, nnode, connec, connecSOD2D);
        }
    }

    // Close the CGNS file
    cgp_close( cgns_file );

    // Writer:

    // Create a HDF5 file with parallel I/O
    hid_t fileOut;
    herr_t status_hdf;
    fileOut = H5Fcreate( output_h5.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    // Create a dataset for the connectivity table
    hsize_t dims_coord[2];
    dims_coord[0] = nelem;
    dims_coord[1] = (porder+1)*(porder+1)*(porder+1);

    hid_t dataspace_id = H5Screate_simple( 2, dims_coord, NULL );
    hid_t dataset_id = H5Dcreate( fileOut, "/connec", H5T_STD_I64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the connectivity table
    status_hdf = H5Dwrite( dataset_id, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, connecSOD2D );
    if ( status_hdf < 0 )
    {
        if ( mpi_rank == 0 )
        {
            std::cerr << "Error: Cannot write the connectivity table" << std::endl;
        }
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }

    // Destroy connec and close the dataset and dataspace
    delete[] connec;
    status_hdf = H5Dclose( dataset_id );
    status_hdf = H5Sclose( dataspace_id );

    dims_coord[0] = npoin;
    dims_coord[1] = 3;

    dataspace_id = H5Screate_simple( 2, dims_coord, NULL );
    dataset_id = H5Dcreate( fileOut, "/coord", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the coordinates
    status_hdf = H5Dwrite( dataset_id, H5T_IEEE_F64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz );

    // Destroy coord and close the dataset and dataspace
    status_hdf = H5Dclose( dataset_id );
    status_hdf = H5Sclose( dataspace_id );

    // Close the HDF5 file
    status_hdf = H5Fclose( fileOut );

    // Finalize MPI
    MPI_Finalize();
    return 0;
}