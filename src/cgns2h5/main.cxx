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
#include "H5Writer.h"

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
    double *tmp = new double[npoin];
    double *xyz = new double[npoin*3];
    cgsize_t irmin, irmax, istart, iend;

    // Lower and upper range indexes
    irmin = 1;
    irmax = npoin;

    // Read the coordinates using serial CGNS
    cg_coord_read( cgns_file, idx_Base, idx_Zone, "CoordinateX",
                   CGNS_ENUMV(RealDouble), &irmin, &irmax, tmp );
    #pragma acc parallel loop
    for ( uint64_t i = 0; i < npoin; i++ )
    {
        xyz[i*3 + 0] = tmp[i];
    }

    cg_coord_read( cgns_file, idx_Base, idx_Zone, "CoordinateY",
                   CGNS_ENUMV(RealDouble), &irmin, &irmax, tmp );
    #pragma acc parallel loop
    for ( uint64_t i = 0; i < npoin; i++ )
    {
        xyz[i*3 + 1] = tmp[i];
    }
    cg_coord_read( cgns_file, idx_Base, idx_Zone, "CoordinateZ",
                   CGNS_ENUMV(RealDouble), &irmin, &irmax, tmp );
    #pragma acc parallel loop
    for ( uint64_t i = 0; i < npoin; i++ )
    {
        xyz[i*3 + 2] = tmp[i];
    }

    // Free the temporary buffer
    delete[] tmp;

    // Get the number of sections
    int nsections;
    cg_nsections( cgns_file, idx_Base, idx_Zone, &nsections );
    if ( mpi_rank == 0 ) printf("Number of sections: %d\n", nsections);

    // Extract basic info from each section
    int nnode_sec[nsections];             // Nodes/elem. in each section
    uint64_t nelem_sec[nsections];        // Number of elements in each section
    uint64_t nbelem_sec[nsections];       // Boundary elem. iin each section
    cgsize_t elemStartEnd[nsections][2];  // Element start and end in each section
    cgsize_t belemStartEnd[nsections][2]; // Boundary elem. start and end in each section
    int ibound;                           // Boundary index
    int snode, nnode, nbnode, porder, eldim;     // Element info
    int nbound = 0;                       // Total number of physical boundaries
    uint64_t nbelem = 0;                  // Total number of boundary elements
    CGNS_ENUMT(ElementType_t) itype;
    int iparent_flag;

    for ( int indx_sec = 1; indx_sec <= nsections; indx_sec++ )
    {
        // Read the section
        char section_name[33];
        cg_section_read( cgns_file, idx_Base, idx_Zone, indx_sec, section_name,
                         &itype, &istart, &iend, &ibound, &iparent_flag );
        if ( mpi_rank == 0 )
        {
            printf("Section %d: %s\n", indx_sec, section_name);
            printf("  Element type: %s\n", ElementTypeName[itype]);
            printf("  Number of bounds: %d\n", ibound);
        }

        // Get elem base type
        std::string elemTypeStr = ElementTypeName[itype];
        elemTypeStr = elemTypeStr.substr(0, elemTypeStr.find("_"));

        // Get element info
        getElemInfo(itype, snode, porder, eldim);
        nnode_sec[indx_sec-1] = snode;

        // Based on the element type,  fill section info
        if ( elemTypeStr == "QUAD" )
        {
            belemStartEnd[indx_sec-1][0] = istart;
            belemStartEnd[indx_sec-1][1] = iend;
            nbelem += iend - istart + 1;
            nbelem_sec[indx_sec-1] = iend - istart + 1;
            nbound++;
            nbnode = snode;
            nnode_sec[indx_sec-1] = nbnode;
            elemStartEnd[indx_sec-1][0] = 0;
            elemStartEnd[indx_sec-1][1] = 0;
            nelem_sec[indx_sec-1] = 0;
        }
        else if ( elemTypeStr == "HEXA" )
        {
            belemStartEnd[indx_sec-1][0] = 0;
            belemStartEnd[indx_sec-1][1] = 0;
            nbelem_sec[indx_sec-1] = 0;
            elemStartEnd[indx_sec-1][0] = istart;
            elemStartEnd[indx_sec-1][1] = iend;
            nelem_sec[indx_sec-1] = iend - istart + 1;
            nnode = snode;
            nnode_sec[indx_sec-1] = nnode;
        }

        // Print section info
        if ( mpi_rank == 0 )
        {
            printf("  Number of nodes: %d\n", nnode_sec[indx_sec-1]);
            printf("  Polynomial order: %d\n", porder);
            printf("  Element dimension: %d\n", eldim);
            printf("  Boundary info:\n");
            printf("    Start: %d, End: %d\n", belemStartEnd[indx_sec-1][0], belemStartEnd[indx_sec-1][1]);
            printf("    Number of boundary elements: %ld\n", nbelem_sec[indx_sec-1]);
            printf("  Element info:\n");
            printf("    Start: %d, End: %d\n", elemStartEnd[indx_sec-1][0], elemStartEnd[indx_sec-1][1]);
            printf("    Number of elements: %ld\n", nelem_sec[indx_sec-1]);
        }
    }

    printf("Total number of physical boundaries: %d\n", nbound);
    printf("Total number of boundary elements: %ld\n", nbelem);

    // Get the HEXA connectivity (only works for Xevi cases, asssumes 1st section is HEXA)
    std::cout << "Reading HEXA connectivity..." << std::endl;

    // Allocate data for the connectivity table
    cgsize_t iparent_data;
    cgsize_t *connecHEXA = new cgsize_t[nelem*nnode];

    // Read the connectivity
    cg_elements_read( cgns_file, idx_Base, idx_Zone, 1, connecHEXA, &iparent_data );

    // Convert to SOD format
    Conversor conv;
    cgsize_t *connecHEXA_SOD = new cgsize_t[nelem*nnode];
    memset( connecHEXA_SOD, 0, nelem*nnode*sizeof(cgsize_t) );
    std::cout << "Converting HEXA connectivity to SOD format..." << std::endl;
    conv.convert2sod_HEXA( porder, nelem, nnode, connecHEXA, connecHEXA_SOD );

    // Free the original connectivity
    delete[] connecHEXA;

    // Remaining sections are BCs, read them. Also crreate belemId array
    std::cout << "Reading boundary element connectivity..." << std::endl;

    istart = 0;
    cgsize_t *connecQUAD = new cgsize_t[nbelem*nbnode]; // Total BC connectivity
    uint64_t *belemId = new uint64_t[nbelem];           // Boundary element ID

    for ( int idx_sec = 2; idx_sec <= nsections; idx_sec++ )
    {
        printf("Reading section %d ...\n", idx_sec);
        // belemId is simply idxSec-1
        for ( int i = istart; i < istart+nbelem_sec[idx_sec-1]; i++ )
        {
            belemId[i] = idx_sec-1;
        }

        // Alloc. a buffer for the current BC connectivity
        cgsize_t *connecQUAD_sec = new cgsize_t[nbelem_sec[idx_sec-1]*nbnode];

        // Read the connectivity into the buffer
        cg_elements_read( cgns_file, idx_Base, idx_Zone, idx_sec, connecQUAD_sec, &iparent_data );

        // Copy to the total BC connectivity
        for ( int i = 0; i < nbelem_sec[idx_sec-1]; i++ )
        {
            for ( int j = 0; j < nbnode; j++ )
            {
                connecQUAD[(istart+i)*nbnode+j] = connecQUAD_sec[i*nbnode+j];
            }
        }
        istart += nbelem_sec[idx_sec-1];

        // Free the buffer
        delete[] connecQUAD_sec;
    }

    // Convert to SOD format
    cgsize_t *connecQUAD_SOD = new cgsize_t[nbelem*nbnode];
    if ( nbound > 0 )
    {
        memset( connecQUAD_SOD, 0, nbelem*nbnode*sizeof(cgsize_t) );
        std::cout << "Converting QUAD connectivity to SOD format..." << std::endl;
        conv.convert2sod_QUAD( porder, nbelem, nbnode, connecQUAD, connecQUAD_SOD );
    }

    // Free the original connectivity
    delete[] connecQUAD;

    // Close the CGNS file
    cgp_close( cgns_file );

    // Writer:

    // Create a HDF5 file with parallel I/O
    hid_t fileOut;
    fileOut = H5Fcreate( output_h5.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    // Basic hdf5 variables
    hid_t dataspace, dataset, group;
    herr_t status_hdf;
    hsize_t data1d[1];
    hsize_t data2d[2];

    // Create a dataset for the bounndary connectivity table
    data2d[0] = nbelem;
    data2d[1] = nbnode;

    dataspace = H5Screate_simple( 2, data2d, NULL );
    dataset = H5Dcreate( fileOut, "/boundFaces", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the boundary connectivity table
    status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, connecQUAD_SOD );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // Create a dataset for belemId
    data1d[0] = nbelem;

    dataspace = H5Screate_simple( 1, data1d, NULL );
    dataset = H5Dcreate( fileOut, "/boundFacesId", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the boundary element ID
    status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, belemId );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // Create a dataset for the HEXA connectivity table
    data2d[0] = nelem;
    data2d[1] = nnode;

    dataspace = H5Screate_simple( 2, data2d, NULL );
    dataset = H5Dcreate( fileOut, "/connec", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the HEXA connectivity table
    status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, connecHEXA_SOD );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // Create a dataset for the coordinates
    data2d[0] = npoin;
    data2d[1] = 3;

    dataspace = H5Screate_simple( 2, data2d, NULL );
    dataset = H5Dcreate( fileOut, "/coords", H5T_IEEE_F64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the coordinates
    status_hdf = H5Dwrite( dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xyz );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // From "/" create the dims group
    group = H5Gcreate( fileOut, "/dims", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Create a dataset for the number of boundary elements
    data1d[0] = 1;

    dataspace = H5Screate_simple( 1, data1d, NULL );
    dataset = H5Dcreate( fileOut, "/dims/numBoundaryFaces", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the number of boundary elements
    status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nbelem );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // Create a dataset for the number of elements
    dataspace = H5Screate_simple( 1, data1d, NULL );
    dataset = H5Dcreate( fileOut, "/dims/numElements", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the number of elements
    status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nelem );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // Create a dataset for the number of mapped faces
    dataspace = H5Screate_simple( 1, data1d, NULL );
    dataset = H5Dcreate( fileOut, "/dims/numMappedFaces", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the number of mapped faces (0 for now)
    uint64_t numMappedFaces = 0;
    status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &numMappedFaces );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // Create a dataset for the nodes per element
    data1d[0] = 1;

    dataspace = H5Screate_simple( 1, data1d, NULL );
    dataset = H5Dcreate( fileOut, "/dims/numNodes", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the nodes per element
    uint64_t numNodes = npoin;
    status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &numNodes );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // Create a dataset for the number of periodic faces
    dataspace = H5Screate_simple( 1, data1d, NULL );
    dataset = H5Dcreate( fileOut, "/dims/numPeriodicFaces", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the number of periodic faces (0 for now)
    uint64_t numPeriodicFaces = 0;
    status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &numPeriodicFaces );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // Create a dataset for the number of periodic links
    dataspace = H5Screate_simple( 1, data1d, NULL );
    dataset = H5Dcreate( fileOut, "/dims/numPeriodicLinks", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the number of periodic links (0 for now)
    uint64_t numPeriodicLinks = 0;
    status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, &numPeriodicLinks );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // Create a dataset for the element order
    dataspace = H5Screate_simple( 1, data1d, NULL );
    dataset = H5Dcreate( fileOut, "/dims/order", H5T_STD_I32LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the element order
    int order = porder;
    status_hdf = H5Dwrite( dataset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &order );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // Close the dims group
    status_hdf = H5Gclose( group );

    // On the root group, create a dataset for the mapped faces
    data2d[0] = numMappedFaces;
    data2d[1] = nbnode;

    dataspace = H5Screate_simple( 2, data2d, NULL );
    dataset = H5Dcreate( fileOut, "/mappedFaces", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the mapped faces (empty for now)
    uint64_t *mappedFaces;
    if ( numMappedFaces > 0 )
    {
        status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, mappedFaces );
    }

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // On the root group, create a dataset for the periodic faces
    data2d[0] = numPeriodicFaces;
    data2d[1] = nbnode;

    dataspace = H5Screate_simple( 2, data2d, NULL );
    dataset = H5Dcreate( fileOut, "/periodicFaces", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the periodic faces (empty for now)
    uint64_t *periodicFaces;
    if ( numPeriodicFaces > 0 )
    {
        status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, periodicFaces );
    }

    // Create a dataset for the periodic links
    data2d[0] = numPeriodicLinks;
    data2d[1] = 2;

    dataspace = H5Screate_simple( 2, data2d, NULL );
    dataset = H5Dcreate( fileOut, "/periodicLinks", H5T_STD_I64LE, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the periodic links (empty for now)
    uint64_t *periodicLinks;
    if ( numPeriodicLinks > 0 )
    {
        status_hdf = H5Dwrite( dataset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, H5P_DEFAULT, periodicLinks );
    }

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset );
    status_hdf = H5Sclose( dataspace );

    // Close the HDF5 file
    status_hdf = H5Fclose( fileOut );

    // Finalize MPI
    MPI_Finalize();
    return 0;
}