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

    // Read the connectivities
    std::vector<cgsize_t> connecQUAD;
    std::vector<cgsize_t> connecHEXA;

    // Loop through the sections
    for ( int idx_sec = 1; idx_sec <= nsections; idx_sec++ )
    {
        // Read the section
        char section_name[33];
        cg_section_read( cgns_file, idx_Base, idx_Zone, idx_sec, section_name,
                         &itype, &istart, &iend, &nbound, &iparent_flag );

        // Get elem base type
        std::string elemTypeStr = ElementTypeName[itype];
        elemTypeStr = elemTypeStr.substr(0, elemTypeStr.find("_"));

        // Read the connectivity
        if ( elemTypeStr == "QUAD" )
        {
            // Allocate memory for the boundary connectivity
            cgsize_t* connecBound = new cgsize_t[nbelem_sec[idx_sec-1]*nnode_sec[idx_sec-1]];

            // Read the boundary elements
            cg_elements_read( cgns_file, idx_Base, idx_Zone, idx_sec, connecBound, &iparent_flag );

            // Append to the vector
            for ( int idx = 0; idx < nbelem_sec[idx_sec-1]*nnode_sec[idx_sec-1]; idx++ )
            {
                connecQUAD.push_back(connecBound[idx]);
            }

            // Free memory
            delete[] connecBound;
        }
        else if ( elemTypeStr == "HEXA" )
        {
            // Allocate memory for the boundary connectivity
            cgsize_t* connecElem = new cgsize_t[nelem_sec[idx_sec-1]*nnode_sec[idx_sec-1]];

            // Read the boundary elements
            cg_elements_read( cgns_file, idx_Base, idx_Zone, idx_sec, connecElem, &iparent_flag );

            // Append to the vector
            for ( int idx = 0; idx < nelem_sec[idx_sec-1]*nnode_sec[idx_sec-1]; idx

    /*
    // Read the element connectivity
    int nnode, nbound, porder, eldim;
    int sec_dims_nnode[nsections];
    int iparent_flag;
    char section_name[33];
    uint64_t nbelem;
    cgsize_t iparent_data;
    CGNS_ENUMT(ElementType_t) itype;
    cgsize_t* connecBound;
    cgsize_t* connecBoundSOD2D;
    cgsize_t* connec = new cgsize_t[nelem];
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
        printf("  Number of nodes: %d\n", nnode);
        printf("  Polynomial order: %d\n", porder);
        printf("  Element dimension: %d\n", eldim);
        connec = (cgsize_t*)realloc(connec, nelem*nnode*sizeof(cgsize_t));
        connecSOD2D = (cgsize_t*)realloc(connecSOD2D, nelem*nnode*sizeof(cgsize_t));

        // Extract the element base type, without number of nodes, to a string
        std::string elemTypeStr = ElementTypeName[itype];
        elemTypeStr = elemTypeStr.substr(0, elemTypeStr.find("_"));

        // If the base is a HEXA, read the connectivity and convert to SOD2D
        if ( elemTypeStr == "HEXA" )
        {
            cg_elements_read( cgns_file, idx_Base, idx_Zone, idx_sec, connec, &iparent_data );

            // Conversion to SO2D format:
            Conversor conv;
            conv.convert2sod_HEXA(porder, nelem, nnode, connec, connecSOD2D);
        }

        // If the base is a QUAD, read the boundary connectivity and convert to SOD2D
        if ( elemTypeStr == "QUAD" )
        {
            // From start and end, get the number of elements
            nbelem = iend - istart + 1;
            if ( mpi_rank == 0 ) printf("  Number of b_elements: %d\n", nbelem);

            // Allocate memory for the boundary connectivity
            connecBound = new cgsize_t[nbelem*nnode];
            connecBoundSOD2D = new cgsize_t[nbelem*nnode];

            // Read the boundary elements
            cg_elements_read( cgns_file, idx_Base, idx_Zone, idx_sec, connecBound, &iparent_data );

            // Conversion to SOD2D format:
            Conversor conv;
            conv.convert2sod_QUAD(porder, nbelem, nnode, connecBound, connecBoundSOD2D);
        }
    }
    */

    // Close the CGNS file
    cgp_close( cgns_file );

    // Writer:

    // Create a HDF5 file with parallel I/O
    hid_t fileOut;
    herr_t status_hdf;
    fileOut = H5Fcreate( output_h5.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    /*
    // Create an empty dataset for boundFaces
    hsize_t dims_bound[2];
    dims_bound[0] = 0;
    dims_bound[1] = 0;

    hid_t dataspace_bound = H5Screate_simple( 2, dims_bound, NULL );
    hid_t dataset_bound = H5Dcreate( fileOut, "/boundFaces", H5T_STD_I64LE, dataspace_bound, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Create an empty dataset for boundFacesId
    hsize_t dims_boundId[1];
    dims_boundId[0] = 0;

    hid_t dataspace_boundId = H5Screate_simple( 1, dims_boundId, NULL );
    hid_t dataset_boundId = H5Dcreate( fileOut, "/boundFacesId", H5T_STD_I64LE, dataspace_boundId, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset_bound );
    status_hdf = H5Sclose( dataspace_bound );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset_bound );
    status_hdf = H5Sclose( dataspace_bound );

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
    if ( status_hdf < 0 )
    {
        if ( mpi_rank == 0 )
        {
            std::cerr << "Error: Cannot write the coordinates" << std::endl;
        }
        MPI_Abort( MPI_COMM_WORLD, 1 );
        return 1;
    }

    // Destroy coord and close the dataset and dataspace
    status_hdf = H5Dclose( dataset_id );
    status_hdf = H5Sclose( dataspace_id );

    // From "/", create the group "/dims"
    hid_t group_dims = H5Gcreate( fileOut, "/dims", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Create the dataset "/dims/numBoundaryFaces"
    hsize_t dims_numBound[1];
    dims_numBound[0] = 1;

    hid_t dataspace_numBound = H5Screate_simple( 1, dims_numBound, NULL );
    hid_t dataset_numBound = H5Dcreate( fileOut, "/dims/numBoundaryFaces", H5T_STD_I64LE, dataspace_numBound, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the number of boundary faces
    int numBound = 0;
    status_hdf = H5Dwrite( dataset_numBound, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &numBound );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset_numBound );
    status_hdf = H5Sclose( dataspace_numBound );

    // Create the dataset "/dims/numElements"
    hsize_t dims_numElem[1];
    dims_numElem[0] = 1;

    hid_t dataspace_numElem = H5Screate_simple( 1, dims_numElem, NULL );
    hid_t dataset_numElem = H5Dcreate( fileOut, "/dims/numElements", H5T_STD_I64LE, dataspace_numElem, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the number of elements
    status_hdf = H5Dwrite( dataset_numElem, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &nelem );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset_numElem );
    status_hdf = H5Sclose( dataspace_numElem );

    // Create the dataset "/dims/numMappedFaces"
    hsize_t dims_numMapFace[1];
    dims_numMapFace[0] = 1;

    hid_t dataspace_numMapFace = H5Screate_simple( 1, dims_numMapFace, NULL );
    hid_t dataset_numMapFace = H5Dcreate( fileOut, "/dims/numMappedFaces", H5T_STD_I64LE, dataspace_numMapFace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the number of mapped faces
    int numMapFace = 0;
    status_hdf = H5Dwrite( dataset_numMapFace, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &numMapFace );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset_numMapFace );
    status_hdf = H5Sclose( dataspace_numMapFace );

    // Create the dataset "/dims/numNodes"
    hsize_t dims_numNodes[1];
    dims_numNodes[0] = 1;

    hid_t dataspace_numNodes = H5Screate_simple( 1, dims_numNodes, NULL );
    hid_t dataset_numNodes = H5Dcreate( fileOut, "/dims/numNodes", H5T_STD_I64LE, dataspace_numNodes, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );

    // Write the number of nodes
    status_hdf = H5Dwrite( dataset_numNodes, H5T_STD_I64LE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &npoin );

    // Close the dataset and dataspace
    status_hdf = H5Dclose( dataset_numNodes );
    status_hdf = H5Sclose( dataspace_numNodes );

    /// Cloose the group "/dims
    status_hdf = H5Gclose( group_dims )
    */

    // Close the HDF5 file
    status_hdf = H5Fclose( fileOut );

    // Finalize MPI
    MPI_Finalize();
    return 0;
}