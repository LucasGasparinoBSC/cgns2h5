#include "Conversor.h"

// Empty constructor
Conversor::Conversor()
{
}

// Destructor
Conversor::~Conversor()
{
}

/**
 * @brief Generates the index table for a CGNS quad element
 *
 * @description Recursive routine for creating the index table for a HEXA element.
 * If the object 2D array is NULL, allocates memory for it, then proceeds to compute the
 * index table. If the order is greater than 1, the routine is called recursively for the
 * edges and teh faces of the element with mOrder-2.
 *
 * @param mOrder Order of the element
 */
int** Conversor::createHexaIndexTable( int &mOrder, int &opt )
{
    // Allocate memory for the index table
    int nNodes = (mOrder+1) * (mOrder+1) * (mOrder+1);
    int **indexTable = new int*[nNodes];
    for ( int i = 0; i < nNodes; i++ )
    {
        indexTable[i] = new int[3];
    }

    // Node 0
    int iNode = 0;
    indexTable[iNode][0] = 0;
    indexTable[iNode][1] = 0;
    indexTable[iNode][2] = 0;

    // If mOrder > 0, compute corner node indexes
    if ( mOrder > 0 )
    {
        // Node 1
        iNode++;
        indexTable[iNode][0] = mOrder;
        indexTable[iNode][1] = 0;
        indexTable[iNode][2] = 0;

        // Node 2
        iNode++;
        indexTable[iNode][0] = mOrder;
        indexTable[iNode][1] = mOrder;
        indexTable[iNode][2] = 0;

        // Node 3
        iNode++;
        indexTable[iNode][0] = 0;
        indexTable[iNode][1] = mOrder;
        indexTable[iNode][2] = 0;

        // Node 4
        iNode++;
        indexTable[iNode][0] = 0;
        indexTable[iNode][1] = 0;
        indexTable[iNode][2] = mOrder;

        // Node 5
        iNode++;
        indexTable[iNode][0] = mOrder;
        indexTable[iNode][1] = 0;
        indexTable[iNode][2] = mOrder;

        // Node 6
        iNode++;
        indexTable[iNode][0] = mOrder;
        indexTable[iNode][1] = mOrder;
        indexTable[iNode][2] = mOrder;

        // Node 7
        iNode++;
        indexTable[iNode][0] = 0;
        indexTable[iNode][1] = mOrder;
        indexTable[iNode][2] = mOrder;

        // If mOrder > 1, compute high-order nodes
        if ( mOrder > 1)
        {
            int i0, i1, i3;
            // Compute edge indexes
            iNode++;
            for ( int iedge = 0; iedge < 12; iedge++ )
            {
                // Edge end nodes
                switch ( opt )
                {
                    case 1:
                        i0 = cgns_HexaEdges[iedge][0];
                        i1 = cgns_HexaEdges[iedge][1];
                        break;
                    case 2:
                        i0 = sod2d_HexaEdges[iedge][0];
                        i1 = sod2d_HexaEdges[iedge][1];
                        break;
                }

                // Subtract the indexTable of the end nodes, normalize with mOrder
                int u[3];
                u[0] = (indexTable[i1][0] - indexTable[i0][0])/mOrder;
                u[1] = (indexTable[i1][1] - indexTable[i0][1])/mOrder;
                u[2] = (indexTable[i1][2] - indexTable[i0][2])/mOrder;

                // Add all the high-order nodes for this edge
                for ( int i = 1; i < mOrder; i++ )
                {
                    indexTable[iNode][0] = indexTable[i0][0] + i*u[0];
                    indexTable[iNode][1] = indexTable[i0][1] + i*u[1];
                    indexTable[iNode][2] = indexTable[i0][2] + i*u[2];
                    iNode++;
                }
            }

            // Faces: Generate a generic HO face with p=mOrder-2
            int pf = mOrder - 2;
            int nFaceNodes = (pf+1) * (pf+1);
            int **tableFace = createQuadIndexTable( pf, opt );

            // Add 1 to the values in tableFace
            for ( int i = 0; i < nFaceNodes; i++ )
            {
                tableFace[i][0]++;
                tableFace[i][1]++;
            }

            // Create face nodes
            for ( int iFace = 0; iFace < 6; iFace++ )
            {
                // Get corners 0, 1, 3
                switch ( opt )
                {
                    case 1:
                        i0 = cgns_HexaFaces[iFace][0];
                        i1 = cgns_HexaFaces[iFace][1];
                        i3 = cgns_HexaFaces[iFace][3];
                        break;
                    case 2:
                        i0 = sod2d_HexaFaces[iFace][0];
                        i1 = sod2d_HexaFaces[iFace][1];
                        i3 = sod2d_HexaFaces[iFace][3];
                        break;
                }

                // Subtract the indexTable of the end nodes, normalize with mOrder
                int u[3], v[3];
                u[0] = (indexTable[i1][0] - indexTable[i0][0])/mOrder;
                u[1] = (indexTable[i1][1] - indexTable[i0][1])/mOrder;
                u[2] = (indexTable[i1][2] - indexTable[i0][2])/mOrder;

                v[0] = (indexTable[i3][0] - indexTable[i0][0])/mOrder;
                v[1] = (indexTable[i3][1] - indexTable[i0][1])/mOrder;
                v[2] = (indexTable[i3][2] - indexTable[i0][2])/mOrder;

                // Populate the interior of the faces
                for ( int i = 0; i < nFaceNodes; i++ )
                {
                    indexTable[iNode][0] = indexTable[i0][0] + u[0]*tableFace[i][0] + v[0]*tableFace[i][1];
                    indexTable[iNode][1] = indexTable[i0][1] + u[1]*tableFace[i][0] + v[1]*tableFace[i][1];
                    indexTable[iNode][2] = indexTable[i0][2] + u[2]*tableFace[i][0] + v[2]*tableFace[i][1];
                    iNode++;
                }
            }

            // Volume nodes: generate by recursive call
            int p = mOrder - 2;
            int nVolNodes = (p+1) * (p+1) * (p+1);
            int sizeVol[2] = {nVolNodes, 3};
            int **tableVol = createHexaIndexTable( p, opt );

            // Increment the values in tableVol by 1
            for ( int i = 0; i < nVolNodes; i++ )
            {
                tableVol[i][0]++;
                tableVol[i][1]++;
                tableVol[i][2]++;
            }

            // Join tableVol to indexTable
            joinTables( iNode, sizeVol, tableVol, indexTable );
        }
    }
    return indexTable;
}

/**
 * @brief Generates the index table for a CGNS quad element
 *
 * @description Recursive routine for creating the index table for a QUAD element.
 * If the object 2D array is NULL, allocates memory for it, then proceeds to compute the
 * index table. If the order is greater than 1, the routine is called recursively for the
 * edges and teh faces of the element with mOrder-2.
 *
 * @param mOrder Order of the element
 */
int** Conversor::createQuadIndexTable( int &mOrder, int &opt )
{
    // Allocate memory for the index table
    int nNodes = (mOrder+1)*(mOrder+1);
    int **indexTable = new int*[nNodes];
    for ( int i = 0; i < nNodes; i++ )
    {
        indexTable[i] = new int[2];
    }

    // Node 0
    int iNode = 0;
    indexTable[iNode][0] = 0;
    indexTable[iNode][1] = 0;

    // If mOrder > 0, compute corner node indexes
    if ( mOrder > 0 )
    {
        // Node 1
        iNode++;
        indexTable[iNode][0] = mOrder;
        indexTable[iNode][1] = 0;

        // Node 2
        iNode++;
        indexTable[iNode][0] = mOrder;
        indexTable[iNode][1] = mOrder;

        // Node 3
        iNode++;
        indexTable[iNode][0] = 0;
        indexTable[iNode][1] = mOrder;

        // If mOrder > 1, compute high-order nodes
        if ( mOrder > 1)
        {
            // Compute edge indexes
            iNode++;
            for ( int iedge = 0; iedge < 4; iedge++ )
            {
                // Edge end nodes
                int i0, i1;
                switch ( opt )
                {
                    case 1:
                        i0 = cgns_QuadEdges[iedge][0];
                        i1 = cgns_QuadEdges[iedge][1];
                        break;
                    case 2:
                        i0 = sod2d_QuadEdges[iedge][0];
                        i1 = sod2d_QuadEdges[iedge][1];
                        break;
                }

                // Subtract the indexTable of the end nodes, normalize with mOrder
                int u[2];
                u[0] = (indexTable[i1][0] - indexTable[i0][0])/mOrder;
                u[1] = (indexTable[i1][1] - indexTable[i0][1])/mOrder;

                // Add all the high-order nodes for this edge
                for ( int i = 1; i < mOrder; i++ )
                {
                    indexTable[iNode][0] = indexTable[i0][0] + i*u[0];
                    indexTable[iNode][1] = indexTable[i0][1] + i*u[1];
                    iNode++;
                }
            }

            // Compute face nodes by a recursive call with mOrder-2, use a tableFace to store the face nodes
            int p = mOrder - 2;
            int nFaceNodes = (p+1) * (p+1);
            int sizeFace[2] = {nFaceNodes, 2};
            int **tableFace = createQuadIndexTable( p, opt );

            // Increment the values in tableFace by 1
            for ( int i = 0; i < nFaceNodes; i++ )
            {
                tableFace[i][0]++;
                tableFace[i][1]++;
            }

            // Join tableFace to indexTable
            joinTables( iNode, sizeFace, tableFace, indexTable );
        }
    }
    return indexTable;
}

void Conversor::joinTables( int &indexDesti, int* size1, int** &table1, int** &table2 )
{
    int j = indexDesti;
    for ( int i = 0; i < size1[0]; i++ )
    {
        for ( int k = 0; k < size1[1]; k++ )
        {
            table2[j][k] = table1[i][k];
        }
        j++;
    }
}

int** Conversor::createQuadIJ( int &mOrder, int** &indexTable)
{
    // Allocate memory for the tensor
    int **ij = new int*[mOrder+1];
    for ( int i = 0; i < mOrder+1; i++ )
    {
        ij[i] = new int[mOrder+1];
    }

    // Initialize the tensor to zeros
    for ( int i = 0; i < mOrder+1; i++ )
    {
        for ( int j = 0; j < mOrder+1; j++ )
        {
            ij[i][j] = 0;
        }
    }

    // Fill the appropriate values
    for ( int inode = 0; inode < (mOrder+1)*(mOrder+1); inode++ )
    {
        ij[indexTable[inode][0]][indexTable[inode][1]] = inode;
    }

    return ij;
}

int*** Conversor::createHexaIJK( int &mOrder, int** &indexTable )
{
    // Allocate memory for the tensor
    int ***ijk = new int**[mOrder+1];
    for ( int i = 0; i < mOrder+1; i++ )
    {
        ijk[i] = new int*[mOrder+1];
        for ( int j = 0; j < mOrder+1; j++ )
        {
            ijk[i][j] = new int[mOrder+1];
        }
    }

    // Initialize the tensor to zeros
    for ( int i = 0; i < mOrder+1; i++ )
    {
        for ( int j = 0; j < mOrder+1; j++ )
        {
            for ( int k = 0; k < mOrder+1; k++ )
            {
                ijk[i][j][k] = 0;
            }
        }
    }

    // Fill the appropriate values
    for ( int inode = 0; inode < (mOrder+1)*(mOrder+1)*(mOrder+1); inode++ )
    {
        ijk[indexTable[inode][0]][indexTable[inode][1]][indexTable[inode][2]] = inode;
    }

    return ijk;
}

void Conversor::convert2sod_QUAD( int &pOrder, uint64_t &nElem, int &nNode, cgsize_t* connecBoundCGNS, cgsize_t* connecBoundSOD2D )
{
    int opt = 1;
    cgns_QuadIndexTable = createQuadIndexTable( pOrder, opt );
    cgns_QuadIJ = createQuadIJ( pOrder, this->cgns_QuadIndexTable );

    opt = 2;
    sod2d_QuadIndexTable = createQuadIndexTable( pOrder, opt );
    sod2d_QuadIJ = createQuadIJ( pOrder, this->sod2d_QuadIndexTable );

    // Loop over all elements
    cgsize_t *connec = new cgsize_t[nNode];
    cgsize_t *connecConv = new cgsize_t[nNode];
    int i;
    int j;
    int iNodeSOD2D;
    #pragma acc parallel loop gang private(connec, connecConv) present(connecBoundCGNS, connecBoundSOD2D)
    for ( uint64_t iElem = 0; iElem < nElem; iElem++ )
    {
        // Ensure connecConv is zero
        #pragma acc loop vector
        for ( int i = 0; i < nNode; i++ )
        {
            connecConv[i] = 0;
        }
        #pragma acc loop vector
        for ( int iNode = 0; iNode < nNode; iNode++ )
        {
            // Extract local connectivity
            connec[iNode] = connecBoundCGNS[iElem*nNode+iNode];

            // Get the CGNS i,j indices for iNode
            i = cgns_QuadIndexTable[iNode][0];
            j = cgns_QuadIndexTable[iNode][1];

            // Get the corrresponding SOD2D index
            iNodeSOD2D = sod2d_QuadIJ[i][j];

            // Put the node on the new position
            connecConv[iNodeSOD2D] = connec[iNode];
        }
        // Add local converted connectivity to global SOD table
        #pragma acc loop vector
        for ( int iNode = 0; iNode < nNode; iNode++ )
        {
            connecBoundSOD2D[iElem*nNode+iNode] = connecConv[iNode];
        }
    }
    delete[] connec;
    delete[] connecConv;
}

void Conversor::convert2sod_HEXA( int &pOrder, uint64_t &nElem, int &nNode, cgsize_t* connecCGNS, cgsize_t* connecSOD2D )
{
    // Generate the ijk tables for CGNS
    int opt = 1;
    cgns_HexaIndexTable = createHexaIndexTable( pOrder, opt );
    cgns_HexaIJK = createHexaIJK( pOrder, this->cgns_HexaIndexTable );

    // Convert the ijk tables to SOD2D
    opt = 2;
    sod2d_HexaIndexTable =  createHexaIndexTable( pOrder, opt );
    sod2d_HexaIJK = createHexaIJK( pOrder, this->sod2d_HexaIndexTable );

    // Loop over all elements
    cgsize_t *connec = new cgsize_t[nNode];
    cgsize_t *connecConv = new cgsize_t[nNode];
    //int i;
    //int j;
    //int k;
    //int iNodeSOD2D;
    #pragma acc parallel loop gang private(connec, connecConv) present(connecCGNS, connecSOD2D)
    for ( uint64_t iElem = 0; iElem < nElem; iElem++ )
    {
        // Ensure connecConv is zero
        #pragma acc loop vector
        for ( int i = 0; i < nNode; i++ )
        {
            connecConv[i] = 0;
        }

        #pragma acc loop vector
        for ( int iNode = 0; iNode < nNode; iNode++ )
        {
            int kc = testCGNS[iNode]-1;
            int kg = testGMSH[iNode];
            connecConv[kg] = connecCGNS[iElem*nNode+kc];
        }

        // Add local converted connectivity to global SOD table
        #pragma acc loop vector
        for ( int iNode = 0; iNode < nNode; iNode++ )
        {
            connecSOD2D[iElem*nNode+iNode] = connecConv[iNode];
        }
        /*
        #pragma acc loop vector
        for ( int iNode = 0; iNode < nNode; iNode++ )
        {
            // Extract local connectivity
            connec[iNode] = connecCGNS[iElem*nNode+iNode];

            // Get the CGNS i,j,k indices for iNode
            i = cgns_HexaIndexTable[iNode][0];
            j = cgns_HexaIndexTable[iNode][1];
            k = cgns_HexaIndexTable[iNode][2];

            // Get the corrresponding SOD2D index
            iNodeSOD2D = sod2d_HexaIJK[i][j][k];

            // Put the node on the new position
            connecConv[iNodeSOD2D] = connec[iNode];
        }
        */

        // Add local converted connectivity to global SOD table
        #pragma acc loop vector
        for ( int iNode = 0; iNode < nNode; iNode++ )
        {
            connecSOD2D[iElem*nNode+iNode] = connecConv[iNode];
        }
    }
    delete[] connec;
    delete[] connecConv;
}