#include "Conversor.h"

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
    int tableSize[2] = {nNodes, 3};
    int **indexTable = new int*[nNodes];
    for ( int i = 0; i < nNodes; i++ )
    {
        indexTable[i] = new int[2];
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
                    case 2:
                        i0 = sod2d_HexaEdges[iedge][0];
                        i1 = sod2d_HexaEdges[iedge][1];
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
                    case 2:
                        i0 = sod2d_HexaFaces[iFace][0];
                        i1 = sod2d_HexaFaces[iFace][1];
                        i3 = sod2d_HexaFaces[iFace][3];
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
                for ( int i = 1; i < nFaceNodes+1; i++ )
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
            joinTables( iNode, sizeVol, tableVol, tableSize, indexTable );
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
    int tableSize[2] = {nNodes, 2};
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
                    case 2:
                        i0 = sod2d_QuadEdges[iedge][0];
                        i1 = sod2d_QuadEdges[iedge][1];
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
            joinTables( iNode, sizeFace, tableFace, tableSize, indexTable );
        }
    }
    return indexTable;
}

void Conversor::joinTables( int &indexDesti, int* size1, int** &table1, int* size2, int** &table2 )
{
    // Print size1
    for ( int i = 0; i < 2; i++ )
    {
        std::cout << size1[i] << " ";
    }

    // Print size2
    for ( int i = 0; i < 2; i++ )
    {
        std::cout << size2[i] << " ";
    }
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