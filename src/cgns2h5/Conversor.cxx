#include "Conversor.h"

/**
 * @brief Generates the index table for a CGNS quad element
 *
 * @description Recursive routine for creating the index table for a CGNS quad element.
 * If the object 2D array is NULL, allocates memory for it, then proceeds to compute the
 * index table. If the order is greater than 1, the routine is called recursively for the
 * edges and teh faces of the element with mOrder-2.
 *
 * @param mOrder Order of the element
 */
int** Conversor::cgns_createQuadIndexTable( int &mOrder )
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
                int i0 = cgns_QuadEdges[iedge][0];
                int i1 = cgns_QuadEdges[iedge][1];

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
            int **tableFace = cgns_createQuadIndexTable(p);

            // Increment the values in tableFace by 1
            for ( int i = 0; i < p*p; i++ )
            {
                tableFace[i][0]++;
                tableFace[i][1]++;
            }

            // Join tableFace to indexTable
            // TODO: Implement this1
        }
    }
    return indexTable;
}