#ifndef CONVERSOR_H_
#define CONVERSOR_H_

#include <iostream>
#include <vector>
#include <cstdlib>
#include <cstdint>
#include <cstdio>
#include <pcgnslib.h>
#include "cgnsElemInfo.h"

class Conversor
{
    private:
        const int cgns_QuadEdges[4][2] = {
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 0}
        };

        const int cgns_HexaEdges[12][2] = {
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 0},
            {0, 4},
            {1, 5},
            {2, 6},
            {3, 7},
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 4}
        };

        const int cgns_HexaFaces[6][4] = {
            {0, 3, 2, 1},
            {0, 1, 5, 4},
            {1, 2, 6, 5},
            {2, 3, 7, 6},
            {0, 4, 7, 3},
            {4, 5, 6, 7}
        };

        const int sod2d_QuadEdges[4][2] = {
            {0, 1},
            {1, 2},
            {2, 3},
            {3, 0}
        };

        const int sod2d_HexaEdges[12][2] = {
            {0, 1},
            {0, 3},
            {0, 4},
            {1, 2},
            {1, 5},
            {2, 3},
            {2, 6},
            {3, 7},
            {4, 5},
            {4, 7},
            {5, 6},
            {6, 7}
        };

        const int sod2d_HexaFaces[6][4] = {
            {0, 3, 2, 1},
            {0, 1, 5, 4},
            {0, 4, 7, 3},
            {1, 2, 6, 5},
            {2, 3, 7, 6},
            {4, 5, 6, 7}
        };

/*
        const int testGMSH[27] = {
            0,8,1,9,20,11,3,13,2,
            10,21,12,22,26,23,15,24,14,
            4,16,5,17,25,18,7,19,6
        };
*/

        const int testGMSH[64] = {
            0,8,9,1,10,32,35,14,11,33,34,15,3,19,18,2,
            12,36,37,16,40,56,57,44,43,59,58,45,22,49,48,20,
            13,39,38,17,41,60,61,47,42,63,62,46,23,50,51,21,
            4,24,25,5,26,52,53,28,27,55,54,29,7,31,30,6
        };

/*
        const int testCGNS[27] = {
            1,9,2,12,21,10,4,11,3,
            13,22,14,25,27,23,16,24,15,
            5,17,6,20,26,18,8,19,7
        };
*/
        const int testCGNS[64] = {
            1,9,10,2,16,33,34,11,15,36,35,12,4,14,13,3,
            17,37,38,19,50,57,58,41,49,60,59,42,23,46,45,21,
            18,40,39,20,51,61,62,44,52,64,63,43,24,47,48,22,
            5,25,26,6,32,53,54,27,31,56,55,28,8,30,29,7
        };

        int **cgns_QuadIndexTable;
        int **cgns_HexaIndexTable;
        int **cgns_QuadIJ;
        int ***cgns_HexaIJK;
        int **sod2d_QuadIndexTable;
        int **sod2d_HexaIndexTable;
        int **sod2d_QuadIJ;
        int ***sod2d_HexaIJK;
        int** createQuadIndexTable( int &mOrder, int &opt );
        int** createHexaIndexTable( int &mOrder, int &opt );
        int** createQuadIJ( int &mOrder, int** &indexTable);
        int*** createHexaIJK( int &mOrder, int** &indexTable );
        void joinTables( int &indexDesti, int* size1, int** &table1, int** &table2 );
    public:
        // Empty constructor
        Conversor();
        // Destructor
        ~Conversor();
        void convert2sod_QUAD( int &pOrder, uint64_t &nElem, int &nNode, cgsize_t* connecBoundCGNS, cgsize_t* connecBoundSOD2D );
        void convert2sod_HEXA( int &pOrder, uint64_t &nElem, int &nNode, cgsize_t* connecCGNS, cgsize_t* connecSOD2D );
};

#endif // !CONVERSOR_H_