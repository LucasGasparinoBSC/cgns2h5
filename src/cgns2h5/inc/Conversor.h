#ifndef CONVERSOR_H_
#define CONVERSOR_H_

#include <vector>
#include <cstdlib>
#include <cstdint>
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
            {4, 5},
            {5, 6},
            {6, 7},
            {7, 4},
            {0, 4},
            {1, 5},
            {2, 6},
            {3, 7}
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
        const int sod3d_HexaEdges[12][2] = {
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
        const int sod3d_HexaFaces[6][4] = {
            {0, 3, 2, 1},
            {0, 1, 5, 4},
            {0, 4, 7, 3},
            {1, 2, 6, 5},
            {2, 3, 7, 6},
            {4, 5, 6, 7}
        };
        int **cgns_QuadIndexTable;
        int **cgns_HexaIndexTable;
        int** cgns_createQuadIndexTable(int &mOrder);
    public:
};

#endif // !CONVERSOR_H_