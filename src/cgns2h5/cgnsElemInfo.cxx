#include "cgnsElemInfo.h"

int* genQuadEdges()
{
    static int quadEdges[4][2] = {
        {0, 1},
        {1, 2},
        {2, 3},
        {3, 0}
    };

    return *quadEdges;
}

int* genHexaEdges()
{
    static int hexaEdges[12][2] = {
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

    return *hexaEdges;
}


/**
 * @brief Get information on a specific CGNS element type
 * 
 */
void getElemInfo( const CGNS_ENUMT(ElementType_t) &eltype, int &numElemNodes, int &elemOrder, int &elemDim )
{
    // Get the element type name into a string
    std::string elname = ElementTypeName[eltype];

    // Strip the base geometrical type from the string
    std::string baseType = elname.substr(0, elname.find("_"));

    // Set element dimension based on the base geometrical type
    // If TETRA, PYRA, PENTA, or HEXA, the element dimension is 3
    if ( baseType == "TETRA" || baseType == "PYRA" || baseType == "PENTA" || baseType == "HEXA" )
    {
        elemDim = 3;
    }
    // If TRI or QUAD, the element dimension is 2
    else if ( baseType == "TRI" || baseType == "QUAD" )
    {
        elemDim = 2;
    }
    // If BAR, the element dimension is 1
    else if ( baseType == "BAR" )
    {
        elemDim = 1;
    }
    // If NODE, the element dimension is 0
    else if ( baseType == "NODE" )
    {
        elemDim = 0;
    }
    // If none of the above, abort
    else
    {
        std::cerr << "Error: Unknown element type: " << baseType << std::endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }

    // Strip the number of nodes from the string
    numElemNodes = std::stoi(elname.substr(elname.find("_") + 1, elname.find("_") + 2));

    // Compute the order accordding to the base geometrical type andd number of nodes
    if ( baseType == "QUAD" ) {
        elemOrder = sqrt(numElemNodes) - 1;
    } else if ( baseType == "HEXA") {
        elemOrder = cbrt(numElemNodes) - 1;
    } else {
        std::cerr << "Error: Element type not supported yet!: " << baseType << std::endl;
        MPI_Abort( MPI_COMM_WORLD, 1 );
    }
}