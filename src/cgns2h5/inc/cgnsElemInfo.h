#ifndef CGNS_ELEM_INFO_H_
#define CGNS_ELEM_INFO_H_

#include <iostream>
#include <string>
#include <cmath>
#include <mpi.h>
#include <pcgnslib.h>

int* genQuadEdges();
int* genHexaEdges();
int* genHexaFaces();
void getElemInfo( const CGNS_ENUMT(ElementType_t) &eltype, int &numElemNodes, int &elemOrder, int &elemDim );

#endif // !CGNS_ELEM_INFO_H_