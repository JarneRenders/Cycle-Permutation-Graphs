#ifndef READ_GRAPH6
#define READ_GRAPH6

#include "bitset.h"

//	Returns the number of vertices of a graph in graph6 format.
int getNumberOfVertices(const char * graphString);

//	Loads a graph in graph6 format into an adjacencylist representation
//	consisting of a list of bitsets.
int loadGraph(const char * graphString, int numberOfVertices, bitset
adjacencyList[]);

#endif