#include <stdio.h>
#include <string.h>
#include "readGraph6.h"
#include "bitset.h"

// Return the previous bit of a char. Unsafe because no defined behaviour if
// character = 0. ctz and clz work with 32 bit numbers.
#define unsafePrev(character, current) (__builtin_ctz(character) - (current) >= 0 ? -1 : (current) -__builtin_clz((character) << (32 - current)) - 1)
#define prev(character,current) ((character) ? unsafePrev(character,current) : -1)

int getNumberOfVertices(const char * graphString) {
	if(strlen(graphString) == 0){
        printf("Error: String is empty.\n");
        return -1;
    }
    else if((graphString[0] < 63 || graphString[0] > 126) &&
     graphString[0] != '>') {
    	printf("Error: Invalid start of graphstring.\n");
    	return -1;
    }

	int index = 0;

	// Skip >>graph6<< header.
	if (graphString[index] == '>') { 
		index += 10;
	}

	//	If first character is not 126 its value - 63 is the number of vertices.
	if(graphString[index] < 126) { // 0 <= n <= 62
		return (int) graphString[index] - 63;
	}

	//	If first character is 126 and second not 126. Subtract 63 from the
	//	second, third and fourth characters and concatenate them to get
	//	the 18-bit binary form of the number of vertices.
	else if(graphString[++index] < 126) { // 63 <= n <= 258047
		int number = 0;
		for(int i = 2; i >= 0; i--) {
			number |= (graphString[index++] - 63) << i*6;
		}
		return number;
	}

	//	If both first and second character are 126. Subtract 63 from the
	//	third to eight characters and concatenate them to get the 36-bit
	//	binary form of the number of vertices.
	else if (graphString[++index] < 126) { // 258048 <= n <= 68719476735
		int number = 0;
		for (int i = 5; i >= 0; i--) {
			number |= (graphString[index++] - 63) << i*6;
		}
		return number;
	}

	else {
		fprintf(stderr,
		 "Error: Format only works for graphs up to 68719476735 vertices.\n");
		return -1;
	}
}

int loadGraph(const char * graphString, int numberOfVertices, bitset
adjacencyList[]) {

	//	First position after the information relating to the number of vertices.
	int startIndex = 0;
	if (graphString[startIndex] == '>') { // Skip >>graph6<< header.
		startIndex += 10;
	}
	if (numberOfVertices <= 62) {
		startIndex += 1;
	}
	else if (numberOfVertices <= MAXBITSETSIZE) {
		startIndex += 4;
	}
	//	MAXBITSETSIZE will never get close to 258047.
	else {
		fprintf(stderr,
		 "Error: Program can only handle graphs with %d vertices or fewer.\n",
		 MAXBITSETSIZE);
		return -1;
	}

	// Initialize adjacencyList.
	for (int vertex = 0; vertex < numberOfVertices; vertex++) { 
		adjacencyList[vertex] = EMPTY;
	}

	//	Taking the remaining characters, subtracting by 63 and concatenating
	//	them represents in binary the concatenation of the upper
	//	triangle of the adjacency matrix of the graph. I.e. the first
	//	bit of this concatenation represents (0,1), the second (0,2), the
	//	third (1,2), the fourth (0,3), etc.
	int currentVertex = 1;
	int sum = 0;
	char finalChar = '\n';
	for (int index = startIndex; graphString[index] != '\n' && 
	 (finalChar = graphString[index]) != '\0'; index++) {
		int i;
		for (i = prev(graphString[index] - 63, 6); i != -1;
		 i = prev(graphString[index] - 63, i)) {
			while(5-i+(index-startIndex)*6 - sum >= 0) {
				sum += currentVertex;
				currentVertex++;
			}
			sum -= --currentVertex;
			int neighbour = 5-i+(index - startIndex)*6 - sum;
			add(adjacencyList[currentVertex], neighbour);
			add(adjacencyList[neighbour], currentVertex);
		}
	}
	if(finalChar == '\0') {
		fprintf(stderr,
		 "Error: The g6 string should end with a newline character.\n");
		return -1;
	}
	return 0;
}