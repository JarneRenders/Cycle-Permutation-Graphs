/**
 * generateCDC.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 *
 */

#define USAGE "Usage: ./generateCDC [-a|-j] [-hpv]"
#define HELPTEXT "Helptext:\n\
Check various properties involving CDCs (cycle double covers) for cubic\n\
graphs, i.e. whether the graph contains a CDC, a given set of cycles can\n\
be extended to a CDC, any set of pairwise disjoint cycles can be\n\
extended to a CDC (test Jackson's conjecture if the graph is cyclically\n\
5-edge-connected) or simply generate all CDCs.\n\
\n\
Graphs are read from stdin in graph6 format. Graphs are sent to stdout\n\
in graph6 format. If the input graph had a graph6 header, so will the\n\
output graph (if it passes through the filter).\n\
\n\
For checking whether a given set of cycles can be extended to a CDC,\n\
these cycles must also be sent to stdin **before** the graph6 string of\n\
the graph is passed. Their vertices must be separated by a space and\n\
the cycles must be separated by newlines. For example, to check if\n\
`0,2,9,8,3` and `1,4,6,7,5` belong to a CDC in a given labelling of the\n\
Petersen graph, send the following to stdin.\n\
\n\
0 2 9 8 3\n\
1 4 6 7 5\n\
IsP@OkWHG\n\
\n\
If a pair of cycles cannot be extended to a CDC the graph will be sent\n\
to stdout.\n\
\n\
    -a, --all       compute all CDCs of the graph\n\
    -h, --help      print help message\n\
    -j, --jackson   tests whether any set of pairwise disjoint cycles\n\
                     extends to a CDC; graphs containing a set which\n\
                     cannot be extended are sent to stdout; if the graph\n\
                     is cyclically 5-edge-connected and not the Petersen\n\
                     graph this tests a conjecture by Jackson\n\
    -p, --print     prints found CDCs to stderr\n\
    -v, --verbose   sends verbose output to stderr\n\
 "

#include <ctype.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include <string.h>
#include "utilities/readGraph6.h"
#include "utilities/bitset.h"

struct graph {
    int numberOfVertices;
    bitset *adjacencyList;
    int *edgeIndices;
    bitset incidentEdges[MAXBITSETSIZE];
    int indexToEdge[2*MAXBITSETSIZE];
};

struct options {
    bool allFlag;
    bool jacksonFlag;
    bool printFlag;
    bool verboseFlag;
};

struct inputCycles {
    int numCycles;
    int lenCycle[MAXBITSETSIZE];
    int inputCycles[MAXBITSETSIZE][MAXBITSETSIZE];
};

struct counters {
    long long unsigned int CDCs;
    long long unsigned int skippedGraphs;
};

void printBitset(bitset set) {
    forEach(element, set) {
        fprintf(stderr, "%d ", element);
    }
    fprintf(stderr, "\n");
}

void printGraph(bitset adjacencyList[], int numberOfVertices) {
    for(int i = 0; i < numberOfVertices; i++) {
        fprintf(stderr, "%d: ", i);
        forEach(nbr, adjacencyList[i]) {
            fprintf(stderr, "%d ", nbr);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

enum WithNewline {
    WITHNEWLINE,
    WITHOUTNEWLINE
};

void printCycle(struct graph *g, bitset cycle,
 enum WithNewline newline) {
    forEach(el, cycle) {
        fprintf(stderr, "%d--%d ", g->indexToEdge[2*el],
         g->indexToEdge[2*el+1]);
    }
    if(newline == WITHNEWLINE) {
        fprintf(stderr, "\n");
    }
}

void printCycleStdout(struct graph *g, bitset cycle,
 enum WithNewline newline) {
    forEach(el, cycle) {
        fprintf(stdout, "%d--%d ", g->indexToEdge[2*el],
         g->indexToEdge[2*el+1]);
    }
    if(newline == WITHNEWLINE) {
        fprintf(stdout, "\n");
    }
}

int readGraph(const char *graphString, struct graph *g,
 struct options *options, struct counters *counters) {
    int n = getNumberOfVertices(graphString);
    g->numberOfVertices = n; 
    if(n == -1 || n > MAXBITSETSIZE) {
        if(options->verboseFlag){
            fprintf(stderr, "Skipping invalid graph!\n");
        }
        counters->skippedGraphs++;
        return 1;
    }
    g->adjacencyList = malloc(sizeof(bitset)*n);
    if(loadGraph(graphString, n , g->adjacencyList) == -1) {
        if(options->verboseFlag){
            fprintf(stderr, "Skipping invalid graph!\n");
        }
        counters->skippedGraphs++;
        return 1;
    }

    g->edgeIndices = malloc(sizeof(int) * n * n);

    int idx = 0;
    for(int i = 0; i < g->numberOfVertices; i++) {
        forEachAfterIndex(j, g->adjacencyList[i], i) {
            g->edgeIndices[g->numberOfVertices * i + j] = idx;
            g->edgeIndices[g->numberOfVertices * j + i] = idx;
            g->indexToEdge[2*idx] = i;
            g->indexToEdge[2*idx+1] = j;
            add(g->incidentEdges[i], idx);
            add(g->incidentEdges[j], idx);
            idx++;
        }
    }

    return 0;
}

// Assume cycles are a line from stdin of the form
// a1 a2 a3 ... ak (with ai numbers of vertices)
void readCycle(char *str, struct inputCycles *iC) {
    char *p = str;
    while(*p) {
        if(isdigit(*p)) {
            int u = strtol(p, &p, 10);
            int n = iC->numCycles; // Add nth cycle
            if(n >= MAXBITSETSIZE) {
                fprintf(stderr,
                 "Error: overflow reading cycle. numCycles = %d\n", n);
                exit(1);
            }
            int l = iC->lenCycle[n]; // l is its current length
            iC->inputCycles[n][l] = u;
            iC->lenCycle[n]++;
        } else {
            p++;
        }
    }
    iC->numCycles++;
}

void resetInputCycles(struct inputCycles *iC) {
    for(int i = 0; i < iC->numCycles; i++) {
        iC->lenCycle[i] = 0;
    }
    iC->numCycles = 0;
}

void freeGraph(struct graph *g) {
    free(g->edgeIndices);
    free(g->adjacencyList);
}

#define getIndex(g, ep1, ep2)\
 (g)->edgeIndices[(g)->numberOfVertices * (ep1) + (ep2)]

//*****************************************************************
//
//                  Dynamic bitset array
//
//*****************************************************************

typedef struct {
  bitset *array;
  size_t used;
  size_t size;
} Array;

void initArray(Array *a, size_t initialSize) {
  a->array = malloc(initialSize * sizeof(bitset));
  if(a->array == NULL) {
    fprintf(stderr, "%lu\n", initialSize);
    fprintf(stderr, "Error: out of memory\n");
    exit(1);
  }
  a->used = 0;
  a->size = initialSize;
}

void insertArray(Array *a, bitset element) {

  // a->used is the number of used entries, because a->array
  // [a->used++] updates a->used only *after* the array has been
  // accessed. Therefore a->used can go up to a->size 
  if (a->used == a->size) {
    a->size *= 2;
    a->array = realloc(a->array, a->size * sizeof(bitset));
    if(a->array == NULL) {
        fprintf(stderr, "%lu\n", a->size);
        fprintf(stderr, "Error: out of memory at insert\n");
        exit(1);
    }
  }
  a->array[a->used++] = element;
}

bitset popArray(Array *a) {
    return a->array[--a->used];
}

void clearArray(Array *a) {
    a->used = 0;
}

void printCycleArray(struct graph *g, Array *a, char *s) {
    fprintf(stderr, "Printing %s: \n", s);
    for(size_t i = 0; i < a->used; i++) {
        printCycle(g, a->array[i], WITHNEWLINE);
    }
    fprintf(stderr, "__________________________________\n");
}

void printCycleArrayStdout(struct graph *g, Array *a, char *s) {
    fprintf(stdout, "Printing %s: \n", s);
    for(size_t i = 0; i < a->used; i++) {
        printCycleStdout(g, a->array[i], WITHNEWLINE);
    }
    fprintf(stdout, "__________________________________\n");
}

void freeArray(Array *a) {
  free(a->array);
  a->array = NULL;
  a->used = a->size = 0;
}

//**********************************************************************
//
//              Generation of cycles
//
//**********************************************************************

long long unsigned int cyclesRecursion(struct graph *g,
 struct options *options, bitset cycle, bitset cycleEdges,
 bitset *remainingVertices, int start, int current, Array *cycles, 
 bitset excludedEdges) {
    long long unsigned int nCycles = 0;
    removeElement(*remainingVertices, current);
    bool canBeClosed = contains(g->adjacencyList[current], start);
    bool finalEdgeIsExcluded = contains(excludedEdges,
     getIndex(g, current, start));
    if(canBeClosed && !finalEdgeIsExcluded) {
        int closingEdge = getIndex(g, current, start);
        add(cycleEdges, closingEdge);
        // if(options->verboseFlag) {
        //     // print(cycle);
        //     fprintf(stderr, "Cycle_edges: ");
        //     printBitset(cycleEdges);
        //     // fprintf(stderr, "\n");
        // }
        insertArray(cycles, cycleEdges);
        nCycles++;
        removeElement(cycleEdges, closingEdge);
        // fprintf(stderr, "s:%d c:%d ",start, current);
    }
    forEach(nbr, intersection(g->adjacencyList[current], *remainingVertices)) {
        if(contains(excludedEdges, getIndex(g, nbr, current))) {
            continue;
        }
        bitset newCycle = union(cycle, singleton(nbr));
        bitset newCycleEdges = union(cycleEdges,
         singleton(getIndex(g, current, nbr)));
        nCycles += cyclesRecursion(g, options, newCycle, newCycleEdges,
         remainingVertices, start, nbr, cycles, excludedEdges);
    }
    add(*remainingVertices, current);
    return nCycles;
}


// Return number of cycles.
long long unsigned int generateCycles(struct graph *g,
 struct options *options, Array *cycles, bitset excludedEdges) {
    long long unsigned int nCycles = 0;
    bitset remainingVertices = complement(EMPTY, g->numberOfVertices);
    for(int i = 0; i < g->numberOfVertices; i++) {
        removeElement(remainingVertices, i);
        
        // Generate cycles of length > 2 containing i
        forEach(nbr1,
         intersection(g->adjacencyList[i], remainingVertices)) {
            if(contains(excludedEdges, getIndex(g, i, nbr1))) {
                continue;
            }
            removeElement(remainingVertices, nbr1);
            forEachAfterIndex(nbr2,
             intersection(g->adjacencyList[i], remainingVertices), 
             nbr1) {
                if(contains(excludedEdges, getIndex(g, i, nbr2))) {
                    continue;
                }
                bitset cycle = union(union(singleton(i),
                 singleton(nbr1)), singleton(nbr2));
                bitset cycleEdges = union(singleton(
                 getIndex(g, nbr1, i)), singleton(getIndex(g, i, nbr2)));
                nCycles += cyclesRecursion(g, options, cycle,
                 cycleEdges, &remainingVertices, nbr1, nbr2, cycles, 
                 excludedEdges);
            }
            add(remainingVertices, nbr1);
        } 
    }
    return nCycles;
}

// Finds all cycles in the graph and stores them in an array.
// Faster method might exist, but probably not bottleneck.
Array *findAllCycles(struct graph *g, struct options *options) {
    Array *cycles = malloc(sizeof(Array));
    initArray(cycles, 100*g->numberOfVertices);
    generateCycles(g, options, cycles, EMPTY);
    // for(size_t i = 0; i < cycles.used; i++) {
    //     printBitset(cycles.array[i]);
    // }
    if(options->verboseFlag) {
        fprintf(stderr, "Cycles: %zu\n", cycles->used);
    }
    return cycles;
}

bool isCubic(struct graph *g) {
    for(int i = 0; i < g-> numberOfVertices; i++) {
        if(size(g->adjacencyList[i]) != 3) {
            return false;
        }
    }
    return true;
}

//**********************************************************************
//
//              Generation of CDCs
//
//**********************************************************************

// Some checks for pruning immediately
bool canAddCycle(struct graph *g, bitset cycle, bitset edgesInCDC,
 bitset edgesInCDCTwice) {

    // Adding cycle yields edge in CDC thrice
    if(!isEmpty(intersection(cycle, edgesInCDCTwice))) {
        return false;
    }

    // Adding cycle yields two cycles in CDC that pass a vertex in the
    // same way, meaning the third incident edge to that vertex cannot
    // be covered.
    for(int i = 0; i < g->numberOfVertices; i++) {
        bitset consecutiveEdges =
         intersection(g->incidentEdges[i], cycle);

         if(size(intersection(edgesInCDC, consecutiveEdges)) < 2)
          continue; 

        // Can only occur if third edge is in the CDC twice
        int e3 = next(difference(g->incidentEdges[i], consecutiveEdges),
         -1);
        if(!contains(edgesInCDCTwice, e3)) return false;
    }

    return true;
}

bool isCycleDoubleCover(struct graph *g, Array *CDC) {
    bitset once = EMPTY;
    bitset twice = EMPTY;

    // Check if every edge in CDC at most twice.
    for(size_t i = 0; i < CDC->used; i++) {
        bitset cycle = CDC->array[i];
        forEach(el, cycle) {
            if(contains(twice, el)) return false;
            if(contains(once, el)) add(twice, el);
            add(once, el);
        }
    }

    // Check if every edge of the graph is in the CDC twice.
    if(!isEmpty(complement(twice, g->numberOfVertices * 3 / 2))) {
        return false;
    }
    return true;
}

bool extendCDC(struct graph *g, struct options *options,
 struct counters *counters, bitset remainingVertices, Array *cycles,
 Array *CDC, bitset edgesInCDC, bitset edgesInCDCTwice);

// Add all possible cycles going through u to the CDC under the
// condition that exactly one edge incident to u is already in the CDC
// twice. Ignore all cycles up to index s.
bool oneEdgeTwice(struct graph *g, struct options *options,
 struct counters *counters, bitset remainingVertices, Array *cycles, 
 Array *CDC, bitset edgesInCDC, bitset edgesInCDCTwice, int u, 
 size_t s) {

    // printCycleArray(g, CDC, "partial CDC oneEdgeTwice");

    bitset edgesOnce = difference(g->incidentEdges[u], edgesInCDCTwice);
    int e1 = next(edgesOnce, -1);
    int e2 = next(edgesOnce, e1);

    // Cycle must contain both edges which are not yet twice in the CDC.
    for(size_t i = s+1; i < cycles->used; i++) {
        bitset cycle = cycles->array[i];
        if(!contains(cycle, e1) || !contains(cycle, e2)) continue;
        if(!canAddCycle(g, cycle, edgesInCDC, edgesInCDCTwice))
            continue;

        insertArray(CDC, cycle);
        bitset newEdges = union(edgesInCDC, cycle);
        bitset newEdgesTwice = union(edgesInCDCTwice,
         intersection(edgesInCDC, cycle));

        if(extendCDC(g, options, counters, 
         difference(remainingVertices, singleton(u)), cycles, CDC,
         newEdges, newEdgesTwice)) {
            if(!options->allFlag) return true;
        }

        popArray(CDC);
    }
    return false;
}

// Add all possible cycles going through u to the CDC under the
// condition that exactly two edges incident to u are already in the
// CDC and no edges incident with u are in the CDC twice. Ignore all
// cycles up to index s.
bool twoEdgesOnce(struct graph *g, struct options *options,
 struct counters *counters, bitset remainingVertices, Array *cycles, 
 Array *CDC, bitset edgesInCDC, bitset edgesInCDCTwice, int u, 
 size_t s) {

    // printCycleArray(g, CDC, "partial CDC twoEdgesOnce");

    bitset edgesOnce = intersection(g->incidentEdges[u], edgesInCDC);
    int e1 = next(edgesOnce, -1);
    int e2 = next(edgesOnce, e1);

    // Cycle must contain at least one of the edges incident with u but
    // not both which are already in the CDC.
    for(size_t i = s+1; i < cycles->used; i++) {
        bitset cycle = cycles->array[i];
        if(!contains(cycle, e1) && !contains(cycle, e2)) continue;
        if(contains(cycle, e1) && contains(cycle, e2)) continue;
        if(!canAddCycle(g, cycle, edgesInCDC, edgesInCDCTwice))
            continue;


        insertArray(CDC, cycle);
        bitset newEdges = union(edgesInCDC, cycle);
        bitset newEdgesTwice = union(edgesInCDCTwice,
         intersection(edgesInCDC, cycle));

        // Now one edge will be present twice in the CDC
        if(oneEdgeTwice(g, options, counters, remainingVertices, cycles,
         CDC, newEdges, newEdgesTwice, u, i)) {
            if(!options->allFlag) return true;
        }

        popArray(CDC);
    }

    return false;
}


// Add all possible cycles going through u to the CDC under the
// condition that no edges incident to u are already in the CDC. Ignore
// all cycles up to index s.
bool noEdgesOnce(struct graph *g, struct options *options,
 struct counters *counters, bitset remainingVertices, Array *cycles,
 Array *CDC, bitset edgesInCDC, bitset edgesInCDCTwice, int u, 
 size_t s) {

    // printCycleArray(g, CDC, "partial CDC noEdgesOnce");

    // Cycle must contain at least one of the edges incident with u but
    // not both which are already in the CDC.
    for(size_t i = s+1; i < cycles->used; i++) {
        bitset cycle = cycles->array[i];
        if(isEmpty(intersection(cycle, g->incidentEdges[u]))) continue;
        if(!canAddCycle(g, cycle, edgesInCDC, edgesInCDCTwice))
            continue;

        insertArray(CDC, cycle);
        bitset newEdges = union(edgesInCDC, cycle);
        bitset newEdgesTwice = union(edgesInCDCTwice,
         intersection(edgesInCDC, cycle));

        // Now two edges incident with u will be present exactly once in
        // the CDC.
        if(twoEdgesOnce(g, options, counters, remainingVertices, cycles,
         CDC, newEdges, newEdgesTwice, u, i)) {
            if(!options->allFlag) return true;
        }

        popArray(CDC);
    }

    return false;
}


// Only works for cubic graphs. Vertices which are not remaining have
// already three cycles passing through it.
bool extendCDC(struct graph *g, struct options *options,
 struct counters *counters, bitset remainingVertices, Array *cycles,
 Array *CDC, bitset edgesInCDC, bitset edgesInCDCTwice) {

    int u = next(remainingVertices, -1);
    if(u == -1) {

        if(!isCycleDoubleCover(g, CDC)) {
            fprintf(stderr, "Error: not a CDC.\n");
            exit(1);
        }
        if(options->printFlag) {
            printCycleArray(g, CDC, "Complete CDC");
        }

        counters->CDCs++;
        return true;
    }

    // printCycleArray(g, CDC, "partial CDC main recursion");
    // fprintf(stderr, "u: %d\n", u);
    // printCycle(g, g->incidentEdges[u], WITHNEWLINE);

    int edgesTwice = size(
     intersection(g->incidentEdges[u], edgesInCDCTwice));

    // If exactly two incident edges are in the CDC twice, we cannot
    // finish it.
    if(edgesTwice == 2) return false;

    // If all three are already present, we can ignore u.
    if(edgesTwice == 3) {
        removeElement(remainingVertices, u);
        if(extendCDC(g, options, counters, remainingVertices, cycles, CDC,
         edgesInCDC, edgesInCDCTwice)) {
            if(!options->allFlag) return true;
        }
        return false;
    }

    // Other edges are precisely once in CDC
    if(edgesTwice == 1) {

        // Try to add cycles containing both edges incident with u which
        // are only present once and continue recursion.
        if(oneEdgeTwice(g, options, counters, remainingVertices, cycles,
         CDC, edgesInCDC, edgesInCDCTwice, u, -1)) {
            if(!options->allFlag) return true;
        }

        return false;
    }

    // No incident edges of u are present twice already
    int edgesOnce = size(intersection(g->incidentEdges[u], edgesInCDC));

    if(edgesOnce == 3 || edgesOnce == 1) {
        fprintf(stderr, "Error: should not happen.\n");
        exit(1);
    }

    if(edgesOnce == 2) {

        // Try to add cycles containing exactly one edge incident with u
        // and present in the CDC and continue recursion.
        if(twoEdgesOnce(g, options, counters, remainingVertices, cycles,
         CDC, edgesInCDC, edgesInCDCTwice, u, -1)) {
            if(!options->allFlag) return true;
        }

        return false;
    }

    // Here no edges incident with u are present in the CDC.

    // Try to add cycles containing passing through u to the CDC.
    if(noEdgesOnce(g, options, counters, remainingVertices, cycles, CDC,
     edgesInCDC, edgesInCDCTwice, u, -1)) {
        if(!options->allFlag) return true;
    }

    return false;

}

// Find all CDCs given the cycles of the graph.
void findAllCDCs(struct graph *g, struct options *options,
 struct counters *counters, Array *cycles) {

    long long unsigned int CDCsOld = counters->CDCs;

    Array CDC;
    initArray(&CDC, g->numberOfVertices);

    extendCDC(g, options, counters, 
     complement(EMPTY, g->numberOfVertices), cycles, &CDC, EMPTY,
     EMPTY);

    if(options->allFlag && options->verboseFlag) {
        fprintf(stderr, "CDCs: %llu\n", counters->CDCs - CDCsOld);
    }

    freeArray(&CDC);
}

bool hasCDC(struct graph *g, struct options *options,
 struct counters *counters, Array *cycles) {

    Array CDC;
    initArray(&CDC, g->numberOfVertices);

    bool hasOne = extendCDC(g, options, counters, 
     complement(EMPTY, g->numberOfVertices), cycles, &CDC, EMPTY,
     EMPTY);

    if(options->verboseFlag) {
        if(hasOne) {
            printCycleArray(g, &CDC, "CDC");
        }
        else {
            fprintf(stderr, "No CDC found.\n");
        }
    }

    freeArray(&CDC);

    return hasOne;

}

//**********************************************************************
//
//                    Jackson's conjecture
//
//**********************************************************************

bool extendCycles(struct graph *g, struct options *options,
 struct counters *counters, Array *cycles, Array *cycleCollection, 
 size_t s) {

    bool extended = false;

    // Try to add the cycles to the collection pairwise disjoint from
    // the previous ones.
    for(size_t i = s+1; i < cycles->used; i++) {
        bitset C = cycles->array[i];
        bool disjoint = true;
        for(size_t j = 0; j < cycleCollection->used; j++) {
            if(!isEmpty(intersection(C, cycleCollection->array[j]))) {
                disjoint = false;
                break;
            }
        }
        if(!disjoint) continue;

        insertArray(cycleCollection, C);
        if(!extendCycles(g, options, counters, cycles, cycleCollection,
         i)) {
            return false;
        }
        extended = true;
        popArray(cycleCollection);
    }

    // This means we have checked a supercollection of disjoint cycles.
    if(extended) return true;

    // Get all edges in the cycle collection (if bottleneck do
    // dynamically) and make copy of collection
    Array CDC;
    initArray(&CDC, 100);
    bitset once = EMPTY;
    for(size_t i = 0; i < cycleCollection->used; i++) {
        bitset cycle = cycleCollection->array[i];
        once = union(once, cycle);
        insertArray(&CDC, cycle);
    }


    // Check if the collection can be extended to a CDC.
    if(!extendCDC(g, options, counters,
     complement(EMPTY, g->numberOfVertices), cycles, &CDC, once,
     EMPTY)) {
        freeArray(&CDC);
        return false;
    }

    freeArray(&CDC);
    return true;
} 

// Assumes the graph is cyclically 5-edge-connected. Any pair of
// pairwise disjoint cycles of G is a subset of a CDC.
bool satisfiesJacksonsConjecture(struct graph *g,
 struct options *options, struct counters *counters, Array *cycles) {

    Array cycleCollection;
    initArray(&cycleCollection, 100);

    bool satisfies = extendCycles(g, options, counters, cycles,
     &cycleCollection, -1);

    if(!satisfies && options->verboseFlag) {
        printCycleArray(g, &cycleCollection,
         "Cycles that cannot be extended to CDC");
    }

    freeArray(&cycleCollection);

    return satisfies;
}

// Assume that input cycles are permutation 2-factors
bool canExtendInputCyclesToCDC(struct graph *g, struct options *options,
 struct counters *counters, struct inputCycles *iC, Array *cycles) {

    for(int i = 0; i < iC->numCycles; i+=2) {

        // Start partial CDC
        Array CDC;
        initArray(&CDC, 100);
        bitset once = EMPTY;
        bitset twice = EMPTY;

        for(int k = 0; k < 2; k++) {
            int *C = iC->inputCycles[i+k];
            int len = iC->lenCycle[i+k];
            bitset cycleEdges = EMPTY;
            for(int j = 0; j < len; j++) {
                int u = C[j];
                int v = C[(j+1)%len]; 
                if(u >= g->numberOfVertices || // u is vertex
                 v >= g->numberOfVertices || // v is vertex
                 !contains(g->adjacencyList[u], v)) { // uv is edge
                    fprintf(stderr, "Error: incorrect cycle.\n");
                    exit(1); 
                }
                add(cycleEdges, getIndex(g, u, v));
            }
            twice = union(twice, intersection(once, cycleEdges));
            once = union(once, cycleEdges);
            insertArray(&CDC, cycleEdges);
        }

        // extend CDC here
        bool canExtend = extendCDC(g, options, counters,
         complement(EMPTY, g->numberOfVertices), cycles, &CDC, once,
          twice);

        if(!canExtend && options->verboseFlag) {
            printCycleArray(g, &CDC,
             "Cycles that cannot be extended to CDC");
        }

        freeArray(&CDC);
        if(!canExtend) return false;
    }

    return true;

}

int main(int argc, char ** argv) {
    struct counters counters = {0};
    struct options options = {0};
    int opt;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = 
        {
            {"all", no_argument, NULL, 'a'},
            {"help", no_argument, NULL, 'h'},
            {"jackson", no_argument, NULL, 'j'},
            {"print", no_argument, NULL, 'p'},
            {"verbose", no_argument, NULL, 'v'}
        };

        opt = getopt_long(argc, argv, "ahpjv", long_options,
         &option_index);
        if (opt == -1) break;
        switch(opt) {
            case 'a':
                options.allFlag = true; 
                break;
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'j':
                options.jacksonFlag = true;
                fprintf(stderr,
                 "Warning: assumes graphs are cyclically"
                 " 5-edge-connected.\n");
                break;
            case 'p':
                options.printFlag = true;
                break;
            case 'v':
                options.verboseFlag = true;
                break;
            case '?':
                fprintf(stderr,"Error: Unknown option: %c\n", optopt);
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr,
                 "Use ./generateCDC --help for more detailed"
                 " instructions.\n");
                return 1;
        }
    }

    if(options.jacksonFlag && options.allFlag) {
        fprintf(stderr, "Error: do not combine -a and -j.\n");
        exit(1);
    }


    unsigned long long int counter = 0;
    unsigned long long int passedGraphs = 0;
    clock_t start = clock();

    struct inputCycles inputCycles = {0};

    //  Start looping over lines of stdin.
    char * graphString = NULL;
    size_t size;
    while(getline(&graphString, &size, stdin) != -1) {

        // Reads sequence of integers on line before graph6 string as
        // cycle. Cycles are separated by newlines.
        if(graphString[0] >= '0' && graphString[0] <= '9') {
            if(options.jacksonFlag || options.allFlag) {
                fprintf(stderr,
                 "Error: do not send cycles to stdin, when -a or -j is"
                 " present.\n");
                exit(1);
            }

            // Store cycle
            readCycle(graphString, &inputCycles);
            continue;
        }

        struct graph g = {0};
        if(readGraph(graphString, &g, &options, &counters) == 1) {
            fprintf(stderr, "Error: problem loading graph.\n");
            exit(1);
        }

        counter++;
        if(!isCubic(&g)) {
            fprintf(stderr, "Error: graph not cubic\n");
            exit(1);
        }

        if(options.verboseFlag) {
            fprintf(stderr, "Looking at: %s", graphString);
        }

        // Store all cycles of g
        Array *cycles = findAllCycles(&g, &options);

        if(inputCycles.numCycles != 0) {

            // If cycles were sent to stdin, try to extend these to
            // CDC. 
            if(!canExtendInputCyclesToCDC(&g, &options, &counters,
             &inputCycles, cycles)) {
                passedGraphs++;
                printf("%s", graphString);
            }
        }
        else if(options.jacksonFlag) {

            // Build the maximal sets of disjoint cycles and try to
            // extend them
            if(!satisfiesJacksonsConjecture(&g, &options, &counters,
             cycles)) {
                passedGraphs++;
                printf("%s", graphString);
            }
        }
        else if(options.allFlag) {
            findAllCDCs(&g, &options, &counters, cycles);
        }
        else {
            if(hasCDC(&g, &options, &counters, cycles)) {
                passedGraphs++;
                printf("%s", graphString);
            }
        }

        resetInputCycles(&inputCycles);
        freeArray(cycles);
        free(cycles);
        freeGraph(&g);
    }
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    free(graphString);

    fprintf(stderr,"\rChecked %lld graphs in %f seconds: %llu passed.\n",
     counter, time_spent, passedGraphs);
    if(counters.skippedGraphs > 0) {
        fprintf(stderr, "Warning: %lld graphs were skipped.\n",
         counters.skippedGraphs);
    }

    return 0;
}