/**
 * genPermutationGraphs.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 * 
 * Generate all pairwise non-isomorphic cycle permutation graphs of a given
 * order n. With -n only non-hamiltonian cycle permutation graphs are
 * generated.
 *
 */


#define USAGE \
"Usage: ./genPermutationGraphs [-hv] [-ng#] n [res/mod]\n"

#define HELPTEXT "\
Helptext: Generate all pairwise non-isomorphic cycle permutation graphs of a\n\
given order n. With -n only non-hamiltonian cycle permutation graphs are\n\
generated. With -n# only cycle permutation graphs of girth at least # are\n\
generated. Generated graphs are sent to stdout in graph6 format. For more\n\
information on the format, see\n\
http://users.cecs.anu.edu.au/~bdm/data/formats.txt.\n\
\n\
The `res/mod` argument, should always appear after the specified order `n`.\n\
Otherwise, the order in which the arguments appear does not matter. Be\n\
careful not to put an argument immediately after one with an option. E.g.\n\
-g#n will not recognize the -n argument.\n\
\n\
Mandatory arguments to long options are mandatory for short options too.\n\
    -h, --help             print help message\n\
    -p, --permutations     use a different generation method; this method WILL\n\
                            print out isomorphic copies, but is faster\n\
    -g, --girth=#          only generate cycle permutation graphs of girth at\n\
                            least #\n\
    -n, --non-hamiltonian  only generate non-hamiltonian cycle permutation\n\
                            graphs\n\
    -v, --verbose          print out extra statistics\n\
    res/mod                split the generation in mod (not necessarily\n\
                            equally large) parts. Here part res will be\n\
                            executed. Splitting will cause extra overhead.\n\
"

#include <stdio.h>
#include <stdbool.h>
#include <getopt.h>
#include <time.h>
#include "utilities/bitset.h" 
#include "utilities/gtools.h" // Includes nauty.h

// Macro for looping over nauty sets.
#define FOREACH(element,nautySet, m)\
    for(int element = nextelement((nautySet),(m),-1); (element) >= 0;\
     (element) = nextelement((nautySet),(m),(element)))

// Struct for keeping a labelling of the edges in a graph. Edges are labelled
// from 0 to ne-1. The endpoints of edge with label/index i are found in
// indexToEdge[2*i] and indexToEdge[2*i+1].
struct edgeLabelling {
    int totalEdges;
    int edgeToIndex[BITSETSIZE][BITSETSIZE];
    int indexToEdge[2*BITSETSIZE*BITSETSIZE];
};

// Struct for keeping information related to the graph. nv is number of
// vertices, nde is the number of directed edges, numberOfWords is the number
// of setwords needed for a nauty set which can hold nv elements, it should be
// ceil(n/WORDSIZE), adjacencyList is the nauty format for storing graphs,
// bitsetList uses the bitset format defined by bitset.h, both should be
// updated simultaneously, calledNauty is a boolean which keeps track of
// whether or not the generators of the automorphism group are up to date, if
// nauty has been called after the last change to the graph, then it will be up
// to date, degree2 and degree3 are bitsets containing the vertices of degree 2
// and of degree 3, respectively. 

// inducedCyclePairs is an array of bitsets. Every entry represents a
// consecutive permutation 2-factor of the graph(see manuscript for its
// definition), but only the induced cycle containing 0 is stored. nCyclePairs
// is the number of consecutive permutation 2-factors. Both are only up to date
// after calling removeOldCyclePairs, updateEnds and addNewCyclePairs in this
// order. SpokesCycle is an array where every entry on index i contains two
// bitsets and corresponds to the consecutive permutation 2-factor at
// inducedCyclePairs[i]. The first bitset contains all degree 3 vertices which
// lie on the induced cycle containing 0, the second bitset contains the
// remaining degree 3 vertices. Similarly, endsCycle has a first and second
// bitset containing the degree 2 vertices which are adjacent to a degree 3
// vertex on respectively, the induced cycle containing 0 and the other induced
// cycle.

// eL contains a labelling of the edges, dynamically updated. stats, options,
// lab, ptn and orbits are variables needed for and giving information about
// the calls to nauty.
typedef struct {
    int nv;
    int nde;
    int numberOfWords;
    graph *adjacencyList;
    bitset *bitsetList;
    bool calledNauty;
    bitset degree2;
    bitset degree3;
    // Pair of induced cycles of length n/2 with the spokes between them
    // consecutive in one of the cycles, vertices of cycle with 0 are
    // stored in a bitset.
    bitset *inducedCyclePairs; 
    // The end vertices of these consecutive spokes. If consecutive on one
    // cycle there will be 2 vertices if consecutive on both cycle there
    // will be 4. If there is only one cycle there will be two vertices in
    // the bitset.
    int nCyclePairs;
    bitset *spokesCycle[2];
    bitset *endsCycle[2];

    // Edge labelling
    struct edgeLabelling *eL; 

    // For calling nauty
    statsblk *stats;
    optionblk *options;
    int *lab;
    int *ptn;
    int *orbits;
} nautygraph;

// Method for freeing the allocated variables of the above struct.
void freeNautyGraph(nautygraph *g) {
    free(g->adjacencyList);
    if(g->bitsetList != NULL) {
        free(g->bitsetList);
        free(g->inducedCyclePairs);
        for(int k = 0; k < 2; k++) {
            free(g->endsCycle[k]); 
            free(g->spokesCycle[k]);
        }
    }
}

// Macro determining whether two vertices are neighbours using nauty format.
#define areNbrs(g, u, v)\
    ISELEMENT(GRAPHROW((g)->adjacencyList, (u), (g)->numberOfWords), (v))

// Macro for looping over all nbrs if a vertex u using nauty format. nbr is the
// variable which contains the neighbour in each iteration.
#define forNbrOf(g, u, nbr)\
    FOREACH((nbr), GRAPHROW((g)->adjacencyList, (u), (g)->numberOfWords),\
    (g)->numberOfWords)

// Macro for adding an edge to the graph, we assume that g is a subcubic graph
// containing a consecutive permutation 2-factor and that the endpoints lie on
// different induced cycles of such a 2-factor. This dynamically updates, the
// edgeLabelling, and the degree3 and degree2 sets.
#define addSpoke(g, i, j) {\
    ADDONEEDGE((g)->adjacencyList, (i), (j),\
    (g)->numberOfWords); (g)->nde += 2;\
    add((g)->bitsetList[i], j); add((g)->bitsetList[j], i);\
    labelEdge((i), (j), (g)->eL);\
    removeElement((g)->degree2, (i));\
    removeElement((g)->degree2, (j));\
    add((g)->degree3, (i));\
    add((g)->degree3, (j));\
}

// Same as previous but for removal.
#define removeSpoke(g, i, j) {\
    DELONEEDGE((g)->adjacencyList, (i), (j),\
    (g)->numberOfWords); (g)->nde -= 2;\
    removeElement((g)->bitsetList[i], j); removeElement((g)->bitsetList[j], i);\
    (g)->eL->totalEdges--;\
    (g)->eL->edgeToIndex[i][j] = 0; (g)->eL->edgeToIndex[j][i] = 0;\
    removeElement((g)->degree3, (i));\
    removeElement((g)->degree3, (j));\
    add((g)->degree2, (i));\
    add((g)->degree2, (j));\
}

// Macro used for building the original graph consisting of two disjoint cycles.
#define addCycleEdges(g, i, lenCycle, offset) {\
    ADDONEARC((g)->adjacencyList, (offset) + (i),\
     (offset) + (((lenCycle) + (i) - 1)%(lenCycle)),\
     SETWORDSNEEDED(2*(lenCycle)));\
    ADDONEARC((g)->adjacencyList, (offset) + (i),\
     (offset) + (((i) + 1)%(lenCycle)),\
     SETWORDSNEEDED(2*(lenCycle)));\
    add((g)->bitsetList[(offset) + (i)],\
     (offset) + (((lenCycle) + (i) - 1)%(lenCycle)));\
    add((g)->bitsetList[(offset) + (((lenCycle) + (i) - 1)%(lenCycle))],\
     (offset) + (i));\
    add((g)->bitsetList[(offset) + (i)], (offset) + (((i) + 1)%(lenCycle)));\
    add((g)->bitsetList[(offset) + (((i) + 1)%(lenCycle))], (offset) + (i));\
}

// Struct containing the options with which the program was called.
struct options {
    int remainder;
    int modulo;
    int splitCounter;
    int splitLevel;
    int minimalGirth;
    bool haveModResPair;
    bool nonHamFlag;
    bool verboseFlag;
    bool permutationMethodFlag;
};

// Many counters for keeping statistics.
#define NUMBEROFINVARIANTS 5
struct counters {
    long long unsigned int recursionSteps;
    long long unsigned int timesPrunedBecauseIsomorphism;
    long long unsigned int callsToCompareMethod;
    long long unsigned int generatedGraphs;
    long long unsigned int numberOfVertexOrbits;
    long long unsigned int sumOfGrpSizes;
    long long unsigned int grpSizeGrTh1;
    long long unsigned int nonIsoIntermediateGraph;
    long long unsigned int pruneGirth;
    long long unsigned int checkHam;
    long long unsigned int checkHamSpokes[BITSETSIZE];
    long long unsigned int pruneHam;
    long long unsigned int pruneHamSpokes[BITSETSIZE];
    long long unsigned int childrenWithoutOrbitRemoval;
    long long unsigned int removedSameOrbit;
    long long unsigned int totalNautyCalls;
    // 0 -> deg2Nbrs, 1-> deg2AtDistAtMost2,
    // 2 -> verticesAtDistAtMost2 + AtMost3 , 3 -> 4Cycles + 5Cycle
     // NINVARIANTS -1 -> highest canonical labelling
    long long unsigned int returnAtInvariant[NUMBEROFINVARIANTS];
    long long unsigned int returnAtInvariantKthSpoke[BITSETSIZE]
     [NUMBEROFINVARIANTS];
    long long unsigned int returnAtInvariantHeuristic[NUMBEROFINVARIANTS];
    long long unsigned int addedSpoke;
    long long unsigned int wasCanonical;
    long long unsigned int addedKthSpoke[BITSETSIZE];
    long long unsigned int calledNauty[BITSETSIZE];
    long long unsigned int canonical[BITSETSIZE];
    long long unsigned int returnBecauseOneEdgeLeft[BITSETSIZE];
    long long unsigned int canonDeterminedByHeuristicTotal;
    long long unsigned int canonDeterminedByHeuristic[BITSETSIZE];
    long long unsigned int needToCheckAllEdgesTotal;
    long long unsigned int needToCheckAllEdges[BITSETSIZE];
    long long unsigned int canonHeurSufficientTotal;
    long long unsigned int canonHeurSufficient[BITSETSIZE];
    long long unsigned int canonHeurNotSufficientTotal;
    long long unsigned int canonHeurNotSufficient[BITSETSIZE];
    long long unsigned int subsetWasAllTotal;
    long long unsigned int subsetWasAll[BITSETSIZE];
    long long unsigned int calledNautyEdgeCanonicalTotal;
    long long unsigned int calledNautyEdgeCanonical[BITSETSIZE];
    long long unsigned int calledNautyEdgeNotCanonicalTotal;
    long long unsigned int calledNautyEdgeNotCanonical[BITSETSIZE];
    long long unsigned int calledNautyForGeneratorsTotal;
    long long unsigned int calledNautyForGenerators[BITSETSIZE];
    long long unsigned int calledAddNewCycles;
    long long unsigned int addNewCyclesAndMoreThanOneDeg2Nbr;
    long long unsigned int didNotComputeGenerators;
    long long unsigned int notRepresentative;
    long long unsigned int expensiveHam;
    long long unsigned int expensiveHamSpokes[BITSETSIZE];
};

// Global variables needed for storing the generators of the automorphism group
// using nauty.
int *generators = NULL;
int numberOfGenerators = 0;




//******************************************************************************
//
//              Various methods for outputting information
//
//******************************************************************************

void printBitsetList(nautygraph *g) {

    fprintf(stderr, "BitsetList on %d vertices: \n", g->nv);
    for(int i = 0; i < g->nv; i++) {
        fprintf(stderr, "%d: ", i);
        forEach(nbr, g->bitsetList[i]) {
            fprintf(stderr, "%d ", nbr);
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

void printGraph(nautygraph *g) {

    fprintf(stderr, "Graph on %d vertices: \n", g->nv);

    for(int i = 0; i < g->nv; i++) {
        fprintf(stderr, "%d: ", i);
        FOREACH(neighbour, GRAPHROW(g->adjacencyList, i, g->numberOfWords),
         g->numberOfWords) {
            fprintf(stderr, " %d", neighbour);
        }
        fprintf(stderr, "\n");
    }
}

void printBitset(bitset set) {
    forEach(element, set) {
        fprintf(stderr, "%d ", element);
    }
    fprintf(stderr, "\n");
}

void writeToG6(nautygraph *g) {
    writeg6(stdout,g->adjacencyList,g->numberOfWords, g->nv);
}

//******************************************************************************
//
//                Initializing program and checking flags
//
//******************************************************************************

bool hasNonOptionArgument(int argc, int *optind) {

    if (*optind >= argc) {
        fprintf(stderr, "Error: add number of vertices.\n");
        fprintf(stderr, "%s", USAGE);
        return false;
    }
    return true;
}

// Only execute once after parsing option arguments.
bool isValidNumberOfVertices(int argc, char** argv, int *optind, int *n) {

    if(!hasNonOptionArgument(argc, optind)) {
        return false;
    }

    char* endptr;
    int numberOfVertices = strtol(argv[*optind], &endptr, 10);

    if(numberOfVertices <= 5 || numberOfVertices > BITSETSIZE) {
        fprintf(stderr, "Error: n needs to be a number between 6 and %d.\n",
         BITSETSIZE);
        fprintf(stderr, "%s",USAGE);
        return false;
    }

    if(numberOfVertices%2) {
        fprintf(stderr, "Error: n needs to be even.\n");
        fprintf(stderr, "%s",USAGE);
        return false;
    }

    *n = numberOfVertices;
    (*optind)++;

    return true;
}


// Only execute once after parsing isValidNumberOfVertices. Checks if the next
// non-option argument is of the form #1/#2. And let res = #1, mod=#2. If no
// non-option arguments were found res = 0, mod = 1.
bool isValidResModPair(int argc, char ** argv, int *optind, 
 struct options *options) {

    // Check if there are still non-option arguments
    while (*optind < argc) {

        char* endptr;

        //  Return false if we already encountered a res/mod pair before and we
        //  are looking at a new non-option argument.
        if(options->haveModResPair) {
            return false;
        }

        options->remainder = strtol(argv[*optind], &endptr, 10);

        // Check if strol succeeded (endptr != 0), i.e. strol found an int, the
        // next character is '/' and that the string containing the argument
        // has a character after that.
        if( !endptr || *endptr != '/' || *(endptr+1) == '\0') {
            return false;
        }

        // Returns false if no integer was found, i.e. if endptr is same as
        // string pointer.
        if(endptr == argv[*optind]) {
            return false;
        }

        options->modulo = strtol(endptr+1, &endptr, 10);

        // Check if strol succeeded and the string containing the argument has
        // no characters left.
        if( !endptr || *endptr != '\0') {
            return false;
        }

        if(options->modulo <= options->remainder) {
            return false;
        }

        options->haveModResPair = true;
        (*optind)++;
    }
    return true;
}

//******************************************************************************
//
//                      Methods dealing with nauty
//
//******************************************************************************

//  Replaces nauty's userautomproc method which gets called for every generator
//  of the automorphism group. generators and numberOfGenerators are global
//  variables.
void saveGenerators(int count, permutation perm[], nvector orbits[],
 int numorbits, int stabvertex, int n) {

    memcpy(generators + n*numberOfGenerators, perm, sizeof (permutation) * n);
    numberOfGenerators++;

}

// Do not forget to free canonForm->adjacencyList and canonForm.  
// Does not initialize the bitsetList.
nautygraph* initCanonForm(nautygraph *g) {

    // Initialize graph struct for storing canonical form.
    graph *gc = malloc(sizeof(graph)*g->nv*g->numberOfWords);
    if(gc == NULL) {
        fprintf(stderr, "Error: out of memory\n");
        exit(EXIT_FAILURE);
    }

    nautygraph *canonGraph = malloc(sizeof *canonGraph);
    if(canonGraph == NULL) {
        fprintf(stderr, "Error: failed allocating canonGraph.\n");
        exit(EXIT_FAILURE);
    }

    canonGraph->nv = g->nv;
    canonGraph->bitsetList = NULL;
    canonGraph->adjacencyList = gc;
    canonGraph->numberOfWords = g->numberOfWords;
    canonGraph->inducedCyclePairs = NULL;

    return canonGraph;
}

// Second option can be NULL if getcanon == false.
void callNauty(nautygraph *g, nautygraph *canonGraph,
 struct counters *counters) {

    if(generators == NULL) {
        generators = malloc(sizeof(int)*g->nv*(g->nv-1));
    }
    numberOfGenerators = 0;

    if(g->options->getcanon == false) {
        densenauty(g->adjacencyList, g->lab, g->ptn, g->orbits, g->options,
         g->stats, g->numberOfWords, g->nv, NULL);
    }
    else {
        densenauty(g->adjacencyList, g->lab, g->ptn, g->orbits, g->options,
         g->stats, g->numberOfWords, g->nv, canonGraph->adjacencyList);
    }

    counters->totalNautyCalls++;
    counters->calledNauty[g->nde/2 - g->nv]++;

    // printGraph(g);
    // fprintf(stderr, "Num gen: %d\n", numberOfGenerators);
    // for(int j = 0; j < numberOfGenerators; j++) {
    //     for(int i = 0; i < g->nv; i++) {
    //         fprintf(stderr, "%d ", generators[g->nv*j + i]);
    //     }
    //     fprintf(stderr, "\n");
    // }
}

// Returns the nth power of 10 up to n = 9.
int getPowersOf10(int n)
{
    static int pow10[10] = {
        1, 10, 100, 1000, 10000, 
        100000, 1000000, 10000000, 100000000, 1000000000
    };

    return pow10[n]; 
}

void updateStats(struct counters *counters, statsblk *stats) {
    counters->numberOfVertexOrbits += stats->numorbits;
    long long unsigned int grpsize =
        stats->grpsize1 * getPowersOf10(stats->grpsize2);
    counters->sumOfGrpSizes += grpsize;
    if (grpsize > 1) {
      counters->grpSizeGrTh1++;
    }
}

//******************************************************************************
//
//                  Dealing with edge labellings
//
//******************************************************************************

// Gives edge uv the next available label/index.
void labelEdge(int u, int v, struct edgeLabelling *eL) {

    //  If edge is not used edgeToIndex should be 0. 
    if(eL->edgeToIndex[u][v] != 0) {
        return;
    } 

    eL->edgeToIndex[u][v] = eL->totalEdges;
    eL->edgeToIndex[v][u] = eL->totalEdges;
    eL->indexToEdge[2*eL->totalEdges] = u;
    eL->indexToEdge[2*eL->totalEdges + 1] = v;
    eL->totalEdges++;

    // The label of the first edge should be 0, but since non-labeled edges have
    // value 0 in edgeToIndex, we temporarily set its value to -1. Call
    // relabelFirstEdge once we are done labelling.
    if(eL->totalEdges == 1) {
        eL->edgeToIndex[u][v] = -1;
        eL->edgeToIndex[v][u] = -1;
    }

}

// Sets the label of the first edge equal to 0.
void relabelFirstEdge(struct edgeLabelling *eL) {

    int u1 = eL->indexToEdge[0];
    int u2 = eL->indexToEdge[1];
    eL->edgeToIndex[u1][u2] = 0;
    eL->edgeToIndex[u2][u1] = 0;
}

//******************************************************************************
//
//                  Union-find by size with path compression
//
//******************************************************************************

// Find with simple path compression (path halving) Make every other node on the
// path between i and the root point to its grandparent. This reduces the
// worst-case time for find.
int find(int sets[], int i) {

    while(i != sets[i]) {
        sets[i] = sets[sets[i]];
        i = sets[i];
    }

    return i;
}

// Union by size
void unionOfSets(int sets[], int setSizes[], int* numberOfSets, int a, int b) {
    int repA = find(sets, a);
    int repB = find(sets, b);

    if(repA != repB) {
        if(setSizes[repA] < setSizes[repB]) {
            sets[repA] = repB;
            setSizes[repB] += setSizes[repA];
        }
        else {
            sets[repB] = repA;
            setSizes[repA] += setSizes[repB];
        }
        (*numberOfSets)--;
    }
}

//******************************************************************************
//
//             Union-find by canonical labelling with path compression
//
//******************************************************************************

// If (u1, v1) < (u2, v2) returns -1, if equal 0 and otherwise 1.
int compareLabelling(int u1, int v1 ,int u2, int v2) {
    if(v1 < u1) {
        int temp = u1;
        u1 = v1;
        v1 = temp;
    }
    if(v2 < u2) {
        int temp = u2;
        u2 = v2;
        v2 = temp;
    }
    if(u1 < u2) {
        return -1;
    }
    else if (u2 < u1) {
        return 1;
    }
    else {
        if(v1 < v2) {
            return -1;
        }
        else if (v2 < v1) {
            return  1;
        }
        else return 0;
    }
}

// Union by size
void unionOfSetsWithLab(int sets[], int* numberOfSets, int a, int b,
 int labelling[], struct edgeLabelling *eL) {
    int repA = find(sets, a);
    int repB = find(sets, b);

    if(repA != repB) {
        int a1 = eL->indexToEdge[2*repA];
        int a2 = eL->indexToEdge[2*repA+1];
        int b1 = eL->indexToEdge[2*repB];
        int b2 = eL->indexToEdge[2*repB+1];
        if(compareLabelling(labelling[a1], labelling[a2], labelling[b1],
         labelling[b2]) < 0) {
            sets[repA] = repB;
        }
        else {
            sets[repB] = repA;
        }
        (*numberOfSets)--;
    }
}

//******************************************************************************
//
//                   Computing edge orbits
//
//******************************************************************************



// Uses union by root with highest labelling in canonical form.
void findEdgeOrbitsWithLab(int nv, struct edgeLabelling *eL, int edgeOrbits[],
 int *numberOfOrbits, int labelling[]) {

    // Initialize array of edge orbits and of orbit sizes.
    for(int i = 0; i < eL->totalEdges; i++) {
        edgeOrbits[i] = i;
    }
    (*numberOfOrbits) = eL->totalEdges;

    // Loop over all generators of the automorphism group. (I.e., permutations
    // of the vertex set.)
    int *perm;

    for(int i = 0; i < numberOfGenerators; i++) {

        // fprintf(stderr, "Generator %d\n", i);

        //  generators is a global variable, used to store all generators given
        //  by nauty.
        perm = &generators[i*nv];

        // Look at the image of every (labeled) edge e under the current
        // generator.
        for(int e = 0; e < eL->totalEdges; e++) {

            // If all edges are already in same orbit we are done.
            if((*numberOfOrbits) == 1) {
                break;
            }

            // Get the image of edge e 
            int u = eL->indexToEdge[2*e];
            int v = eL->indexToEdge[2*e+1];

            int uNew = perm[u];
            int vNew = perm[v];

            if(u == uNew && v == vNew) {
                continue;
            }


            // Merge the orbits of e and its image.
            unionOfSetsWithLab(edgeOrbits, numberOfOrbits, e,
             eL->edgeToIndex[uNew][vNew], labelling, eL);
        }
    }

    //  Point every edge to root of orbit. Is this sometimes not the case?
    for(int e = 0; e < eL->totalEdges; e++) {
        edgeOrbits[e] = find(edgeOrbits, e);
    }
}

// Uses union by size.
void findEdgeOrbits(int nv, struct edgeLabelling *eL, int edgeOrbits[],
 int *numberOfOrbits) {

    // Initialize array of edge orbits and of orbit sizes.
    int orbitSizes[eL->totalEdges];
    for(int i = 0; i < eL->totalEdges; i++) {
        edgeOrbits[i] = i;
        orbitSizes[i] = 1;
    }
    (*numberOfOrbits) = eL->totalEdges;

    // Loop over all generators of the automorphism group. (I.e., permutations
    // of the vertex set.)
    int *perm;

    for(int i = 0; i < numberOfGenerators; i++) {
        // fprintf(stderr, "Generator %d\n", i);

        //  generators is a global variable, used to store all generators given
        //  by nauty.
        perm = &generators[i*nv];

        // Look at the image of every (labeled) edge e under the current
        // generator.
        for(int e = 0; e < eL->totalEdges; e++) {

            // If all edges are already in same orbit we are done.
            if((*numberOfOrbits) == 1) {
                break;
            }

            // Get the image of edge e 
            int u = eL->indexToEdge[2*e];
            int v = eL->indexToEdge[2*e+1];

            int uNew = perm[u];
            int vNew = perm[v];
            if(u == uNew && v == vNew) {
                continue;
            }

            // Merge the orbits of e and its image.
            unionOfSets(edgeOrbits, orbitSizes, numberOfOrbits, e,
             eL->edgeToIndex[uNew][vNew]);
        }
    }

    //  Point every edge to root of orbit. Is this sometimes not the case?
    for(int e = 0; e < eL->totalEdges; e++) {
        edgeOrbits[e] = find(edgeOrbits, e);
    }
}

//******************************************************************************
//
//                          Pruning girth
//
//******************************************************************************

//  Check recursively whether a path can be extended to a cycle of length
//  smaller than minimalGirth.
bool canBeForbiddenCycle(nautygraph *g, int start, int minimalGirth,
 int prev, int curr, int pathLength) {

    if(areNbrs(g, curr, start) &&
        pathLength < minimalGirth) {
        return true;
    }

    if(pathLength + 1 >= minimalGirth) {
        return false;
    }

    forNbrOf(g, curr, el) { 
        if(el == prev) {
            continue;
        }
        if(canBeForbiddenCycle(g, start, minimalGirth, curr, el,
         pathLength + 1)) {
            return true;
        }
    }

    return false;
}

// Check if there exists a cycle of length smaller than minimalGirth
// containing the specified edge.
bool containsForbiddenCycleWithEdge(nautygraph *g, int minimalGirth,
 int u, int v) {

    forNbrOf(g, v, el) {
        if(el == u) {
            continue;
        }
        if(canBeForbiddenCycle(g, u, minimalGirth, v, el, 3)) {
            return true;
        }

    }

    return false;
}

//******************************************************************************
//
//                        Pruning hamiltonian cycles 
//
//******************************************************************************

/**************************************************************************/

// The following code was adapted from cubhamg included in nauty. (This used
// version 2.8.8).

/* cubham.h */

#define MAXNE ((3 * BITSETSIZE) / 2)
#define YES 1
#define DUNNO 0
#define NO (-1)

typedef int cubgraph[BITSETSIZE][4];
typedef int vertvec[BITSETSIZE];
typedef int edgevec[MAXNE+1];

#define POP(x) (onstack[x = *(--stackptr)] = 0)
#define PUSH(x) if(onstack[x]!=stacklev){onstack[x]=stacklev;*(stackptr++)=x;}
#define RESETSTACK {stacklev++; stackptr = stack;}

typedef struct
{
    edgevec class;
    vertvec din,dout,farend;
} nodedata;

static vertvec stack;       /* stack contains vertex numbers */
static int *stackptr,stacklev;      /* stackptr points above top */
static int classstack[4*MAXNE];      /* stack of classifications */
   /* x >= 0        : edge number
     (x < 0 above y : farend[-x-1] = y  */
static int *classstackptr;       /* points above top of classstack */

static int
classout(nodedata *nodat, int v, int w, int en) 
/* classify edge en = vw out */
{
    nodedata *np;

    np = nodat;
    ++np->dout[v];
    ++np->dout[w];
    np->class[en] = NO;
    *classstackptr++ = en; // Put the edge number on the classification stack

    return DUNNO;
}

static int
classin(nautygraph *g, nodedata *nodat,
                                int v, int w, int en, int *nin, int nv, vertvec onstack)
/* classify edge en = vw in */
{
    nodedata *np;
    int *farend,fv,fw,i;

    np = nodat;
    ++np->din[v];
    ++np->din[w];
    np->class[en] = YES;
    *classstackptr++ = en; // Put the edge number on the classification stack

    ++*nin;
    if (*nin == nv) // We have a hamiltonian cycle
    {
        return DUNNO;
    }

    farend = np->farend;
    fv = farend[v]; // fv, fw are old farends of v and w. 
    fw = farend[w];
    *classstackptr++ = farend[fv]; // Put vertex fv on classification stack
    // Put negative index on top to indicate next is vertex
    *classstackptr++ = -fv-1; 
    *classstackptr++ = farend[fw];
    *classstackptr++ = -fw-1;

    farend[fv] = fw;
    farend[fw] = fv;

    if(!contains(g->bitsetList[fv], fw)) {
        return DUNNO;
    }
    i = g->eL->edgeToIndex[fv][fw];
    if (np->class[i] == DUNNO)
    {
        PUSH(fv); // Put vertex fv on vertex number stack
        PUSH(fw);
        if (*nin == nv - 1)
            return classin(g,np,fv,fw,i,nin,nv,onstack);
        else
            return classout(np,fv,fw,i);
    }

    return DUNNO;
}

// long long unsigned int freq[BITSETSIZE] = {0};

// Only speedup in 64 bit case
#ifdef USE_64_BIT
    #define unsafeFirst(set) __builtin_ctzll((set))
    #define safeFirst(set) isEmpty(set) ? -1 : unsafeFirst(set)
    // Set will be changed after calling the forEach
    #define forEachFast(element, set) for (int element = unsafeFirst(set);\
     (set) = ((set) & ((set) - 1)), (element) != -1; (element) = safeFirst(set))
#else
    #define unsafeFirst(set) unsafeNext(set, -1)
    #define safeFirst(set) next(set, -1)
    #define forEachFast(element, set) forEach(element, set)
#endif


static int
propagate(nautygraph *g, nodedata *ndptr, int *nin, int nv, vertvec onstack)
/* propagate classifications: */
/*   ans = YES, NO or DUNNO */
{
    int v,w,status;
    nodedata *np;
    int *class,*din,*dout;

    status = DUNNO;
    np = ndptr;
    class = np->class;
    din = np->din;
    dout = np->dout;

    int (*edgeToIndex)[BITSETSIZE] = g->eL->edgeToIndex;

    w = -1;

    // While not decided and vertex stack is non-empty. Stack contains the
    // vertices which are the endpoints of changed edges?
    while (status == DUNNO && stackptr > stack) 
    {
        POP(v); // Get top of vertex stack

        // If v has no NO edges 
        if (dout[v] == 0)
        {
            if (din[v] == 2) // but v has 2 YES edges.
            {
                // Find the DUNNO edge and make it NO
                bitset vNbrs = g->bitsetList[v];
                forEachFast(el, vNbrs) {
                    if(class[edgeToIndex[v][el]] == DUNNO) {
                        w = el;
                        break;
                    }
                }
                int en = edgeToIndex[v][w];
                status = classout(np,v,w,en);
                // edge of w has changed. v is saturated so no need to put on
                // stack.
                PUSH(w); 
            }
            else if (din[v] == 3) // but v has 3 YES edges.
                status = NO; // Cannot become hamiltonian
        }
        else if (dout[v] == 1) // If v has 1 NO edge:
        {
            // Find the DUNNO edges and make them YES
            bitset vNbrs = g->bitsetList[v];
            forEachFast(w, vNbrs) {
                int en = edgeToIndex[v][w];
                if (class[en] == DUNNO) {
                    if ((status = classin(g,np,v,w,en,nin,nv,onstack)) != DUNNO)
                        break;
                    else
                        PUSH(w); // edge of w has changed and status is DUNNO
                        // This for loop ends with decision or with saturated v.
                }
            }
        }
        else // v has two NO edges. 
            status = NO;
    }

    // Still undecided, but there are nv YES edges, this implies a hamiltonian
    // cycle.
    if (status != NO && *nin == nv) { 
        return YES;
    }
    else {
        return status;
    }
}

static int
hamnode(nautygraph *g, nodedata *nodat, int level, int nin, int nv, vertvec onstack)
/* main node for recursion */
{


    int p,q,status;
    int v,en;
    int *csptr;

    status = propagate(g,nodat,&nin,nv,onstack);

    if (status != DUNNO) {

        return status;
    }

    // For every vertex check the amount of YES neighbours. Take one with one
    // YES neighbour or vertex 0.
    for (v = nv; --v >= 0;)
        if (nodat->din[v] == 1) break;

    // No v with 1 YES nbr, take zero.
    if (v < 0) v = 0;

    struct edgeLabelling *eL = g->eL;

    // Check the neighbours of v.
    bitset vNbrs = g->bitsetList[v];
    forEachFast(w, vNbrs) {
        en = eL->edgeToIndex[v][w]; 

        // v will have one YES edge and two DUNNO.
        if (nodat->class[en] == DUNNO)
        {
            // Classify vw as out
            csptr = classstackptr;
            status = classout(nodat,v,w,en); // Always returns DUNNO.

            RESETSTACK;
            PUSH(v); // Put v on vertex stack
            PUSH(w);

            // Recurse with vw as NO
            status = hamnode(g,nodat,level+1,nin,nv,onstack);
            if (status == YES) break;

            // Reset the stack to before
            while (classstackptr > csptr)
            {
                p = *--classstackptr; // Get top of classification stack
                if (p >= 0) // It indicated the edge with index p
                {
                    // If p was YES, remove YES degree from its endpoints.
                    if (nodat->class[p] == YES) 
                    {
                        --nodat->din[eL->indexToEdge[2*p]];
                        --nodat->din[eL->indexToEdge[2*p+1]];
                    }
                    else // If p was NO, remove NO degree from its endpoints.
                    {
                        --nodat->dout[eL->indexToEdge[2*p]];
                        --nodat->dout[eL->indexToEdge[2*p+1]];
                    }
                    nodat->class[p] = DUNNO; // Mark it DUNNO again.
                }
                else // Points to an index of farend.
                {
                    q = *--classstackptr; // q is the vertex farend used to be.
                    nodat->farend[-p-1] = q;
                }
            }
        }
    }

    if (status == DUNNO)
        fprintf(stderr,"hamnode returning DUNNO, this can't happen\n");

    return status;
}

// Makes this method 64 bit specific
static int farendTemplate[64] = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,
    19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,
    44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63};

// Highly optimized for this specific algorithm version of cubhamg
static int
cubham(nautygraph *gn, edgevec initclass, int setEdges[], int nSetEdges, int en, int v, int w)
/* external interface */
{
    int status,nin;
    int nv = gn->nv;

    // Puts all entries of hcnodat.class to DUNNO.
    // and all entries of din, dout to 0;
    nodedata hcnodat = {0};
    vertvec onstack = {0};

    memcpy(hcnodat.farend, farendTemplate, nv*sizeof(int));

    nin = 0;
    stacklev = 0;
    RESETSTACK;

    bitset deg2 = gn->degree2;
    forEach(el, deg2) {
        hcnodat.dout[el] = 1;
        PUSH(el);
    }

    status = DUNNO;
    classstackptr = classstack;

    // We require presence of vw in hamiltonian cycle.
    // Initializing other edges to YES or NO has no effect or makes it slower.
    status = classin(gn,&hcnodat,v,w,en,&nin,nv,onstack);
    PUSH(v);
    PUSH(w);

    if (status == DUNNO)
        status = hamnode(gn,&hcnodat,0,nin,nv,onstack);

    // if(status == YES) {
    //     freq[0]++;
    //     int crossings = 0;
    //     for(int i = 0; i < nv/2; i++) {
    //         forEachAfterIndex(j, gn->bitsetList[i], nv/2-1) {
    //             int en = gn->eL->edgeToIndex[i][j];
    //             if(hcnodat.class[en] == YES) crossings++;
    //         }
    //     }
    //     ++freq[crossings];
    //     if((crossings == 4 || crossings == 6) ) {
    //      // && gn->nde != 3*nv) {
    //         for(int i = 0; i < nv; i++) {
    //             fprintf(stderr, "%d: ", i);
    //             forEach(j, gn->bitsetList[i]) {
    //                 fprintf(stderr, "%d (%d)", j, hcnodat.class[gn->eL->edgeToIndex[i][j]]);
    //             }
    //             fprintf(stderr, "\n");
    //         }
    //         // writeToG6(gn);
    //         // fprintf(stderr, "Should not\n");
    //     }
    // }

    return status;
}

// Let P be a path such that the degree 3 vertices of one of the cycles of a
// consecutive permutation 2-factor induce P. We assume that u is the endpoint
// of such a path.
bool addingEdgeGivesHamiltonianCycleCubHam(nautygraph *g,
 struct counters *counters, int u, int v) {

    bitset remainingVertices = setComplement(singleton(u), g->nv);
    removeElement(remainingVertices, v);

    // Heuristic for whether hamiltonicity needs to be checked.
    bitset vNbrs = g->bitsetList[v];
    bitset degree2NbrsOfV = intersection(vNbrs, g->degree2);
    if(size(degree2NbrsOfV) == 2) {
        return false;
    }

    counters->expensiveHam++;
    counters->expensiveHamSpokes[g->nde/2 - g->nv]++;

    edgevec initclass = {0};
    int setEdges[g->nde];
    int nSetEdges = 0;

    int en = g->eL->edgeToIndex[u][v];
    initclass[en] = YES;
    setEdges[nSetEdges++] = en; 

    // Setting edges to YES/NO here does not yield any speedup.

    int status = cubham(g, initclass, setEdges, nSetEdges, en, u, v);
    if(status == YES)
        counters->calledNautyEdgeCanonicalTotal++;
    return status == YES;

}

//******************************************************************************
//
//               Finding reducible and eligible edges
//
//******************************************************************************

// Check for consecutive permutation 2-factor i and induced cycle k (0 or 1)
// whether the degree 3 vertices on this cycle induce a path.
bool degree3VerticesInducePath(nautygraph *g, int k, int i) {

    // Ends are the degree 2 vertices adjacent to a degree 3 vertex on the
    // cycle.
    bitset ends = g->endsCycle[k][i];
    int nEnds= size(ends);

    // A cycle in which the degree 3 vertices induce a path has at most two
    // ends.
    if(nEnds > 2) return false;

    // If 0 or 1 ends the degree 3 vertices must induce a path on this cycle.
    if(nEnds < 2) return true;

    // If 2 ends, need to check they are both only adjacent to one degree 3
    // vertex.
    forEach(end, ends) {
        int nAdjSpokes = size(
         intersection(g->bitsetList[end], g->spokesCycle[k][i]));
        if(nAdjSpokes == 2) {
            return false;
        }
    }
    return true;
}

// An edge is reducible if it does not lie on a consecutive permutation 2-factor
// and one of its endpoints is the end vertex of a degree 3 vertex path on a
// cycle. 
void getReducibleEdges(nautygraph *g, struct edgeLabelling *eL,
 int reducibleEdges[], int *nEdges) {

    bitset edges[BITSETSIZE] = {EMPTY};

    // We assume that all nonEmpty permutation 2-factors are consecutive 
    // in either cycle 0 or cycle 1.
    for(int i = 0; i < g->nCyclePairs; i++) {
        if(isEmpty(g->inducedCyclePairs[i])) continue;

        // Check both cycle 0 and cycle 1. If both have consecutive spokes, the
        // first and last spoke need not necessarily be the same edges.
        for(int k = 0; k < 2; k++) {

            if(!degree3VerticesInducePath(g, k, i)) continue; 

            // Ends are the deg 2 vertices on the cycle adjacent to a degree 3
            // vertex.
            bitset ends = g->endsCycle[k][i];
            bitset spokes = g->spokesCycle[k][i];
            bitset spokesOtherCycle = g->spokesCycle[1-k][i];

            // If no ends, all spokes are reducible. (Graph is cubic.)
            if(isEmpty(ends)) {
                forEach(u, spokes) {
                    int v = next(
                     intersection(g->bitsetList[u], spokesOtherCycle), -1);
                    if(contains(edges[u], v)) continue;
                    add(edges[u], v);
                    add(edges[v], u);
                    reducibleEdges[*nEdges] = eL->edgeToIndex[u][v];
                    (*nEdges)++;
                }
            }

            // If there are ends, there is either 1 or 2. Find the degree 3
            // vertices adjacent to them and add these spokes to the edge
            // labelling and the reducible edges (stored in canonical form).
            forEach(end, ends) {
                bitset spokeNbrs = intersection(g->bitsetList[end], spokes);
                forEach(u, spokeNbrs) {
                    int v = next(
                     intersection(g->bitsetList[u], spokesOtherCycle), -1);
                    if(contains(edges[u], v)) continue;
                    add(edges[u], v);
                    add(edges[v], u);
                    reducibleEdges[*nEdges] = eL->edgeToIndex[u][v];
                    (*nEdges)++;
                }
            }
        }
    }
}

// Let P be a path induced by the degree 3 vertices on an induced cycle of a
// consecutive permutation 2-factor. We label all the nonedges of g which have
// an endpoint that is adjacent to an end of such a path P.
void getEligibleEdges(nautygraph *g, struct counters *counters,
 struct options *options, struct edgeLabelling *eL) {

    for(int i = 0; i < g->nCyclePairs; i++) {

        bitset cycle0 = g->inducedCyclePairs[i];
        if(isEmpty(cycle0)) continue;

        bitset cycle1 = setComplement(cycle0, g->nv); 
        bitset cycles[2] = {cycle0, cycle1};

        // Check for cycle 0 and for cycle 1
        for(int k = 0; k < 2; k++) {

            if(!degree3VerticesInducePath(g, k, i)) continue;

            // Ends belong to consecutive spokes
            bitset ends = g->endsCycle[k][i];
            forEach(u, ends) {
                bitset degTwoOppCycle = 
                 difference(cycles[1-k], g->spokesCycle[1-k][i]);
                forEach(v, degTwoOppCycle) {

                    if(options->minimalGirth > 4) {
                        if(containsForbiddenCycleWithEdge(g,
                         options->minimalGirth, u, v)) {
                            counters->pruneGirth++;
                            continue;
                        }
                    }

                    labelEdge(u, v, eL);
                }
            }
        }
    }
    relabelFirstEdge(eL);
}

//******************************************************************************
//
//                      Determining canonical edge
//
//******************************************************************************

// Note that a neighbour cannot occur twice (unless when n is really small)
// Sums the number of degree 2 neighbours of u and of v
int numberOfDeg2Nbrs(nautygraph *g, int u, int v) {
    int uvSum = 0;
    bitset nbrs = g->bitsetList[u];
    forEach(nbr, nbrs) {
        if(size(g->bitsetList[nbr]) == 2) {
            uvSum++;
        }
    } 
    nbrs = g->bitsetList[v];
    forEach(nbr, nbrs) {
        if(size(g->bitsetList[nbr]) == 2) {
            uvSum++;
        }
    } 
    return uvSum;
}


// Returns true if uv is the reducible edge which has the maximum number of
// degree 2 neighbours out of all reducible edges.
bool mostDegree2Nbrs(nautygraph *g, int u, int v, struct edgeLabelling *eL,
 int reducibleEdges[], int *nEdges) {

    int newNumberOfEdges = 0;
    int uvSum = numberOfDeg2Nbrs(g, u, v);

    for(int i = 0; i < *nEdges; i++) {

        int x = eL->indexToEdge[2*reducibleEdges[i]];
        int y = eL->indexToEdge[2*reducibleEdges[i] + 1];
        int xySum = numberOfDeg2Nbrs(g, x, y);

        if(uvSum < xySum) return false;
        if(uvSum > xySum) continue;
        reducibleEdges[newNumberOfEdges] = reducibleEdges[i];
        newNumberOfEdges++;

    }
    *nEdges = newNumberOfEdges;
    return true;
}

// Not that a degree 2 vertex cannot occur twice (unless when is really small)
// Same as above, but for degree 2 vertices at distance at most 2 of u and v.
int numberOfDeg2AtDistanceAtMostTwo(nautygraph *g, int u, int v) {
    int uvSum = 0;
    bitset uNbrs = g->bitsetList[u];
    forEach(nbr, uNbrs) {
        bitset nbrNbrs = g->bitsetList[nbr];
        forEach(nbr2, nbrNbrs) {
            if(size(g->bitsetList[nbr2]) == 2) {
                uvSum++;
            }
        }
    }
    bitset vNbrs = g->bitsetList[v];
    forEach(nbr, vNbrs) {
        bitset nbrNbrs = g->bitsetList[nbr];
        forEach(nbr2, nbrNbrs) {
            if(size(g->bitsetList[nbr2]) == 2) {
                uvSum++;
            }
        }
    }
    return uvSum;
}

bool fewestDegree2AtDistance2(nautygraph *g, int u, int v,
 struct edgeLabelling *eL, int reducibleEdges[], int *nEdges) {

    int newNumberOfEdges = 0;
    int uvSum = numberOfDeg2AtDistanceAtMostTwo(g, u, v);

    for(int i = 0; i < *nEdges; i++) {

        int x = eL->indexToEdge[2*reducibleEdges[i]];
        int y = eL->indexToEdge[2*reducibleEdges[i] + 1];
        int xySum = numberOfDeg2AtDistanceAtMostTwo(g, x, y);

        if(uvSum > xySum) return false;
        if(uvSum < xySum) continue;
        reducibleEdges[newNumberOfEdges] = reducibleEdges[i];
        newNumberOfEdges++;

    }
    *nEdges = newNumberOfEdges;
    return true;
}

// Sums the number of vertices at distance at most 2 of u and v.
int numberOfNbrsAtDistanceAtMost2(nautygraph *g, int u, int v) {
    bitset vertices = EMPTY;
    bitset nbrs = g->bitsetList[u];
    forEach(nbr, nbrs) {
        vertices = union(vertices, g->bitsetList[nbr]);
    }
    bitset vNbrs = g->bitsetList[v];
    forEach(nbr, vNbrs) {
        vertices = union(vertices, g->bitsetList[nbr]);
    }
    return size(vertices);
}

bool fewestNbrsAtDistanceAtMost2(nautygraph *g, int u, int v,
 struct edgeLabelling *eL, int reducibleEdges[], int *nEdges) {

    int newNumberOfEdges = 0;
    int uvSum = numberOfNbrsAtDistanceAtMost2(g, u, v);

    for(int i = 0; i < *nEdges; i++) {

        int x = eL->indexToEdge[2*reducibleEdges[i]];
        int y = eL->indexToEdge[2*reducibleEdges[i] + 1];
        int xySum = numberOfNbrsAtDistanceAtMost2(g, x, y);

        if(uvSum > xySum) return false;
        if(uvSum < xySum) continue;
        reducibleEdges[newNumberOfEdges] = reducibleEdges[i];
        newNumberOfEdges++;

    }
    *nEdges = newNumberOfEdges;
    return true;
}

int numberOfNbrsAtDistanceAtMost3(nautygraph *g, int u, int v) {
    bitset vertices = EMPTY;
    bitset uNbrsNoV = difference(g->bitsetList[u], singleton(v));
    forEach(nbr, uNbrsNoV) {
        bitset nbrNbrs = g->bitsetList[nbr];
        vertices = union(vertices, nbrNbrs);
        forEach(nbr2, nbrNbrs) {
            vertices = union(vertices, g->bitsetList[nbr2]);
        }
    }
    bitset vNbrsNoU = difference(g->bitsetList[v], singleton(u));
    forEach(nbr, vNbrsNoU) {
        bitset nbrNbrs = g->bitsetList[nbr];
        vertices = union(vertices, nbrNbrs);
        forEach(nbr2, nbrNbrs) {
            vertices = union(vertices, g->bitsetList[nbr2]);
        }
    }
    return size(vertices);
}

bool mostNbrsAtDistanceAtMost3(nautygraph *g, int u, int v,
 struct edgeLabelling *eL, int reducibleEdges[], int *nEdges) {

    int newNumberOfEdges = 0;
    int uvSum = numberOfNbrsAtDistanceAtMost3(g, u, v);

    for(int i = 0; i < *nEdges; i++) {

        int x = eL->indexToEdge[2*reducibleEdges[i]];
        int y = eL->indexToEdge[2*reducibleEdges[i] + 1];
        int xySum = numberOfNbrsAtDistanceAtMost3(g, x, y);

        if(uvSum < xySum) return false;
        if(uvSum > xySum) continue;
        reducibleEdges[newNumberOfEdges] = reducibleEdges[i];
        newNumberOfEdges++;

    }
    *nEdges = newNumberOfEdges;
    return true;
}

int numberOfNbrsAtDistanceAtMost4(nautygraph *g, int u, int v) {
    bitset vertices = EMPTY;
    bitset uNbrsNoV = difference(g->bitsetList[u], singleton(v));
    forEach(nbr, uNbrsNoV) {
        bitset nbrNbrs = g->bitsetList[nbr];
        vertices = union(vertices, nbrNbrs);
        forEach(nbr2, nbrNbrs) {
            bitset nbrNbrNbrs = g->bitsetList[nbr2];
            vertices = union(vertices, nbrNbrNbrs);
            forEach(nbr3, nbrNbrNbrs) {
                vertices = union(vertices, g->bitsetList[nbr3]);
            }
        }
    }
    bitset vNbrsNoU = difference(g->bitsetList[v], singleton(u));
    forEach(nbr, vNbrsNoU) {
        bitset nbrNbrs = g->bitsetList[nbr];
        vertices = union(vertices, nbrNbrs);
        forEach(nbr2, nbrNbrs) {
            bitset nbrNbrNbrs = g->bitsetList[nbr2];
            vertices = union(vertices, nbrNbrNbrs);
            forEach(nbr3, nbrNbrNbrs) {
                vertices = union(vertices, g->bitsetList[nbr3]);
            }
        }
    }
    return size(vertices);
}

bool mostNbrsAtDistanceAtMost4(nautygraph *g, int u, int v,
 struct edgeLabelling *eL, int reducibleEdges[], int *nEdges) {

    int newNumberOfEdges = 0;
    int uvSum = numberOfNbrsAtDistanceAtMost4(g, u, v);

    for(int i = 0; i < *nEdges; i++) {

        int x = eL->indexToEdge[2*reducibleEdges[i]];
        int y = eL->indexToEdge[2*reducibleEdges[i] + 1];
        int xySum = numberOfNbrsAtDistanceAtMost4(g, x, y);

        if(uvSum < xySum) return false;
        if(uvSum > xySum) continue;
        reducibleEdges[newNumberOfEdges] = reducibleEdges[i];
        newNumberOfEdges++;

    }
    *nEdges = newNumberOfEdges;
    return true;
}

// Counts the number of 4-cycles containing uv.
int numberOf4CyclesContainingEdge(nautygraph *g, int u, int v) {
    int sum = 0;
    bitset uNbrs = g->bitsetList[u];
    bitset vNbrsNoU = difference(g->bitsetList[v], singleton(u));
    forEach(nbr, uNbrs) {
        if(nbr == v) continue;
        bitset nbrNbrs = g->bitsetList[nbr];
        if(!isEmpty(intersection(nbrNbrs, vNbrsNoU))) {
            sum++;
        }
    }

    return sum;
}

// Returns true if uv has the fewest number of 4-cycles containing it out of all
// reducible edges.
bool fewest4Cycles(nautygraph *g, int u, int v,
 struct edgeLabelling *eL, int reducibleEdges[], int *nEdges) {

    int newNumberOfEdges = 0;
    int uvSum = numberOf4CyclesContainingEdge(g, u, v);

    for(int i = 0; i < *nEdges; i++) {

        int x = eL->indexToEdge[2*reducibleEdges[i]];
        int y = eL->indexToEdge[2*reducibleEdges[i] + 1];
        int xySum = numberOf4CyclesContainingEdge(g, x, y);

        if(uvSum > xySum) return false;
        if(uvSum < xySum) continue;
        reducibleEdges[newNumberOfEdges] = reducibleEdges[i];
        newNumberOfEdges++;

    }
    *nEdges = newNumberOfEdges;
    return true;
}

int numberOf5CyclesContainingEdge(nautygraph *g, int u, int v) {
    int sum = 0;
    bitset uNbrsNoV = difference(g->bitsetList[u], singleton(v));
    bitset vNbrsNoU = difference(g->bitsetList[v], singleton(u));
    forEach(uNbr, uNbrsNoV) {
        bitset nbrsUNbr = g->bitsetList[uNbr];
        forEach(vNbr, vNbrsNoU) {
            if(!isEmpty(intersection(nbrsUNbr, g->bitsetList[vNbr]))) {
                sum++;
            }
        }
    }
    return sum;
}

bool fewest5Cycles(nautygraph *g, int u, int v,
 struct edgeLabelling *eL, int reducibleEdges[], int *nEdges) {

    int newNumberOfEdges = 0;
    int uvSum = numberOf5CyclesContainingEdge(g, u, v);

    for(int i = 0; i < *nEdges; i++) {

        int x = eL->indexToEdge[2*reducibleEdges[i]];
        int y = eL->indexToEdge[2*reducibleEdges[i] + 1];
        int xySum = numberOf5CyclesContainingEdge(g, x, y);

        if(uvSum > xySum) return false;
        if(uvSum < xySum) continue;
        reducibleEdges[newNumberOfEdges] = reducibleEdges[i];
        newNumberOfEdges++;

    }
    *nEdges = newNumberOfEdges;
    return true;
}

// Recursive method for finding cycles of a given length.
void countGivenLengthRecursive(nautygraph *g, int requiredLength, int start,
 int last, bitset remainingVertices, int length, int *sum) {

    bitset nbrs = g->bitsetList[last];
    if(length == requiredLength && contains(nbrs, start)) {
        (*sum)++;
        return;
    }

    if(length >= requiredLength) return;

    bitset remainingNbrs = intersection(nbrs, remainingVertices);
    forEach(nbr, remainingNbrs) {
        countGivenLengthRecursive(g, requiredLength, start, nbr,
        difference(remainingVertices, singleton(nbr)), length + 1, sum);
    }
}

// Counts the number of cycles containing uv of a specific given length.
int numberOfCycles(nautygraph *g, int length, int u, int v) {
    int sum = 0;

    bitset uNbrsNoV = difference(g->bitsetList[u], singleton(v));
    bitset vNbrsNoU = difference(g->bitsetList[v], singleton(u));

    bitset remainingVertices = setComplement(union(singleton(u), singleton(v)),
     g->nv);

    forEach(uNbr, uNbrsNoV) {
        removeElement(remainingVertices, uNbr);
        forEach(vNbr, vNbrsNoU) {
            removeElement(remainingVertices, vNbr);
            countGivenLengthRecursive(g, length, uNbr, vNbr, remainingVertices,
             4, &sum);
            add(remainingVertices, vNbr);
        }
        add(remainingVertices, uNbr);
    }

    return sum;
}

// Returns true if uv is the reducible edge with fewest cycles of a given length
// containing it.
bool fewestCyclesOfLength(nautygraph *g, int length, int u, int v,
 struct edgeLabelling *eL, int reducibleEdges[], int *nEdges) {

    int newNumberOfEdges = 0;
    int uvSum = numberOfCycles(g, length, u, v);

    for(int i = 0; i < *nEdges; i++) {

        int x = eL->indexToEdge[2*reducibleEdges[i]];
        int y = eL->indexToEdge[2*reducibleEdges[i] + 1];
        int xySum = numberOfCycles(g, length, x, y);

        if(uvSum > xySum) return false;
        if(uvSum < xySum) continue;
        reducibleEdges[newNumberOfEdges] = reducibleEdges[i];
        newNumberOfEdges++;

    }
    *nEdges = newNumberOfEdges;
    return true;
}

// Checks all of the heuristics for the reducible edge uv in order to determine
// whether or not it is canonical. Reducible edges which do not attain the
// maximum or minimum for the checked invariant are removed from the list as
// they cannot be canonical. The order of the heuristics was empirically
// determined to be (one of) the most efficient.
bool passesHeuristics(nautygraph *g, struct counters *counters,
 struct options *options, int u, int v, struct edgeLabelling *eL,
 int reducibleEdges[], int *nEdges) {


    if(!fewestNbrsAtDistanceAtMost2(g, u, v, eL, reducibleEdges, nEdges)) {
        counters->returnAtInvariantKthSpoke[g->nde/2-g->nv][2]++;
        counters->returnAtInvariant[2]++;
        return false;
    }

    // If one edge left, it is uv and it must be canonical.
    if(*nEdges == 1) {
     return true;
    } 

    if(options->minimalGirth <= 4 && 
     !fewest4Cycles(g, u, v, eL, reducibleEdges, nEdges)) {
        counters->returnAtInvariantKthSpoke[g->nde/2-g->nv][3]++;
        counters->returnAtInvariant[3]++;
        return false;
    }

    // If one edge left, it is uv and it must be canonical.
    if(*nEdges == 1) {
     return true;
    } 

    if(options->minimalGirth <= 5 && 
     !fewest5Cycles(g, u, v, eL, reducibleEdges, nEdges)) {
        counters->returnAtInvariantKthSpoke[g->nde/2-g->nv][3]++;
        counters->returnAtInvariant[3]++;
        return false;
    }

    // If one edge left, it is uv and it must be canonical.
    if(*nEdges == 1) {
     return true;
    } 

    if(!mostNbrsAtDistanceAtMost3(g, u, v, eL, reducibleEdges, nEdges)) {
        counters->returnAtInvariantKthSpoke[g->nde/2-g->nv][2]++;
        counters->returnAtInvariant[2]++;
        return false;
    }

    // If one edge left, it is uv and it must be canonical.
    if(*nEdges == 1) {
     return true;
    } 

    // Checking 7-cycles Slows down code for girth 7
    if(options->minimalGirth <= 6) { 
        if(!fewestCyclesOfLength(g, 6, u, v, eL, reducibleEdges, nEdges)) {
            counters->returnAtInvariant[3]++;
            counters->returnAtInvariantKthSpoke[g->nde/2-g->nv][3]++;
            return false;
        }
    }

    // The canonical edge has the most nbrs of degree 2. Check if uv is in this
    // category and remove edges with fewer from reducibleEdges.
    if(options->minimalGirth < 5) { // Always 0 otherwise
        if(!mostDegree2Nbrs(g, u, v, eL, reducibleEdges, nEdges)) {
            counters->returnAtInvariant[0]++;
            counters->returnAtInvariantKthSpoke[g->nde/2-g->nv][0]++;
            return false;
        }
    }

    // If one edge left, it is uv and it must be canonical.
    if(*nEdges == 1) {
     return true;
    }

    if(!fewestDegree2AtDistance2(g, u, v, eL, reducibleEdges, nEdges)) {
        counters->returnAtInvariantKthSpoke[g->nde/2-g->nv][1]++;
        counters->returnAtInvariant[1]++;
        return false;
    }


    // If one edge left, it is uv and it must be canonical.
    if(*nEdges == 1) {
     return true;
    } 

    if(!mostNbrsAtDistanceAtMost4(g, u, v, eL, reducibleEdges, nEdges)) {
        counters->returnAtInvariantKthSpoke[g->nde/2-g->nv][2]++;
        counters->returnAtInvariant[2]++;
        return false;
    }

    return true;
}

// Checks whether uv lies in the orbit of the edge with the lexicographically
// highest labelling in canonical form.
bool hasHighestCanonicalLabelling(nautygraph *g, struct counters *counters,
 struct edgeLabelling *eL, int u, int v, int reducibleEdges[], int nEdges) {

    /**
     * lab is the order in which the vertices of the original graph should be
     * relabeled in order to obtain the canonical graph. So lab[k] is the
     * original label of vertex k
     */
    int labelling[g->nv];
    for(int i = 0; i < g->nv; i++) {
        labelling[g->lab[i]] = i;
    }

    // Compute the edge orbits of the reducible edges. The representative of
    // each orbit is the edge with the lexicographically highest labelling in
    // canonical form.
    int edgeOrbits[eL->totalEdges];
    int numberOfOrbits;
    findEdgeOrbitsWithLab(g->nv, eL, edgeOrbits, &numberOfOrbits, labelling); 

    // The edge in uv's orbit with highest labelling in canon form.
    int e1 = edgeOrbits[eL->edgeToIndex[u][v]];
    int x1 = labelling[eL->indexToEdge[2*e1]];
    int y1 = labelling[eL->indexToEdge[2*e1 + 1]];

    // Compare all representatives of reducible edges orbits to representative
    // of uv orbit.
    for(int e = 0; e < nEdges; e++) {

        int e2 = edgeOrbits[reducibleEdges[e]];
        if(e1 == e2) continue;

        int x2 = labelling[eL->indexToEdge[2*e2]];
        int y2 = labelling[eL->indexToEdge[2*e2 + 1]];

        if(compareLabelling(x2, y2, x1, y1) > 0) {
            counters->returnAtInvariantKthSpoke[g->nde/2-g->nv]
             [NUMBEROFINVARIANTS - 1]++;
            counters->returnAtInvariant[NUMBEROFINVARIANTS - 1]++;
            return false; 
        }
    }
    return true;
}

// uv is canonical if it lies in the same orbit as the edge with the
// lexicographically highest labelling in canonical form which passes all
// heuristics. 
bool isCanonical(nautygraph *g, struct counters *counters,
 struct options *options, struct edgeLabelling *eL, int u, int v) {

    counters->needToCheckAllEdgesTotal++;
    counters->needToCheckAllEdges[g->nde/2-g->nv]++;

    // Store all reducible edges only works correctly if permutation 2-factors
    // are up to date.
    int nEdges = 0;
    int reducibleEdges[g->nde/2];
    getReducibleEdges(g, eL, reducibleEdges, &nEdges);

    if(!passesHeuristics(g, counters, options, u, v, eL, reducibleEdges,
     &nEdges)) {
        return false;
    }

    // If there is only one reducible edge left, then it must be uv and hence it
    // is canonical.
    if(nEdges == 1) {
        counters->returnBecauseOneEdgeLeft[g->nde/2-g->nv]++;
        counters->canonDeterminedByHeuristicTotal++;
        counters->canonDeterminedByHeuristic[g->nde/2-g->nv]++;
        return true;
    }

    // If heuristics do not discriminate canonical edge, we need to obtain the
    // canonical form using nauty, unless we called it earlier. The mapping
    // from the canonical form will be/is already stored in g->lab.
    g->options->getcanon = true;
    nautygraph *canonGraph = initCanonForm(g);
    if(!g->calledNauty) {
        callNauty(g, canonGraph, counters);
        updateStats(counters, g->stats);
        g->calledNauty = true;
    }

    if(!hasHighestCanonicalLabelling(g, counters, eL, u, v, reducibleEdges,
     nEdges)) {
        freeNautyGraph(canonGraph);
        free(canonGraph);
        return false;
    }

    freeNautyGraph(canonGraph);
    free(canonGraph);
    return true;
}

// Same as above, but we call this when the consecutive permutation 2-factors
// are not up to date and hence, the set of reducible edges we find is a subset
// of the complete set of reducible edges. Only works if we removed all non-
// (consecutive permutation 2-factors) from the list.
bool isCanonicalHeur(nautygraph *g, struct counters *counters,
 struct options *options, struct edgeLabelling *eL, int u, int v) {

    int nEdges = 0;
    int reducibleEdges[g->nde/2];
    getReducibleEdges(g, eL, reducibleEdges, &nEdges);
    
    if(!passesHeuristics(g, counters, options, u, v, eL, reducibleEdges,
     &nEdges)) {
        return false;
    }

    if(nEdges == 1) {
        counters->returnBecauseOneEdgeLeft[g->nde/2-g->nv]++;
        return true;
    }

    // Get canonical labelling of extension
    g->options->getcanon = true;
    nautygraph *canonGraph = initCanonForm(g);
    if(!g->calledNauty) {
        callNauty(g, canonGraph, counters);
        updateStats(counters, g->stats);
        g->calledNauty = true;
    }

    if(!hasHighestCanonicalLabelling(g, counters, eL, u, v, reducibleEdges,
     nEdges)) {
        freeNautyGraph(canonGraph);
        free(canonGraph);
        return false;
    }

    freeNautyGraph(canonGraph);
    free(canonGraph);
    return true;
}

//******************************************************************************
//
//                  Dealing with consecutive permutation 2-factors
//
//******************************************************************************

// Given the induced cycle containing 0 (pair) of a consecutive permutation
// 2-factor, together with the spokevertices and the ends, store it in
// g->inducedCyclePairs[], g->spokesCycle[], g->endsCycle[], respectively, if
// it is not already present. 
// Spokescycle is an array of two bitset, the first one contains the degree 3
// vertices of g which lie on the cycle containing 0, the second bitset
// contains the remaining degree 3 vertices. EndsCycle is also an array of two
// bitsets, the first one containing the degree 2 vertices adjacent to a degree
// 3 vertex which lie on the induced cycle containing 0, the second bitset
// contains the remaining such degree 2 vertices.
void insertInducedPair(nautygraph *g, bitset pair, bitset spokesCycle[],
 bitset endsCycle[]) {

    if(isEmpty(pair)) {
        return;
    }

    // Find out whether the permutation 2-factor is already stored, or if we
    // have an empty entry in our list.
    int emptyIdx = -1;
    bitset *inducedCyclePairs = g->inducedCyclePairs;
    for(int i = 0; i < g->nCyclePairs; i++) {
        if(emptyIdx == -1 && isEmpty(inducedCyclePairs[i])) {
            emptyIdx = i;
        }
        if(equals(pair, inducedCyclePairs[i])) {
            return;
        }
    }

    // If an entry was empty, replace it with the new cons. permutation
    // 2-factor.
    if(emptyIdx != -1) {
        inducedCyclePairs[emptyIdx] = pair;
        for(int k = 0; k < 2; k++) {
            g->endsCycle[k][emptyIdx] = endsCycle[k];
            g->spokesCycle[k][emptyIdx] = spokesCycle[k];
        }
    }
    else {

        // Some bound we use on the max number of allowed consecutive
        // permutation 2-factors. Not proven to be an actual upper bound, but
        // is not exceeded in any of the computations we performed.
        int cycleListSize = 100 > (1 << (g->nv/4)) ? 100 : (1 << (g->nv/4)); 
        if(g->nCyclePairs >= cycleListSize) {
            fprintf(stderr,
             "Error: out of memory when storing induced cycle pair.\n");
            exit(EXIT_FAILURE);
        }

        inducedCyclePairs[g->nCyclePairs] = pair;
        for(int k = 0; k < 2; k++) {
            g->endsCycle[k][g->nCyclePairs] = endsCycle[k];
            g->spokesCycle[k][g->nCyclePairs] = spokesCycle[k];
        }
        g->nCyclePairs++;
    }
}

// Assume that cycle0 is an induced n/2 cycle such that the degree 3 vertices of
// g induce a path on cycle0. This function checks whether cycle0 and its
// complement form a permutation 2-factor, i.e. we check whether the complement
// contains a cycle containing all remaining vertices. If so, we return the
// induced cycle containing 0. We also compute the degree 3 vertices of the
// second cycle and the degree 2 vertices adjacent to degree 3 vertices
// (ends) for each cycle and store them in spokesCycle[] and endsCycle[].
bitset inducedPairOrEmptySet(nautygraph *g, bitset cycle0, bitset spokesCycle[],
 bitset endsCycle[]) {

    bitset cycle1 = setComplement(cycle0, g->nv);
    bitset remaining = cycle1;
    int start = next(cycle1, -1);
    int currentVertex = start;
    int nextVertex = next(intersection(g->bitsetList[currentVertex], cycle1),
     -1);

    // Need to check if all vertices in complement of cycle0 form a cycle. Then
    // it will be induced (since complement has same number of degree 3
    // vertices and is induced).
    while(nextVertex != -1) {

        removeElement(remaining, currentVertex);
        if(contains(g->degree3, currentVertex)) {
            if(contains(g->degree2, nextVertex)) {
                add(endsCycle[1], nextVertex);
            }
        }
        else if(contains(g->degree3, nextVertex)) {
            add(endsCycle[1], currentVertex);
        }

        currentVertex = nextVertex;

        bitset nbrs = intersection(g->bitsetList[currentVertex], remaining);
        nextVertex = next(nbrs, -1);
    } 

    // Did not use all vertices in complement or did not form cycle
    if(size(remaining) != 1) return EMPTY;
    if(!contains(g->bitsetList[currentVertex], start)) return EMPTY;

    // Compute spokes and ends for final vertex.
    if(contains(g->degree2, currentVertex) && contains(g->degree3, start)) {
        add(endsCycle[1], currentVertex);
    }

    spokesCycle[1] = intersection(cycle1, g->degree3);

    // Store appropriate cycle (one containing 0).
    if(!contains(cycle0, 0)) {
        cycle0 = cycle1;
        bitset temp = spokesCycle[0];
        spokesCycle[0] = spokesCycle[1];
        spokesCycle[1] = temp;
        temp = endsCycle[0];
        endsCycle[0] = endsCycle[1];
        endsCycle[1] = temp;
    }

    return cycle0;
}

// Recursively build an induced cycle which will be part of a consecutive
// permutation 2-factor. We also keep track of the degree 3 vertices of this
// cycle and the degree 2 vertices adjacent to a degree 3 vertex.
bool partOfConsecutiveCyclePair(nautygraph *g, int start, int last, 
 bitset cycle, bitset remainingVertices, int nConsecutive, int nSpokes,
 bitset ends) {

    bitset nbrs = g->bitsetList[last];

    // If the number of consecutive degree 3 vertices is the same as the number
    // of spokes a permutation 2-factor should have, the first degree 2 vertex
    // added to the cycle will be an end.
    if(nConsecutive == nSpokes) {
        if(size(ends) < 2 && size(nbrs) == 2) {
            add(ends, last);
        }
    }

    // If we have an induced cycle such that the degree 3 vertices of g induce a
    // path in the cycle, check whether it is the cycle of a consecutive
    // permutation 2-factor.
    if(size(cycle) == g->nv/2 && contains(nbrs, start)) {
        bitset spokesCycle[2] = {intersection(cycle, g->degree3), EMPTY};
        bitset endsCycle[2] = {ends, EMPTY};
        bitset pair = inducedPairOrEmptySet(g, cycle, spokesCycle, endsCycle);
        insertInducedPair(g, pair, spokesCycle, endsCycle);
        return !isEmpty(pair);
    }

    // If not full cycle, but last adjacent to start, we will not make induced
    // cycle. (Unless in the first iteration when our cycle to be consists of
    // only u and v).
    if(contains(nbrs, start) && size(cycle) > 2) {
        return false;
    }

    bool found = false;

    // Add nbrs of last to the path. If already have enough consecutive degree 3
    // vertices, then only add degree 2 vertices to the path.
    bitset remainingNbrs = intersection(nbrs, remainingVertices);
    bitset newRemainingVertices = difference(remainingVertices, remainingNbrs);
    if(nConsecutive == nSpokes) {
        remainingNbrs = intersection(remainingNbrs, g->degree2);
        forEach(nbr, remainingNbrs) {
            if(partOfConsecutiveCyclePair(g, start, nbr,
             union(cycle, singleton(nbr)), newRemainingVertices, nConsecutive, 
             nSpokes, ends)) {
                found = true;
            }
        }
    }
    else {
        remainingNbrs = intersection(remainingNbrs, g->degree3);
        forEach(nbr, remainingNbrs) {
            if(partOfConsecutiveCyclePair(g, start, nbr,
             union(cycle, singleton(nbr)), newRemainingVertices, nConsecutive+1, 
             nSpokes, ends)) {
                found = true;
            }
        }
    }

    return found;
}

// Add all the consecutive permutation 2-factors containing uv to the list.
bool addNewCyclePairs(nautygraph *g, struct counters *counters, int u, int v) {

    counters->calledAddNewCycles++;

    // Such a new pair cannot exist if the newly added edge had 2 or more degree
    // 2 neighbours. (Proof in manuscript.)
    bitset deg2UNbrs = intersection(g->degree2, g->bitsetList[u]);
    int deg2NbrsU = size(deg2UNbrs);
    int deg2NbrsV = size(intersection(g->degree2, g->bitsetList[v]));

    if(deg2NbrsU+deg2NbrsV > 1) {
        counters->addNewCyclesAndMoreThanOneDeg2Nbr++;
        return false;
    }

    // Now v has two deg 3 nbrs and u at most 1, so we have to do the full
    // check.
    bitset cycle0 = union(singleton(u), singleton(v));
    bitset remainingVertices = setComplement(cycle0, g->nv);

    // We now build induced cycles using backtracking starting with uv, building
    // a path of degree 3 vertices and finishing with degree 2 vertices
    // (unless g is cubic already). One can show that all consecutive
    // permutation 2-factors must contain an induced cycle of this form.
    return partOfConsecutiveCyclePair(g, u, v, cycle0, remainingVertices, 2, 
     g->nde/2 - g->nv, deg2UNbrs);
}

// We assume that uv is not an edge of any 2-factor in the list and remove
// consecutive permutation 2-factors of g-uv in which u and v lie on the same
// induced cycle. 
void removeOldCyclePairs(nautygraph *g, int u, int v) {

    bitset *inducedCyclePairs = g->inducedCyclePairs;
    bitset uv = union(singleton(u), singleton(v));
    for(int i = 0; i < g->nCyclePairs; i++) {
        if(isEmpty(inducedCyclePairs[i])) {
            continue;
        }
        int s = size(intersection(uv, inducedCyclePairs[i]));

        // uv connects two elements of induced cycle or of its complement 
        if(s == 2 || s == 0) {
            inducedCyclePairs[i] = EMPTY;
        }
    } 

}

//******************************************************************************
//
//                      Canonical construction path method
//
//******************************************************************************


// The first time we call this for a given set of edges in an edge orbit, we
// mark e as the representative of that edge orbit.
bool isRepresentative(int e, int edgeOrbits[]) {

    int root = edgeOrbits[e];

    // Expansion already chosen. Root can be -1 if e was the root of its orbit.
    if(root == -1 || edgeOrbits[root] == -1) {
        return false; 
    }

    //  Point root to -1 if we choose the expansion.
    edgeOrbits[root] = -1;

    return true;

}

// When adding uv to g, the ends and spokes of the stored consecutive
// permutation 2-factors change. Update them.
void updateEnds(nautygraph *g, int u, int v) {
    int uv[2] = {u,v};
    for(int i = 0; i < g->nCyclePairs; i++) {

        bitset cycle1 = g->inducedCyclePairs[i];
        if(isEmpty(cycle1)) continue;

        bitset cycle2 = setComplement(cycle1, g->nv);
        bitset cycles[2] = {cycle1, cycle2};
        for(int j = 0; j < 2; j++) {

            int k = contains(cycle1, uv[j]) ? 0 : 1;

            // Update endsCycle and SpokesCycle
            removeElement(g->endsCycle[k][i], uv[j]);
            add(g->spokesCycle[k][i], uv[j]);
            bitset cycleNbrsOfDeg2 =
             difference(intersection(g->bitsetList[uv[j]], cycles[k]),
             g->spokesCycle[k][i]);
            g->endsCycle[k][i] = union(g->endsCycle[k][i], cycleNbrsOfDeg2);
        }
    }
}

// Macro for cleaning up main recursive method. Makes a copy of all information
// surrounding the stored consecutive permutation 2-factors.
#define COPYGRAPHSTATE \
    int nCyclePairsCopy = g->nCyclePairs;\
    bitset inducedCyclePairsCopy[g->nCyclePairs];\
    memcpy(inducedCyclePairsCopy, g->inducedCyclePairs,\
     g->nCyclePairs*sizeof(bitset));\
    bitset endsCycleCopy[2][g->nCyclePairs];\
    bitset spokesCycleCopy[2][g->nCyclePairs];\
    for(int k = 0; k < 2; k++) {\
        memcpy(endsCycleCopy[k], g->endsCycle[k],\
         g->nCyclePairs*sizeof(bitset));\
        memcpy(spokesCycleCopy[k], g->spokesCycle[k],\
         g->nCyclePairs*sizeof(bitset));\
    }\

// Resets the stored consecutive permutation 2-factors to the state in which it
// was copied.
#define RESETGRAPHSTATE \
    memcpy(g->inducedCyclePairs, inducedCyclePairsCopy,\
     nCyclePairsCopy*sizeof(bitset));\
    for(int k = 0; k < 2; k++) {\
        memcpy(g->endsCycle[k], endsCycleCopy[k],\
         nCyclePairsCopy*sizeof(bitset));\
        memcpy(g->spokesCycle[k], spokesCycleCopy[k],\
         nCyclePairsCopy*sizeof(bitset));\
    }\
    g->nCyclePairs = nCyclePairsCopy;\


// Main recursive method. We search eligible edges, compute their edge orbits,
// add the representative of each orbit, and check if this addition was
// canonical, if so we continue with this new graph, otherwise we prune. Spokes
// is the number of spokes the permutation 2-factors have.
void recursionCCPM(nautygraph *g, struct options *options,
 struct counters *counters, int spokes) {

    // Used for res/mod split.
    if(spokes == options->splitLevel) {

        //  Count how long the common part of each split takes. Algorithm
        //  stops as soon as the common part is finished.
#ifdef COUNT_SPLIT_TIME 
        return;
#endif

        options->splitCounter++;
        if(options->splitCounter == options->modulo) {
            options->splitCounter = 0;
        }
        if(options->splitCounter != options->remainder) {
            return;
        }
    }

    counters->recursionSteps++;
    counters->nonIsoIntermediateGraph++;

    // Graph is a permutation graph.
    if(spokes == g->nv/2) {

        writeToG6(g);
        counters->generatedGraphs++;

        return;
    }

    // 1. Find expansions/eligible edges: Given a consecutive permutation
    // 2-factor with cycles C1, C2 and the degree 3 vertices of G induce a path
    // P on C1, then eligible edges are (u,v) such that u has degree 2, lies on
    // C1, v has degree 2, lies on C2 and u is adjacent to (an end of) P.

    // 2. Compute classes of equivalent expansions. We do so by computing edge
    // orbits.

    // Find edges which are eligible for adding to the graph and label them. We
    // also check if the vertices satisfy girth requirements mentioned in
    // options->minimalGirth.
    struct edgeLabelling eLNonEdges = {0};
    getEligibleEdges(g, counters, options, &eLNonEdges);
    if(!eLNonEdges.totalEdges) {
        if(!g->calledNauty) {
            counters->didNotComputeGenerators++;
        }
        return;
    } 

    // If we computed canonical form in previous step we will already have
    // generators of the automorphism group. Otherwise call nauty to obtain
    // them.
    if(!g->calledNauty) {

        //  Set options and call nauty.
        g->calledNauty = true;
        g->options->getcanon = false;
        callNauty(g, NULL, counters); 
        counters->calledNautyForGenerators[g->nde/2-g->nv]++;
        counters->calledNautyForGeneratorsTotal++;
    }

    //  With the edge labels and the generators, we can compute edge orbits.
    int numberOfOrbits;
    int edgeOrbits[eLNonEdges.totalEdges];
    findEdgeOrbits(g->nv, &eLNonEdges, edgeOrbits, &numberOfOrbits); 

    // Copy the currently stored consecutive permutation 2-factors.
    COPYGRAPHSTATE;

    // 3. For each equivalence class do:

    // Loop over each eligible edge. If an edge belongs to an equivalence class
    // (orbit) already encountered, we continue.
    for(int e = 0; e < eLNonEdges.totalEdges; e++) {

        int u = eLNonEdges.indexToEdge[2*e];
        int v = eLNonEdges.indexToEdge[2*e + 1];

        // 4. Choose one expansion X

        //  Check if e is the representative of the equivalence class of
        //  extensions. (Not necessarily the root!)
        if(!isRepresentative(e, edgeOrbits)) {
            counters->notRepresentative++;
            continue;
        }

        // Girth restrictions are checked at the getEligibleEdges stage.


        // 5. Apply expansion X
        addSpoke(g, u, v);
        g->calledNauty = false;

        counters->addedSpoke++;
        counters->addedKthSpoke[g->nde/2-g->nv]++;

        // If adding uv destroys a consecutive permutation 2-factor, we remove
        // it and we update the ends of consecutive permutation 2-factors which
        // were not removed. 
        removeOldCyclePairs(g, u, v); 
        updateEnds(g, u, v);

        // We are now working with a subset of the reducible edges, but might
        // already deduce that uv is not canonical without calling expensive
        // check for new cycle pairs.
        if(!isCanonicalHeur(g, counters, options, g->eL, u, v)) {

            counters->canonHeurSufficientTotal++;
            counters->canonHeurSufficient[g->nde/2-g->nv]++;
            if(g->calledNauty) {
                counters->calledNautyEdgeNotCanonical[g->nde/2 - g->nv]++;
                counters-> calledNautyEdgeNotCanonicalTotal++;
            }

            removeSpoke(g, u, v);
            RESETGRAPHSTATE;
            
            continue; 
        }

        counters->canonHeurNotSufficientTotal++;
        counters->canonHeurNotSufficient[g->nde/2-g->nv]++;

        // Heuristics did not work, so we need to compute all new induced cycle
        // pairs to determine all reducible edges. Addednew will be true if we
        // added any. If we did not, we have already determined that uv is
        // canonical.
        bool addedNew = addNewCyclePairs(g, counters, u, v);

        if(!addedNew) {
            counters->canonHeurSufficientTotal++;
            counters->canonHeurSufficient[g->nde/2-g->nv]++;
            counters->subsetWasAllTotal++;
            counters->subsetWasAll[g->nde/2-g->nv]++;
            if(!g->calledNauty) {
                counters->canonDeterminedByHeuristicTotal++;
                counters->canonDeterminedByHeuristic[g->nde/2-g->nv]++;
            }
        }

        // 6. If expansion is canonical: recurse
        if(!addedNew || isCanonical(g, counters, options, g->eL, u, v)) {
            counters->wasCanonical++;
            counters->canonical[g->nde/2 - g->nv]++;
            if(g->calledNauty) {
                counters->calledNautyEdgeCanonical[g->nde/2 - g->nv]++;
                counters-> calledNautyEdgeCanonicalTotal++;
            }

            // Check if this edges adds a hamiltonian cycle if flag is present.
            // Seems to be bottleneck so perform this as little as possible.
            if(options->nonHamFlag) {
                counters->checkHam++;
                counters->checkHamSpokes[g->nde/2-g->nv]++;
                if(addingEdgeGivesHamiltonianCycleCubHam(g, counters, u, v)) {

                    counters->pruneHam++;
                    counters->pruneHamSpokes[g->nde/2-g->nv]++;

                    removeSpoke(g, u, v);
                    RESETGRAPHSTATE;

                    continue;
                }
            }

            // Recurse
            recursionCCPM(g, options, counters, spokes + 1);

        } 
        else {
            if(g->calledNauty) {
                counters->calledNautyEdgeNotCanonical[g->nde/2 - g->nv]++;
                counters-> calledNautyEdgeNotCanonicalTotal++;
            }
        }

        // 7. Perform reduction X^-1
        removeSpoke(g, u, v);
        RESETGRAPHSTATE;
    }
}

bool mirroringIsGreater(int list[], int list2[], int spokes, int n){
    for(int i = 1; i <= spokes; i++) {
        int normal = list[i];
        int mirrored = n/2 - list2[i];
        if(normal > mirrored) {
            // fprintf(stderr, "Mirror\n");
            return false; 
        }
        if(normal < mirrored) break;
    }
    return true;
}

bool reversingIsGreater(int list[], int list2[], int spokes, int n) {
    for(int i = 0; i <= spokes; i++) {
        int normal = list[i];
        int reversed = (n/2 - list2[spokes] + list2[spokes - i]) % (n/2);
        if(normal > reversed) {
            // fprintf(stderr, "Reverse\n");
            return false;
        }
        if(normal < reversed) break;
    }
    return true;
}

// Compares the cyclic lists list1 and list2 starting at idx1 and idx2 resp. We
// compare the lexicographically cumulative lists.
int compareLists(int list1[], int idx1, int list2[], int idx2, int len, int n) {
    int sum1 = 0;
    int sum2 = 0;
    for(int i = 0; i < len; i++) {
        sum1 = (sum1 + list1[(idx1+i)%len])%(n/2);
        sum2 = (sum2 + list2[(idx2+i)%len])%(n/2);
        if(sum1 < sum2) return -1;
        if(sum1 > sum2) return 1;
    }
    return 0;
}

// Let L be a list with L[0] = 0 and distinct elements between 0 and n/2. Let
// distList contain at index i value (L[i+1] - L[i])%n/2. This method returns
// the index k such that if we cyclically rotate L such that k is at position 0
// and subtract L[k] modulo n/2 we get the lexicographically smallest list.
int leastRotation(int distList[], int length, int n) {
    int smallest = 0;
    for(int i = 1; i < length; i++) {
        if(compareLists(distList, i, distList, smallest, length, n) == -1) {
            smallest = i;
        }
    }
    return smallest;
}

bool distListIsLexMinimalRot(int list[], int n) {
    // int correctList[] = {0,2,4,6,8,5,7,1,3};
    // int isCorrect = (compareLists(list, 0, correctList, 0, 6, n) == 0);
    // if(isCorrect) {
    //     fprintf(stderr, "Correct\n");
    // }

    int distList[n/2];
    int reverseList[n/2];
    for(int i = 0; i < n/2-1; i++) {
        distList[i] = (n/2 + list[i+1] - list[i]) % (n/2);
        reverseList[n/2-1-i] = distList[i];
    }
    distList[n/2-1] = (n/2 + list[0] - list[n/2-1]) % (n/2);
    reverseList[0] = distList[n/2-1];

    // if(isCorrect) {
    //     fprintf(stderr, "distList: ");
    //     for(int i = 0; i < n/2; i++) {
    //         fprintf(stderr, "%d ", distList[i]);
    //     }
    //     fprintf(stderr, "\n");
    //             fprintf(stderr, "revList: ");
    //     for(int i = 0; i < n/2; i++) {
    //         fprintf(stderr, "%d ", reverseList[i]);
    //     }
    //     fprintf(stderr, "\n");
    //     fprintf(stderr, "%d\n",leastRotation(distList, n/2, n));
    // }

    if(compareLists(distList, 0, distList, leastRotation(distList, n/2, n),
     n/2, n) == 1) {
        // if(isCorrect) {
        //     // fprintf(stderr, "Normal rotation\n");
        // }
        return false;
    }
    if(compareLists(distList, 0, reverseList,
     leastRotation(reverseList, n/2, n), n/2, n) == 1) {
        // if(isCorrect) fprintf(stderr, "Reverse\n");
        return false;
    }

    int revReverseList[n/2];
    for(int i = 0; i < n/2-1; i++) {
        reverseList[n/2 - 1 - i] = (n/2 - list[i+1] + list[i]) % (n/2);
        revReverseList[i] = reverseList[n/2 - 1 -i]; 
    }
    reverseList[0] = (n/2 - list[0] + list[n/2-1]) % (n/2);
    revReverseList[n/2-1] = reverseList[0];
        // fprintf(stderr, "revList: ");
        // for(int i = 0; i < n/2; i++) {
        //     fprintf(stderr, "%d ", reverseList[i]);
        // }
        // fprintf(stderr, "\n");

    // if(!compareSuffixes(distList, reverseList, n/2)) return false;
    if(compareLists(distList, 0, reverseList,
     leastRotation(reverseList, n/2, n), n/2, n) == 1)
    {
        // if(isCorrect) fprintf(stderr, "Reverse+Compl\n");
        return false;
    }
    if(compareLists(distList, 0, revReverseList,
     leastRotation(revReverseList, n/2, n), n/2, n) == 1) {
        // if(isCorrect) fprintf(stderr, "Reverse+Compl\n");
        return false;
    }

    int swappedList[n/2];
    for(int i = 0; i < n/2; i++) {
        swappedList[list[i]] = i;
    }
    int distListSwapped[n/2];
    for(int i = 0; i < n/2-1; i++) {
        distListSwapped[i] = (n/2 + swappedList[i+1] - swappedList[i]) % (n/2);
        reverseList[n/2-1-i] = distListSwapped[i];
    }
    distListSwapped[n/2-1] = (n/2 + swappedList[0] - swappedList[n/2-1]) % (n/2);
    reverseList[0] = distListSwapped[n/2-1];
    if(compareLists(distList, 0, distListSwapped,
     leastRotation(distListSwapped, n/2, n), n/2, n) == 1)
    {
        // if(isCorrect) {
        //     for(int i = 0; i < n/2; i++) {
        //         fprintf(stderr, "%d ", swappedList[i]);
        //     }
        //     fprintf(stderr, "\n %d, ", leastRotation(distListSwapped, n/2));
        //     fprintf(stderr, "Swapped\n");
        // }
        return false;
    }
    if(compareLists(distList, 0, reverseList,
     leastRotation(reverseList, n/2, n), n/2, n) == 1) {
        // if(isCorrect) fprintf(stderr, "Swapped+rev\n");
        return false;
    }


    for(int i = 0; i < n/2-1; i++) {
        reverseList[n/2 - 1 - i] = (n/2 - swappedList[i+1] + swappedList[i]) % (n/2);
        revReverseList[i] = reverseList[n/2 - 1 -i]; 
    }
    reverseList[0] = (n/2 - swappedList[0] + swappedList[n/2-1]) % (n/2);
    revReverseList[n/2-1] = reverseList[0];

    if(compareLists(distList, 0, reverseList, leastRotation(reverseList, n/2, n),
     n/2, n) == 1) {
        // if(isCorrect) fprintf(stderr, "Swapped + rev + compl\n");
        return false;
    }
    if(compareLists(distList, 0, revReverseList,
     leastRotation(revReverseList, n/2, n), n/2, n) == 1) {
        // if(isCorrect) fprintf(stderr, "Swapped+Reverse+Compl\n");
        return false;
    }

    return true;   
}

void addForbiddenSteps(int list[], int inverseList[], int spokes, bitset remainingOptions, bitset forbiddenNextSteps[], int n) {
    int fj1 = list[spokes - 1];
    int fj2 = list[spokes];

    int k = n/2;
    int inc[] = {(fj1+1)%k, (fj2+1)%k};
    int dec[] = {(k+fj1-1)%k, (k+fj2-1)%k};


    // We have i1 != j1-1, because otherwise fi2 = fi1 - 1

    // fi1 = fj1 + 1 or fi2 = fj1 + 1
    if(!contains(remainingOptions, inc[0])) {
        int fi1 = inc[0];
        int i1 = inverseList[fi1];
        int fi2 = list[(i1+1)%k];

        // na fi1, fj2 dan fi2
        if((fj2 > fi1 && fi2 > fj2) || (fj2 > fi1 && fi2 < fi1-1) || (fi2 < fj1 && fj2 < fi2)) {
            add(forbiddenNextSteps[dec[1]], (fi2+1)%k);
            add(forbiddenNextSteps[inc[1]], (fi2+1)%k);
            add(forbiddenNextSteps[(fi2 + 1)%k], inc[1]);
        }
        // na fi1, fi2 dan fj2
        else {
            add(forbiddenNextSteps[(k + fi2-1)%k], dec[1]);
            add(forbiddenNextSteps[(k + fi2-1)%k], inc[1]);
            add(forbiddenNextSteps[(fi2+1)%k], dec[1]);
            add(forbiddenNextSteps[dec[1]], (fi2 + 1)%k);
            add(forbiddenNextSteps[inc[1]], (fi2 + 1)%k);
        }

        fi2 = inc[0];
        int i2 = inverseList[fi2];
        if(i2 != 0) {
            fi1 = list[(k + i2 -1)%k];

            // na fi2, fi1 dan fj2
            if((fi1 > fi2 && fj2 > fi1) || (fi1 > fi2 && fj2 < fi2-1) || (fj2 < fj1 && fi1 < fj2)) {
                add(forbiddenNextSteps[(k + fi1-1)%k], inc[1]);
                add(forbiddenNextSteps[(fi1+1)%k], dec[1]);
                add(forbiddenNextSteps[dec[1]], (k + fi1-1)%k);
                add(forbiddenNextSteps[dec[1]], (fi1+1)%k);
                add(forbiddenNextSteps[inc[1]], (k + fi1-1)%k);
                add(forbiddenNextSteps[inc[1]], (fi1+1)%k);
            }
            else {
                add(forbiddenNextSteps[dec[1]], (k+fi1-1)%k);
                add(forbiddenNextSteps[dec[1]], (fi1+1)%k);
                add(forbiddenNextSteps[inc[1]], (k+fi1-1)%k);
                add(forbiddenNextSteps[inc[1]], (fi1+1)%k);
                add(forbiddenNextSteps[(k + fi1-1)%k], inc[1]);
                add(forbiddenNextSteps[(fi1+1)%k], dec[1]);
            }
        }
    }

    // fi1 = fj1 - 1 or fi2 = fj1 - 1
    if(!contains(remainingOptions, dec[0])) {
        int fi1 = dec[0];
        int i1 = inverseList[fi1];
        int fi2 = list[(i1+1)%k];

        // We have i1 != j1-1, because otherwise fi2 = fi1 + 1

        // na fj1, fj2 dan fi2
        if((fj2 > fj1 && fi2 > fj2) || (fj2 > fj1 && fi2 < fj1-1) || (fi2 < fi1 && fj2 < fi2)) {
            add(forbiddenNextSteps[dec[1]], (k + fi2 - 1)%k);
            add(forbiddenNextSteps[inc[1]], (k+fi2-1)%k);
            add(forbiddenNextSteps[(k + fi2-1)%k], inc[1]);
            add(forbiddenNextSteps[(fi2 + 1)%k], dec[1]);
            add(forbiddenNextSteps[(fi2 + 1)%k], inc[1]);
        }

        // na fj1, fi2 dan fj2
        else {
            add(forbiddenNextSteps[(k + fi2-1)%k], dec[1]);
            add(forbiddenNextSteps[dec[1]], (k + fi2 - 1)%k);
            add(forbiddenNextSteps[inc[1]], (k + fi2 - 1)%k);
        }

        fi2 = dec[0];
        int i2 = inverseList[fi2];
        if(i2 != 0) {
            fi1 = list[(k + i2 -1)%k];

            // na fj1, fi1 dan fj2
            if((fi1 > fj1 && fj2 > fi1) || (fi1 > fj1 && fj2 < fj1-1) || (fj2 < fi2 && fi1 < fj2)) {
                add(forbiddenNextSteps[(k + fi1-1)%k], inc[1]);
                add(forbiddenNextSteps[(fi1+1)%k], dec[1]);
                add(forbiddenNextSteps[dec[1]], (k + fi1-1)%k);
                add(forbiddenNextSteps[dec[1]], (fi1+1)%k);
                add(forbiddenNextSteps[inc[1]], (k + fi1-1)%k);
                add(forbiddenNextSteps[inc[1]], (fi1+1)%k);
            }
            else {
                add(forbiddenNextSteps[dec[1]], (k+fi1-1)%k);
                add(forbiddenNextSteps[dec[1]], (fi1+1)%k);
                add(forbiddenNextSteps[inc[1]], (k+fi1-1)%k);
                add(forbiddenNextSteps[inc[1]], (fi1+1)%k);
                add(forbiddenNextSteps[(k + fi1-1)%k], inc[1]);
                add(forbiddenNextSteps[(fi1+1)%k], dec[1]);
            }
        }
    }

    // fi1 = fj2+1 or fi2 = fj2+1
    if(!contains(remainingOptions, inc[1])) {
        int fi1 = inc[1];
        int i1 = inverseList[fi1];
        int fi2 = list[(i1+1)%k];

        if(fi2 != fj1) {
            // na fi1, fi2 dan fj1
            if((fi2 > fi1 && fj1 > fi2) || (fi2 > fi1 && fj1 < fi1-1) || (fj1 < fj2 && fi2 < fj1)) {
                add(forbiddenNextSteps[(k + fi2 - 1)%k], dec[0]);
                add(forbiddenNextSteps[(fi2 + 1)%k], dec[0]);
                add(forbiddenNextSteps[(fi2 + 1)%k], inc[0]);
                add(forbiddenNextSteps[dec[0]], (fi2+1)%k);
            }
            else {
                add(forbiddenNextSteps[dec[0]], (k+fi2-1)%k);
                add(forbiddenNextSteps[inc[0]], (k+fi2-1)%k);
                add(forbiddenNextSteps[inc[0]], (fi2+1)%k);
                add(forbiddenNextSteps[(k + fi2 - 1)%k], inc[0]);

            }
        }

        fi2 = inc[1];
        int i2 = inverseList[fi2];
        fi1 = list[(k+i2-1)%k];

        if(i2 != 0) {
            // na fi2, fi1 dan fj1
            if((fi1 > fi2 && fj1 > fi1) || (fi1 > fi2 && fj1 < fi2-1) || (fj1 < fj2 && fi1 < fj1)) {
                add(forbiddenNextSteps[(k + fi1 - 1)%k], inc[0]);
                add(forbiddenNextSteps[(fi1 + 1)%k], dec[0]);
                add(forbiddenNextSteps[(fi1 + 1)%k], inc[0]);
                add(forbiddenNextSteps[dec[0]], (k + fi1-1)%k);
                add(forbiddenNextSteps[dec[0]], (fi1+1)%k);
            }
            else {
                add(forbiddenNextSteps[dec[0]], (k+fi1-1)%k);
                add(forbiddenNextSteps[dec[0]], (fi1+1)%k);
                add(forbiddenNextSteps[(k + fi1 - 1)%k], dec[0]);

            }
        }
    }

    // fi1 = fj2 - 1 or fi2 = fj2 - 1
    if(!contains(remainingOptions, dec[1])) {
        int fi1 = dec[1];
        int i1 = inverseList[fi1];
        int fi2 = list[(i1+1)%k];

        if(fi2 != fj1) {
            // na fj2, fi2 dan fj1
            if((fi2 > fj2 && fj1 > fi2) || (fi2 > fj2 && fj1 < fj2-1) || (fj1 < fi1 && fi2 < fj1)) {
                add(forbiddenNextSteps[(fi2 + 1)%k], dec[0]);
                add(forbiddenNextSteps[dec[0]], (k + fi2-1)%k);
                add(forbiddenNextSteps[dec[0]], (fi2+1)%k);
                add(forbiddenNextSteps[inc[0]], (fi2+1)%k);
            }
            else {
                add(forbiddenNextSteps[inc[0]], (k + fi2-1)%k);
                add(forbiddenNextSteps[(k + fi2 - 1)%k], dec[0]);
                add(forbiddenNextSteps[(k + fi2 - 1)%k], inc[0]);
                add(forbiddenNextSteps[(fi2 + 1)%k], inc[0]);

            }
        }

        fi2 = dec[1];
        int i2 = inverseList[fi2];
        fi1 = list[(k+i2-1)%k];

        if(i2 != 0) {
            // na fj2, fi1 dan fj1
            if((fi1 > fj2 && fj1 > fi1) || (fi1 > fj2 && fj1 < fj2-1) || (fj1 < fi2 && fi1 < fj1)) {
                add(forbiddenNextSteps[(fi1 + 1)%k], inc[0]);
                add(forbiddenNextSteps[inc[0]], (k+fi1-1)%k);
                add(forbiddenNextSteps[inc[0]], (fi1+1)%k);
            }
            else {
                add(forbiddenNextSteps[inc[0]], (k+fi1-1)%k);
                add(forbiddenNextSteps[inc[0]], (fi1+1)%k);
                add(forbiddenNextSteps[(k+fi1 - 1)%k], dec[0]);
                add(forbiddenNextSteps[(k+fi1 - 1)%k], inc[0]);
                add(forbiddenNextSteps[(fi1 + 1)%k], dec[0]);
            }
        }
    }
}


// Another method which builds permutations and only accepts the
// lexicographically smallest after performing several operations. This method
// does not eliminate all isomorphisms.
void recursionLexSmallestPerm(nautygraph *g, int list[], int inverseList[],
 bitset remainingOptions, bitset forbiddenNextSteps[], struct options *options,
  struct counters *counters, int spokes) {

    // Used for res/mod split.
    if(spokes == options->splitLevel) {

        //  Count how long the common part of each split takes. Algorithm
        //  stops as soon as the common part is finished.
#ifdef COUNT_SPLIT_TIME 
        return;
#endif

        options->splitCounter++;
        if(options->splitCounter == options->modulo) {
            options->splitCounter = 0;
        }
        if(options->splitCounter != options->remainder) {
            return;
        }
    }

    counters->recursionSteps++;
    counters->nonIsoIntermediateGraph++;

    // Graph is a permutation graph.
    if(spokes == g->nv/2) {

        writeToG6(g);
        counters->generatedGraphs++;

        return;
    }

    // if(!isEmpty(intersection(forbiddenNextSteps[list[spokes - 1]], remainingOptions))) {
    // //     printGraph(g);
    // //     printBitset(badOpts);
    //     counter++;
    // }

    bitset goodOpts = difference(remainingOptions, forbiddenNextSteps[list[spokes - 1]]);

    bitset forbiddenNextStepsCopy[BITSETSIZE];
    memcpy(forbiddenNextStepsCopy, forbiddenNextSteps, sizeof(bitset) * BITSETSIZE);

    forEach(el, goodOpts) {

        if(options->minimalGirth > 4) {
            if(containsForbiddenCycleWithEdge(g, options->minimalGirth, spokes, el+g->nv/2)) continue;
        }

        list[spokes] = el;
        inverseList[el] = spokes;
        bool isSmallest = true;

        // fprintf(stderr, "List: ");
        // for(int i = 0; i < spokes+1; i++) {
        //     fprintf(stderr, "%d ", list[i]);
        // }
        // fprintf(stderr, "\n");

        if(spokes == g->nv/2 -  1 && !distListIsLexMinimalRot(list, g->nv)) {
            // fprintf(stderr, "Not lex minimal rotation\n");
            continue;
        }

        // Check if mirroring inner cycle is lexicographically smaller:
        if(!mirroringIsGreater(list, list, spokes, g->nv)) continue;

        // Check if reversing permutation is smaller:
        if(!reversingIsGreater(list, list, spokes, g->nv)) continue;

        // Check if both cycles are consecutive
        bitset opts = difference(remainingOptions, singleton(el));
        int firstOption = next(opts, -1);
        int secondNonOption = next(setComplement(opts, g->nv/2), 0);
        if(firstOption <= spokes && secondNonOption < g->nv/2-spokes) {

            addSpoke(g, spokes, el+g->nv/2);

            if(options->nonHamFlag) {
                add(forbiddenNextSteps[el - 1], (g->nv/2 + list[spokes - 1] - 1)%(g->nv/2));
                add(forbiddenNextSteps[el - 1], (list[spokes - 1] + 1) % (g->nv/2));
                add(forbiddenNextSteps[(el + 1)%(g->nv/2)], (g->nv/2 + list[spokes - 1] - 1)%(g->nv/2));
                add(forbiddenNextSteps[(el + 1)%(g->nv/2)], (list[spokes - 1] + 1) % (g->nv/2));
                add(forbiddenNextSteps[(g->nv/2 + list[spokes - 1] - 1)%(g->nv/2)], (el+1)%(g->nv/2));
                add(forbiddenNextSteps[(list[spokes - 1] + 1)%(g->nv/2)], (g->nv/2+el-1)%(g->nv/2));
                addForbiddenSteps(list, inverseList, spokes, remainingOptions, forbiddenNextSteps, g->nv);

                counters->checkHam++;
                counters->checkHamSpokes[spokes+1]++;
                if(
                 addingEdgeGivesHamiltonianCycleCubHam(g, counters, spokes, el+g->nv/2)) {

                    counters->pruneHamSpokes[spokes+1]++;
                    counters->pruneHam++;
                    removeSpoke(g, spokes, el+g->nv/2);
                    memcpy(forbiddenNextSteps, forbiddenNextStepsCopy, sizeof(bitset) * BITSETSIZE);


                    continue;
                }
            }

            recursionLexSmallestPerm(g, list, inverseList,
             difference(remainingOptions, singleton(el)), forbiddenNextSteps,
              options, counters, spokes+1);

            removeSpoke(g, spokes, el+g->nv/2);
            memcpy(forbiddenNextSteps, forbiddenNextStepsCopy, sizeof(bitset) * BITSETSIZE);

            continue;
        }

        // If both are consecutive we can swap.

        // Swapping cycles is smaller:
        int swappedList[spokes+1];
        if(firstOption > spokes) {
            for(int i = 0; i <= spokes; i++) {
                swappedList[list[i]] = i;
            }
        }
        else {
            for(int i = 0; i <= spokes; i++) {
                swappedList[(g->nv/2 - list[i]) % (g->nv/2)] = i;
            }
        }
        for(int i = 0; i <= spokes; i++) {
            if(list[i] > swappedList[i]) {
                isSmallest = false;
                break;
            }
            if(list[i] < swappedList[i]) break;
        }
        if(!isSmallest) continue;


        // Check if mirroring inner cycle after swap is lexicographically
        // smaller:
        if(!mirroringIsGreater(list, swappedList, spokes, g->nv)) continue;

        // Check if reversing permutation after swap is smaller:
        if(!reversingIsGreater(list, swappedList, spokes, g->nv)) continue;

        // Is this the best place for this check?
        if(spokes == g->nv/2 - 1 && contains(forbiddenNextSteps[el], list[0])) {
            continue;
        }

        addSpoke(g, spokes, el+g->nv/2);

       
        if(options->nonHamFlag) {

            // Seems to have no effect, but do not see why it should not for very
            // large n.
            add(forbiddenNextSteps[el - 1], (g->nv/2 + list[spokes - 1] - 1)%(g->nv/2));
            add(forbiddenNextSteps[el - 1], (list[spokes - 1] + 1) % (g->nv/2));
            add(forbiddenNextSteps[(el + 1)%(g->nv/2)], (g->nv/2 + list[spokes - 1] - 1)%(g->nv/2));
            add(forbiddenNextSteps[(el + 1)%(g->nv/2)], (list[spokes - 1] + 1) % (g->nv/2));
            add(forbiddenNextSteps[(g->nv/2 + list[spokes - 1] - 1)%(g->nv/2)], (el+1)%(g->nv/2));
            add(forbiddenNextSteps[(list[spokes - 1] + 1)%(g->nv/2)], (g->nv/2+el-1)%(g->nv/2));
            addForbiddenSteps(list, inverseList, spokes, remainingOptions, forbiddenNextSteps, g->nv);
            counters->checkHam++;
            counters->checkHamSpokes[spokes+1]++;
            if(
             addingEdgeGivesHamiltonianCycleCubHam(g, counters, spokes, el+g->nv/2)) {

                counters->pruneHamSpokes[spokes+1]++;
                counters->pruneHam++;
                removeSpoke(g, spokes, el+g->nv/2);
                memcpy(forbiddenNextSteps, forbiddenNextStepsCopy, sizeof(bitset) * BITSETSIZE);

                continue;
            }
        }

        // Seems to be lexicographically smallest, so recurse.
        recursionLexSmallestPerm(g, list, inverseList,
         difference(remainingOptions, singleton(el)), forbiddenNextSteps,
          options, counters, spokes+1);

        removeSpoke(g, spokes, el+g->nv/2);
        memcpy(forbiddenNextSteps, forbiddenNextStepsCopy, sizeof(bitset) * BITSETSIZE);
    }
}

//******************************************************************************
//
//                  Initialise algorithm
//
//******************************************************************************

// Add two disjoint cycles of length n/2 to our empty graph and add all vertices
// to degree 2 set.
void addTwoCycles(nautygraph *g) {
    int n = g->nv;
    int lenCycle = n / 2;

    for(int i = 0; i < lenCycle; i++) {

        addCycleEdges(g, i, lenCycle, 0);
        addCycleEdges(g, i, lenCycle, lenCycle);

    }
    g->degree2 = setComplement(EMPTY, g->nv);
    g->degree3 = EMPTY;
}

// Initialize graph struct.
void initializeGraphWithTwoCycles(nautygraph *g, int n) {

    int m = SETWORDSNEEDED(n);
    nauty_check(WORDSIZE,m,n,NAUTYVERSIONID);

    graph* adjacencyList = malloc(sizeof(graph)*n*m);
    if(adjacencyList == NULL) {
        fprintf(stderr, "Error: out of memory\n");
        exit(1);
    }
    EMPTYGRAPH(adjacencyList, m, n);

    g->nv = n;
    g->nde = 2*n; // Start with 2n directed edges from the two induced cycles.
    g->nCyclePairs = 1;
    g->numberOfWords = m;
    g->adjacencyList = adjacencyList;
    g->calledNauty = false;
    g->bitsetList = calloc(n, sizeof(bitset));

    // An estimate of how many induced cycle pairs we might encounter at most.
    int cycleListSize = 100 > (1 << (n/4)) ? 100 : (1 << (n/4)); 

    // This array will contain a bitset for every consecutive permutation
    // 2-factor. The bitset will contain the vertices of its induced cycle
    // containing 0. Might contain the EMPTY set, but these entries should be
    // ignored and replaced when new such 2-factors are found. Will refer to
    // the cycle containing 0 as cycle 0 and the complement as cycle 1.
    g->inducedCyclePairs = malloc(cycleListSize*sizeof(bitset));

    // For cycle 0 and 1 keep track of the ends of the spoke vertices, i.e. keep
    // track of the degree 2 vertices on each cycle which are adjacent to a
    // degree 3 vertex.
    g->endsCycle[0] = malloc(cycleListSize*sizeof(bitset));
    g->endsCycle[1] = malloc(cycleListSize*sizeof(bitset));

    // For cycle 0 and 1 keep track of the degree 3 vertices, i.e. the vertices
    // which are the endpoint of a spoke between the cycles.
    g->spokesCycle[0] = malloc(cycleListSize*sizeof(bitset));
    g->spokesCycle[1] = malloc(cycleListSize*sizeof(bitset));

    // Since we start with two cycles and edge (0,n/2) update the arrays.
    g->spokesCycle[0][0] = singleton(0);
    g->spokesCycle[1][0] = singleton(g->nv/2);
    g->inducedCyclePairs[0] = setComplement(EMPTY, g->nv/2);
    g->endsCycle[0][0] = union(singleton(g->nv/2-1), singleton(1));
    g->endsCycle[1][0] = union(singleton(g->nv/2+1), singleton(g->nv-1));

    // Fix the two cycles and connect in one place to avoid symmetries.
    addTwoCycles(g);
}

void initializeNauty(nautygraph *g, statsblk *stats, optionblk *options, 
 int *lab, int *ptn, int *orbits) {
    g->lab = lab;
    g->ptn = ptn;
    g->stats = stats;
    g->orbits = orbits;
    g->options = options;
    g->options->userautomproc = saveGenerators;
}

void generatePermutationGraphs(struct options *options,
 struct counters *counters, int n) {

    nautygraph g;
    initializeGraphWithTwoCycles(&g, n);

    // Label all edges in the graph. 
    struct edgeLabelling eL = {0};
    for(int i = 0; i < n; i++) { 
        forEachAfterIndex(j, g.bitsetList[i], i) {
            labelEdge(i, j, &eL);
        }
    }
    relabelFirstEdge(&eL);
    g.eL = &eL;

    // Initialize nauty variables
    int lab[g.nv];
    int ptn[g.nv];
    statsblk stats;
    int orbits[g.nv];
    DEFAULTOPTIONS_GRAPH(nautyOptions);

    if(!options->permutationMethodFlag) {
        nauty_check(WORDSIZE,g.numberOfWords,n,NAUTYVERSIONID);
        initializeNauty(&g, &stats, &nautyOptions, lab, ptn, orbits);
    }

    //  Add first spoke, wlog it can always be this one.
    addSpoke(&g, 0, n/2); 

    if(!options->permutationMethodFlag) {
        recursionCCPM(&g, options, counters, 1);
    }
    else {
        int list[n/2];
        int inverseList[n/2];
        list[0] = 0;
        inverseList[0] = 0;

        bitset forbiddenNextSteps[BITSETSIZE] = {0};
        bitset remainingOptions = setComplement(singleton(0), n/2);
        recursionLexSmallestPerm(&g, list, inverseList, remainingOptions, forbiddenNextSteps,
         options, counters, 1);
    }

    freeNautyGraph(&g);
    free(generators);
}

void printStatistics(struct options * options, struct counters * counters,
 int n) {

    fprintf(stderr, "Recursion steps: %llu\n",counters->recursionSteps);
    fprintf(stderr, "Times isomorphic intermediate graph generated: %llu\n",
     counters->timesPrunedBecauseIsomorphism);
    fprintf(stderr, "Number of non-isomorphic intermediate graphs: %llu\n",
     counters->nonIsoIntermediateGraph);
    fprintf(stderr, "Calls to compare method: %llu\n",
     counters->callsToCompareMethod);

    if(counters->nonIsoIntermediateGraph == 0)  {
        fprintf(stderr, "Error: divide by zero.");
        exit(EXIT_FAILURE);
    }

    fprintf(stderr, "Average number of vertex orbits: %llu\n",
     counters->numberOfVertexOrbits/counters->nonIsoIntermediateGraph);
    fprintf(stderr, "Average group sizes: %llu\n", 
     counters->sumOfGrpSizes/counters->nonIsoIntermediateGraph);
    fprintf(stderr, "Times group size > 1: %llu\n", 
     counters->grpSizeGrTh1);

    if(options->minimalGirth >= 3) {
        fprintf(stderr,
         "Times pruned girth: %llu (%.2f%% counts edges in same orbit)\n",
         counters->pruneGirth, 
         100.0*counters->pruneGirth / 
          (counters->pruneGirth + 
           counters->pruneHam + 
           counters->notRepresentative + 
           counters->addedSpoke));
    }

    if(options->nonHamFlag) {
        fprintf(stderr, 
         "Times check hamiltonian: %llu, " 
         "Times expensive check: %llu, "
         "Times pruned hamiltonian: %llu (%.2f%% of checks, %.2f%% of exp. checks)\n",
         counters->checkHam, counters->expensiveHam, counters->pruneHam, 
         100.0 * counters->pruneHam / counters->checkHam,
         100.0 * counters->pruneHam / counters->expensiveHam);

        for(int i = 1; i < n/2+1; i++) {
            fprintf(stderr,
             "Spokes: %2d, \t check ham: %7llu "
             "\t expensive check: %7llu, "
             "\t pruned because ham: %7llu (%.2f%% of checks, %.2f%% of exp. checks)\n",
             i, counters->checkHamSpokes[i], 
             counters->expensiveHamSpokes[i], counters->pruneHamSpokes[i],
             100.0 * counters->pruneHamSpokes[i] / counters->checkHamSpokes[i],
             100.0 * counters->pruneHamSpokes[i] / counters->expensiveHamSpokes[i]);
        }
    }

    fprintf(stderr, "Total nauty calls: %llu\n", counters->totalNautyCalls);
    fprintf(stderr, 
     "Times added spoke: %llu, Times spoke was canonical: %llu (%.2f%%)\n",
     counters->addedSpoke, counters->wasCanonical, 
     100.0 * counters->wasCanonical / counters->addedSpoke);

    fprintf(stderr, 
     "Edge determined canonical by heuristics (ie 1 left): "
     "%llu (%.2f%% of canonical edges)\n", 
     counters->canonDeterminedByHeuristicTotal, 
     100.0*counters->canonDeterminedByHeuristicTotal / counters->wasCanonical);
    for(int i = 1; i < n/2+1; i++) {
        fprintf(stderr, 
         "Spokes: %2d, \t Added spoke: %7llu, \tCalled nauty: %7llu, "
         "\tEdge was canonical: %7llu (%.2f%% of nauty) (%.2f%% of addedspoke)"
         "\t Edge det by heur: %7llu (%.2f%% of canon edges)\n",
         i, counters->addedKthSpoke[i], counters->calledNauty[i],
         counters->canonical[i], 
         100.0 * counters->canonical[i] / counters->calledNauty[i],
         100.0 * counters->canonical[i] / counters->addedKthSpoke[i], 
         counters->canonDeterminedByHeuristic[i], 
         100.0*counters->canonDeterminedByHeuristic[i]/counters->canonical[i]);
    }

    fprintf(stderr, "Return because of invariant:\n\t");
    int totalReturnAtInv = 0;
    for(int i = 0; i < NUMBEROFINVARIANTS; i++) {
        totalReturnAtInv += counters->returnAtInvariant[i];
    }
    for(int i = 0; i < NUMBEROFINVARIANTS; i++) {
        fprintf(stderr, "%d: %llu (%.2f%%) ", i, counters->returnAtInvariant[i],
         100.0 * counters->returnAtInvariant[i] / totalReturnAtInv);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, 
     "\nReturn because of invariant for each number of spokes:\n");
    for(int i = 1; i < n/2+1; i++) {
        fprintf(stderr, "Spokes: %2d     ", i);
        totalReturnAtInv = 0;
        for(int j = 0; j < NUMBEROFINVARIANTS; j++) {
            totalReturnAtInv += counters->returnAtInvariantKthSpoke[i][j];
        }
        for(int j = 0; j < NUMBEROFINVARIANTS; j++) {
            fprintf(stderr, "Inv %d: %7llu (%5.2f%%)    ", j, 
             counters->returnAtInvariantKthSpoke[i][j], 
             100.0*counters->returnAtInvariantKthSpoke[i][j]/totalReturnAtInv);
        }
        fprintf(stderr, "\n");
    }

    fprintf(stderr, 
     "Times need to check all reducible edges: %llu (%.2f%%) "
     "Times checking subset sufficient: %llu (%.2f%%)"
     "\t Subset was all reducible: %7llu "
     "(%.2f%% of times subset not sufficient)\n",
     counters->needToCheckAllEdgesTotal, 
     100.0 * counters->needToCheckAllEdgesTotal / 
     (counters->needToCheckAllEdgesTotal + counters->canonHeurSufficientTotal),
     counters->canonHeurSufficientTotal, 
     100.0 * counters->canonHeurSufficientTotal / 
     (counters->needToCheckAllEdgesTotal + counters->canonHeurSufficientTotal),
     counters->subsetWasAllTotal, 
     100.0 * counters->subsetWasAllTotal / 
     (counters->canonHeurNotSufficientTotal));

    for(int i = 1; i < n/2+1; i++) {
        fprintf(stderr, "Spokes: %2d, \t Times check all: %7llu (%.2f%%) \t "
         "Times subset sufficient: %7llu (%.2f%%)"
         "\t Subset was all reducible: %7llu "
         "(%.2f%% of times subset not sufficient)\n",
         i, counters->needToCheckAllEdges[i], 
         100.0 * counters->needToCheckAllEdges[i] / 
         (counters->needToCheckAllEdges[i] + counters->canonHeurSufficient[i]),
         counters->canonHeurSufficient[i],
         100.0 * counters->canonHeurSufficient[i] / 
         (counters->needToCheckAllEdges[i] + counters->canonHeurSufficient[i]),
         counters->subsetWasAll[i], 
         100.0 * counters->subsetWasAll[i] / 
         counters->canonHeurNotSufficient[i]);
    }

    fprintf(stderr, 
     "Called nauty: %llu,  "
     "Called nauty edge canonical: %llu (%.2f%% of nauty calls)"
     "Called nauty edge not canonical: %7llu (%.2f%% of nauty calls), "
     "Called nauty only for generators: %7llu (%.2f%% of nauty calls)\n",
     counters->totalNautyCalls, 
     counters->calledNautyEdgeCanonicalTotal, 
     100.0*counters->calledNautyEdgeCanonicalTotal / counters->totalNautyCalls,
     counters->calledNautyEdgeNotCanonicalTotal, 
     100.0*counters->calledNautyEdgeNotCanonicalTotal/counters->totalNautyCalls,
     counters->calledNautyForGeneratorsTotal, 
     100.0*counters->calledNautyForGeneratorsTotal / counters->totalNautyCalls);

    for(int i = 1; i < n/2+1; i++) {
        fprintf(stderr, 
         "Spokes: %2d, \t Called nauty: %7llu \t "
         "Edge canon: %7llu (%.2f%% of nauty calls)"
         "\t Edge not canon: %7llu (%.2f%% of nauty calls) "
         "\t For generators: %7llu (%.2f%% of nauty calls)\n",
         i, counters->calledNauty[i], counters->calledNautyEdgeCanonical[i], 
         100.0*counters->calledNautyEdgeCanonical[i] / counters->calledNauty[i],
         counters->calledNautyEdgeNotCanonical[i], 
         100.0 * counters->calledNautyEdgeNotCanonical[i] / 
          counters->calledNauty[i],
         counters->calledNautyForGenerators[i], 
         100.0 * counters->calledNautyForGenerators[i] /
          counters->calledNauty[i]);
    }
    fprintf(stderr, "Total calls to addNewCycle %llu,"
     " Calls in which we have two or more deg 2 neighbours: %7llu (%.2f%%) \n",
      counters->calledAddNewCycles, 
      counters->addNewCyclesAndMoreThanOneDeg2Nbr,
      100.0 * counters->addNewCyclesAndMoreThanOneDeg2Nbr / 
       counters->calledAddNewCycles);

    fprintf(stderr,
     "Times did not need generators because of girth: %llu (%.2f%%)\n",
     counters->didNotComputeGenerators, 
     100.0 * counters->didNotComputeGenerators / 
      (counters->didNotComputeGenerators + 
       counters->calledNautyForGeneratorsTotal));
}

int main(int argc, char ** argv) {

    struct counters counters = {0};
    struct options options = {0};
    int opt;
    char* ptr;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = 
        {
            {"girth", required_argument, NULL, 'g'},
            {"help", no_argument, NULL, 'h'},
            {"non-hamiltonian", no_argument, NULL, 'n'},
            {"permutation-method", no_argument, NULL, 'p'},
            {"verbose", no_argument, NULL, 'v'}
        };

        opt = getopt_long(argc, argv, "g:hnpv", long_options, &option_index);
        if (opt == -1) break;
        switch(opt) {
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'g': 
                options.minimalGirth = strtol(optarg, &ptr, 10);
                if(options.minimalGirth < 3) {
                    fprintf(stderr,
                     "Error: the minimal girth needs to be at least 3.\n");
                    fprintf(stderr, "%s", USAGE);
                    return 1;
                }
                break;
            case 'n':
                options.nonHamFlag = true;
                break;
            case 'p':
                options.permutationMethodFlag = true;
                break;
            case 'v':
                options.verboseFlag = true;
                fprintf(stderr,
                 "Outputting non-hamiltonian permutation graphs.\n");
                break;
            case '?':
                fprintf(stderr,"Error: Unknown option: %c\n", optopt);
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr,
                   "Use ./genPermutationGraphs --help" 
                   "for more detailed instructions.\n");
                return 1;
        }
    }

    //  First non-option argument should be number of vertices and is mandatory.
    int n;
    if(!isValidNumberOfVertices(argc, argv, &optind, &n)) {
        return 1;
    }

    //  Check for res/mod pair.
    options.modulo = 1; // Will change if command line argument is given.
    if(!isValidResModPair(argc, argv, &optind, &options)) {
        if(options.haveModResPair) {
            fprintf(stderr,
             "Error: You can only add one res/mod pair as an argument.\n");
        }
        else {
            fprintf(stderr,
             "Error: Invalid res/mod pair: '%s'.\n", argv[optind]);
        }
        fprintf(stderr, "%s", USAGE);
        fprintf(stderr,
         "Use ./genPermutationGraphs --help for more detailed instructions.\n");
        return 1;
    }

    // Cycle of length 4 gives hamiltonian cycle.
    if(options.nonHamFlag && options.minimalGirth < 5) {
        options.minimalGirth = 5;
    }

    if(options.minimalGirth > n/2) {
        fprintf(stderr, "Error: for n=%d the girth is at most %d.\n",n, n/2);
        fprintf(stderr,
         "Use ./genPermutationGraphs --help for more detailed instructions.\n");
        return 1;
    }

    // Some values which seem good based on a couple of experiments. For the
    // heavy computations a value which makes sure the common part is 5-10
    // minutes seems ideal.
    options.splitLevel = n/2 - 7 > 6 ? n/2 - 7 : 6;
    options.splitLevel = options.splitLevel < 9 ? options.splitLevel : 9;

    // Some heuristic for higher splitlevel with higher girth
    if(options.minimalGirth > 4) {
        options.splitLevel += options.minimalGirth-5;
    }

    if(options.permutationMethodFlag && BITSETSIZE > 64) {
        fprintf(stderr, "Error: -p only works for 64 bit version.\n");
        exit(1);
    } 

    fprintf(stderr, "Class=%d/%d. Splitlevel = %d.\n",
     options.remainder, options.modulo, options.splitLevel);

    //************************************************************************** 
    // Main part

    clock_t start = clock();

    // Initialize the generation.
    generatePermutationGraphs(&options, &counters, n);

    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;

    //**************************************************************************
    // Printing output

    if(options.verboseFlag) {
        printStatistics(&options, &counters, n);
    }

    fprintf(stderr,"Generated %lld graphs in %f seconds.\n",
     counters.generatedGraphs, time_spent);

    return 0;
}