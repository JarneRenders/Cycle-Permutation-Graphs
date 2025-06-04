/**
 * isPermutationGraph.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 *
 */

#define USAGE "Usage: ./isPermutationGraph [-acgo#p] [-hv]"
#define HELPTEXT "\
Filter cycle permutation graphs. We assume the input graphs are cubic.\n\
\n\
Graphs are read from stdin in graph6 format. Graphs are sent to\n\
stdout in graph6 format. If the input graph had a graph6 header, so\n\
will the output graph (if it passes through the filter).\n\
\n\
    -a, --all       compute all permutation 2-factors\n\
    -c, --count     tabulate number of permutation 2-factors of the\n\
                     input graphs; with -g also tabulates number of\n\
                     removable cycles\n\
    -g, --goddyn    output graphs which have removable cycles, i.e.\n\
                     satisfy Goddyn's conjecture\n\
    -h, --help      print help message\n\
    -o#, --output=# output graphs admitting exactly # permutation\n\
                     2-factors if -g is not present, or admitting \n\
                     exactly # removable cycles if -g is present\n\
    -p, --print     send computed permutation 2-factors to stdout\n\
    -v, --verbose   make output more verbose\n\
    res/mod         only check the ith graph if its remainder after\n\
                     dividing by mod is res; ignore the other graphs\n\
"

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <getopt.h>
#include <limits.h>
#include <time.h>
#include "utilities/readGraph6.h"
#include "utilities/bitset.h"

#define MAXCYCLES 10000 // Increase if overflow

struct graph {
    int numberOfVertices;
    bitset *adjacencyList;
    struct cycle *uniqueCycles[MAXBITSETSIZE][MAXCYCLES];
    int numUniqueCycles[MAXCYCLES];
    bool foundRemovableCycle;
};

struct options {
    bool allFlag;
    bool printFlag;
    bool countFlag;
    bool goddynFlag;
    bool outputFlag;
    bool verboseFlag;
    int modulo;
    int remainder;
    bool haveModResPair;
    int output;
};

struct counters {
    unsigned long long int frequencies[MAXBITSETSIZE];
    long long unsigned int skippedGraphs;
    long long unsigned int removableCycles;
    long long unsigned int remCycFreq[MAXCYCLES];
};

//**********************************************************************
//
//                  Input or output methods
//
//**********************************************************************

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

void printBitset(bitset set) {
    forEach(el, set) {
        fprintf(stderr, "%d ", el);
    }
    fprintf(stderr, "\n");
}

int readGraph(const char *graphString, struct graph *g,
 struct options *options, struct counters *counters) {

    g->numberOfVertices = getNumberOfVertices(graphString);

    if(g->numberOfVertices == -1 || 
     g->numberOfVertices > MAXBITSETSIZE) {

        if(options->verboseFlag){
            fprintf(stderr, "Skipping invalid graph!\n");
        }

        counters->skippedGraphs++;
        return 1;
    }

    g->adjacencyList = calloc(g->numberOfVertices, sizeof(bitset));

    // Returns -1 if there was an error.
    if(loadGraph(graphString, g->numberOfVertices, g->adjacencyList)
     == -1) {

        if(options->verboseFlag){
            fprintf(stderr, "Skipping invalid graph!\n");
        }
        counters->skippedGraphs++;
        return 1;
    }

    return 0;
}

void freeUniqueCycles(struct graph *g) {
    for(int i = 0; i < MAXBITSETSIZE; i++) {
        for(int j = 0; j < g->numUniqueCycles[i]; j++) {
            free(g->uniqueCycles[i][j]);
        }
    }
}

void freeGraph(struct graph *g) {
    freeUniqueCycles(g);
    free(g->adjacencyList);
}

//**********************************************************************
//
//              Methods for finding permutation 2-factors
//
//**********************************************************************

// Return a bitset with the vertices of the cycle containing u in the
// 2-factor G - M 
bitset uCycleInTwoFactor(struct graph *g, int u, int M[]) {

    bitset cycle = singleton(u);

    int currentVertex = u;
    int nextVertex = next(difference(g->adjacencyList[currentVertex],
     singleton(M[currentVertex])), -1);

    while(nextVertex != u) {

        add(cycle, nextVertex);

        int previousVertex = currentVertex;
        currentVertex = nextVertex;
        nextVertex = next(difference(g->adjacencyList[currentVertex],
         union(singleton(previousVertex), singleton(M[currentVertex]))),
         -1);
    }

    return cycle;
}


// Return a bitset with the vertices of the cycle containing u in the
// 2-factor G - M and print the cycle.
bitset uCycleInTwoFactorPrintCycle(struct graph *g, int u, int M[]) {

    bitset cycle = singleton(u);

    int currentVertex = u;
    int nextVertex = next(difference(g->adjacencyList[currentVertex],
     singleton(M[currentVertex])), -1);

    while(nextVertex != u) {
        printf("%d ", currentVertex);

        add(cycle, nextVertex);

        int previousVertex = currentVertex;
        currentVertex = nextVertex;
        nextVertex = next(difference(g->adjacencyList[currentVertex],
         union(singleton(previousVertex), singleton(M[currentVertex]))),
         -1);
    }
    printf("%d\n", currentVertex);

    return cycle;
}


// Check if the cycle C in G - M containing 0 is of length n/2 and
// induced. If so, then we need to check that G - M - C is a cycle of
// length n/2. (It will automatically be induced.)
bool complementConsistsOfTwoInducedCycles(struct graph *g, 
 struct options *options, struct counters *counters, int M[]) {

    bitset cycleWith0 = uCycleInTwoFactor(g, 0, M);

    // Check cycle length
    if(size(cycleWith0) != g->numberOfVertices / 2) {
        return false;
    }

    // Check if induced
    forEach(el, cycleWith0) {
        if(contains(cycleWith0, M[el])) {
            return false;
        }
    }

    // Construct a cycle in the complement
    bitset complementOfCycleWith0 = complement(cycleWith0,
     g->numberOfVertices);
    bitset cycleInComplement = uCycleInTwoFactor(g,
     next(complementOfCycleWith0, -1), M);

    // Check length
    if(size(cycleInComplement) != g->numberOfVertices / 2) {
        return false;
    }

    if(options->printFlag) {
        if(options->verboseFlag) fprintf(stderr, "Induced cycles: \n");
        uCycleInTwoFactorPrintCycle(g, 0, M);
        uCycleInTwoFactorPrintCycle(g, next(complementOfCycleWith0, -1),
         M);
    } 

    return true;
}

bool hasNewRemovableCycle(struct graph *g, struct options *options,
 struct counters *counters, int M[]);

int countEdgesPassingInPetersenWay(struct graph *g,
 struct options *options, int M[]);

//  Good perfect matching is one where the complement is a permutation
//  2-factor.
bool isPartOfGoodPerfectMatching(struct graph *g,
 struct options *options, struct counters *counters,
 bitset remainingVertices, int M[], 
 long long unsigned int *nPerm2Factors) {

    //  If this holds, M is a perfect matching. Check if its complement
    //  is a permutation 2-factor.
    int nextVertex = next(remainingVertices, -1);
    if(nextVertex == -1) {
        if(complementConsistsOfTwoInducedCycles(g, options, counters,
         M)) {

            (*nPerm2Factors)++;

            if(options->goddynFlag) {
                if(hasNewRemovableCycle(g, options, counters, M)) {
                    return true;
                };
                return false;
            }


            return true;
        }

        return false;
    }

    // Otherwise M is not yet a perfect matching. 
    forEach(neighbor, intersection(g->adjacencyList[nextVertex],
     remainingVertices)) {

        M[neighbor] = nextVertex;
        M[nextVertex] = neighbor;

        bitset newRemainingVertices = difference(remainingVertices,
         union(singleton(nextVertex), singleton(neighbor)));

        if(isPartOfGoodPerfectMatching(g, options, counters,
         newRemainingVertices, M, nPerm2Factors)) {
            if(!options->allFlag && !options->countFlag) return true;
            if(!options->countFlag && options->goddynFlag && 
             g->foundRemovableCycle) return true;
        }
    }

    return (*nPerm2Factors);
}

//**********************************************************************
//
//              Connecivity
//
//**********************************************************************

// Hopcroft and Tarjan's linear time algorithm
bool hasCutVertex(struct graph *g, bitset *checked, int depths[],
 int lowpoints[], int parent[], int u, int depth) {
    add(*checked, u);
    depths[u] = depth;
    lowpoints[u] = depth; //min depth of nbrs of all descendants of v
    int nChildren = 0;
    bool isCutVertex = false;

    bitset nbrs = g->adjacencyList[u];
    forEach(v, nbrs) {

        // If already visited update lowpoints
        if(contains(*checked, v)) {
            if(depths[v] < lowpoints[u]) {
                lowpoints[u] = depths[v];
            }
            continue;
        }

        // New child
        parent[v] = u;
        if(hasCutVertex(g, checked, depths, lowpoints, parent,
         v, depth+1)) {
            return true;
        }
        nChildren++;

        // u is cut-vertex if it has a child v with lowpoint greater
        // than u's depth.
        if(lowpoints[v] >= depths[u]) {
            isCutVertex = true;
        }

        // Update lowpoint
        if(lowpoints[v] < lowpoints[u]) {
            lowpoints[u] = lowpoints[v];
        }
    }

    // Check if root or not.
    if(parent[u] != -1) { 
        return isCutVertex;
    }
    else { // If root has at least two children it is cut vertex.
        return nChildren > 1;
    } 
}

// Check if the graph is disconnected or has a cut vertex.
bool is2Conn(struct graph *g) {
    bitset checked = EMPTY;
    int depths[MAXBITSETSIZE] = {0};
    int lowpoints[MAXBITSETSIZE] = {0};
    int parent[MAXBITSETSIZE] = {0};
    parent[0] = -1;
    bool cutVertex = hasCutVertex(g, &checked, depths, lowpoints,
     parent, 0, 0);
    bool connected = size(checked) == g->numberOfVertices; 
    return connected && !cutVertex;

}

//**********************************************************************
//
//                  Goddyn's conjecture
//
//**********************************************************************

struct cycle {
    int list[MAXBITSETSIZE];
    int len;
};

// Check if removal of cycle edges lying on perm 2-factor is
// 2-connected.
bool isRemovableCycle(struct graph *g, struct options *options,
 struct cycle *C) {

    // Remove all edges of C lying on perm 2-factor
    for(int i = 0; i < C->len; i+=2) {
        int u = C->list[i];
        int v = C->list[i+1];
        removeElement(g->adjacencyList[u], v);
        removeElement(g->adjacencyList[v], u);
    }

    bool is2Connected = is2Conn(g); 

    for(int i = 0; i < C->len; i+=2) {
        int u = C->list[i];
        int v = C->list[i+1];
        add(g->adjacencyList[u], v);
        add(g->adjacencyList[v], u);
    }

    if(is2Connected) {

        if(options->verboseFlag) {
            fprintf(stderr, "Removable:");
            for(int i = 0; i < C->len; i++) {
                fprintf(stderr, "%d ", C->list[i]);
            }
            fprintf(stderr, "\n");
        }
        g->foundRemovableCycle = true;
    }

    return is2Connected;
}

struct cycle *canoniseCycle(struct cycle *C) {

    struct cycle *canC = malloc(sizeof(struct cycle));

    canC->len = C->len;

    int min = INT_MAX;
    int minIdx = -1;
    for(int i = 0; i < C->len; i++) {
        if(C->list[i] < min) {
            minIdx = i;
            min = C->list[i];
        }
    }

    // If C[i-1] > C[i+1], loop forwards else reverse
    if(C->list[(C->len + minIdx - 1) % C->len] >
     C->list[(minIdx + 1) % C->len]) {
        for(int i = minIdx, j = 0; j < C->len; i = (i+1)%C->len, j++) {
            canC->list[j] = C->list[i];
        }
    }
    else {
        for(int i = minIdx, j = 0; j < C->len;
         i = (C->len+i-1)%C->len, j++) {
            canC->list[j] = C->list[i];
        }

    }
    return canC;
}

bool areEqual(struct cycle *C1, struct cycle *C2) {
    if(C1->len !=
     C2->len) return false;
    for(int i = 0; i < C1->len; i++) {
        if(C1->list[i] != C2->list[i]) return false;
    }
    return true;
}

bool storeCycle(struct graph *g, struct cycle *C) {
    struct cycle *canC = canoniseCycle(C);

    for(int i = 0; i < g->numUniqueCycles[C->len]; i++) {
        if(areEqual(g->uniqueCycles[C->len][i], canC)) {
            free(canC);
            return false;
        }
    }
    if(g->numUniqueCycles[C->len] >= MAXCYCLES) {
        fprintf(stderr, "Error: overflow: %d\n",
         g->numUniqueCycles[C->len]);
        exit(1);
    }
    g->uniqueCycles[C->len][g->numUniqueCycles[C->len]++] = canC;
    return true;
}

void printUniqueCycles(struct graph *g) {
    for(int i = 0; i < MAXBITSETSIZE; i++) {
        if(g->numUniqueCycles[i] == 0) continue;
        fprintf(stderr, "Length %d:\n", i);
        for(int j = 0; j < g->numUniqueCycles[i]; j++) {
            fprintf(stderr, "\t ");
            struct cycle *C = g->uniqueCycles[i][j];
            for(int k = 0; k < C->len; k++) {
                fprintf(stderr, "%d ", C->list[k]);
            }
            fprintf(stderr, "\n");
        }
    }
    fprintf(stderr, "___\n");
}


long long unsigned int extendRemovableCycle(struct graph *g,
 struct options *options, struct counters *counters, int M[], 
 bitset cycle, int start, int end, bitset checked, struct cycle *C) {

    long long unsigned int numCycles = 0;

    if(M[start] == end) {
        // We have cycle, check if removable
        if(isRemovableCycle(g, options, C)) {
            if(options->countFlag) {
                if(!storeCycle(g, C)) return 0;
            }
            counters->removableCycles++;
            return 1;
        }
        return 0;
    }

    int u = M[end];
    add(cycle, u);
    C->list[C->len++] = u;

    bitset nbrs = g->adjacencyList[u];
    forEach(nbr, nbrs) {
        if(contains(cycle, nbr)) continue;
        if(contains(checked, nbr)) continue;
        add(cycle, nbr);
        C->list[C->len++] = nbr;

        numCycles += extendRemovableCycle(g, options, counters, M,
         cycle, start, nbr, checked, C);
        if(!options->countFlag && numCycles > 0) {
            return 1;
        }

        removeElement(cycle, nbr);
        C->len--;
    }

    C->len--;
    return numCycles;
}

// M is perfect matching complement of the current perm 2-factor
bool hasNewRemovableCycle(struct graph *g, struct options *options,
 struct counters *counters, int M[]) {

    if(!options->countFlag && g->foundRemovableCycle) return false;

    long long unsigned int removableCyclesOld =
     counters->removableCycles;

    struct cycle C = {0}; // Struct keeping ordered cycle
    bitset checked = EMPTY;

    for(int i = 0; i < g->numberOfVertices; i++) {
        if(i > M[i]) continue;
        bitset cycle = union(singleton(i), singleton(M[i]));
        C.list[2] = i;
        C.list[1] = M[i];
        bitset iNbrs = g->adjacencyList[i];
        bitset MiNbrs = g->adjacencyList[M[i]];
        forEach(nbr, iNbrs) {
            if(nbr == M[i]) continue;
            if(contains(checked, nbr)) continue;
            add(cycle, nbr);
            C.list[3] = nbr;

            forEach(nbr2, MiNbrs) {
                if(nbr2 == i) continue;
                if(contains(checked, nbr2)) continue;
                add(cycle, nbr2);
                C.list[0] = nbr2;
                C.len = 4;
                if(extendRemovableCycle(g, options, counters, M, cycle,
                 nbr2, nbr, checked, &C) && !options->countFlag) {
                    return true;
                }

                removeElement(cycle, nbr2);
            }

            removeElement(cycle, nbr);
        } 

        // All cycles with (i,M[i]) have been found already. 
        add(checked,i);
        add(checked,M[i]);
    }

    long long unsigned int newRemovableCycles = 
     counters->removableCycles - removableCyclesOld;
    return newRemovableCycles;
}

//**********************************************************************
//
//                      Initialisation
//
//**********************************************************************

// Checks if the next non-option argument is of the form #1/#2. And lets
// res = #1, mod=#2. If no non-option arguments were found res = 0, mod
// = 1.
bool isValidResModPair(int argc, char ** argv, int *optind, 
 struct options *options) {

    // Check if there are still non-option arguments
    while (*optind < argc) {

        char* endptr;

        //  Return false if we already encountered a res/mod pair before
        //  and we are looking at a new non-option argument.
        if(options->haveModResPair) {
            return false;
        }

        options->remainder = strtol(argv[*optind], &endptr, 10);

        // Check if strol succeeded (endptr != 0), i.e. strol found an
        // int, the next character is '/' and that the string
        // containing the argument has a character after that.
        if( !endptr || *endptr != '/' || *(endptr+1) == '\0') {
            return false;
        }

        // Returns false if no integer was found, i.e. if endptr is same
        // as string pointer.
        if(endptr == argv[*optind]) {
            return false;
        }

        options->modulo = strtol(endptr+1, &endptr, 10);

        // Check if strol succeeded and the string containing the
        // argument has no characters left.
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

int main(int argc, char ** argv) {
    struct counters counters = {0};
    struct options options = {0};
    int opt;
    while (1) {
        int option_index = 0;
        static struct option long_options[] = 
        {
            {"all", no_argument, NULL, 'a'},
            {"goddyn", no_argument, NULL, 'g'},
            {"output", required_argument, NULL, 'o'},
            {"print", no_argument, NULL, 'p'},
            {"count", no_argument, NULL, 'c'},
            {"help", no_argument, NULL, 'h'},
            {"verbose", no_argument, NULL, 'v'}
        };

        opt = getopt_long(argc, argv, "acgho:pv", long_options,
         &option_index);
        if (opt == -1) break;
        switch(opt) {
            case 'a':
                options.allFlag = true;
                break;
            case 'c':
                options.countFlag = true;
                break;
            case 'g':
                options.goddynFlag = true;
                options.allFlag = true; // Need to check all perm 2-f.
                break;
            case 'h':
                fprintf(stderr, "%s\n", USAGE);
                fprintf(stderr, "%s", HELPTEXT);
                return 0;
            case 'o':
                options.output = 
                 (int) strtol(optarg, (char **)NULL, 10);
                options.outputFlag = true;
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
                 "Use ./isPermutationGraph --help "
                 "for more detailed instructions.\n");
                return 1;
        }
    }

    //  Check for res/mod pair.
    options.modulo = 1; 
    if(!isValidResModPair(argc, argv, &optind, &options)) {
        if(options.haveModResPair) {
            fprintf(stderr,
             "Error: You can only add one res/mod pair "
             "as an argument.\n");
        }
        else {
            fprintf(stderr,
             "Error: Invalid res/mod pair: '%s'.\n", argv[optind]);
        }
        fprintf(stderr, "%s", USAGE);
        fprintf(stderr,
         "Use ./isPermutationGraphs --help "
         "for more detailed instructions.\n");
        return 1;
    }

    if(options.goddynFlag && options.printFlag) {
        fprintf(stderr,
         "Warning: only the first permutation 2-factor with a removable"
         " cycle will be written.\n");
    }

    unsigned long long int counter = 0;
    unsigned long long int totalGraphs = 0;
    unsigned long long int passedGraphs = 0;

    clock_t start = clock();

    //  Start looping over lines of stdin.
    char * graphString = NULL;
    size_t size;
    while(getline(&graphString, &size, stdin) != -1) {

        // Read graph
        struct graph g = {0};
        if(readGraph(graphString, &g, &options, &counters) == 1) {
            continue;
        }

        // Continue if not in the right part.
        if(totalGraphs++ % options.modulo != options.remainder){
            continue; 
        } 

        counter++;

        int M[g.numberOfVertices];
        long long unsigned int nPerm2Factors = 0;
        long long unsigned int remCyclesOld = counters.removableCycles;
        if(isPartOfGoodPerfectMatching(&g, &options, &counters,
         complement(EMPTY, g.numberOfVertices), M, &nPerm2Factors)) {
            if(!options.goddynFlag && !options.countFlag) {
                passedGraphs++;
                printf("%s", graphString);
            }
        }
        if(options.countFlag) {
            if(nPerm2Factors >= MAXBITSETSIZE) {
                fprintf(stderr, "Error: overflow!\n");
                exit(1);
            }
            counters.frequencies[nPerm2Factors]++;
            if(!options.goddynFlag) {
                if(!options.outputFlag) {
                    if(nPerm2Factors > 0) {
                        passedGraphs++;
                        printf("%s", graphString);
                    }
                }
                else if(options.output == nPerm2Factors) {
                    passedGraphs++;
                    printf("%s", graphString);
                }
            }

        }
        if(options.goddynFlag) {
            long long unsigned int newRemovableCycles =
             counters.removableCycles - remCyclesOld;
             if(newRemovableCycles >= MAXCYCLES) {
                fprintf(stderr, "Error: overflow\n");
                exit(1);
            }
            counters.remCycFreq[newRemovableCycles]++;
            if(newRemovableCycles > 0) {
                if(!options.outputFlag) {
                    passedGraphs++;
                    printf("%s", graphString);
                }
                else if (options.countFlag &&
                 newRemovableCycles == options.output) {
                    passedGraphs++;
                    printf("%s", graphString);
                }
            }

        }

        freeGraph(&g);
    }
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    free(graphString);

    if(options.countFlag) {
        fprintf(stderr, "Permutation 2-factors\n");
        for(int i = 0; i < MAXBITSETSIZE; ++i) {
            if(counters.frequencies[i] != 0) {
                fprintf(stderr,
                 "\n \t\t%2d perm. 2-factors: %8lld graphs.",
                 i, counters.frequencies[i]);
            }
        }
        fprintf(stderr, "\n");
        if(options.goddynFlag) {
            fprintf(stderr, "Removable cycles:\n");
            for(int i = 0; i < MAXCYCLES; ++i) {
                if(counters.remCycFreq[i] != 0) {
                    fprintf(stderr,
                     "\n \t\t%2d removable cycles: %8lld graphs.",
                     i, counters.remCycFreq[i]);
                }
            }
            fprintf(stderr, "\n");
        }
    } 

    fprintf(stderr,
        "\rChecked %lld graphs in %f seconds: %llu passed.\n",
        counter, time_spent, passedGraphs);

    if(counters.skippedGraphs > 0) {
        fprintf(stderr, "Warning: %lld graphs were skipped.\n",
         counters.skippedGraphs);
    }

    return 0;
}