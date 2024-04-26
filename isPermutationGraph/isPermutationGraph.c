/**
 * isPermutationGraph.c
 * 
 * Author: Jarne Renders (jarne.renders@kuleuven.be)
 *
 */

#define USAGE "Usage: ./isPermutationGraph [-a] [-hv]"
#define HELPTEXT "\
Filter cycle permutation graphs. We assume the input graphs are cubic.\n\
\n\
Graphs are read from stdin in graph6 format. Graphs are sent to\n\
stdout in graph6 format. If the input graph had a graph6 header, so\n\
will the output graph (if it passes through the filter).\n\
\n\
    -a, --all       compute all permutation 2-factors and output the\n\
                     induced cycles line per line to stdout before \n\
                     outputting the graph\n\
    -h, --help      print help message\n\
    -v, --verbose   make output more verbose\n\
    res/mod         only check the ith graph if its remainder after\n\
                     dividing by mod is res; ignore the other graphs\n\
"

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <getopt.h>
#include <time.h>
#include "utilities/readGraph6.h"
#include "utilities/bitset.h"

struct graph {
    int numberOfVertices;
    bitset *adjacencyList;
};

struct options {
    bool allFlag;
    bool verboseFlag;
    int modulo;
    int remainder;
    bool haveModResPair;
};

struct counters {
    long long unsigned int skippedGraphs;
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

    g->adjacencyList = malloc(sizeof(bitset)*g->numberOfVertices);

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

void freeGraph(struct graph *g) {
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

    if(options->allFlag) {
        if(options->verboseFlag) fprintf(stderr, "Induced cycles: \n");
        uCycleInTwoFactorPrintCycle(g, 0, M);
        uCycleInTwoFactorPrintCycle(g, next(complementOfCycleWith0, -1),
         M);
    } 

    return true;
}

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
            if(!options->allFlag) return true;
        }
    }

    return (*nPerm2Factors);
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
            {"help", no_argument, NULL, 'h'},
            {"verbose", no_argument, NULL, 'v'}
        };

        opt = getopt_long(argc, argv, "ahv", long_options,
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

    unsigned long long int counter = 0;
    unsigned long long int totalGraphs = 0;
    unsigned long long int passedGraphs = 0;

    clock_t start = clock();

    //  Start looping over lines of stdin.
    char * graphString = NULL;
    size_t size;
    while(getline(&graphString, &size, stdin) != -1) {

        // Read graph
        struct graph g;
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
        if(isPartOfGoodPerfectMatching(&g, &options, &counters,
         complement(EMPTY, g.numberOfVertices), M, &nPerm2Factors)) {
            passedGraphs++;
            printf("%s", graphString);
        }

        freeGraph(&g);
    }
    clock_t end = clock();
    double time_spent = (double)(end - start) / CLOCKS_PER_SEC;
    free(graphString);

    fprintf(stderr,
        "\rChecked %lld graphs in %f seconds: %llu passed.\n",
        counter, time_spent, passedGraphs);

    if(counters.skippedGraphs > 0) {
        fprintf(stderr, "Warning: %lld graphs were skipped.\n",
         counters.skippedGraphs);
    }

    return 0;
}