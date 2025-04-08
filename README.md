# Cycle Permutation Graphs
This repository contains two programs which were created for the article "J. Goedgebeur, J. Renders, Generation of Cycle Permutation Graphs and Permutation Snarks, manuscript".

In the folder `genPermutationGraphs` one can find the program `genPermutationGraphs`, which allows for the generation of all pairwise non-isomorphic (non-hamiltonian) cycle permutation graphs of a given order `n` (with at least a given girth).

In the folder `isPermutationGraph`, one can find the program `isPermutationGraph`, which can be used to filter cycle permutation graphs.

Each of these subfolders contain their own readme as well as a makefile for compiling the programs.

These programs all use Brendan McKay's graph6 format to read and write graphs. See <http://users.cecs.anu.edu.au/~bdm/data/formats.txt>.

The latest version of the programs can be obtained from: <https://github.com/JarneRenders/Cycle-Permutation-Graphs/>.

## genPermutationGraphs
### Short manual

The program supports the generation of graphs up to and including 128 vertices.

### Installation

This requires a working shell and `make`. 

- Download, extract and configure the latest version of [`nauty`](https://pallini.di.uniroma1.it/).
- Copy the following files to the `genPermutationGraphs/utilities`: 
	* gtools.h, 
	* naurng.h, 
	* nausparse.h
	* naututil.h, 
	* nauty.a, 
	* nauty.h
- Navigate to the folder `genPermutationGraphs` containing `genPermutationGraphs.c` and compile using: 
	* `make` to create a binary for the 64 bit version
	* `make 128bit` to create a binary for the 128 bit version 
	* `make 128bitarray` to create a binary for an alternative 128 bit version 
	* `make all` to create all of the above

The 64-bit version supports graphs only up to 64 vertices, the 128-bit versions up to 128 vertices. For graphs containing up to 64 vertices, the 64-bit version performs significantly faster than the 128-bit versions. Use `make clean` to remove all binaries created in this way.

### Usage of genPermutationGraphs 

All options can be found by executing `./genPermutationGraphs -h`.

Usage: `./genPermutationGraphs [-hv] [-nopsSg#] n [res/mod]`

Generate all pairwise non-isomorphic cycle permutation graphs of a given order n. Allows to restrict to graphs with a given lower bound on the girth, non-hamiltonian graphs and permutation snarks generated. Generated graphs are sent to stdout in graph6 format. For more information on the format, see <http://users.cecs.anu.edu.au/~bdm/data/formats.txt>.

The `res/mod` argument, should always appear after the specified order `n`. Otherwise, the order in which the arguments appear does not matter. Be careful not to put an argument immediately after one with an option. E.g. -g#n will not recognize the -n argument.

Without optional arguments, the program will output cycle permutation graphs of order n using the canonical construction path method to generate them. The same result can be obtained by adding -o, but in this case the faster orderly generation approach will be used for the generation.

Mandatory arguments to long options are mandatory for short options too.
```
    -h, --help             print help message
    -g, --girth=#          only generate cycle permutation graphs of girth at
                            least #
    -n, --non-hamiltonian  only generate non-hamiltonian cycle permutation
                            graphs
    -o, --strong-orderly-generation
                           uses the strong orderly generation approach, which
                            is faster than the canonical construction path method
                            and outputs no isomorphic copies; slightly slower
                            than using only -p, which does output isomorphic
                            copies
    -p, --permutations     uses the weak orderly generation approach; this
                            method WILL print out isomorphic copies unless
                            paired with -o, but is faster
    -s, --snark-heuristic  can only be used with -p or -o; uses a heuristic
                            to remove non-permutation snarks from the output;
                            using -ps and filtering the output afterwards is the
                            fastest way to obtain all permutation snarks
    -S, --snark            can only be used with -p or -o; removes all
                            non-permutation snarks from the output, but is slower
                            than only using -s
    -v, --verbose          print out extra statistics
    res/mod                split the generation in mod (not necessarily
                            equally large) parts. Here part res will be
                            executed. Splitting will cause extra overhead.
```

### Examples

`./genPermutationGraphs 10`
Generate all cycle permutation graphs of order 10 and send them to stdout.

`./genPermutationGraphs -o 10`
Generate all cycle permutation graphs of order 10 and send them to stdout. This method is slightly faster.

`./genPermutationGraphs -n 18`
Generate all non-hamiltonian cycle permutation graphs of order 18 and send them to stdout.

`./genPermutationGraphs -g8 10`
Generate all cycle permutation graphs of order 10 and girth at least 8 and send them to stdout.

`./genPermutationGraphs -n -g7 18`
Generate all non-hamiltonian cycle permutation graphs of order 18 and girth at least 7 and send them to stdout.

`./genPermutationGraphs -pn 18`
Generate all non-hamiltonian cycle permutation graphs of order 18 and send them to stdout. This is slightly faster than without -p, but graphs which are isomorphic will most likely be output.

`./genPermutationGraphs -on 18`
Generate all non-hamiltonian cycle permutation graphs of order 18 and send them to stdout. This is slightly slower than the previous method but faster than without -o. No isomorphic copies are output.

`./genPermutationGraphs -ps 18`
Generate all permutation snarks of order 18 and send them to stdout. Isomorphic copies will be output. In theory, this method can output graphs which are not permutation snarks, but in practice this does not happen up to order at least 46. This is the fastest way to generate permutation snarks.

`./genPermutationGraphs -pS 18`
Generate all permutation snarks of order 18 and send them to stdout. Isomorphic copies will be output. This will only output permutation snarks, but is slightly slower than using -s.

`./genPermutationGraphs -oS 18`
Generate all permutation snarks of order 18 and send them to stdout. No isomorphic copies will be output, but is slower than using -p. 

## isPermutationGraph
### Short manual

This program supports graphs up to and including 128 vertices.

### Installation

This requires a working shell and `make`. 

- Navigate to the folder `isPermutationGraph` containing `isPermutationGraph.c` and compile using: 
	* `make` to create a binary for the 64 bit version
	* `make 128bit` to create a binary for the 128 bit version 
	* `make 128bitarray` to create a binary for an alternative 128 bit version 
	* `make all` to create all of the above

The 64-bit version supports graphs only up to 64 vertices, the 128-bit versions up to 128 vertices. For graphs containing up to 64 vertices, the 64-bit version performs significantly faster than the 128-bit versions. Use `make clean` to remove all binaries created in this way.

### Usage of isPermutationGraph 

All options can be found by executing `./isPermutationGraph -h`.

Usage: `./isPermutationGraph [-a] [-hv]`

Filter cycle permutation graphs. We assume the input graphs are cubic.

Graphs are read from stdin in graph6 format. Graphs are sent to stdout in graph6 format. If the input graph had a graph6 header, so will the output graph (if it passes through the filter).

```
    -a, --all       compute all permutation 2-factors and output the
                     induced cycles line per line to stdout before 
                     outputting the graph
    -h, --help      print help message
    -v, --verbose   make output more verbose
    res/mod         only check the ith graph if its remainder after
                     dividing by mod is res; ignore the other graphs
```
### Examples

Let graphs.g6 be a file containing cubic graphs in graph6 format in the `isPermutationGraph` folder.

`./genPermutationGraphs < graphs.g6`
Sends all graphs in `graphs.g6` which are cycle permutation graphs to stdout.

`./genPermutationGraphs -a < graphs.g6`
Sends all graphs in `graphs.g6` which are cycle permutation graphs to stdout and for each cycle permutation graph the permutation 2-factors are sent to stdout.
