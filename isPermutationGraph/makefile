compiler=gcc
flags=-std=gnu11 -march=native -Wall -Wno-missing-braces -O3
profileflags=-std=gnu11 -march=native -Wall -fsanitize=address -g -pg

# The 64-bit version of this program is faster but only supports graphs up to 64 vertices.
64bit: isPermutationGraph.c utilities/readGraph6.c utilities/bitset.h 
	$(compiler) -DUSE_64_BIT -o isPermutationGraph isPermutationGraph.c utilities/readGraph6.c $(flags)

# There are two different implementations of the 128-bit version. The array version generally performs faster.
128bit: isPermutationGraph.c utilities/readGraph6.c utilities/bitset.h 
	$(compiler) -DUSE_128_BIT -o isPermutationGraph-128 isPermutationGraph.c utilities/readGraph6.c $(flags)

128bitarray: isPermutationGraph.c utilities/readGraph6.c utilities/bitset.h 
	$(compiler) -DUSE_128_BIT_ARRAY -o isPermutationGraph-128a isPermutationGraph.c utilities/readGraph6.c $(flags)	

profile: isPermutationGraph.c utilities/readGraph6.c utilities/bitset.h 
	$(compiler) -DUSE_64_BIT -o isPermutationGraph-pr isPermutationGraph.c utilities/readGraph6.c $(profileflags)

all: 64bit 128bit 128bitarray

.PHONY: clean
clean:
	rm -f isPermutationGraph isPermutationGraph-128 isPermutationGraph-128a isPermutationGraph-pr

