compiler=gcc
flags=-std=gnu11 -march=native -Wall -Wno-missing-braces -O3
profileflags=-std=gnu11 -march=native -Wall -fsanitize=address -g -pg

# The 64-bit version of this program is faster but only supports graphs up to 64 vertices.
64bit: genPermutationGraphs.c utilities/bitset.h 
	$(compiler) -DUSE_64_BIT -o genPermutationGraphs genPermutationGraphs.c utilities/nauty.a $(flags)

# sparse-64bit: genPermutationGraphs.c utilities/bitset.h  
# 	$(compiler) -DUSE_64_BIT -DSPARSE_FORMAT -o genPermutationGraphs-sp genPermutationGraphs.c utilities/nausparse.c utilities/nauty.a $(flags)

# There are two different implementations of the 128-bit version. The array version generally performs faster.
128bit: genPermutationGraphs.c utilities/bitset.h 
	$(compiler) -DUSE_128_BIT -o genPermutationGraphs-128 genPermutationGraphs.c utilities/nauty.a $(flags)

128bitarray: genPermutationGraphs.c utilities/bitset.h
	$(compiler) -DUSE_128_BIT_ARRAY -o genPermutationGraphs-128a genPermutationGraphs.c utilities/nauty.a $(flags)	

profile: genPermutationGraphs.c utilities/bitset.h  
	$(compiler) -DUSE_64_BIT -o genPermutationGraphs-pr genPermutationGraphs.c utilities/nauty.a $(profileflags)

split: genPermutationGraphs.c utilities/bitset.h 
	$(compiler) -DUSE_64_BIT -DCOUNT_SPLIT_TIME -o genPermutationGraphs-sp genPermutationGraphs.c utilities/nauty.a $(flags)

all: 64bit 128bit 128bitarray

.PHONY: clean
clean:
	rm -f genPermutationGraphs genPermutationGraphs-128 genPermutationGraphs-128a genPermutationGraphs-pr genPermutationGraphs-sp

