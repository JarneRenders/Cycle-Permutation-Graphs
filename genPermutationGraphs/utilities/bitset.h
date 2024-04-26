#ifndef BITSETCHOOSER
#define BITSETCHOOSER

#ifdef USE_64_BIT
	#include "bitset64Vertices.h" // IWYU pragma: export
	#define BITSETSIZE 64

#elif defined(USE_128_BIT)
	#include "bitset128Vertices.h" // IWYU pragma: export
	#define BITSETSIZE 128

#elif defined(USE_128_BIT_ARRAY)
	#include "bitset128VerticesArray.h" // IWYU pragma: export
	#define BITSETSIZE 128

#endif

#endif
