#ifndef BITSETCHOOSER
#define BITSETCHOOSER

#ifdef USE_64_BIT
	#include "bitset64Vertices.h"
	#define MAXBITSETSIZE 64

#elif defined(USE_128_BIT)
	#include "bitset128Vertices.h"
	#define MAXBITSETSIZE 128

#elif defined(USE_128_BIT_ARRAY)
	#include "bitset128VerticesArray.h"
	#define MAXBITSETSIZE 128

#endif

#endif
