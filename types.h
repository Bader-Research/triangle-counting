#ifndef _TYPES_H
#define _TYPES_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <strings.h>
#include <stdbool.h>
#include <sys/time.h>
#include <stdbool.h>

#ifdef ICX
#define UINT_t uint
#define INT_t int
#endif

#ifndef UINT_t
#define UINT_t uint32_t
#endif
#ifndef INT_t
#define INT_t int32_t
#endif

typedef struct {
    UINT_t numVertices;
    UINT_t numEdges;
    UINT_t* rowPtr;
    UINT_t* colInd;
} GRAPH_TYPE;

typedef struct {
  UINT_t src;
  UINT_t dst;
} edge_t;


static struct timeval  tp;
static struct timezone tzp;
#define get_seconds()   (gettimeofday(&tp, &tzp),			\
                        (double)tp.tv_sec + (double)tp.tv_usec / 1000000.0)

#define ODD(n) ((n)&1)==1
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

static void assert_malloc(const void *ptr) {
  if (ptr==NULL) {
    fprintf(stderr,"ERROR: Null pointer\n");
    exit(1);
  }
}

extern bool QUIET;
extern bool PRINT;
extern bool NCUBED;

#ifdef PARALLEL
extern bool BENCHMARK_BFS;
extern bool PARALLEL_MAX;
extern int  PARALLEL_PROCS;
#endif

#ifdef GCC
#define INLINE inline
/* #define INLINE */
#else
#define INLINE
#endif


#endif
