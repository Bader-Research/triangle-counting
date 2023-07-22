#ifdef PARALLEL

#include "types.h"
#include "graph.h"
#include "tc_parallel.h"
#include <omp.h>


#define PBODY(foo) \
  int numThreads; \
  UINT_t *mycount; \
  _Pragma("omp parallel") \
  { \
    _Pragma("omp master") \
    numThreads = omp_get_num_threads(); \
  } \
  mycount = (UINT_t *)calloc(numThreads, sizeof(UINT_t)); \
  assert_malloc(mycount); \
  _Pragma("omp parallel") \
  { \
    int myID = omp_get_thread_num(); \
    _Pragma("omp for schedule(dynamic)") \
      foo \
  } \
  _Pragma("omp parallel for reduction(+:count)") \
  for (int i = 0; i < numThreads ; i++) \
    count += mycount[i];

#define myCount mycount[myID]


UINT_t tc_triples_P(const GRAPH_TYPE *graph) {
  /* Algorithm: for each triple (i, j, k), determine if the three triangle edges exist. */
  
  UINT_t count = 0;

  const UINT_t n = graph->numVertices;

  PBODY(
	for (UINT_t i = 0; i < n; i++)
	  for (UINT_t j = 0; j < n; j++)
	    for (UINT_t k = 0; k < n; k++)
	      if (check_edge(graph, i, j) && check_edge(graph, j, k) && check_edge(graph, k, i))
		myCount++;
	);

  return (count/6);
}

UINT_t tc_triples_DO_P(const GRAPH_TYPE *graph) {
  /* Algorithm: for each triple (i, j, k), determine if the three triangle edges exist. */
  /* Direction oriented. */
  
  UINT_t count = 0;

  const UINT_t n = graph->numVertices;

  PBODY(
	  for (UINT_t i = 0; i < n; i++)
	    for (UINT_t j = i; j < n; j++)
	      for (UINT_t k = j; k < n; k++)
		if (check_edge(graph, i, j) && check_edge(graph, j, k) && check_edge(graph, k, i))
		  myCount++;
	);

  return count;
}

UINT_t tc_wedge_P(const GRAPH_TYPE *graph) {
  /* Algorithm: For each vertex i, for each open wedge (j, i, k), determine if there's a closing edge (j, k) */
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  PBODY(
	for (UINT_t i = 0; i < n; i++) {
	  UINT_t s = Ap[i];
	  UINT_t e = Ap[i + 1];
    
	  for (UINT_t j = s; j < e; j++) {
	    UINT_t neighbor1 = Ai[j];

	    for (UINT_t k = s; k < e; k++) {
	      UINT_t neighbor2 = Ai[k];
	      
	      if (neighbor1 != neighbor2) {
		UINT_t s_n1 = Ap[neighbor1];
		UINT_t e_n1 = Ap[neighbor1 + 1];
		
		for (UINT_t l = s_n1; l < e_n1; l++) {
		  if (Ai[l] == neighbor2) {
		    myCount++;
		    break;
		  }
		}
	      }
	    }
	  }
	}
	);

  return (count/6);
}


UINT_t tc_wedge_DO_P(const GRAPH_TYPE *graph) {
  /* Algorithm: For each vertex i, for each open wedge (j, i, k), determine if there's a closing edge (j, k) */
  /* Direction oriented. */
  
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  PBODY(
	for (UINT_t i = 0; i < n; i++) {
	  UINT_t s = Ap[i];
	  UINT_t e = Ap[i + 1];

	  for (UINT_t j = s; j < e; j++) {
	    UINT_t neighbor1 = Ai[j];
	    if (neighbor1 > i) {

	      for (UINT_t k = s; k < e; k++) {
		UINT_t neighbor2 = Ai[k];

		if ((neighbor1 != neighbor2) && (neighbor2 > neighbor1)) {
		  UINT_t s_n1 = Ap[neighbor1];
		  UINT_t e_n1 = Ap[neighbor1 + 1];
	  
		  for (UINT_t l = s_n1; l < e_n1; l++) {
		    if (Ai[l] == neighbor2) {
		      myCount++;
		      break;
		    }
		  }
		}
	      }
	    }
	  }
	}
	);

  return count;
}


UINT_t tc_intersectMergePath_P(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a linear scan. */
  
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  PBODY(
	for (UINT_t v = 0; v < n; v++) {
	  UINT_t b = Ap[v  ];
	  UINT_t e = Ap[v+1];
	  for (UINT_t i=b ; i<e ; i++) {
	    UINT_t w = Ai[i];
	    myCount += intersectSizeMergePath(graph, v, w);
	  }
	}
	);

  return (count/6);
}

UINT_t tc_intersectMergePath_DO_P(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a linear scan. */
  /* Direction oriented. */
  
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  PBODY(
	for (UINT_t v = 0; v < n; v++) {
	  UINT_t b = Ap[v  ];
	  UINT_t e = Ap[v+1];
	  for (UINT_t i=b ; i<e ; i++) {
	    UINT_t w = Ai[i];
	    if (v < w)
	      myCount += intersectSizeMergePath(graph, v, w);
	  }
	}
	);

  return (count/3);
}




#endif

