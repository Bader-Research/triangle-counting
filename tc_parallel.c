#ifdef PARALLEL

#include "types.h"
#include "graph.h"
#include "tc_parallel.h"
#include <omp.h>

UINT_t tc_triples_P(const GRAPH_TYPE *graph) {
  /* Algorithm: for each triple (i, j, k), determine if the three triangle edges exist. */
  
  UINT_t count = 0;
  register UINT_t i, j, k;

  const UINT_t n = graph->numVertices;

  int numThreads;
  UINT_t *mycount;
  
#pragma omp parallel
  {
#pragma omp master
    numThreads = omp_get_num_threads();
  }
  mycount = (UINT_t *)calloc(numThreads, sizeof(UINT_t));
  assert_malloc(mycount);
  
#pragma omp parallel private(i, j, k)
  {
    int myID = omp_get_thread_num();
#pragma omp for
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
	for (k = 0; k < n; k++)
	  if (check_edge(graph, i, j) && check_edge(graph, j, k) && check_edge(graph, k, i))
	    mycount[myID]++;
  }

#pragma omp parallel for reduction(+:count)
  for (i = 0; i < numThreads ; i++)
    count += mycount[i];
  
  return (count/6);
}

UINT_t tc_triples_DO_P(const GRAPH_TYPE *graph) {
  /* Algorithm: for each triple (i, j, k), determine if the three triangle edges exist. */
  /* Direction oriented. */
  
  UINT_t count = 0;
  register UINT_t i, j, k;

  const UINT_t n = graph->numVertices;

  int numThreads;
  UINT_t *mycount;

#pragma omp parallel
  {
#pragma omp master
    numThreads = omp_get_num_threads();
  }
  mycount = (UINT_t *)calloc(numThreads, sizeof(UINT_t));
  assert_malloc(mycount);

#pragma omp parallel private(i, j, k)
  {
    int myID = omp_get_thread_num();
#pragma omp for
    for (i = 0; i < n; i++)
      for (j = i; j < n; j++)
	for (k = j; k < n; k++)
	  if (check_edge(graph, i, j) && check_edge(graph, j, k) && check_edge(graph, k, i))
	    mycount[myID]++;
  }

#pragma omp parallel for reduction(+:count)
  for (i = 0; i < numThreads ; i++)
    count += mycount[i];
  
  return count;
}



#endif
