#ifdef PARALLEL

#include "types.h"
#include "graph.h"
#include "bfs.h"
#include "tc_parallel.h"
#include <unistd.h>
#include <omp.h>

#define PBODY(foo)				\
  int numThreads;				\
  UINT_t *mycount;				\
  _Pragma("omp parallel")			\
  {							      \
  int myID = omp_get_thread_num();			      \
  if (myID == 0) {					      \
    numThreads = omp_get_num_threads();			      \
    mycount = (UINT_t *)calloc(numThreads, sizeof(UINT_t));   \
    assert_malloc(mycount);				      \
  }							      \
  _Pragma("omp barrier")				      \
  _Pragma("omp for schedule(dynamic)")			      \
  foo							      \
  _Pragma("omp for reduction(+:count)")			      \
  for (int i = 0; i < numThreads ; i++)			      \
    count += mycount[i];				      \
  }


#define PBODY2(foo, bar)				\
  int numThreads;					\
  UINT_t *mycount;					\
  _Pragma("omp parallel")				\
  {							\
  int myID = omp_get_thread_num();			      \
  if (myID == 0) {					      \
      numThreads = omp_get_num_threads();		      \
      mycount = (UINT_t *)calloc(numThreads, sizeof(UINT_t)); \
      assert_malloc(mycount);				      \
      bar						      \
  }							      \
  _Pragma("omp barrier")				      \
  _Pragma("omp for schedule(dynamic)")			      \
  foo							      \
  _Pragma("omp for reduction(+:count)")			      \
  for (int i = 0; i < numThreads ; i++)			      \
    count += mycount[i];				      \
  }

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


UINT_t tc_intersectBinarySearch_P(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a binary search. */

  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  PBODY(
	for (UINT_t v = 0; v < n; v++) {
	  UINT_t b = Ap[v  ];
	  UINT_t e = Ap[v+1];
	  for (UINT_t i=b ; i<e ; i++) {
	    UINT_t w  = Ai[i];
	    myCount += intersectSizeBinarySearch(graph, v, w);
	  }
	}
	);

  return (count/6);
}

UINT_t tc_intersectBinarySearch_DO_P(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a binary search. */
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
	    UINT_t w  = Ai[i];
	    if (v < w)
	      myCount += intersectSizeBinarySearch(graph, v, w);
	  }
	}
	);

  return (count/3);
}


UINT_t tc_intersectPartition_P(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a binary search-based partition. */
  
  UINT_t count = 0;

  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  PBODY(
	for (UINT_t v = 0; v < n ; v++) {
	  UINT_t b = Ap[v  ];
	  UINT_t e = Ap[v+1];
	  for (UINT_t i=b ; i<e ; i++) {
	    UINT_t w  = Ai[i];
	    myCount += searchLists_with_partitioning((UINT_t *)Ai, (INT_t) Ap[v], (INT_t)Ap[v+1]-1, (UINT_t *)Ai, (INT_t)Ap[w], (INT_t)Ap[w+1]-1);
	  }
	}
	);

  return (count/6);
}

UINT_t tc_intersectPartition_DO_P(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a binary search-based partition. */
  /* Direction oriented. */
  
  UINT_t count = 0;

  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  PBODY(
	for (UINT_t v = 0; v < n; v++) {
	  UINT_t b = Ap[v  ];
	  UINT_t e = Ap[v+1];
	  for (UINT_t i=b ; i<e ; i++) {
	    UINT_t w  = Ai[i];
	    if (v < w)
	      myCount += searchLists_with_partitioning((UINT_t *)Ai, (INT_t) Ap[v], (INT_t)Ap[v+1]-1, (UINT_t *)Ai, (INT_t)Ap[w], (INT_t)Ap[w+1]-1);
	  }
	}
	);

  return (count/3);
}


UINT_t tc_intersectHash_P(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a hash. */

  UINT_t count = 0;

  bool *Hash;
  
  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  PBODY2(
	for (UINT_t v = 0; v < n ; v++) {
	  UINT_t b = Ap[v  ];
	  UINT_t e = Ap[v+1];
	  for (UINT_t i=b ; i<e ; i++) {
	    UINT_t w  = Ai[i];
	    myCount += intersectSizeHash(graph, Hash + (myID * n), v, w);
	  }
	}
	,
	 Hash = (bool *)calloc(n * omp_get_num_threads(), sizeof(bool));
	 assert_malloc(Hash);
	);

  free(Hash);

  return (count/6);
}




UINT_t tc_intersectHash_DO_P(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a hash. */
  /* Direction oriented. */

  UINT_t count = 0;

  bool *Hash;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

  PBODY2(
	 for (UINT_t v = 0; v < n ; v++) {
	   UINT_t b = Ap[v  ];
	   UINT_t e = Ap[v+1];
	   for (UINT_t i=b ; i<e ; i++) {
	     UINT_t w  = Ai[i];
	     if (v < w)
	       myCount += intersectSizeHash(graph, Hash + (myID * n), v, w);
	   }
	   }
	 ,
	 Hash = (bool *)calloc(n * omp_get_num_threads(), sizeof(bool));
	 assert_malloc(Hash);
	 );

  free(Hash);

  return (count/3);
}



void bfs_mark_horizontal_edges_P(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* restrict level, Queue* queue, bool* visited, bool* horiz) {
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  UINT_t listsize;

  UINT_t *list = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(list);

  visited[startVertex] = true;
#if 0
  enqueue(queue, startVertex);
#else
  UINT_t liststart = 0;
  list[liststart] = startVertex;
  UINT_t listend = 1;
  listsize = listend - liststart;

  UINT_t curlistend = listend;
#endif
  level[startVertex] = 1;

  while (
#if 0
  !isEmpty(queue)
#else
  listsize > 0
#endif
	 ) {
#pragma omp parallel for schedule(dynamic)
    for (UINT_t i = liststart ; i<listend ; i++) {
#if 0
      UINT_t v = dequeue(queue);
#else
      listsize = curlistend - liststart;
      UINT_t v = list[i];
#endif
      for (UINT_t i = Ap[v]; i < Ap[v + 1]; i++) {
	UINT_t w = Ai[i];
	if (!visited[w])  {
	  horiz[i] = false;
	  visited[w] = true;
#if 0
	  enqueue(queue, w);
#else
#pragma omp critical
	  {
	    list[listend++] = w;
	  }
#endif
	  level[w] = level[v] + 1;
	}
	else {
	  horiz[i] = (level[w] == 0) || (level[w] == level[v]);
	}
      }
    }
#pragma omp single
    {
      curlistend = listend;
      liststart += listsize;
    }
#pragma omp barrier
  }
  
  free(list);
}



UINT_t tc_bader_bfs1_P(const GRAPH_TYPE *graph) {
  /* Bader's new algorithm for triangle counting based on BFS */
  /* Uses Hash array to detect triangles (v, w, x) if x is adjacent to v */
  /* For level[], 0 == unvisited. Needs a modified BFS starting from level 1 */
  /* Mark horizontal edges during BFS */
  /* Direction orientied. */
  UINT_t* restrict level;
  UINT_t c1, c2;
  static bool *Hash;
  bool *horiz;
  bool *visited;
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

  level = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(level);
  
  visited = (bool *)calloc(n, sizeof(bool));
  assert_malloc(visited);

  horiz = (bool *)malloc(m * sizeof(bool));
  assert_malloc(horiz);

  Queue *queue = createQueue(n);

  for (UINT_t v = 0 ; v < n ; v++) {
    if (!level[v])
      bfs_mark_horizontal_edges(graph, v, level, queue, visited, horiz);
  }

  static int numThreads;
  static UINT_t *myc1;
  static UINT_t *myc2;

  c1 = 0; c2 = 0;
#pragma omp parallel \
  shared(numThreads, myc1, myc2, Hash, Ap, Ai, n, m, level)
  {
    int myID = omp_get_thread_num();
    if (myID==0) {
      numThreads = omp_get_num_threads();
      
      myc1 = (UINT_t *)calloc(numThreads, sizeof(UINT_t));
      assert_malloc(myc1);
      myc2 = (UINT_t *)calloc(numThreads, sizeof(UINT_t));
      assert_malloc(myc2);
  
      Hash = (bool *)calloc(n * numThreads, sizeof(bool));
      assert_malloc(Hash);
    }
#pragma omp barrier
#pragma omp for schedule(dynamic)
    for (UINT_t v = 0 ; v < n ; v++) {
      bool *myHash = Hash + (n * myID);
    
      const UINT_t s = Ap[v  ];
      const UINT_t e = Ap[v+1];
      const UINT_t l = level[v];

      for (UINT_t p = s ; p<e ; p++)
	myHash[Ai[p]] = true;

      for (UINT_t j = s ; j<e ; j++) {
	if (horiz[j]) {
	  const UINT_t w = Ai[j];
	  if (v < w) {
	    for (UINT_t k = Ap[w]; k < Ap[w+1] ; k++) {
	      UINT_t x = Ai[k];
	      if (myHash[x]) {
		if (level[x] != l) {
		  myc1[myID]++;
		}
		else {
		  myc2[myID]++;
		}
	      }
	    }
	  }
	}
      }

      for (UINT_t p = s ; p<e ; p++)
	myHash[Ai[p]] = false;
    }

#pragma omp for reduction(+:c1,c2)
    for (int i = 0; i < numThreads ; i++) {
      c1 += myc1[i];
      c2 += myc2[i];
    }
  }

  free_queue(queue);

  free(Hash);
  free(visited);
  free(level);
  free(horiz);

  return c1 + (c2/3);
}


static UINT_t tc_bader_bfs_core_P(const GRAPH_TYPE* graph, void (*f)(const GRAPH_TYPE*, const UINT_t, UINT_t *, bool *)) {
  /* Bader's new algorithm for triangle counting based on BFS */
  /* Uses Hash array to detect triangles (v, w, x) if x is adjacent to v */
  /* For level[], 0 == unvisited. Needs a modified BFS starting from level 1 */
  /* Mark horizontal edges during BFS */
  /* Direction orientied. */
  UINT_t* restrict level;
  UINT_t c1, c2;
  static bool *Hash;
  bool *visited;

  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

  level = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(level);
  
  visited = (bool *)calloc(n, sizeof(bool));
  assert_malloc(visited);

  Queue *queue = createQueue(n);

  for (UINT_t v = 0 ; v < n ; v++) {
    if (!visited[v])
      (*f)(graph, v, level, visited);
  }

  static int numThreads;
  static UINT_t *myc1;
  static UINT_t *myc2;

  c1 = 0; c2 = 0;
#pragma omp parallel \
  shared(numThreads, myc1, myc2, Hash, Ap, Ai, n, m, level)
  {
    int myID = omp_get_thread_num();
    if (myID==0) {
      numThreads = omp_get_num_threads();
      
      myc1 = (UINT_t *)calloc(numThreads, sizeof(UINT_t));
      assert_malloc(myc1);
      myc2 = (UINT_t *)calloc(numThreads, sizeof(UINT_t));
      assert_malloc(myc2);
  
      Hash = (bool *)calloc(n * numThreads, sizeof(bool));
      assert_malloc(Hash);
    }
#pragma omp barrier
#pragma omp for schedule(dynamic)
    for (UINT_t v = 0 ; v < n ; v++) {
      bool *myHash = Hash + (n * myID);
    
      const UINT_t s = Ap[v  ];
      const UINT_t e = Ap[v+1];
      const UINT_t l = level[v];

      for (UINT_t p = s ; p<e ; p++)
	myHash[Ai[p]] = true;

      for (UINT_t j = s ; j<e ; j++) {
	const UINT_t w = Ai[j];
	if ((v<w) && (l == level[w])) {
	  for (UINT_t k = Ap[w]; k < Ap[w+1] ; k++) {
	    UINT_t x = Ai[k];
	    if (myHash[x]) {
	      if (level[x] != l) {
		myc1[myID]++;
	      }
	      else {
		myc2[myID]++;
	      }
	    }
	  }
	}
      }

      for (UINT_t p = s ; p<e ; p++)
	myHash[Ai[p]] = false;
    }

#pragma omp for reduction(+:c1,c2)
    for (int i = 0; i < numThreads ; i++) {
      c1 += myc1[i];
      c2 += myc2[i];
    }
  }
  
  free_queue(queue);

  free(visited);
  free(Hash);
  free(level);

  return c1 + (c2/3);
}


UINT_t tc_bader_bfs3_P(const GRAPH_TYPE *graph) {
  return tc_bader_bfs_core_P(graph, bfs_visited);
}

UINT_t tc_bader_bfs_visited_P(const GRAPH_TYPE *graph) {
  return tc_bader_bfs_core_P(graph, bfs_visited_P);
}

UINT_t tc_bader_bfs_hybrid_P(const GRAPH_TYPE *graph) {
  return tc_bader_bfs_core_P(graph, bfs_hybrid_visited);
}

UINT_t tc_bader_bfs_hybrid2_P(const GRAPH_TYPE *graph) {
  return tc_bader_bfs_core_P(graph, bfs_hybrid_visited_P);
}

UINT_t tc_bader_bfs_chatgpt_P(const GRAPH_TYPE *graph) {
  return tc_bader_bfs_core_P(graph, bfs_chatgpt_P);
}

UINT_t tc_bader_bfs_locks_P(const GRAPH_TYPE *graph) {
  return tc_bader_bfs_core_P(graph, bfs_locks_P);
}


UINT_t tc_forward_hash_P(const GRAPH_TYPE *graph) {
  
/* Schank, T., Wagner, D. (2005). Finding, Counting and Listing All Triangles in Large Graphs, an Experimental Study. In: Nikoletseas, S.E. (eds) Experimental and Efficient Algorithms. WEA 2005. Lecture Notes in Computer Science, vol 3503. Springer, Berlin, Heidelberg. https://doi.org/10.1007/11427186_54 */

  register UINT_t s, t;
  register UINT_t b, e;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

  bool* Hash = (bool *)calloc(m, sizeof(bool));
  assert_malloc(Hash);

  UINT_t* Size = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(Size);
  
  UINT_t* A = (UINT_t *)calloc(m, sizeof(UINT_t));
  assert_malloc(A);

  for (s = 0; s < n ; s++) {
    b = Ap[s  ];
    e = Ap[s+1];
    for (UINT_t i=b ; i<e ; i++) {
      t  = Ai[i];
      if (s<t) {
	count += intersectSizeHash_forward_P(graph, Hash, s, t, A, Size);
	A[Ap[t] + Size[t]] = s;
	Size[t]++;
      }
    }
  }

  free(A);
  free(Size);
  free(Hash);
  
  return count;
}



int intCompare(const void * a, const void * b) {
   return ( *(int*)a - *(int*)b );
}


/* MapJIK is the triangle counting algorithm in this paper:

   "Exploring Optimizations on Shared-memory Platforms for Parallel Triangle Counting Algorithms." Ancy Sarah Tom, Narayanan Sundaram, Nesreen K. Ahmed, Shaden Smith, Stijn Eyerman, Midhunchandra Kodiyath, Ibrahim Hur, Fabrizio Petrini, and George Karypis. IEEE High Performance Extreme Computing Conference (HPEC), 2017.

   The source code is adapted from: https://github.com/KarypisLab/TriangleCounting
*/


#define GK_SBSIZE  64
#define GK_DBSIZE  8

/*************************************************************************/
/*! Reorders the vertices in the graph in inc degree order and returns
    the re-ordered graph in which the adjancency lists are sorted in 
    increasing order. In addition, a diagonal entry is added to each
    row to make the code that follows tighter.
*/
/*************************************************************************/
GRAPH_TYPE *ptc_Preprocess(const GRAPH_TYPE *originalGraph)
{
  int32_t vi, nvtxs, nthreads, maxdegree=0, csrange=0;
#if 0
  ssize_t *xadj, *nxadj, *psums;
#else
  int32_t *xadj, *nxadj, *psums;
#endif
  int32_t *adjncy, *nadjncy, *perm=NULL, *iperm=NULL, *chunkptr=NULL;
  int32_t *gcounts;
  GRAPH_TYPE *graph;

#pragma omp parallel
  {
    #pragma omp single
    nthreads = omp_get_num_threads();
  }

  nvtxs  = (int32_t)originalGraph->numVertices;
#if 0
  xadj   = (ssize_t *)originalGraph->rowPtr;
#else
  xadj   = (int32_t *)originalGraph->rowPtr;
#endif
  adjncy = (int32_t *)originalGraph->colInd;

  graph = (GRAPH_TYPE *)malloc(sizeof(GRAPH_TYPE));
  assert_malloc(graph);
  graph->numVertices = nvtxs;
  graph->numEdges = /* nvtxs+xadj[nvtxs] */ xadj[nvtxs];
  /* graph->rowPtr   = nxadj = (ssize_t *)malloc((nvtxs+1) * sizeof(ssize_t)); */
  nxadj = (int32_t *)malloc((nvtxs+1) * sizeof(int32_t));
  assert_malloc(nxadj);
  graph->rowPtr = (UINT_t *)nxadj;
  nadjncy = (int32_t *)malloc((nvtxs+xadj[nvtxs]) * sizeof(int32_t));
  assert_malloc(nadjncy);
  graph->colInd = (UINT_t *)nadjncy;
  
  perm  = (int32_t *)malloc(nvtxs * sizeof(int32_t));
  assert_malloc(perm);
  iperm = (int32_t *)malloc(nvtxs * sizeof(int32_t));
  assert_malloc(iperm);

  /* Determine maxdegree/csrange */
  #pragma omp parallel for schedule(static,4096) default(none) \
     shared(nvtxs, xadj) \
     reduction(max: maxdegree)
  for (vi=0; vi<nvtxs; vi++) 
    maxdegree = max(maxdegree, (int32_t)(xadj[vi+1]-xadj[vi]));

  csrange = maxdegree+1; 
  csrange = 16*((csrange+15)/16); /* get the per thread arrays to be alligned 
                                     at the start of the cache line */

  gcounts = (int32_t *)malloc((nthreads*csrange) * sizeof(int32_t));
  assert_malloc(gcounts);
#if 0
  psums   = (ssize_t *)malloc(nthreads * sizeof(ssize_t));
#else
  psums   = (int32_t *)malloc(nthreads * sizeof(int32_t));
#endif
  assert_malloc(psums);

  #pragma omp parallel default(none) \
    shared(nvtxs, nthreads, maxdegree, csrange, xadj, adjncy, nxadj, nadjncy, perm, iperm, gcounts, psums, chunkptr) 
  {
    int32_t vi, vistart, viend, vj, nedges, nchunks;
    int32_t ti, di, ci, dstart, dend;
    int32_t *counts, *buffer;
#if 0
    ssize_t ej, ejend, psum, chunksize;
#else
    int32_t ej, ejend, psum, chunksize;
#endif
    int mytid = omp_get_thread_num();

    vistart = mytid*((nvtxs+nthreads-1)/nthreads);
    viend   = min(nvtxs, (mytid+1)*((nvtxs+nthreads-1)/nthreads));

    dstart = mytid*((csrange+nthreads-1)/nthreads);
    dend   = min(csrange, (mytid+1)*((csrange+nthreads-1)/nthreads));

    /* Compute the local counts */
    counts = gcounts + mytid*csrange;
    for (int i = 0 ; i < csrange ; i++)
      counts[i] = 0;

    for (vi=vistart; vi<viend; vi++) 
      counts[xadj[vi+1]-xadj[vi]]++;
    #pragma omp barrier

    /* Compute the partial sum of the range assigned to each thread */
    for (psum=0, ti=0; ti<nthreads; ti++) {
      counts = gcounts + ti*csrange;
      for (di=dstart; di<dend; di++) 
        psum += counts[di];
    }
    psums[mytid] = psum;
    #pragma omp barrier

    #pragma omp single
    for (ti=1; ti<nthreads; ti++)
      psums[ti] += psums[ti-1];
    #pragma omp barrier

    /* Compute the actual prefix sums of the range assigned to each thread.
       This is done from right to left to get it into the desired exclusive
       prefix sum stage. */
    psum = psums[mytid];
    for (di=dend-1; di>=dstart; di--) { 
      counts = gcounts + (nthreads-1)*csrange;
      for (ti=nthreads-1; ti>=0; ti--) {
        psum -= counts[di];
        counts[di] = psum;
        counts -= csrange;
      }
    }
    #pragma omp barrier

    /* Create the perm/iperm arrays and the nxadj array of the re-ordered graph */
    counts = gcounts + mytid*csrange;

    /* TODO: This can be optimized by pre-sorting the per-thread vertices according 
             to their degree and processing them in increasing degree order */
    for (vi=vistart; vi<viend; vi++) {
      perm[vi] = counts[xadj[vi+1]-xadj[vi]]++;
      nxadj[perm[vi]] = xadj[vi+1]-xadj[vi]+1; /* the +1 is for the diagonal */
      iperm[perm[vi]] = vi;
    }
    #pragma omp barrier

    #pragma omp barrier
    /* compute the local sums and their prefix sums */
    for (psum=0, vi=vistart; vi<viend; vi++)
      psum += nxadj[vi];
    psums[mytid] = psum;
    #pragma omp barrier

    #pragma omp single
    for (ti=1; ti<nthreads; ti++)
      psums[ti] += psums[ti-1];
    #pragma omp barrier

    /* Compute the actual prefix sums of the nxadj[] array assigned to each thread.
       This is done from right to left to get it into the desired exclusive
       prefix sum stage. */
    psum = psums[mytid];
    if (mytid == nthreads-1)
      nxadj[nvtxs] = psum;
    for (vi=viend-1; vi>=vistart; vi--) { 
      psum -= nxadj[vi];
      nxadj[vi] = psum;
    }
    #pragma omp barrier
      
    /* Compute the chunk-based partitioning of the work for the reordered/sorted graph */
    chunksize = 1+psums[nthreads-1]/(100*nthreads);
    for (nchunks=0, psum=0, vi=vistart; vi<viend; vi++) {
      if ((psum += nxadj[vi+1]-nxadj[vi]) >= chunksize) {
        nchunks++;
        psum = 0;
      }
    }
    psums[mytid] = nchunks+1;

    #pragma omp barrier
    #pragma omp single
    for (ti=1; ti<nthreads; ti++)
      psums[ti] += psums[ti-1];
    #pragma omp barrier

    #pragma omp single
    chunkptr = (int32_t *)malloc((psums[nthreads-1]+1) * sizeof(int32_t));
    assert_malloc(chunkptr);
    #pragma omp barrier

    nchunks = psums[mytid];
    chunkptr[nchunks] = viend;
    for (psum=0, vi=viend-1; vi>=vistart; vi--) {
      if ((psum += nxadj[vi+1]-nxadj[vi]) >= chunksize) {
        chunkptr[--nchunks] = vi;
        psum = 0;
      }
    }
    if (mytid == 0)
      chunkptr[0] = 0;
    #pragma omp barrier

    nchunks = psums[nthreads-1]; /* this is the total # of chunks */
    /*
    #pragma omp single
    {
      for (vi=0; vi<nchunks; vi++) {
        printf("%4d: %6d - %6d [%5d: %zd]\n", vi, chunkptr[vi], chunkptr[vi+1], 
            chunkptr[vi+1]-chunkptr[vi], nxadj[chunkptr[vi+1]]-nxadj[chunkptr[vi]]);
      }
    }
    #pragma omp barrier
    */

    /* create the reordered/sorted graph by processing the chunks in parallel */
    #pragma omp for schedule(dynamic, 1) nowait
    for (ci=nchunks-1; ci>=0; ci--) {
      for (vi=chunkptr[ci]; vi<chunkptr[ci+1]; vi++) {
        vj = iperm[vi];
        buffer = nadjncy+nxadj[vi];
        for (nedges=0, ej=xadj[vj], ejend=xadj[vj+1]; ej<ejend; ej++, nedges++)  {
          buffer[nedges] = perm[adjncy[ej]];
	  
	}
        buffer[nedges++] = vi; /* put the diagonal */

        if (nedges > 1)
	  qsort(buffer, nedges, sizeof(int32_t), intCompare);
      }
    }

  }

  free(perm);
  free(iperm);
  free(gcounts);
  free(psums);
  free(chunkptr);

  return graph;
}


/*************************************************************************/
/*! The hash-map-based triangle-counting routine that uses the JIK
    triangle enumeration scheme.

    This version implements the following:
     - It does not store location information in L
     - Adds diagonal entries in U to make loops tighter
     - Reverts the order within U's adjancency lists to allow ++ traversal
*/
/*************************************************************************/
UINT_t tc_MapJIK_P(const GRAPH_TYPE *originalGraph)
{
  int32_t vi, vj, nvtxs, startv;
#if 0
  ssize_t ei, ej;
#else
  int32_t ei, ej;
#endif
  int64_t ntriangles=0;
#if 0
  ssize_t *xadj, *uxadj;
#else
  int32_t *xadj, *uxadj;
#endif
  int32_t *adjncy;
  int32_t l2, maxhmsize=0;
  GRAPH_TYPE *graph;

  graph = ptc_Preprocess(originalGraph);

  nvtxs  = (int32_t)graph->numVertices;
#if 0
  xadj   = (ssize_t *)graph->rowPtr;
#else
  xadj   = (int32_t *)graph->rowPtr;
#endif
  adjncy = (int32_t *)graph->colInd;

#if 0
  uxadj = (ssize_t *)malloc(nvtxs * sizeof(ssize_t)); /* the locations of the upper triangular part */
#else
  uxadj = (int32_t *)malloc(nvtxs * sizeof(int32_t)); /* the locations of the upper triangular part */
#endif
  assert_malloc(uxadj);

  /* populate uxadj[] and determine the size of the hash-map */
  startv = nvtxs;
  #pragma omp parallel for schedule(dynamic,1024) \
     default(none) \
     shared(nvtxs, xadj, adjncy, uxadj) \
     private(vj, ei, ej) \
     reduction(max: maxhmsize) \
     reduction(min: startv)
  for (vi=nvtxs-1; vi>=0; vi--) {
    for (ei=xadj[vi+1]-1; adjncy[ei]>vi; ei--);
    uxadj[vi] = ei;
    maxhmsize = max(maxhmsize, (int32_t)(xadj[vi+1]-uxadj[vi]));
    startv = (uxadj[vi] != xadj[vi] ? vi : startv);

    /* flip the order of Adj(vi)'s upper triangular adjacency list */
    for (ej=xadj[vi+1]-1; ei<ej; ei++, ej--) {
      vj = adjncy[ei];
      adjncy[ei] = adjncy[ej];
      adjncy[ej] = vj;
    }
  }

  /* convert the hash-map is converted into a format that is compatible with a 
     bitwise AND operation */
  for (l2=1; maxhmsize>(1<<l2); l2++);
  maxhmsize = (1<<(l2+4))-1;

#if 0
  printf("& compatible maxhmsize: %"PRId32", startv: %d\n", maxhmsize, startv);
#endif

#if 1
  if (nvtxs < maxhmsize) {
    fprintf(stderr,"ERROR: MapJIK does not work when nvtxs < maxhmsize\n");
    return 0;
  }
#endif

  #pragma omp parallel default(none) \
    shared(nvtxs, xadj, adjncy, uxadj, maxhmsize, startv) \
    reduction(+: ntriangles)
  {
    int32_t vi, vj, vk, vl, nlocal;
#if 0
    ssize_t ei, eiend, eistart, ej, ejend, ejstart;
#else
    int32_t ei, eiend, eistart, ej, ejend, ejstart;
#endif
    int32_t l, nc;
    int32_t l2=1, hmsize=(1<<(l2+4))-1, *hmap;
    int mytid = omp_get_thread_num();

    hmap = (int32_t *)calloc(maxhmsize+1, sizeof(int32_t));
    assert_malloc(hmap);

    /* Phase 1: Count triangles for vj < nvtxs-maxhmsize */
    #pragma omp for schedule(dynamic,GK_SBSIZE) nowait
    for (vj=startv; vj<nvtxs-maxhmsize; vj++) {
      if (xadj[vj+1]-uxadj[vj] == 1 || xadj[vj] == uxadj[vj])
        continue;
  
      /* adjust hmsize if needed */
      if (xadj[vj+1]-uxadj[vj] > (1<<l2)) {
        for (++l2; (xadj[vj+1]-uxadj[vj])>(1<<l2); l2++);
        hmsize = (1<<(l2+4))-1;
      }

      /* hash Adj(vj) */
      for (nc=0, ej=uxadj[vj], ejend=xadj[vj+1]-1; ej<ejend; ej++) {
        vk = adjncy[ej];
        for (l=(vk&hmsize); hmap[l]!=0; l=((l+1)&hmsize), nc++);
        hmap[l] = vk;
      }
  
      nlocal = 0;

      /* find intersections */
      if (nc > 0) { /* we had collisions */
        for (ej=xadj[vj], ejend=uxadj[vj]; ej<ejend; ej++) {
          vi = adjncy[ej];
          for (ei=uxadj[vi]; adjncy[ei]>vj; ei++) {
            vk = adjncy[ei];
            for (l=vk&hmsize; hmap[l]!=0 && hmap[l]!=vk; l=((l+1)&hmsize));
            if (hmap[l] == vk) 
              nlocal++;
          }
        }
  
        /* reset hash */
        for (ej=uxadj[vj], ejend=xadj[vj+1]-1; ej<ejend; ej++) {
          vk = adjncy[ej];
          for (l=(vk&hmsize); hmap[l]!=vk; l=((l+1)&hmsize));
          hmap[l] = 0;
        }
      }
      else { /* there were no collisons */
        for (ej=xadj[vj], ejend=uxadj[vj]; ej<ejend; ej++) {
          vi = adjncy[ej];
#ifdef TC_VECOPT 
          for (eiend=uxadj[vi]; adjncy[eiend]>vj; eiend++);
          for (ei=uxadj[vi]; ei<eiend; ei++) 
#else
          for (ei=uxadj[vi]; adjncy[ei]>vj; ei++) 
#endif
          {
            vk = adjncy[ei];
            nlocal += (hmap[vk&hmsize] == vk);
          }
        }
  
        /* reset hash */
        for (ej=uxadj[vj], ejend=xadj[vj+1]-1; ej<ejend; ej++) 
          hmap[adjncy[ej]&hmsize] = 0;
      }
      
      if (nlocal > 0)
        ntriangles += nlocal;
    }


    /* Phase 2: Count triangles for the last hmsize vertices, which can be done
                faster by using hmap as a direct map array. */
    hmap -= (nvtxs - maxhmsize);
    #pragma omp for schedule(dynamic,GK_DBSIZE) nowait
    for (vj=nvtxs-1; vj>=nvtxs-maxhmsize; vj--) {
      if (xadj[vj+1]-uxadj[vj] == 1 || xadj[vj] == uxadj[vj])
        continue;
  
      nlocal = 0;

      if (xadj[vj+1]-uxadj[vj] == nvtxs-vj) { /* complete row */
        /* find intersections */
        for (ej=xadj[vj], ejend=uxadj[vj]; ej<ejend; ej++) {
          vi = adjncy[ej];
          for (ei=uxadj[vi]; adjncy[ei]>vj; ei++);
          nlocal  += ei-uxadj[vi];
        }
      }
      else {
        /* hash Adj(vj) */
        for (ej=uxadj[vj], ejend=xadj[vj+1]-1; ej<ejend; ej++) 
          hmap[adjncy[ej]] = 1;
  
        /* find intersections */
        for (ej=xadj[vj], ejend=uxadj[vj]; ej<ejend; ej++) {
          vi = adjncy[ej];
#ifdef TC_VECOPT 
          for (eiend=uxadj[vi]; adjncy[eiend]>vj; eiend++);
          for (ei=uxadj[vi]; ei<eiend; ei++) 
#else
          for (ei=uxadj[vi]; adjncy[ei]>vj; ei++) 
#endif
            nlocal += hmap[adjncy[ei]];
        }
  
        /* reset hash */
        for (ej=uxadj[vj], ejend=xadj[vj+1]-1; ej<ejend; ej++) 
          hmap[adjncy[ej]] = 0;
      }

      if (nlocal > 0)
        ntriangles += nlocal;
    }
    hmap += (nvtxs - maxhmsize);

    free(hmap);
  }

  free_graph(graph);
  free(uxadj);

  return (UINT_t)ntriangles;
}

#endif

