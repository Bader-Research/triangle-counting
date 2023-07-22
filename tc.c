#include "types.h"
#include "queue.h"
#include "graph.h"
#include "tc.h"


static void bfs(const GRAPH_TYPE *, const UINT_t, UINT_t*);
static void bfs_mark_horizontal_edges(const GRAPH_TYPE *, const UINT_t, UINT_t*, Queue*, bool*, bool*);


/* Algorithm from
   T. A. Davis, "Graph algorithms via SuiteSparse: GraphBLAS: triangle counting and K-truss," 2018 IEEE High Performance extreme Computing Conference (HPEC), Waltham, MA, USA, 2018, pp. 1-6, doi: 10.1109/HPEC.2018.8547538.
*/
UINT_t tc_davis // # of triangles
(
#if 1
const GRAPH_TYPE *graph
#else
const UINT_t *restrict Ap, // column pointers, size n+1
const UINT_t *restrict Ai, // row indices
const UINT_t n // A is n-by-n
#endif
 )
{
#if 1
  UINT_t *Ap = graph->rowPtr;
  UINT_t *Ai = graph->colInd;
  UINT_t n = graph->numVertices;
#endif
  bool *restrict Mark = (bool *) calloc (n, sizeof (bool)) ;
  if (Mark == NULL) return (-1) ;
  
  UINT_t ntri = 0 ;
  for (UINT_t j = 0 ; j < n ; j++) {
    // scatter A(:,j) into Mark
    for (UINT_t p = Ap [j] ; p < Ap [j+1] ; p++)
      Mark [Ai [p]] = 1 ;
    // sum(C(:,j)) where C(:,j) = (A * A(:,j)) .* Mark
    for (UINT_t p = Ap [j] ; p < Ap [j+1] ; p++) {
      const UINT_t k = Ai [p] ;
      // C(:,j) += (A(:,k) * A(k,j)) .* Mark
      for (UINT_t pa = Ap [k] ; pa < Ap [k+1] ; pa++)
	// C(i,j) += (A(i,k) * A(k,j)) .* Mark
	ntri += Mark [Ai [pa]] ;
    }
    for (UINT_t p = Ap [j] ; p < Ap [j+1] ; p++)
      Mark [Ai [p]] = 0 ;
  }
  free (Mark) ;
  return (ntri/6) ;
}
 

UINT_t tc_wedge(const GRAPH_TYPE *graph) {
  /* Algorithm: For each vertex i, for each open wedge (j, i, k), determine if there's a closing edge (j, k) */
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

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
	      count++;
	      break;
	    }
	  }
	}
      }
    }
  }

  return (count/6);
}

UINT_t tc_wedge_DO(const GRAPH_TYPE *graph) {
  /* Algorithm: For each vertex i, for each open wedge (j, i, k), determine if there's a closing edge (j, k) */
  /* Direction oriented. */
  
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

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
		count++;
		break;
	      }
	    }
	  }
	}
      }
    }
  }

  return count;
}


UINT_t tc_triples(const GRAPH_TYPE *graph) {
  /* Algorithm: for each triple (i, j, k), determine if the three triangle edges exist. */
  
  register UINT_t i, j, k;
  UINT_t count = 0;

  const UINT_t n = graph->numVertices;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      for (k = 0; k < n; k++) {

	if (check_edge(graph, i, j) && check_edge(graph, j, k) && check_edge(graph, k, i))
	  count++;

      }

  return (count/6);
}

UINT_t tc_triples_DO(const GRAPH_TYPE *graph) {
  /* Algorithm: for each triple (i, j, k), determine if the three triangle edges exist. */
  /* Direction oriented. */
  
  register UINT_t i, j, k;
  UINT_t count = 0;

  const UINT_t n = graph->numVertices;

  for (i = 0; i < n; i++)
    for (j = i; j < n; j++)
      for (k = j; k < n; k++) {

	if (check_edge(graph, i, j) && check_edge(graph, j, k) && check_edge(graph, k, i))
	  count++;

      }

  return count;
}





UINT_t tc_intersectMergePath(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a linear scan. */
  
  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  for (v = 0; v < n; v++) {
    b = Ap[v  ];
    e = Ap[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w = Ai[i];
      count += intersectSizeMergePath(graph, v, w);
    }
  }

  return (count/6);
}

UINT_t tc_intersectMergePath_DO(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a linear scan. */
  /* Direction oriented. */
  
  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  for (v = 0; v < n; v++) {
    b = Ap[v  ];
    e = Ap[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w = Ai[i];
      if (v < w)
	count += intersectSizeMergePath(graph, v, w);
    }
  }

  return (count/3);
}


UINT_t tc_intersectBinarySearch(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a binary search. */

  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  for (v = 0; v < n; v++) {
    b = Ap[v  ];
    e = Ap[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w  = Ai[i];
      count += intersectSizeBinarySearch(graph, v, w);
    }
  }

  return (count/6);
}

UINT_t tc_intersectBinarySearch_DO(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a binary search. */
  /* Direction oriented. */

  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  for (v = 0; v < n; v++) {
    b = Ap[v  ];
    e = Ap[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w  = Ai[i];
      if (v < w)
	count += intersectSizeBinarySearch(graph, v, w);
    }
  }

  return (count/3);
}

static UINT_t binarySearch_partition(const UINT_t* list, const UINT_t start, const UINT_t end, const UINT_t target) {
  register INT_t s=start, e=end, mid;
  while (s < e) {
    mid = s + (e - s) / 2;
#if 0
    if (mid >=e) {
      printf("ERROR, binarypart: mid >= e\n");
      exit(-1);
    }
#endif
    if (list[mid] == target)
      return mid;

    if (list[mid] < target)
      s = mid + 1;
    else
      e = mid;
  }
  return s;
}



static UINT_t searchLists_with_partitioning(const UINT_t* list1, const INT_t s1, const INT_t e1, const UINT_t* list2, const INT_t s2, const INT_t e2) {
  INT_t mid1, loc2;
  UINT_t count = 0;

  if ((s1>e1)||(s2>e2))
    return 0;
  
  mid1 = s1 + (e1 - s1)/2;

#if 0
  if (mid1 > e1) {
    printf("ERROR: searchpart: mid1 > e1\n");
    exit(-1);
  }
#endif

  /* need a binary search that returns the item or the next higher position */
  loc2 = binarySearch_partition(list2, s2, e2, list1[mid1]); 

#if 0
  if (loc2 > e2) {
    printf("ERROR: searchpart loc2 (%d) > e2 (%d)\n",loc2,e2);
    exit(-1);
  }
#endif

  INT_t s11 = s1;
  INT_t e11 = mid1 - 1;
  INT_t s21 = s2;
  INT_t e21 = loc2;

  INT_t s12 = mid1 + 1;
  INT_t e12 = e1;
  INT_t s22 = loc2;
  INT_t e22 = e2;

  if (list1[mid1] == list2[loc2]) {
    count++;
    e21--;
    s22++;
  }
  count += searchLists_with_partitioning(list1, s11, e11, list2, s21, e21);
  count += searchLists_with_partitioning(list1, s12, e12, list2, s22, e22);
  return count;
}


UINT_t tc_intersectPartition(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a binary search-based partition. */
  
  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t count = 0;

  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  for (v = 0; v < n ; v++) {
    b = Ap[v  ];
    e = Ap[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w  = Ai[i];
      count += searchLists_with_partitioning((UINT_t *)Ai, (INT_t) Ap[v], (INT_t)Ap[v+1]-1, (UINT_t *)Ai, (INT_t)Ap[w], (INT_t)Ap[w+1]-1);
    }
  }

  return (count/6);
}

UINT_t tc_intersectPartition_DO(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a binary search-based partition. */
  /* Direction oriented. */
  
  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t count = 0;

  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  for (v = 0; v < n; v++) {
    b = Ap[v  ];
    e = Ap[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w  = Ai[i];
      if (v < w)
	count += searchLists_with_partitioning((UINT_t *)Ai, (INT_t) Ap[v], (INT_t)Ap[v+1]-1, (UINT_t *)Ai, (INT_t)Ap[w], (INT_t)Ap[w+1]-1);
    }
  }

  return (count/3);
}


static UINT_t intersectSizeHash(const GRAPH_TYPE *graph, bool *Hash, const UINT_t v, const UINT_t w) {

  register UINT_t vb, ve, wb, we;
  register UINT_t s1, e1, s2, e2;
  UINT_t count = 0;
  
  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;

  vb = Ap[v  ];
  ve = Ap[v+1];
  wb = Ap[w  ];
  we = Ap[w+1];

  if ((ve-vb) < (we-wb)) {
    s1 = vb;
    e1 = ve;
    s2 = wb;
    e2 = we;
  } else {
    s1 = wb;
    e1 = we;
    s2 = vb;
    e2 = ve;
  }
  
  for (UINT_t i=s1 ; i<e1 ; i++)
    Hash[Ai[i]] = true;

  for (UINT_t i= s2; i<e2 ; i++)
    if (Hash[Ai[i]]) count++;

  for (UINT_t i=s1 ; i<e1 ; i++)
    Hash[Ai[i]] = false;

  return count;
}

UINT_t tc_intersectHash(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a hash. */

  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t count = 0;

  bool *Hash;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  Hash = (bool *)calloc(n, sizeof(bool));
  assert_malloc(Hash);
  
  for (v = 0; v < n ; v++) {
    b = Ap[v  ];
    e = Ap[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w  = Ai[i];
      count += intersectSizeHash(graph, Hash, v, w);
    }
  }

  free(Hash);

  return (count/6);
}


UINT_t tc_intersectHash_DO(const GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a hash. */
  /* Direction oriented. */

  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t count = 0;

  bool *Hash;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

  Hash = (bool *)calloc(m, sizeof(bool));
  assert_malloc(Hash);
  
  for (v = 0; v < n ; v++) {
    b = Ap[v  ];
    e = Ap[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w  = Ai[i];
      if (v < w)
	count += intersectSizeHash(graph, Hash, v, w);
    }
  }

  free(Hash);

  return (count/3);
}

/* T. M. Low, V. N. Rao, M. Lee, D. Popovici, F. Franchetti and S. McMillan,
   "First look: Linear algebra-based triangle counting without matrix multiplication,"
   2017 IEEE High Performance Extreme Computing Conference (HPEC),
   Waltham, MA, USA, 2017, pp. 1-6,
   doi: 10.1109/HPEC.2017.8091046. */
UINT_t tc_low(
#if 1
			const GRAPH_TYPE *graph
#else
			UINT_t *IA , // row indices
			UINT_t *JA  // column indices
#endif
			)
{
#if 1
  UINT_t* IA = graph->rowPtr;
  UINT_t* JA = graph->colInd;
  UINT_t N = graph->numVertices;
#endif
  UINT_t delta = 0 ; // number of triangles

  // for every vertex in V_{BR}
  // # pragma omp parallel for schedule (dynamic) reduction (+ : delta)
  for ( UINT_t i = 1 ; i<N-1; i ++ )
    {
      UINT_t *curr_row_x = IA + i ;
      UINT_t *curr_row_A = IA + i + 1;
      UINT_t num_nnz_curr_row_x = *curr_row_A - *curr_row_x;
      UINT_t *x_col_begin = ( JA + *curr_row_x);
      UINT_t *x_col_end = x_col_begin;
      UINT_t *row_bound = x_col_begin + num_nnz_curr_row_x;
      // UINT_t col_x_min = 0;
      UINT_t col_x_max = i - 1;

      // partition the current row into x and y, where x == a01 ^ T == a10t and y == a12t
      while (x_col_end < row_bound && *x_col_end < col_x_max)
	++x_col_end;
      x_col_end -= (*x_col_end > col_x_max || x_col_end == row_bound);

      UINT_t *y_col_begin = x_col_end + 1;
      UINT_t *y_col_end = row_bound - 1;
      UINT_t num_nnz_y = (y_col_end - y_col_begin) + 1;
      UINT_t num_nnz_x = (x_col_end - x_col_begin) + 1;

      UINT_t y_col_first = i + 1;
      UINT_t x_col_first = 0;
      UINT_t *y_col = y_col_begin;

      // compute y*A20*x ( Equation 5 )
      for (UINT_t j = 0 ; j< num_nnz_y ; ++j ,++ y_col)
	{
	  UINT_t row_index_A = *y_col - y_col_first;
	  UINT_t *x_col = x_col_begin;
	  UINT_t num_nnz_A = *( curr_row_A + row_index_A + 1 ) - *( curr_row_A + row_index_A );
	  UINT_t *A_col = ( JA + *( curr_row_A + row_index_A ) );
	  UINT_t *A_col_max = A_col + num_nnz_A ;

	  for (UINT_t k = 0 ; k < num_nnz_x && *A_col <= col_x_max ; ++k )
	    {
	      UINT_t row_index_x = *x_col - x_col_first;
	      while ( (* A_col < *x_col ) && ( A_col < A_col_max ) )
		++A_col;
	      delta += (* A_col == row_index_x ) ;
	      
	      ++x_col ;
	    }
	}
    }
  return delta ;
 }


void bfs_treelist(const GRAPH_TYPE *graph, const bool* E, UINT_t* parent /* , UINT_t* component */) {

  /* UINT_t c; */
  
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  bool *visited = (bool *)malloc(n * sizeof(bool));
  assert_malloc(visited);
  for (UINT_t i=0 ; i<n ; i++)
    visited[i] = false;

  // c = 0;
  for (UINT_t v=0 ; v<n ; v++) {
    if (!visited[v]) {
      // c++;
  
      Queue *queue = createQueue(n);

      visited[v] = true;
      // component[v] = c;
      enqueue(queue, v);
  
      while (!isEmpty(queue)) {
	UINT_t currentVertex = dequeue(queue);
	UINT_t s = Ap[v];
	UINT_t e = Ap[v+1];
	for (UINT_t i = s; i < e ; i++) {
	  if (E[i]) {
	    UINT_t adjacentVertex = Ai[i];
	    if (!visited[adjacentVertex])  {
	      visited[adjacentVertex] = true;
	      // component[adjacentVertex] = c;
	      parent[adjacentVertex] = currentVertex;
	      enqueue(queue, adjacentVertex);
	    }
	  }
	}
      }
      free_queue(queue);
    }
  }

  free(visited);
}

static bool and_E(const bool* E, const UINT_t m) {
    for (UINT_t i = 0 ; i<m ; i++)
      if (E[i]) return true;
    return false;
}

static bool check_edge_treelist(const GRAPH_TYPE *graph, const bool* E, const UINT_t v, const UINT_t w) {

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  if ((v==n) || (w==n)) return false;
  UINT_t s = Ap[v];
  UINT_t e = Ap[v+1];
  for (UINT_t i = s; i < e; i++)
    if (E[i] && (Ai[i] == w))
      return true;

  return false;
}


static void remove_treelist(const GRAPH_TYPE* graph, bool *E, const UINT_t *parent) {

  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  for (UINT_t v=0 ; v<n ; v++) {
    UINT_t s = Ap[v];
    UINT_t e = Ap[v+1];
    for (UINT_t i = s; i < e; i++) {
      if (E[i]) {
	UINT_t w = Ai[i];
	if (parent[w] == v) {
	  E[i] = false;
	  UINT_t ws = Ap[w];
	  UINT_t we = Ap[w+1];
	  for (UINT_t j=ws ; j<we ; j++) {
	    if (Ai[j] == v) {
	      E[j] = false;
	      break;
	    }
	  }
	}
      }
    }
  }
  
  return;
}

  
UINT_t tc_treelist(const GRAPH_TYPE *graph) {
  /* Itai and Rodeh, SIAM Journal of Computing, 1978 */
  UINT_t count = 0;
  UINT_t *parent;
  // UINT_t* component;
  
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

  // 1. While there remains an edge in E:
  //   1a. compute a covering tree for each connected component of G;
  //   1b. for each edge (u, v) in none of these trees:
  //     1ba. If (father(u), v) in E then output triangle (u, v, father(u))
  //     1bb. else If (father(v), u) in E then output triangle (u, v, father (v)
  //   1c. remove from E all the edges in these trees


  bool* E = (bool *)malloc(m * sizeof(bool));
  assert_malloc(E);

  for (UINT_t i=0 ; i<m ; i++)
    E[i] = true;

  //  component = (UINT_t *)malloc(n * sizeof(UINT_t));
  //  assert_malloc(component);
    
  parent = (UINT_t *)malloc(n * sizeof(UINT_t));
  assert_malloc(parent);
    
  while (and_E(E, m)) {

    for (UINT_t i=0 ; i<n ; i++) {
      // component[i] = 0;
      parent[i] = n;
    }

    bfs_treelist(graph, E, parent /* , component */);

    for (UINT_t u=0 ; u<n ; u++) {
      UINT_t s = Ap[u];
      UINT_t e = Ap[u+1];
      for (UINT_t j=s; j<e ; j++) {
	if (E[j]) {
	  UINT_t v = Ai[j];
	  if (parent[u] != v) {
	    if (check_edge_treelist(graph, E, parent[u], v))
	      count++;
	    else if (check_edge_treelist(graph, E, parent[v], u))
	      count++;
	  }
	}
      }
    }

    remove_treelist(graph, E, parent);
  }
  
  free(parent);
  // free(component);
  free(E);

  return count/2;
}


static void bfs_treelist2(const GRAPH_TYPE *graph, UINT_t* parent /* , UINT_t* component */) {

  /* UINT_t c; */
  
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  bool *visited = (bool *)malloc(n * sizeof(bool));
  assert_malloc(visited);
  for (UINT_t i=0 ; i<n ; i++)
    visited[i] = false;

  // c = 0;
  for (UINT_t v=0 ; v<n ; v++) {
    if (!visited[v]) {
      // c++;
  
      Queue *queue = createQueue(n);

      visited[v] = true;
      // component[v] = c;
      enqueue(queue, v);
  
      while (!isEmpty(queue)) {
	UINT_t currentVertex = dequeue(queue);
	UINT_t s = Ap[v];
	UINT_t e = Ap[v+1];
	for (UINT_t i = s; i < e ; i++) {
	  UINT_t adjacentVertex = Ai[i];
	  if (!visited[adjacentVertex])  {
	    visited[adjacentVertex] = true;
	    // component[adjacentVertex] = c;
	    parent[adjacentVertex] = currentVertex;
	    enqueue(queue, adjacentVertex);
	  }
	}
      }
      free_queue(queue);
    }
  }

  free(visited);
}


static void remove_treelist2(GRAPH_TYPE* graph, const UINT_t *parent) {

  UINT_t *restrict Ap = graph->rowPtr;
  UINT_t *restrict Ai = graph->colInd;
  UINT_t n = graph->numVertices;
  UINT_t m = graph->numEdges;

  bool *E = (bool *)malloc(m * sizeof(bool));
  assert_malloc(E);

  UINT_t *Degree = (UINT_t *)malloc(n * sizeof(UINT_t));
  assert_malloc(Degree);

  for (UINT_t i=0 ; i<m ; i++)
    E[i] = true;

  for (UINT_t v=0 ; v<n ; v++) {
    UINT_t s = Ap[v];
    UINT_t e = Ap[v+1];
    for (UINT_t i = s; i < e; i++) {
      UINT_t w = Ai[i];
      if (parent[w] == v) {
	E[i] = false;
	UINT_t ws = Ap[w];
	UINT_t we = Ap[w+1];
	for (UINT_t j=ws ; j<we ; j++) {
	  if (Ai[j] == v) {
	    E[j] = false;
	    break;
	  }
	}
      }
    }
  }

  UINT_t numEdges_new=0;
  for (UINT_t i=0 ; i<m ; i++) {
    if (E[i]) {
      Ai[numEdges_new] = Ai[i];
      numEdges_new++;
    }
  }

  for (UINT_t v=0 ; v<n ; v++) {
    Degree[v] = 0;
    UINT_t s = Ap[v];
    UINT_t e = Ap[v+1];
    for (UINT_t i = s; i < e; i++) {
      if (E[i]) {
	Degree[v]++;
      }
    }
  }

  Ap[0] = 0;
  for (UINT_t i=1 ; i<=n ; i++) {
    Ap[i] = Ap[i-1] + Degree[i-1];
  }

  graph->numEdges = numEdges_new;
    
  free(Degree);
  free(E);

  return;
}

UINT_t tc_treelist2(const GRAPH_TYPE *graph) {
  /* Itai and Rodeh, SIAM Journal of Computing, 1978 */
  UINT_t edges;
  UINT_t count = 0;
  UINT_t *parent;
  // UINT_t* component;

  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

  // 1. While there remains an edge in E:
  //   1a. compute a covering tree for each connected component of G;
  //   1b. for each edge (u, v) in none of these trees:
  //     1ba. If (father(u), v) in E then output triangle (u, v, father(u))
  //     1bb. else If (father(v), u) in E then output triangle (u, v, father (v)
  //   1c. remove from E all the edges in these trees


  GRAPH_TYPE* graph2;

  graph2 = (GRAPH_TYPE *)malloc(sizeof(GRAPH_TYPE));
  assert_malloc(graph2);

  graph2->numVertices = n;
  graph2->numEdges = m;
  allocate_graph(graph2);
  copy_graph(graph, graph2);
  
  const UINT_t *restrict Ap2 = graph2->rowPtr;
  const UINT_t *restrict Ai2 = graph2->colInd;

  //  component = (UINT_t *)malloc(n * sizeof(UINT_t));
  //  assert_malloc(component);
    
  parent = (UINT_t *)malloc(n * sizeof(UINT_t));
  assert_malloc(parent);

  edges = m;
    
  while (edges>0) {

    for (UINT_t i=0 ; i<n ; i++) {
      // component[i] = 0;
      parent[i] = n;
    }

    bfs_treelist2(graph2, parent /* , component */);

    for (UINT_t u=0 ; u<n ; u++) {
      UINT_t s = Ap2[u];
      UINT_t e = Ap2[u+1];
      for (UINT_t j=s; j<e ; j++) {
	UINT_t v = Ai2[j];
	if (parent[u] != v) {
	  if ((parent[u] < n) && (check_edge(graph2, parent[u], v)))
	    count++;
	  else if ((parent[v] < n) && (check_edge(graph2, parent[v], u)))
	    count++;
	}
      }
    }

    remove_treelist2(graph2, parent);
    edges = graph2->numEdges;
  }
  
  free(parent);
  // free(component);

  free_graph(graph2);

  return count/2;
}


static UINT_t intersectSizeMergePath_forward(const GRAPH_TYPE* graph, const UINT_t v, const UINT_t w, const UINT_t* A, const UINT_t* Size) {
  register UINT_t vb, ve, wb, we;
  register UINT_t ptr_v, ptr_w;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  
  vb = Ap[v  ];
  ve = vb + Size[v];
  wb = Ap[w  ];
  we = wb + Size[w];

  ptr_v = vb;
  ptr_w = wb;
  while ((ptr_v < ve) && (ptr_w < we)) {
    if (A[ptr_v] == A[ptr_w]) {
      count++;
      ptr_v++;
      ptr_w++;
    }
    else
      if (A[ptr_v] < A[ptr_w])
	ptr_v++;
      else
	ptr_w++;
  }
  return count;
}

static UINT_t intersectSizeHash_forward(const GRAPH_TYPE *graph, bool *Hash, const UINT_t v, const UINT_t w, const UINT_t* A, const UINT_t* Size) {

  register UINT_t vb, ve, wb, we;
  register UINT_t s1, e1, s2, e2;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;

  vb = Ap[v  ];
  ve = vb + Size[v];
  wb = Ap[w  ];
  we = wb + Size[w];

  if (Size[v] < Size[w]) {
    s1 = vb;
    e1 = ve;
    s2 = wb;
    e2 = we;
  } else {
    s1 = wb;
    e1 = we;
    s2 = vb;
    e2 = ve;
  }

  for (UINT_t i=s1 ; i<e1 ; i++)
    Hash[A[i]] = true;

  for (UINT_t i= s2; i<e2 ; i++)
    if (Hash[A[i]]) count++;
  
  for (UINT_t i=s1 ; i<e1 ; i++)
    Hash[A[i]] = false;
    
  return count;
}


static UINT_t intersectSizeHashSkip_forward(const GRAPH_TYPE *graph, bool *Hash, const UINT_t v, const UINT_t w, const UINT_t* A, const UINT_t* Size) {

  register UINT_t s1, e1, s2, e2;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;

  if (Size[v] < Size[w]) {
    if (Size[v] == 0) return 0;
    s1 = Ap[v  ];
    e1 = s1 + Size[v];
    s2 = Ap[w  ];
    e2 = s2 + Size[w];
  } else {
    if (Size[w] == 0) return 0;
    s1 = Ap[w  ];
    e1 = s1 + Size[w];
    s2 = Ap[v  ];
    e2 = s2 + Size[v];
  }

  for (UINT_t i=s1 ; i<e1 ; i++)
    Hash[A[i]] = true;

  for (UINT_t i=s2 ; i<e2 ; i++)
    if (Hash[A[i]]) count++;

  for (UINT_t i=s1 ; i<e1 ; i++)
    Hash[A[i]] = false;

  return count;
}


UINT_t tc_forward(const GRAPH_TYPE *graph) {
  
/* Schank, T., Wagner, D. (2005). Finding, Counting and Listing All Triangles in Large Graphs, an Experimental Study. In: Nikoletseas, S.E. (eds) Experimental and Efficient Algorithms. WEA 2005. Lecture Notes in Computer Science, vol 3503. Springer, Berlin, Heidelberg. https://doi.org/10.1007/11427186_54 */

  register UINT_t s, t;
  register UINT_t b, e;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

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
	count += intersectSizeMergePath_forward(graph, s, t, A, Size);
	A[Ap[t] + Size[t]] = s;
	Size[t]++;
      }
    }
  }

  free(A);
  free(Size);
  
  return count;
}

static UINT_t tc_forward_hash_config_size(const GRAPH_TYPE *graph, UINT_t hashSize) {
  
/* Schank, T., Wagner, D. (2005). Finding, Counting and Listing All Triangles in Large Graphs, an Experimental Study. In: Nikoletseas, S.E. (eds) Experimental and Efficient Algorithms. WEA 2005. Lecture Notes in Computer Science, vol 3503. Springer, Berlin, Heidelberg. https://doi.org/10.1007/11427186_54 */

  register UINT_t s, t;
  register UINT_t b, e;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

  bool* Hash = (bool *)calloc( (hashSize == 0)? m: hashSize, sizeof(bool));
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
	count += intersectSizeHash_forward(graph, Hash, s, t, A, Size);
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

UINT_t tc_forward_hash(const GRAPH_TYPE *graph) {
  return tc_forward_hash_config_size(graph, 0);
}

static UINT_t tc_forward_hash_skip_config_size(const GRAPH_TYPE *graph, UINT_t hashSize) {

/* Schank, T., Wagner, D. (2005). Finding, Counting and Listing All Triangles in Large Graphs, an Experimental Study. In: Nikoletseas, S.E. (eds) Experimental and Efficient Algorithms. WEA 2005. Lecture Notes in Computer Science, vol 3503. Springer, Berlin, Heidelberg. https://doi.org/10.1007/11427186_54 */

  register UINT_t s, t;
  register UINT_t b, e;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

  bool* Hash = (bool *)calloc( (hashSize == 0)? m: hashSize, sizeof(bool));
  assert_malloc(Hash);

  UINT_t* Size = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(Size);

  UINT_t* A = (UINT_t *)calloc(m, sizeof(UINT_t));
  assert_malloc(A);

  for (s = 0 ; s < n ; s++) {
    b = Ap[s  ];
    e = Ap[s+1];
    for (UINT_t i=b ; i<e ; i++) {
      t  = Ai[i];
      if (s<t) {
	count += intersectSizeHashSkip_forward(graph, Hash, s, t, A, Size);
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

UINT_t tc_forward_hash_skip(const GRAPH_TYPE *graph) {
  return tc_forward_hash_skip_config_size(graph, 0);
}

UINT_t tc_forward_hash_degreeOrder(const GRAPH_TYPE *graph) {
  
/* Schank, T., Wagner, D. (2005). Finding, Counting and Listing All Triangles in Large Graphs, an Experimental Study. In: Nikoletseas, S.E. (eds) Experimental and Efficient Algorithms. WEA 2005. Lecture Notes in Computer Science, vol 3503. Springer, Berlin, Heidelberg. https://doi.org/10.1007/11427186_54 */

  UINT_t count = 0;

  GRAPH_TYPE *graph2;
  graph2 = reorder_graph_by_degree(graph, REORDER_HIGHEST_DEGREE_FIRST);

  count = tc_forward_hash(graph2);

  free_graph(graph2);
  
  return count;
}


UINT_t tc_forward_hash_degreeOrderReverse(const GRAPH_TYPE *graph) {

/* Schank, T., Wagner, D. (2005). Finding, Counting and Listing All Triangles in Large Graphs, an Experimental Study. In: Nikoletseas, S.E. (eds) Experimental and Efficient Algorithms. WEA 2005. Lecture Notes in Computer Science, vol 3503. Springer, Berlin, Heidelberg. https://doi.org/10.1007/11427186_54 */

  UINT_t count = 0;

  GRAPH_TYPE *graph2;
  graph2 = reorder_graph_by_degree(graph, REORDER_LOWEST_DEGREE_FIRST);

  count = tc_forward_hash(graph2);

  free_graph(graph2);

  return count;
}


// Function to perform breadth-first search
static void bfs(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level) {
  bool *visited = (bool *)calloc(graph->numVertices, sizeof(bool));
  assert_malloc(visited);

  Queue *queue = createQueue(graph->numVertices);

  visited[startVertex] = true;
  enqueue(queue, startVertex);
  level[startVertex] = 0;
  
  while (!isEmpty(queue)) {
    UINT_t v = dequeue(queue);
    for (UINT_t i = graph->rowPtr[v]; i < graph->rowPtr[v + 1]; i++) {
      UINT_t w = graph->colInd[i];
      if (!visited[w])  {
	visited[w] = true;
	enqueue(queue, w);
	level[w] = level[v] + 1;
      }
    }
  }

  free(visited);
  free_queue(queue);
}

static void bfs_bader3(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* restrict level, Queue* queue, bool* visited) {
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;

  visited[startVertex] = true;
  enqueue(queue, startVertex);
  level[startVertex] = 1;
  
  while (!isEmpty(queue)) {
    UINT_t v = dequeue(queue);
    for (UINT_t i = Ap[v]; i < Ap[v + 1]; i++) {
      UINT_t w = Ai[i];
      if (!visited[w])  {
	visited[w] = true;
	enqueue(queue, w);
	level[w] = level[v] + 1;
      }
    }
  }
}

static void bfs_mark_horizontal_edges(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* restrict level, Queue* queue, bool* visited, bool* horiz) {
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;

  visited[startVertex] = true;
  enqueue(queue, startVertex);
  level[startVertex] = 1;
  
  while (!isEmpty(queue)) {
    UINT_t v = dequeue(queue);
    for (UINT_t i = Ap[v]; i < Ap[v + 1]; i++) {
      UINT_t w = Ai[i];
      if (!visited[w])  {
	horiz[i] = false;
	visited[w] = true;
	enqueue(queue, w);
	level[w] = level[v] + 1;
      }
      else {
	horiz[i] = (level[w] == 0) || (level[w] == level[v]);
      }
    }
  }
}


static void bader_intersectSizeMergePath(const GRAPH_TYPE* graph, const UINT_t* level, const UINT_t v, const UINT_t w, UINT_t* restrict c1, UINT_t* restrict c2) {
  register UINT_t vb, ve, wb, we;
  register UINT_t ptr_v, ptr_w;
  UINT_t level_v;
  
  level_v = level[v];

  vb = graph->rowPtr[v ];
  ve = graph->rowPtr[v+1];
  wb = graph->rowPtr[w  ];
  we = graph->rowPtr[w+1];


  ptr_v = vb;
  ptr_w = wb;
  while ((ptr_v < ve) && (ptr_w < we)) {
    if (graph->colInd[ptr_v] == graph->colInd[ptr_w]) {
      if (level_v == level[graph->colInd[ptr_v]]) (*c2)++;
      else (*c1)++;
      ptr_v++;
      ptr_w++;
    }
    else
      if (graph->colInd[ptr_v] < graph->colInd[ptr_w])
	ptr_v++;
      else
	ptr_w++;
  }

  return;
}


double tc_bader_compute_k(const GRAPH_TYPE *graph) {
  /* Direction orientied. */
  UINT_t* restrict level;
  UINT_t s, e, l, w;
  UINT_t c1, c2;
  UINT_t NO_LEVEL;
  UINT_t k;

  level = (UINT_t *)malloc(graph->numVertices * sizeof(UINT_t));
  assert_malloc(level);
  NO_LEVEL = graph->numVertices;
  for (UINT_t i = 0 ; i < graph->numVertices ; i++) 
    level[i] = NO_LEVEL;
  
  for (UINT_t i = 0 ; i < graph->numVertices ; i++) {
    if (level[i] == NO_LEVEL) {
      bfs(graph, i, level);
    }
  }

  k = 0;
  
  c1 = 0; c2 = 0;
  for (UINT_t v = 0 ; v < graph->numVertices ; v++) {
    s = graph->rowPtr[v  ];
    e = graph->rowPtr[v+1];
    l = level[v];
    for (UINT_t j = s ; j<e ; j++) {
      w = graph->colInd[j];
      if ((v < w) && (level[w] == l)) {
	k++;
	bader_intersectSizeMergePath(graph, level, v, w, &c1, &c2);
      }
    }
  }

  free(level);

  return (2.0 * (double)k/(double)graph->numEdges);
}

UINT_t tc_bader(const GRAPH_TYPE *graph) {
  /* Direction orientied. */
  UINT_t* restrict level;
  UINT_t s, e, l, w;
  UINT_t c1, c2;
  UINT_t NO_LEVEL;

  level = (UINT_t *)malloc(graph->numVertices * sizeof(UINT_t));
  assert_malloc(level);
  NO_LEVEL = graph->numVertices;
  for (UINT_t i = 0 ; i < graph->numVertices ; i++) 
    level[i] = NO_LEVEL;
  
  for (UINT_t i = 0 ; i < graph->numVertices ; i++) {
    if (level[i] == NO_LEVEL) {
      bfs(graph, i, level);
    }
  }

  c1 = 0; c2 = 0;
  for (UINT_t v = 0 ; v < graph->numVertices ; v++) {
    s = graph->rowPtr[v  ];
    e = graph->rowPtr[v+1];
    l = level[v];
    for (UINT_t j = s ; j<e ; j++) {
      w = graph->colInd[j];
      if ((v < w) && (level[w] == l))
	bader_intersectSizeMergePath(graph, level, v, w, &c1, &c2);
    }
  }


  free(level);

  return c1 + (c2/3);
}


UINT_t tc_bader3(const GRAPH_TYPE *graph) {
  /* Bader's new algorithm for triangle counting based on BFS */
  /* Uses Hash array to detect triangles (v, w, x) if x is adjacent to v */
  /* For level[], 0 == unvisited. Needs a modified BFS starting from level 1 */
  /* Direction orientied. */
  UINT_t* restrict level;
  UINT_t c1, c2;
  register UINT_t x;
  bool* Hash;
  bool* visited;
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  level = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(level);
  
  visited = (bool *)calloc(n, sizeof(bool));
  assert_malloc(visited);

  Hash = (bool *)calloc(n, sizeof(bool));
  assert_malloc(Hash);

  Queue *queue = createQueue(n);

  c1 = 0; c2 = 0;
  for (UINT_t v = 0 ; v < n ; v++) {
    if (!level[v])
      bfs_bader3(graph, v, level, queue, visited);
    const UINT_t s = Ap[v  ];
    const UINT_t e = Ap[v+1];
    const UINT_t l = level[v];

    for (UINT_t p = s ; p<e ; p++)
      Hash[Ai[p]] = true;
    
    for (UINT_t j = s ; j<e ; j++) {
      const UINT_t w = Ai[j];
      if ((v < w) && (level[w] == l)) {
	/* bader_intersectSizeMergePath(graph, level, v, w, &c1, &c2); */
	for (UINT_t k = Ap[w]; k < Ap[w+1] ; k++) {
	  x = Ai[k];
	  if (Hash[x]) {
	    if (level[x] != l)
	      c1++;
	    else
	      c2++;
	  }
	}
      }
    }

    for (UINT_t p = s ; p<e ; p++)
      Hash[Ai[p]] = false;
  }

  free_queue(queue);

  free(Hash);
  free(visited);
  free(level);

  return c1 + (c2/3);
}



UINT_t tc_bader4(const GRAPH_TYPE *graph) {
  /* Bader's new algorithm for triangle counting based on BFS */
  /* Uses Hash array to detect triangles (v, w, x) if x is adjacent to v */
  /* For level[], 0 == unvisited. Needs a modified BFS starting from level 1 */
  /* Mark horizontal edges during BFS */
  /* Direction orientied. */
  UINT_t* restrict level;
  UINT_t c1, c2;
  register UINT_t x;
  bool *Hash;
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

  Hash = (bool *)calloc(n, sizeof(bool));
  assert_malloc(Hash);

  horiz = (bool *)malloc(m * sizeof(bool));
  assert_malloc(horiz);

  Queue *queue = createQueue(n);

  c1 = 0; c2 = 0;
  for (UINT_t v = 0 ; v < n ; v++) {
    if (!level[v])
      bfs_mark_horizontal_edges(graph, v, level, queue, visited, horiz);
    const UINT_t s = Ap[v  ];
    const UINT_t e = Ap[v+1];
    const UINT_t l = level[v];

    for (UINT_t p = s ; p<e ; p++)
      Hash[Ai[p]] = true;
    
    for (UINT_t j = s ; j<e ; j++) {
      if (horiz[j]) {
	const UINT_t w = Ai[j];
	if (v < w) {
	  for (UINT_t k = Ap[w]; k < Ap[w+1] ; k++) {
	    x = Ai[k];
	    if (Hash[x]) {
	      if (level[x] != l)
		c1++;
	      else
		c2++;
	    }
	  }
	}
      }
    }

    for (UINT_t p = s ; p<e ; p++)
      Hash[Ai[p]] = false;
  }

  free_queue(queue);

  free(Hash);
  free(visited);
  free(level);
  free(horiz);

  return c1 + (c2/3);
}


UINT_t tc_bader4_degreeOrder(const GRAPH_TYPE *graph) {
  /* Bader's new algorithm for triangle counting based on BFS */
  /* Uses Hash array to detect triangles (v, w, x) if x is adjacent to v */
  /* For level[], 0 == unvisited. Needs a modified BFS starting from level 1 */
  /* Mark horizontal edges during BFS */
  /* Direction orientied. */

  GRAPH_TYPE *graph2 = reorder_graph_by_degree(graph, REORDER_HIGHEST_DEGREE_FIRST);
  return tc_bader4(graph2);
}


UINT_t tc_bader5(const GRAPH_TYPE *graph) {
  /* Bader's new algorithm for triangle counting based on BFS */
  /* Uses Hash array to detect triangles (v, w, x) if x is adjacent to v */
  /* For level[], 0 == unvisited. Needs a modified BFS starting from level 1 */
  /* Mark horizontal edges during BFS */
  /* Use directionality to only use one counter for triangles where v < w < x */
  /* Direction orientied. */
  UINT_t* restrict level;
  UINT_t count;
  register UINT_t x;
  bool *Hash;
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

  Hash = (bool *)calloc(n, sizeof(bool));
  assert_malloc(Hash);

  horiz = (bool *)malloc(m * sizeof(bool));
  assert_malloc(horiz);

  Queue *queue = createQueue(n);

  count = 0;
  for (UINT_t v = 0 ; v < n ; v++) {
    if (!level[v])
      bfs_mark_horizontal_edges(graph, v, level, queue, visited, horiz);
    const UINT_t s = Ap[v  ];
    const UINT_t e = Ap[v+1];
    const UINT_t l = level[v];

    for (UINT_t p = s ; p<e ; p++)
      Hash[Ai[p]] = true;
    
    for (UINT_t j = s ; j<e ; j++) {
      if (horiz[j]) {
	const UINT_t w = Ai[j];
	if (v < w) {
	  for (UINT_t k = Ap[w]; k < Ap[w+1] ; k++) {
	    x = Ai[k];
	    if (Hash[x]) {
	      if ( (l != level[x]) || ((l == level[x]) && (w < x)) ) {
		count++;
	      }
	    }
	  }
	}
      }
    }

    for (UINT_t p = s ; p<e ; p++)
      Hash[Ai[p]] = false;
  }

  free_queue(queue);

  free(Hash);
  free(visited);
  free(level);
  free(horiz);

  return count;
}


static UINT_t bader2_intersectSizeMergePath(const GRAPH_TYPE* graph, const UINT_t* restrict level, const UINT_t v, const UINT_t w) {
  register UINT_t vb, ve, wb, we;
  register UINT_t ptr_v, ptr_w;
  UINT_t vlist, wlist, level_v;
  UINT_t count = 0;
  
  level_v = level[v];

  vb = graph->rowPtr[v ];
  ve = graph->rowPtr[v+1];
  wb = graph->rowPtr[w  ];
  we = graph->rowPtr[w+1];


  ptr_v = vb;
  ptr_w = wb;
  while ((ptr_v < ve) && (ptr_w < we)) {
    vlist = graph->colInd[ptr_v];
    wlist = graph->colInd[ptr_w];
    if (vlist == wlist) {
      if (level_v != level[vlist])
	count++;
      else
	if ( (vlist < v) && (vlist < w) ) /* (level_v == level[vlist]) */ 
	  count++;
      ptr_v++;
      ptr_w++;
    }
    else
      if (vlist < wlist)
	ptr_v++;
      else
	ptr_w++;
  }

  return count;
}


UINT_t tc_bader2(const GRAPH_TYPE *graph) {
  /* Instead of c1, c2, use a single counter for triangles */
  /* Direction orientied. */
  UINT_t* restrict level;
  UINT_t s, e, l, w;
  UINT_t count = 0;
  UINT_t NO_LEVEL;

  level = (UINT_t *)malloc(graph->numVertices * sizeof(UINT_t));
  assert_malloc(level);
  NO_LEVEL = graph->numVertices;
  for (UINT_t i = 0 ; i < graph->numVertices ; i++) 
    level[i] = NO_LEVEL;
  
  for (UINT_t i = 0 ; i < graph->numVertices ; i++) {
    if (level[i] == NO_LEVEL) {
      bfs(graph, i, level);
    }
  }

  for (UINT_t v = 0 ; v < graph->numVertices ; v++) {
    s = graph->rowPtr[v  ];
    e = graph->rowPtr[v+1];
    l = level[v];
    for (UINT_t j = s ; j<e ; j++) {
      w = graph->colInd[j];
      if ((v < w) && (level[w] == l)) {
	count += bader2_intersectSizeMergePath(graph, level, v, w);
      }
    }
  }

  free(level);

  return count;
}


UINT_t tc_bader_forward_hash(const GRAPH_TYPE *graph) {
  /* Bader's new algorithm for triangle counting based on BFS */
  /* Uses Hash array to detect triangles (v, w, x) if x is adjacent to v */
  /* For level[], 0 == unvisited. Needs a modified BFS starting from level 1 */
  /* Mark horizontal edges during BFS */
  /* Partition edges into two sets -- horizontal and non-horizontal (spanning two levels). */
  /* Run hash intersections on the non-horizontal edge graph using the horizontal edges. */
  /* Run forward_hash on the graph induced by the horizontal edges. */
  /* Direction oriented. */
  UINT_t* level;
  UINT_t count;
  bool *Hash;
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

  Hash = (bool *)calloc(n, sizeof(bool));
  assert_malloc(Hash);

  horiz = (bool *)malloc(m * sizeof(bool));
  assert_malloc(horiz);

  Queue *queue = createQueue(n);

  for (UINT_t v=0 ; v<n ; v++)
    if (!level[v])
      bfs_mark_horizontal_edges(graph, v, level, queue, visited, horiz);

  free_queue(queue);
  free(visited);

  GRAPH_TYPE *graph0 = (GRAPH_TYPE *)malloc(sizeof(GRAPH_TYPE));
  assert_malloc(graph0);
  graph0->numVertices = n;
  graph0->numEdges = m;
  allocate_graph(graph0);
  UINT_t* restrict Ap0 = graph0->rowPtr;
  UINT_t* restrict Ai0 = graph0->colInd;

  GRAPH_TYPE *graph1 = (GRAPH_TYPE *)malloc(sizeof(GRAPH_TYPE));
  assert_malloc(graph1);
  graph1->numVertices = n;
  graph1->numEdges = m;
  allocate_graph(graph1);
  UINT_t* restrict Ap1 = graph1->rowPtr;
  UINT_t* restrict Ai1 = graph1->colInd;

  UINT_t edgeCountG0 = 0;
  UINT_t edgeCountG1 = 0;
  Ap0[0] = 0;
  Ap1[0] = 0;
  for (UINT_t v=0 ; v<n ; v++) {
    register const UINT_t s = Ap[v];
    register const UINT_t e = Ap[v+1];
    register const UINT_t lv = level[v];
    for (UINT_t j=s ; j<e ; j++) {
      register UINT_t w = Ai[j];
      if (lv == level[w]) {
	/* add (v,w) to G0 */
	Ai0[edgeCountG0] = w;
	edgeCountG0++;
      }
      else {
	/* add (v,w) to G1 */
	Ai1[edgeCountG1] = w;
	edgeCountG1++;
      }
    }
    Ap0[v+1] = edgeCountG0;
    Ap1[v+1] = edgeCountG1;
  }

  graph0->numEdges = edgeCountG0;
  graph1->numEdges = edgeCountG1;

  count = tc_forward_hash_config_size(graph0, m);

  for (UINT_t v=0 ; v<n ; v++) {
    register const UINT_t s0 = Ap0[v  ];
    register const UINT_t e0 = Ap0[v+1];
    register const UINT_t s1 = Ap1[v  ];
    register const UINT_t e1 = Ap1[v+1];

    if (s1<e1) {

      for (UINT_t j=s1 ; j<e1 ; j++)
	Hash[Ai1[j]] = true;
    
      for (UINT_t j=s0 ; j<e0 ; j++) {
	register const UINT_t w = Ai0[j];
	if (v < w) {
	  for (UINT_t k = Ap1[w]; k < Ap1[w+1] ; k++) {
	    if (Hash[Ai1[k]]) {
	      count++;
	    }
	  }
	}
      }

      for (UINT_t j=s1 ; j<e1 ; j++)
	Hash[Ai1[j]] = false;
    }
  }
  
  free_graph(graph1);
  free_graph(graph0);

  free(horiz);
  free(Hash);
  free(level);

  return count;
}


UINT_t tc_bader_forward_hash_degreeOrder(const GRAPH_TYPE *graph) {
  /* Bader's new algorithm for triangle counting based on BFS */
  /* Uses Hash array to detect triangles (v, w, x) if x is adjacent to v */
  /* For level[], 0 == unvisited. Needs a modified BFS starting from level 1 */
  /* Mark horizontal edges during BFS */
  /* Direction orientied. */


  GRAPH_TYPE *graph2 = reorder_graph_by_degree(graph, REORDER_HIGHEST_DEGREE_FIRST);
  UINT_t count = tc_bader_forward_hash(graph2);
  free_graph(graph2);
  return count;
}



