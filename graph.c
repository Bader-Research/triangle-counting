#include "types.h"
#include "graph.h"
#ifdef PARALLEL
#include "omp.h"
#endif


void copy_graph(const GRAPH_TYPE *srcGraph, GRAPH_TYPE *dstGraph) {
  dstGraph->numVertices = srcGraph->numVertices;
  dstGraph->numEdges = srcGraph->numEdges;
  memcpy(dstGraph->rowPtr, srcGraph->rowPtr, (srcGraph->numVertices + 1) * sizeof(UINT_t));
  memcpy(dstGraph->colInd, srcGraph->colInd, srcGraph->numEdges * sizeof(UINT_t));
}

void allocate_graph(GRAPH_TYPE* graph) {
  graph->rowPtr = (UINT_t*)calloc((graph->numVertices + 1), sizeof(UINT_t));
  assert_malloc(graph->rowPtr);
  graph->colInd = (UINT_t*)calloc(graph->numEdges, sizeof(UINT_t));
  assert_malloc(graph->colInd);
}

void free_graph(GRAPH_TYPE* graph) {
    free(graph->rowPtr);
    free(graph->colInd);
    free(graph);
}

void allocate_graph_RMAT(const int scale, const int edgeFactor, GRAPH_TYPE* graph) {
    graph->numVertices = 1 << scale;
    graph->numEdges = 2 * graph->numVertices * edgeFactor; /* Factor of 2 is to store undirected edges (a, b) and (b, a) */

    allocate_graph(graph);
}


static int compareInt_t(const void *a, const void *b) {
    UINT_t arg1 = *(const UINT_t *)a;
    UINT_t arg2 = *(const UINT_t *)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

void convert_edges_to_graph(const edge_t* edges, GRAPH_TYPE* graph) {
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  for (UINT_t i = 0 ; i < n + 1 ; i++)
    Ap[i] = 0;

  // Count the number of edges incident to each vertex
  for (UINT_t i = 0; i < m; i++) {
    UINT_t vertex = edges[i].src;
    Ap[vertex + 1]++;
  }

  // Compute the prefix sum of the rowPtr array
  for (UINT_t i = 1; i <= n; i++) {
    Ap[i] += Ap[i - 1];
  }

  // Populate the col_idx array with the destination vertices
  UINT_t *current_row = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(current_row);
  for (UINT_t i = 0; i < m; i++) {
    UINT_t src_vertex = edges[i].src;
    UINT_t dst_vertex = edges[i].dst;
    UINT_t index = Ap[src_vertex] + current_row[src_vertex];
    Ai[index] = dst_vertex;
    current_row[src_vertex]++;
  }

  // Sort the column indices within each row
  for (UINT_t i = 0; i < graph->numVertices; i++) {
    UINT_t s = Ap[i];
    UINT_t e = Ap[i + 1];
    UINT_t size = e - s;
    UINT_t *row_indices = &Ai[s];
    qsort(row_indices, size, sizeof(UINT_t), compareInt_t);
  }

  free(current_row);

}

void create_graph_RMAT(GRAPH_TYPE* graph, const UINT_t scale) {

    register int good;
    register UINT_t src, dst;

    edge_t* edges = (edge_t*)calloc(graph->numEdges, sizeof(edge_t));
    assert_malloc(edges);

    UINT_t e_start = 0;

    for (UINT_t e = e_start ; e < graph->numEdges ; e+=2) {

      good = 0;

      while (!good) {
	src = 0;
	dst = 0;

	for (UINT_t level = 0; level < scale; level++) {
	  double randNum = (double)rand() / RAND_MAX;

	  double a = 0.57, b = 0.19, c = 0.19; /* d = 1 - a - b - c */

	  if (randNum < a)
	    continue;
	  else if (randNum < a + b)
	    dst |= 1 << level;
	  else if (randNum < a + b + c)
	    src |= 1 << level;
	  else {
	    src |= 1 << level;
	    dst |= 1 << level;
	  }
	}

	good = 1;

	/* Only keep unique edges */
	for (UINT_t i = 0; i<e ; i++)
	  if ((edges[i].src == src) && (edges[i].dst == dst)) good = 0;
	/* Do not keep self-loops */
	if (src == dst) good = 0;
      }

      edges[e  ].src = src;
      edges[e  ].dst = dst;
      edges[e+1].src = dst;
      edges[e+1].dst = src;
#if DEBUG
      fprintf(stdout,"Edge[%5d]: (%5d, %5d)\n",e, edges[e].src, edges[e].dst);
#endif
    }

    convert_edges_to_graph(edges,graph);

    free(edges);
}



void print_graph(const GRAPH_TYPE* graph, FILE *outfile) {
  const UINT_t* Ap = graph->rowPtr;
  const UINT_t* Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  
  fprintf(outfile,"Number of Vertices: %u\n", n);
  fprintf(outfile,"Number of Edges: %u\n", m);
  fprintf(outfile,"RowPtr: ");
  for (UINT_t i = 0; i <= n; i++)
    fprintf(outfile,"%u ", Ap[i]);
  fprintf(outfile,"\n");
  fprintf(outfile,"ColInd: ");
  for (UINT_t i = 0; i < m; i++)
    fprintf(outfile,"%u ", Ai[i]);
  fprintf(outfile,"\n");
}

bool check_edge(const GRAPH_TYPE *graph, const UINT_t v, const UINT_t w) {

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;

  UINT_t s = Ap[v];
  UINT_t e = Ap[v+1];
  for (UINT_t i = s; i < e; i++)
    if (Ai[i] == w)
      return true;

  return false;
}

  

typedef struct {
  UINT_t degree;
  UINT_t index;
} vertexDegree_t;

static int compareVertexDegree_t(const void *a, const void *b) {
    vertexDegree_t v1 = *(const vertexDegree_t *)a;
    vertexDegree_t v2 = *(const vertexDegree_t *)b;
    if (v1.degree > v2.degree) return -1;
    if (v1.degree < v2.degree) return 1;
    if (v1.index < v2.index) return -1;
    if (v1.index > v2.index) return 1;
    return 0;
}

static int compareVertexDegreeLowestFirst_t(const void *a, const void *b) {
    vertexDegree_t v1 = *(const vertexDegree_t *)a;
    vertexDegree_t v2 = *(const vertexDegree_t *)b;
    if (v1.degree < v2.degree) return -1;
    if (v1.degree > v2.degree) return 1;
    if (v1.index < v2.index) return -1;
    if (v1.index > v2.index) return 1;
    return 0;
}

GRAPH_TYPE *reorder_graph_by_degree(const GRAPH_TYPE *graph, enum reorderDegree_t reorderDegree) {
  
  register UINT_t s;
  register UINT_t b, e;

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

  vertexDegree_t* Perm = (vertexDegree_t *)malloc(n * sizeof(vertexDegree_t));
  assert_malloc(Perm);
  
  for (UINT_t i=0; i<n ; i++) {
    Perm[i].degree = Ap[i+1] - Ap[i];
    Perm[i].index  = i;
  }

  qsort(Perm, n, sizeof(vertexDegree_t), ((reorderDegree == REORDER_HIGHEST_DEGREE_FIRST)? compareVertexDegree_t: compareVertexDegreeLowestFirst_t));

  GRAPH_TYPE *graph2;
  graph2 = (GRAPH_TYPE *)malloc(sizeof(GRAPH_TYPE));
  assert_malloc(graph2);

  graph2->numVertices = n;
  graph2->numEdges = m;
  allocate_graph(graph2);
  UINT_t* restrict Ap2 = graph2->rowPtr;
  UINT_t* restrict Ai2 = graph2->colInd;

  Ap2[0] = 0;
  for (UINT_t i=1 ; i<=n ; i++)
    Ap2[i] = Ap2[i-1] + Perm[i-1].degree;

  UINT_t *reverse = (UINT_t *)malloc(n*sizeof(UINT_t));
  assert_malloc(reverse);

  for (UINT_t i=0 ; i<n ; i++)
    reverse[Perm[i].index] = i;

  for (s = 0; s < n ; s++) {
    UINT_t ps = Perm[s].index;
    b = Ap[ps];
    e = Ap[ps+1];
    UINT_t d = 0;
    for (UINT_t i=b ; i<e ; i++) {
      Ai2[Ap2[s] + d] = reverse[Ai[i]];
      d++;
    }
  }

  free(reverse);
  free(Perm);
  free(A);
  free(Size);
  free(Hash);
  
  return graph2;
}

UINT_t intersectSizeMergePath(const GRAPH_TYPE* graph, const UINT_t v, const UINT_t w) {
  register UINT_t vb, ve, wb, we;
  register UINT_t ptr_v, ptr_w;
  UINT_t count = 0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  
  vb = Ap[v ];
  ve = Ap[v+1];
  wb = Ap[w  ];
  we = Ap[w+1];

  ptr_v = vb;
  ptr_w = wb;
  while ((ptr_v < ve) && (ptr_w < we)) {
    if (Ai[ptr_v] == Ai[ptr_w]) {
      count++;
      ptr_v++;
      ptr_w++;
    }
    else
      if (Ai[ptr_v] < Ai[ptr_w])
	ptr_v++;
      else
	ptr_w++;
  }
  return count;
}


static INT_t binarySearch(const UINT_t* list, const UINT_t start, const UINT_t end, const UINT_t target) {
  register INT_t s=start, e=end, mid;
  while (s < e) {
    mid = s + (e - s) / 2;
    if (list[mid] == target)
      return mid;

    if (list[mid] < target)
      s = mid + 1;
    else
      e = mid;
  }
  return -1;
}

UINT_t intersectSizeBinarySearch(const GRAPH_TYPE* graph, const UINT_t v, const UINT_t w) {
  register UINT_t vb, ve, wb, we;
  UINT_t count=0;

  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  if ((v<0) || (v >= n) || (w<0) || (w >= n)) {
    fprintf(stderr,"vertices out of range in intersectSize()\n");
    exit(-1);
  }
  vb = Ap[v  ];
  ve = Ap[v+1];
  wb = Ap[w  ];
  we = Ap[w+1];

  UINT_t size_v = ve-vb;
  UINT_t size_w = we-wb;

  if (size_v <= size_w) {
    for (UINT_t i=vb ; i<ve ; i++)
      if (binarySearch((UINT_t *)Ai, wb, we, Ai[i])>=0) count++;
  } else {
    for (UINT_t i=wb ; i<we ; i++)
      if (binarySearch((UINT_t *)Ai, vb, ve, Ai[i])>=0) count++;
  }

  return count;
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



UINT_t searchLists_with_partitioning(const UINT_t* list1, const INT_t s1, const INT_t e1, const UINT_t* list2, const INT_t s2, const INT_t e2) {
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


UINT_t intersectSizeHash(const GRAPH_TYPE *graph, bool *Hash, const UINT_t v, const UINT_t w) {

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


UINT_t intersectSizeMergePath_forward(const GRAPH_TYPE* graph, const UINT_t v, const UINT_t w, const UINT_t* A, const UINT_t* Size) {
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

UINT_t intersectSizeHash_forward(const GRAPH_TYPE *graph, bool *Hash, const UINT_t v, const UINT_t w, const UINT_t* A, const UINT_t* Size) {

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

#ifdef PARALLEL
UINT_t intersectSizeHash_forward_P(const GRAPH_TYPE *graph, bool *Hash, const UINT_t v, const UINT_t w, const UINT_t* A, const UINT_t* Size) {

  register UINT_t vb, ve, wb, we;
  register UINT_t s1, e1, s2, e2;
  UINT_t count = 0;
  int numThreads;
  UINT_t *mycount;

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

#pragma omp parallel
  {

    int myID = omp_get_thread_num();
    if (myID==0) {
      numThreads = omp_get_num_threads();
      
      mycount = (UINT_t *)calloc(numThreads, sizeof(UINT_t));
      assert_malloc(mycount);
    }
#pragma omp barrier
#pragma omp for
    for (UINT_t i=s1 ; i<e1 ; i++)
      Hash[A[i]] = true;

#pragma omp for
    for (UINT_t i= s2; i<e2 ; i++)
      if (Hash[A[i]]) mycount[myID]++;
  
#pragma omp for
    for (UINT_t i=s1 ; i<e1 ; i++) 
      Hash[A[i]] = false;

#pragma omp for reduction(+:count)
    for (int i = 0; i < numThreads ; i++)
      count += mycount[i];

  }

  free(mycount);
    
  return count;
}
#endif


UINT_t intersectSizeHashSkip_forward(const GRAPH_TYPE *graph, bool *Hash, const UINT_t v, const UINT_t w, const UINT_t* A, const UINT_t* Size) {

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




