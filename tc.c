#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <limits.h>
#include <strings.h>
#include <stdbool.h>

#define DEFAULT_SCALE  10
#define EDGE_FACTOR    16
#define LOOP_CNT       10
#define SCALE_MIN       6
#define DEBUG           0

#ifdef GCC
#define INLINE inline
/* #define INLINE */
#else
#define INLINE
#endif

#define UINT_t uint32_t
#define INT_t int32_t

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

void print_graph(const GRAPH_TYPE*);
void convert_edges_to_graph(const edge_t*, GRAPH_TYPE*);
void copy_graph(const GRAPH_TYPE *, GRAPH_TYPE *);
bool check_triangleCount(const GRAPH_TYPE *, const UINT_t);
void allocate_graph(GRAPH_TYPE*);
void free_graph(GRAPH_TYPE*);
void allocate_graph_RMAT(const int, const int, GRAPH_TYPE*);
void create_graph_RMAT(GRAPH_TYPE*, const UINT_t);
bool check_edge(const GRAPH_TYPE *, const UINT_t, const UINT_t);

void benchmarkTC(UINT_t (*f)(const GRAPH_TYPE*), const GRAPH_TYPE *, GRAPH_TYPE *, const char *);

double tc_bader_compute_k(const GRAPH_TYPE *);

typedef struct {
  UINT_t *items;
  UINT_t front;
  UINT_t rear;
  UINT_t size;
} Queue;

Queue *createQueue(UINT_t);
void free_queue(Queue *);
int isEmpty(Queue *);
int isFull(Queue *);
void enqueue(Queue *, UINT_t);
UINT_t dequeue(Queue *);

void bfs(const GRAPH_TYPE *, const UINT_t, UINT_t*);
void bfs_mark_horizontal_edges(const GRAPH_TYPE *, const UINT_t, UINT_t*, Queue*, UINT_t*, bool*);


#define ODD(n) ((n)&1)==1
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

struct timeval  tp;
struct timezone tzp;

#define get_seconds()   (gettimeofday(&tp, &tzp), \
                        (double)tp.tv_sec + (double)tp.tv_usec / 1000000.0)


void assert_malloc(const void *ptr) {
    if (ptr==NULL) {
	fprintf(stderr,"ERROR: Null pointer\n");
	exit(1);
    }
}

FILE *infile = NULL, *outfile = NULL;
char *INFILENAME = NULL;
int QUIET;
int SCALE = 0;
int PRINT = 0;
int NCUBED = 1;

void usage(void) {

  printf("Triangle Counting\n\n");
  printf("Usage:\n\n");
  printf("Either one of these two must be selected:\n");
  printf(" -f <filename>   [Input Graph in Matrix Market format]\n");
  printf(" -r SCALE        [Use RMAT graph of size SCALE] (SCALE must be >= 5) \n");
  printf("Optional arguments:\n");
  printf(" -o <filename>   [Output File]\n");
  printf(" -p              [Print Input Graph]\n");
  printf(" -q              [Turn on Quiet mode]\n");
  printf(" -x              [Do not run N^3 algorithms]\n");
  exit (8);
}

void parseFlags(int argc, char **argv) {

  if (argc < 1) usage();
  infile = NULL;
  QUIET = 0;

  while ((argc > 1) && (argv[1][0] == '-')) {

    switch (argv[1][1]) {

    case 'f':
      if (!QUIET)
	printf("Input Graph: %s\n",argv[2]);
      infile = fopen(argv[2], "r");
      if (infile == NULL) usage();
      INFILENAME = argv[2];
      argv+=2;
      argc-=2;
      break;

    case 'o':
      if (!QUIET)
	printf("Output file: %s\n",argv[2]);
      outfile = fopen(argv[2], "a");
      if (outfile == NULL) usage();
      argv+=2;
      argc-=2;
      break;

    case 'r':
      SCALE = atoi(argv[2]);
      if (!QUIET)
	printf("RMAT Scale: %d\n",SCALE);
      argv+=2;
      argc-=2;
      break;

    case 'q':
      QUIET = 1;
      argv++;
      argc--;
      break;
	
    case 'p':
      PRINT = 1;
      argv++;
      argc--;
      break;
	
    case 'x':
      NCUBED = 0;
      argv++;
      argc--;
      break;
	
    default:
      fprintf(stderr,"Wrong Argument: %s\n", argv[1]);
      usage();
    }

  }

  if ((INFILENAME == NULL) && (SCALE <= 5)) usage();
  
  return;
}

static UINT_t correctTriangleCount;
/* Check the correctness of the triangle count.
   Return 1 if worked, 0 if failed */
bool check_triangleCount(const GRAPH_TYPE *graph, const UINT_t numTriangles) {
  return (numTriangles==correctTriangleCount);
}


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


int compareEdge_t(const void *a, const void *b) {
    edge_t arg1 = *(const edge_t *)a;
    edge_t arg2 = *(const edge_t *)b;
    if (arg1.src < arg2.src) return -1;
    if (arg1.src > arg2.src) return 1;
    if ((arg1.src == arg2.src) && (arg1.dst < arg2.dst)) return -1;
    if ((arg1.src == arg2.src) && (arg1.dst > arg2.dst)) return 1;
    return 0;
}

int compareInt_t(const void *a, const void *b) {
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



void print_graph(const GRAPH_TYPE* graph) {
  fprintf(outfile,"Number of Vertices: %u\n", graph->numVertices);
  fprintf(outfile,"Number of Edges: %u\n", graph->numEdges);
  fprintf(outfile,"RowPtr: ");
  for (UINT_t i = 0; i <= graph->numVertices; i++)
    fprintf(outfile,"%u ", graph->rowPtr[i]);
  fprintf(outfile,"\n");
  fprintf(outfile,"ColInd: ");
  for (UINT_t i = 0; i < graph->numEdges; i++)
    fprintf(outfile,"%u ", graph->colInd[i]);
  fprintf(outfile,"\n");
}

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

INT_t binarySearch(const UINT_t* list, const UINT_t start, const UINT_t end, const UINT_t target) {
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

UINT_t binarySearch_partition(const UINT_t* list, const UINT_t start, const UINT_t end, const UINT_t target) {
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


UINT_t intersectSizeHash(const GRAPH_TYPE *graph, bool *Hash, const UINT_t v, const UINT_t w) {

  register UINT_t vb, ve, wb, we;
  UINT_t count = 0;
  
  const UINT_t* restrict Ap = graph->rowPtr;
  const UINT_t* restrict Ai = graph->colInd;

  vb = Ap[v  ];
  ve = Ap[v+1];
  wb = Ap[w  ];
  we = Ap[w+1];

  for (UINT_t i=vb ; i<ve ; i++)
    Hash[Ai[i]] = true;

  for (UINT_t i= wb; i<we ; i++)
    if (Hash[Ai[i]]) count++;
  
  for (UINT_t i=vb ; i<ve ; i++)
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
  const UINT_t m = graph->numEdges;

  Hash = (bool *)calloc(m, sizeof(bool));
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

bool and_E(const bool* E, const UINT_t m) {
    for (UINT_t i = 0 ; i<m ; i++)
      if (E[i]) return true;
    return false;
}

bool check_edge_treelist(const GRAPH_TYPE *graph, const bool* E, const UINT_t v, const UINT_t w) {

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


void remove_treelist(const GRAPH_TYPE* graph, bool *E, const UINT_t *parent) {

  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  for (UINT_t v=0 ; v<n ; v++) {
    UINT_t s = Ap[v];
    UINT_t e = Ap[v+1];
    for (UINT_t i = s; i < e; i++)
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
  //   1a. compute a converting tree for each connected component of G;
  //   1b. for each each (u, v) in none of these tree:
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
	  if (parent[v] != n) {
	    if (parent[v] != u) {
	      if (check_edge_treelist(graph, E, parent[u], v))
		count++;
	      else if (check_edge_treelist(graph, E, parent[v], u))
		count++;
	    }
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


UINT_t intersectSizeHash_forward(const GRAPH_TYPE *graph, bool *Hash, const UINT_t v, const UINT_t w, const UINT_t* A, const UINT_t* Size) {

  register UINT_t vb, ve, wb, we;
  UINT_t count = 0;
  
  const UINT_t* restrict Ap = graph->rowPtr;

  vb = Ap[v  ];
  ve = vb + Size[v];
  wb = Ap[w  ];
  we = wb + Size[w];

  for (UINT_t i=vb ; i<ve ; i++)
    Hash[A[i]] = true;

  for (UINT_t i= wb; i<we ; i++)
    if (Hash[A[i]]) count++;
  
  for (UINT_t i=vb ; i<ve ; i++)
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





static UINT_t EMPTY;

// Function to create a new queue
Queue *createQueue(UINT_t size) {
  Queue *queue = (Queue *)malloc(sizeof(Queue));
  assert_malloc(queue);
  queue->items = (UINT_t *)malloc(size * sizeof(UINT_t));
  assert_malloc(queue->items);
  queue->front = EMPTY;
  queue->rear = EMPTY;
  queue->size = size;
  return queue;
}

void free_queue(Queue *queue) {
  free(queue->items);
  free(queue);
}

// Function to check if the queue is empty
int isEmpty(Queue *queue) {
  return queue->rear == EMPTY;
}

// Function to check if the queue is full
int isFull(Queue *queue) {
  return queue->rear == queue->size - 1;
}

// Function to add an element to the queue
void enqueue(Queue *queue, UINT_t value) {
  if (isFull(queue))
    fprintf(stderr,"Queue is full.\n");
  else {
    if (queue->front == EMPTY) {
      queue->front = 0;
      queue->rear  = 0;
      queue->items[queue->rear] = value;
    } else {
      queue->rear++;
      queue->items[queue->rear] = value;
    }
  }
}

// Function to remove an element from the queue
UINT_t dequeue(Queue *queue) {
  UINT_t item;
  if (isEmpty(queue)) {
    printf("Queue is empty.\n");
    item = EMPTY;
  } else {
    item = queue->items[queue->front];
    queue->front++;
    if (queue->front > queue->rear)
      queue->front = queue->rear = EMPTY;
  }
  return item;
}

// Function to perform breadth-first search
void bfs(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level) {
  UINT_t *visited = (UINT_t *)calloc(graph->numVertices, sizeof(UINT_t));
  assert_malloc(visited);
  Queue *queue = createQueue(graph->numVertices);

  visited[startVertex] = 1;
  enqueue(queue, startVertex);
  level[startVertex] = 0;
  
  while (!isEmpty(queue)) {
    UINT_t currentVertex = dequeue(queue);
    for (UINT_t i = graph->rowPtr[currentVertex]; i < graph->rowPtr[currentVertex + 1]; i++) {
      UINT_t adjacentVertex = graph->colInd[i];
      if (!visited[adjacentVertex])  {
	visited[adjacentVertex] = 1;
	enqueue(queue, adjacentVertex);
	level[adjacentVertex] = level[currentVertex] + 1;
      }
    }
  }

  free(visited);
  free_queue(queue);
}

void bfs_bader3(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level, Queue* queue, UINT_t* visited) {
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;

  visited[startVertex] = 1;
  enqueue(queue, startVertex);
  level[startVertex] = 1;
  
  while (!isEmpty(queue)) {
    UINT_t currentVertex = dequeue(queue);
    for (UINT_t i = Ap[currentVertex]; i < Ap[currentVertex + 1]; i++) {
      UINT_t adjacentVertex = Ai[i];
      if (!visited[adjacentVertex])  {
	visited[adjacentVertex] = 1;
	enqueue(queue, adjacentVertex);
	level[adjacentVertex] = level[currentVertex] + 1;
      }
    }
  }
}

void bfs_mark_horizontal_edges(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level, Queue* queue, UINT_t* visited, bool* horiz) {
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;

  visited[startVertex] = 1;
  enqueue(queue, startVertex);
  level[startVertex] = 1;
  
  while (!isEmpty(queue)) {
    UINT_t currentVertex = dequeue(queue);
    for (UINT_t i = Ap[currentVertex]; i < Ap[currentVertex + 1]; i++) {
      UINT_t adjacentVertex = Ai[i];
      if (!visited[adjacentVertex])  {
	horiz[i] = 0;
	visited[adjacentVertex] = 1;
	enqueue(queue, adjacentVertex);
	level[adjacentVertex] = level[currentVertex] + 1;
      }
      else 
	horiz[i] = (level[adjacentVertex] == 0) || (level[adjacentVertex] == level[currentVertex]);
    }
  }
}


void bader_intersectSizeMergePath(const GRAPH_TYPE* graph, const UINT_t* level, const UINT_t v, const UINT_t w, UINT_t* c1, UINT_t* c2) {
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
  UINT_t* level;
  UINT_t s, e, l, w;
  UINT_t c1, c2;
  UINT_t NO_LEVEL;
  UINT_t k;

  level = (UINT_t *)malloc(graph->numVertices * sizeof(UINT_t));
  assert_malloc(level);
  NO_LEVEL = graph->numVertices;
  for (UINT_t i = 0 ; i < graph->numVertices ; i++) 
    level[i] = NO_LEVEL;
  
  EMPTY = graph->numVertices;

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
  UINT_t* level;
  UINT_t s, e, l, w;
  UINT_t c1, c2;
  UINT_t NO_LEVEL;

  level = (UINT_t *)malloc(graph->numVertices * sizeof(UINT_t));
  assert_malloc(level);
  NO_LEVEL = graph->numVertices;
  for (UINT_t i = 0 ; i < graph->numVertices ; i++) 
    level[i] = NO_LEVEL;
  
  EMPTY = graph->numVertices;

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
  /* Uses Mark array to detect triangles (v, w, x) if x is adjacent to v */
  /* For level[], 0 == unvisited. Needs a modified BFS starting from level 1 */
  /* Direction orientied. */
  UINT_t* level;
  UINT_t c1, c2;
  register UINT_t x;
  bool *Mark;
  UINT_t *visited;
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;

  level = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(level);
  
  visited = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(visited);

  EMPTY = n;

  Mark = (bool *)calloc(n, sizeof(bool));
  assert_malloc(Mark);

  Queue *queue = createQueue(n);

  c1 = 0; c2 = 0;
  for (UINT_t v = 0 ; v < n ; v++) {
    if (!level[v])
      bfs_bader3(graph, v, level, queue, visited);
    const UINT_t s = Ap[v  ];
    const UINT_t e = Ap[v+1];
    const UINT_t l = level[v];

    for (UINT_t p = s ; p<e ; p++)
      Mark[Ai[p]] = 1;
    
    for (UINT_t j = s ; j<e ; j++) {
      const UINT_t w = Ai[j];
      if ((v < w) && (level[w] == l)) {
	/* bader_intersectSizeMergePath(graph, level, v, w, &c1, &c2); */
	for (UINT_t k = Ap[w]; k < Ap[w+1] ; k++) {
	  x = Ai[k];
	  if (Mark[x]) {
	    if (level[x] != l)
	      c1++;
	    else
	      c2++;
	  }
	}
      }
    }

    for (UINT_t p = s ; p<e ; p++)
      Mark[Ai[p]] = 0;
  }

  free_queue(queue);

  free(Mark);
  free(visited);
  free(level);

  return c1 + (c2/3);
}



UINT_t tc_bader4(const GRAPH_TYPE *graph) {
  /* Bader's new algorithm for triangle counting based on BFS */
  /* Uses Mark array to detect triangles (v, w, x) if x is adjacent to v */
  /* For level[], 0 == unvisited. Needs a modified BFS starting from level 1 */
  /* Mark horizontal edges during BFS */
  /* Direction orientied. */
  UINT_t* level;
  UINT_t c1, c2;
  register UINT_t x;
  bool *Mark;
  bool *horiz;
  UINT_t *visited;
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

  level = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(level);
  
  visited = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(visited);

  EMPTY = n;

  Mark = (bool *)calloc(n, sizeof(bool));
  assert_malloc(Mark);

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
      Mark[Ai[p]] = 1;
    
    for (UINT_t j = s ; j<e ; j++) {
      if (horiz[j]) {
	const UINT_t w = Ai[j];
	if (v < w) {
	  /* bader_intersectSizeMergePath(graph, level, v, w, &c1, &c2); */
	  for (UINT_t k = Ap[w]; k < Ap[w+1] ; k++) {
	    x = Ai[k];
	    if (Mark[x]) {
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
      Mark[Ai[p]] = 0;
  }

  free_queue(queue);

  free(Mark);
  free(visited);
  free(level);

  return c1 + (c2/3);
}


UINT_t tc_bader5(const GRAPH_TYPE *graph) {
  /* Bader's new algorithm for triangle counting based on BFS */
  /* Uses Mark array to detect triangles (v, w, x) if x is adjacent to v */
  /* For level[], 0 == unvisited. Needs a modified BFS starting from level 1 */
  /* Mark horizontal edges during BFS */
  /* Use directionality to only use one counter for triangles where v < w < x */
  /* Direction orientied. */
  UINT_t* level;
  UINT_t count;
  register UINT_t x;
  bool *Mark;
  bool *horiz;
  UINT_t *visited;
  const UINT_t *restrict Ap = graph->rowPtr;
  const UINT_t *restrict Ai = graph->colInd;
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;

  level = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(level);
  
  visited = (UINT_t *)calloc(n, sizeof(UINT_t));
  assert_malloc(visited);

  EMPTY = n;

  Mark = (bool *)calloc(n, sizeof(bool));
  assert_malloc(Mark);

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
      Mark[Ai[p]] = 1;
    
    for (UINT_t j = s ; j<e ; j++) {
      if (horiz[j]) {
	const UINT_t w = Ai[j];
	if (v < w) {
	  for (UINT_t k = Ap[w]; k < Ap[w+1] ; k++) {
	    x = Ai[k];
	    if (Mark[x]) {
	      if ( (l != level[x]) || ((l == level[x]) && (w < x)) ) {
		count++;
	      }
	    }
	  }
	}
      }
    }

    for (UINT_t p = s ; p<e ; p++)
      Mark[Ai[p]] = 0;
  }

  free_queue(queue);

  free(Mark);
  free(visited);
  free(level);

  return count;
}


UINT_t bader2_intersectSizeMergePath(const GRAPH_TYPE* graph, const UINT_t* level, const UINT_t v, const UINT_t w) {
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



#if 1
UINT_t k;
#endif


UINT_t tc_bader2(const GRAPH_TYPE *graph) {
  /* Instead of c1, c2, use a single counter for triangles */
  /* Direction orientied. */
  UINT_t* level;
  UINT_t s, e, l, w;
  UINT_t count = 0;
  UINT_t NO_LEVEL;
#if 1
  k=0;
#endif

  level = (UINT_t *)malloc(graph->numVertices * sizeof(UINT_t));
  assert_malloc(level);
  NO_LEVEL = graph->numVertices;
  for (UINT_t i = 0 ; i < graph->numVertices ; i++) 
    level[i] = NO_LEVEL;
  
  EMPTY = graph->numVertices;

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
#if 1
	k++;
#endif
	count += bader2_intersectSizeMergePath(graph, level, v, w);
      }
    }
  }

  free(level);

  return count;
}


UINT_t* tc_bader2_bfs(const GRAPH_TYPE *graph) {
  /* Instead of c1, c2, use a single counter for triangles */
  /* Direction orientied. */
  UINT_t* level;
  UINT_t NO_LEVEL;

  level = (UINT_t *)malloc(graph->numVertices * sizeof(UINT_t));
  assert_malloc(level);
  NO_LEVEL = graph->numVertices;
  for (UINT_t i = 0 ; i < graph->numVertices ; i++) 
    level[i] = NO_LEVEL;
  
  EMPTY = graph->numVertices;

  for (UINT_t i = 0 ; i < graph->numVertices ; i++) {
    if (level[i] == NO_LEVEL) {
      bfs(graph, i, level);
    }
  }

  return level;
}

UINT_t tc_bader2_tc(const GRAPH_TYPE *graph, UINT_t* level) {
  /* Instead of c1, c2, use a single counter for triangles */
  /* Direction orientied. */
  UINT_t s, e, l, w;
  UINT_t count = 0;
#if 1
  k=0;
#endif

  for (UINT_t v = 0 ; v < graph->numVertices ; v++) {
    s = graph->rowPtr[v  ];
    e = graph->rowPtr[v+1];
    l = level[v];
    for (UINT_t j = s ; j<e ; j++) {
      w = graph->colInd[j];
      if ((v < w) && (level[w] == l)) {
#if 1
	k++;
#endif
	count += bader2_intersectSizeMergePath(graph, level, v, w);
      }
    }
  }

  return count;
}



void runTC_bader2(UINT_t (*f)(const GRAPH_TYPE*, UINT_t*), const UINT_t scale, const GRAPH_TYPE *originalGraph, GRAPH_TYPE *graph, const char *name) {
  int loop, err;
  double 
    total_time,
    over_time;
  UINT_t numTriangles;
  UINT_t* level;

  level = tc_bader2_bfs(originalGraph);
  
  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    copy_graph(originalGraph, graph);
    numTriangles = (*f)(graph, level);
  }
  total_time = get_seconds() - total_time;
  err = check_triangleCount(graph,numTriangles);
  if (!err) fprintf(stderr,"ERROR with %s\n",name);

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    copy_graph(originalGraph, graph);
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  if (!QUIET) {
    fprintf(outfile," scale: %2d \t tc: %12d \t %s: \t",
	    scale,numTriangles,name);
    if (strlen(name) <= 12) fprintf(outfile,"\t");
    fprintf(outfile, " %9.6f\n",total_time);
  }

  free(level);

}
 
void benchmarkTC(UINT_t (*f)(const GRAPH_TYPE*), const GRAPH_TYPE *originalGraph, GRAPH_TYPE *graph, const char *name) {
  int loop, err;
  double 
    total_time,
    over_time;
  UINT_t numTriangles;
  
  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    copy_graph(originalGraph, graph);
    numTriangles = (*f)(graph);
  }
  total_time = get_seconds() - total_time;
  err = check_triangleCount(graph,numTriangles);
  if (!err) fprintf(stderr,"ERROR with %s\n",name);

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    copy_graph(originalGraph, graph);
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  if (!QUIET)
    fprintf(outfile,"%9.6f \t tc: %12d \t %12d\t %12d\t %s\n", total_time, numTriangles, graph->numVertices, graph->numEdges, name);

}

void readMatrixMarketFile(const char *filename, GRAPH_TYPE* graph) {
  FILE *infile = fopen(filename, "r");
  if (infile == NULL) {
    printf("Error opening file %s.\n", filename);
    exit(1);
  }

  char line[256];
  UINT_t num_rows, num_cols, num_entries;

  // Skip the header lines
  do {
    fgets(line, sizeof(line), infile);
  } while (line[0] == '%');

  sscanf(line, "%d %d %d", &num_rows, &num_cols, &num_entries);

  if (num_rows != num_cols) {
    fprintf(stderr,"ERROR: Matrix Market input file is not square: rows: %d  cols: %d  nnz: %d\n",
	    num_rows, num_cols, num_entries);
    exit(-1);
  }

#if DEBUG
  printf("readMatrixMarketFile: %d %d %d\n",num_rows, num_cols, num_entries);
#endif

  UINT_t edgeCount = 0;
  edge_t* edges = (edge_t*)calloc(2*num_entries, sizeof(edge_t));
  assert_malloc(edges);

  for (UINT_t i = 0; i < num_entries; i++) {
    UINT_t row, col;
    if (fscanf(infile, "%d %d\n", &row, &col) != 2) {
      fprintf(stderr,"Invalid Matrix Market file: bad entry.\n");
      fclose(infile);
      fclose(outfile);
      exit(8);
    }
    if ((row > num_rows) || (col > num_rows)) {
      fprintf(stderr,"Invalid Matrix Market file: entry out of range.\n");
      fclose(infile);
      fclose(outfile);
      exit(8);
    }

    edges[edgeCount].src = row - 1;
    edges[edgeCount].dst = col - 1;
    edgeCount++;
    edges[edgeCount].src = col - 1;
    edges[edgeCount].dst = row - 1;
    edgeCount++;
  }

  qsort(edges, edgeCount, sizeof(edge_t), compareEdge_t);

  edge_t* edgesNoDup = (edge_t*)calloc(2*num_entries, sizeof(edge_t));
  assert_malloc(edgesNoDup);

  UINT_t edgeCountNoDup;
  edge_t lastedge;
  lastedge.src = edges[0].src;
  lastedge.dst = edges[0].dst;
  edgesNoDup[0].src = edges[0].src;
  edgesNoDup[0].dst = edges[0].dst;
  edgeCountNoDup = 1;
  for (UINT_t i=0 ; i<edgeCount ; i++) {
    if (compareEdge_t(&lastedge,&edges[i])!=0) {
      edgesNoDup[edgeCountNoDup].src = edges[i].src;
      edgesNoDup[edgeCountNoDup].dst = edges[i].dst;
      edgeCountNoDup++;
      lastedge.src = edges[i].src;
      lastedge.dst = edges[i].dst;
    }
  }
    
    
  graph->numVertices = num_rows;
  graph->rowPtr = (UINT_t *)calloc((graph->numVertices + 1), sizeof(UINT_t));
  assert_malloc(graph->rowPtr);

  graph->numEdges = edgeCountNoDup;
  graph->colInd = (UINT_t *)calloc(graph->numEdges, sizeof(UINT_t));
  assert_malloc(graph->colInd);

  convert_edges_to_graph(edgesNoDup, graph);

  free(edgesNoDup);
  free(edges);

  fclose(infile);
  return;
}



int
main(int argc, char **argv) {
  UINT_t numTriangles = 0;

  outfile = stdout;

  parseFlags(argc, argv);
  
  GRAPH_TYPE 
    *originalGraph,
    *graph;

  originalGraph = (GRAPH_TYPE *)malloc(sizeof(GRAPH_TYPE));
  assert_malloc(originalGraph);
  graph = (GRAPH_TYPE *)malloc(sizeof(GRAPH_TYPE));
  assert_malloc(graph);


  if (SCALE > 5) {
    allocate_graph_RMAT(SCALE, EDGE_FACTOR, originalGraph);
    create_graph_RMAT(originalGraph, SCALE);
    allocate_graph_RMAT(SCALE, EDGE_FACTOR, graph);
  }
  else {
    if (INFILENAME != NULL) {
      readMatrixMarketFile(INFILENAME, originalGraph);
      graph->numVertices = originalGraph->numVertices;
      graph->numEdges = originalGraph->numEdges;
      allocate_graph(graph);
    }
    else {
      fprintf(stderr,"ERROR: No input graph selected.\n");
      exit(8);
    }
  }

  if (!QUIET)
    fprintf(outfile,"Graph has %d vertices and %d undirected edges.\n",originalGraph->numVertices,originalGraph->numEdges/2);

  if (PRINT)
    print_graph(originalGraph);

  if (!QUIET)
    fprintf(outfile,"%% of horizontal edges from bfs (k): %9.6f\n",tc_bader_compute_k(originalGraph));

  copy_graph(originalGraph, graph);
  numTriangles = tc_wedge(graph);
  correctTriangleCount = numTriangles;

  benchmarkTC(tc_wedge, originalGraph, graph, "tc_wedge");
  benchmarkTC(tc_wedge_DO, originalGraph, graph, "tc_wedge_DO");
  benchmarkTC(tc_intersectMergePath, originalGraph, graph, "tc_intersect_MergePath");
  benchmarkTC(tc_intersectMergePath_DO, originalGraph, graph, "tc_intersect_MergePath_DO");
  benchmarkTC(tc_intersectBinarySearch, originalGraph, graph, "tc_intersect_BinarySearch");
  benchmarkTC(tc_intersectBinarySearch_DO, originalGraph, graph, "tc_intersect_BinarySearch_DO");
  benchmarkTC(tc_intersectPartition, originalGraph, graph, "tc_intersect_Partition");
  benchmarkTC(tc_intersectPartition_DO, originalGraph, graph, "tc_intersect_Partition_DO");
  benchmarkTC(tc_intersectHash, originalGraph, graph, "tc_intersect_Hash");
  benchmarkTC(tc_intersectHash_DO, originalGraph, graph, "tc_intersect_Hash_DO");
  benchmarkTC(tc_treelist, originalGraph, graph, "tc_treelist");
  benchmarkTC(tc_forward, originalGraph, graph, "tc_forward");
  benchmarkTC(tc_davis, originalGraph, graph, "tc_davis");
  benchmarkTC(tc_low, originalGraph, graph, "tc_low");
  benchmarkTC(tc_bader, originalGraph, graph, "tc_bader");
  benchmarkTC(tc_bader2, originalGraph, graph, "tc_bader2");
#if 0
  runTC_bader2(tc_bader2_tc, originalGraph, SCALE, graph, "tc_bader2 (bfs time excluded)");
  if (!QUIET)
    printf("k: %f\n",2.0 * (double)k/(double)graph->numEdges);
#endif
  benchmarkTC(tc_bader3, originalGraph, graph, "tc_bader3");
  benchmarkTC(tc_bader4, originalGraph, graph, "tc_bader4");
  benchmarkTC(tc_bader5, originalGraph, graph, "tc_bader5");
  if (NCUBED)
    benchmarkTC(tc_triples, originalGraph, graph, "tc_triples");
  if (NCUBED)
    benchmarkTC(tc_triples_DO, originalGraph, graph, "tc_triples_DO");
  
  free_graph(originalGraph);
  free_graph(graph);

  if (!QUIET)
    fprintf(outfile,"Number of Triangles: %12d\n",numTriangles);

  fclose(outfile);
  
  return(0);
}






