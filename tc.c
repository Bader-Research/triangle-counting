#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <limits.h>
#include <strings.h>

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

#define ODD(n) ((n)&1)==1
#define max(a,b) ((a)>(b)?(a):(b))
#define min(a,b) ((a)<(b)?(a):(b))

struct timeval  tp;
struct timezone tzp;

#define get_seconds()   (gettimeofday(&tp, &tzp), \
                        (double)tp.tv_sec + (double)tp.tv_usec / 1000000.0)


void assert_malloc(void *ptr) {
    if (ptr==NULL) {
	fprintf(stderr,"ERROR: Null pointer\n");
	exit(1);
    }
}

static UINT_t correctTriangleCount;
/* Check the correctness of the triangle count.
   Return 1 if worked, 0 if failed */
int check_triangleCount(GRAPH_TYPE *graph, UINT_t numTriangles) {
  return (numTriangles==correctTriangleCount);
}


void copy_graph(GRAPH_TYPE *srcGraph, GRAPH_TYPE *dstGraph) {
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

void allocate_graph_RMAT(int scale, int edgeFactor, GRAPH_TYPE* graph) {
    graph->numVertices = 1 << scale;
    graph->numEdges = 2 * graph->numVertices * edgeFactor; /* Factor of 2 is to store undirected edges (a, b) and (b, a) */

    allocate_graph(graph);
}


typedef struct {
  UINT_t src;
  UINT_t dst;
} edge_t;

int compareInt_t(const void *a, const void *b) {
    UINT_t arg1 = *(const UINT_t *)a;
    UINT_t arg2 = *(const UINT_t *)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

void create_graph_RMAT(GRAPH_TYPE* graph, UINT_t scale) {

    register int good;
    register UINT_t src, dst;

    edge_t* edges = (edge_t*)calloc(graph->numEdges, sizeof(edge_t));
    assert_malloc(edges);

    UINT_t e_start = 0;
#if 0
    for (UINT_t i = 0; i < graph->numVertices - 1 ; i++) {
      edges[(2*i)  ].src = i;
      edges[(2*i)  ].dst = i+1;
      edges[(2*i)+1].src = i+1;
      edges[(2*i)+1].dst = i;
    }
    e_start = 2*(graph->numVertices - 1);
#endif
    
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

    // Count the number of edges incident to each vertex
    for (UINT_t i = 0; i < graph->numEdges; i++) {
        UINT_t vertex = edges[i].src;
        graph->rowPtr[vertex + 1]++;
    }

    // Compute the prefix sum of the rowPtr array
    for (UINT_t i = 1; i <= graph->numVertices; i++) {
      graph->rowPtr[i] += graph->rowPtr[i - 1];
    }

    // Populate the col_idx array with the destination vertices
    UINT_t *current_row = (UINT_t *)calloc(graph->numVertices, sizeof(UINT_t));
    assert_malloc(current_row);
    for (UINT_t i = 0; i < graph->numEdges; i++) {
        UINT_t src_vertex = edges[i].src;
        UINT_t dst_vertex = edges[i].dst;
        UINT_t index = graph->rowPtr[src_vertex] + current_row[src_vertex];
        graph->colInd[index] = dst_vertex;
        current_row[src_vertex]++;
    }

    // Sort the column indices within each row
    for (UINT_t i = 0; i < graph->numVertices; i++) {
      UINT_t start = graph->rowPtr[i];
      UINT_t end = graph->rowPtr[i + 1];
      UINT_t size = end - start;
      UINT_t *row_indices = &graph->colInd[start];
      qsort(row_indices, size, sizeof(UINT_t), compareInt_t);
    }

    free(current_row);

    free(edges);

}



void print_graph(const GRAPH_TYPE* graph) {
    printf("Number of Vertices: %u\n", graph->numVertices);
    printf("Number of Edges: %u\n", graph->numEdges);
    printf("RowPtr: ");
    for (UINT_t i = 0; i <= graph->numVertices; i++)
        printf("%u ", graph->rowPtr[i]);
    printf("\n");
    printf("ColInd: ");
    for (UINT_t i = 0; i < graph->numEdges; i++)
        printf("%u ", graph->colInd[i]);
    printf("\n");
}



UINT_t tc_wedge(GRAPH_TYPE *graph) {
  /* Algorithm: For each vertex i, for each open wedge (j, i, k), determine if there's a closing edge (j, k) */
  UINT_t num_triangles = 0;

  UINT_t* row_ptr = graph->rowPtr;
  UINT_t* col_ind = graph->colInd;
  UINT_t num_vertices = graph->numVertices;

  for (UINT_t i = 0; i < num_vertices; i++) {
    UINT_t start = row_ptr[i];
    UINT_t end = row_ptr[i + 1];

    for (UINT_t j = start; j < end; j++) {
      UINT_t neighbor1 = col_ind[j];

      for (UINT_t k = start; k < end; k++) {
	UINT_t neighbor2 = col_ind[k];

	if (neighbor1 != neighbor2) {
	  UINT_t start_n1 = row_ptr[neighbor1];
	  UINT_t end_n1 = row_ptr[neighbor1 + 1];
	  
	  for (UINT_t l = start_n1; l < end_n1; l++) {
	    if (col_ind[l] == neighbor2) {
	      num_triangles++;
	      break;
	    }
	  }
	}
      }
    }
  }

  return (num_triangles/6);
}

UINT_t tc_wedge_DO(GRAPH_TYPE *graph) {
  /* Algorithm: For each vertex i, for each open wedge (j, i, k), determine if there's a closing edge (j, k) */
  /* Direction oriented. */
  
  UINT_t num_triangles = 0;

  UINT_t* row_ptr = graph->rowPtr;
  UINT_t* col_ind = graph->colInd;
  UINT_t num_vertices = graph->numVertices;

  for (UINT_t i = 0; i < num_vertices; i++) {
    UINT_t start = row_ptr[i];
    UINT_t end = row_ptr[i + 1];

    for (UINT_t j = start; j < end; j++) {
      UINT_t neighbor1 = col_ind[j];
      if (neighbor1 > i) {

	for (UINT_t k = start; k < end; k++) {
	  UINT_t neighbor2 = col_ind[k];

	  if ((neighbor1 != neighbor2) && (neighbor2 > neighbor1)) {
	    UINT_t start_n1 = row_ptr[neighbor1];
	    UINT_t end_n1 = row_ptr[neighbor1 + 1];
	  
	    for (UINT_t l = start_n1; l < end_n1; l++) {
	      if (col_ind[l] == neighbor2) {
		num_triangles++;
		break;
	      }
	    }
	  }
	}
      }
    }
  }

  return num_triangles;
}


int check_edge(GRAPH_TYPE *graph, UINT_t src, UINT_t dst) {
  int exists = 0;

  UINT_t start = graph->rowPtr[src];
  UINT_t end = graph->rowPtr[src + 1];
  for (UINT_t i = start; i < end; i++)
    if (graph->colInd[i] == dst) exists = 1;

  return(exists);
}

#if 0
UINT_t tc_wedge_DO2(GRAPH_TYPE *graph) {
  /* Algorithm: For each vertex i, for each open wedge (j, i, k), determine if there's a closing edge (j, k) */
  /* Direction oriented. */
  register UINT_t i, j, k, n1, n2;
  register UINT_t start, end;
  UINT_t num_triangles = 0;

  for (i = 0; i < graph->numVertices; i++) {
    start = graph->rowPtr[i];
    end = graph->rowPtr[i+1];
    for (j = start ; j<end ; j++) {
      n1 = graph->colInd[j];
      if (n1 > i) 
	for (k = start ; k<end ; k++) {
	  n2 = graph->colInd[k];
	  if ((n2 > n1) && check_edge(graph, n1, n2))
	    num_triangles++;
	}
    }
  }

  return num_triangles;
}
#endif
  
UINT_t tc_triples(GRAPH_TYPE *graph) {
  /* Algorithm: for each triple (i, j, k), determine if the three triangle edges exist. */
  
  register UINT_t i, j, k;
  UINT_t num_triangles = 0;


  for (i = 0; i < graph->numVertices; i++)
    for (j = 0; j < graph->numVertices; j++)
      for (k = 0; k < graph->numVertices; k++) {

	if (check_edge(graph, i, j) && check_edge(graph, j, k) && check_edge(graph, k, i))
	  num_triangles++;

      }

  return (num_triangles/6);
}

UINT_t tc_triples_DO(GRAPH_TYPE *graph) {
  /* Algorithm: for each triple (i, j, k), determine if the three triangle edges exist. */
  /* Direction oriented. */
  
  register UINT_t i, j, k;
  UINT_t num_triangles = 0;


  for (i = 0; i < graph->numVertices; i++)
    for (j = i; j < graph->numVertices; j++)
      for (k = j; k < graph->numVertices; k++) {

	if (check_edge(graph, i, j) && check_edge(graph, j, k) && check_edge(graph, k, i))
	  num_triangles++;

      }

  return num_triangles;
}




UINT_t intersectSizeLinear(GRAPH_TYPE* graph, UINT_t v, UINT_t w) {
  register UINT_t vb, ve, wb, we;
  register UINT_t ptr_v, ptr_w;
  UINT_t num_triangles = 0;
  
  vb = graph->rowPtr[v ];
  ve = graph->rowPtr[v+1];
  wb = graph->rowPtr[w  ];
  we = graph->rowPtr[w+1];

  ptr_v = vb;
  ptr_w = wb;
  while ((ptr_v < ve) && (ptr_w < we)) {
    if (graph->colInd[ptr_v] == graph->colInd[ptr_w]) {
      num_triangles++;
      ptr_v++;
      ptr_w++;
    }
    else
      if (graph->colInd[ptr_v] < graph->colInd[ptr_w])
	ptr_v++;
      else
	ptr_w++;
  }
  return num_triangles;
}

UINT_t tc_intersectLin(GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a linear scan. */
  
  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t num_triangles = 0;

  for (v = 0; v < graph->numVertices; v++) {
    b = graph->rowPtr[v  ];
    e = graph->rowPtr[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w = graph->colInd[i];
      num_triangles += intersectSizeLinear(graph, v, w);
    }
  }

  return (num_triangles/6);
}

UINT_t tc_intersectLin_DO(GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a linear scan. */
  /* Direction oriented. */
  
  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t num_triangles = 0;

  for (v = 0; v < graph->numVertices; v++) {
    b = graph->rowPtr[v  ];
    e = graph->rowPtr[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w = graph->colInd[i];
      if (v < w)
	num_triangles += intersectSizeLinear(graph, v, w);
    }
  }

  return (num_triangles/3);
}

INT_t binarySearch(UINT_t* list, UINT_t start, UINT_t end, UINT_t target) {
  INT_t s=start, e=end, mid;
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

UINT_t intersectSizeLog(GRAPH_TYPE* graph, UINT_t v, UINT_t w) {
  register UINT_t vb, ve, wb, we;
  UINT_t count=0;
  if ((v<0) || (v >= graph->numVertices) || (w<0) || (w >= graph->numVertices)) {
    fprintf(stderr,"vertices out of range in intersectSize()\n");
    exit(-1);
  }
  vb = graph->rowPtr[v  ];
  ve = graph->rowPtr[v+1];
  wb = graph->rowPtr[w  ];
  we = graph->rowPtr[w+1];

  for (UINT_t i=vb ; i<ve ; i++)
    if (binarySearch(graph->colInd, wb, we, graph->colInd[i])>=0) count++;

  return count;
}

UINT_t tc_intersectLog(GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a binary search. */

  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t num_triangles = 0;

  for (v = 0; v < graph->numVertices; v++) {
    b = graph->rowPtr[v  ];
    e = graph->rowPtr[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w  = graph->colInd[i];
      num_triangles += intersectSizeLog(graph, v, w);
    }
  }

  return (num_triangles/6);
}

UINT_t tc_intersectLog_DO(GRAPH_TYPE *graph) {
  /* Algorithm: For each edge (i, j), find the size of its intersection using a binary search. */
  /* Direction oriented. */

  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t num_triangles = 0;

  for (v = 0; v < graph->numVertices; v++) {
    b = graph->rowPtr[v  ];
    e = graph->rowPtr[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w  = graph->colInd[i];
      if (v < w)
	num_triangles += intersectSizeLog(graph, v, w);
    }
  }

  return (num_triangles/3);
}

 UINT_t searchLists_with_partitioning(UINT_t* list1, UINT_t start1, UINT_t end1, UINT_t* list2, UINT_t start2, UINT_t end2) {
  UINT_t mid1, loc2;
  INT_t result;
  UINT_t count = 0;
  if ((start1>=end1)||(start2>=end2))
    return 0;
  mid1 = start1 + (end1 - start1)/2;

  result = binarySearch(list2, start2, end2, list1[mid1]); /* need a binary search that returns the item or the next higher position */
  if (result >= 0) loc2 = result;
  
  if (list1[mid1] == list2[loc2]) {
    count++;
  }
  count += searchLists_with_partitioning(list1, start1, mid1-1, list2, start2, loc2-1);
  count += searchLists_with_partitioning(list1, mid1+1, end1, list2, loc2, end2);
  return count;
}

UINT_t intersectSizePartition(GRAPH_TYPE* graph, UINT_t v, UINT_t w) {
  register UINT_t vb, ve, wb, we;

  if ((v<0) || (v >= graph->numVertices) || (w<0) || (w >= graph->numVertices)) {
    fprintf(stderr,"vertices out of range in intersectSize()\n");
    exit(-1);
  }
  vb = graph->rowPtr[v  ];
  ve = graph->rowPtr[v+1];
  wb = graph->rowPtr[w  ];
  we = graph->rowPtr[w+1];

  return searchLists_with_partitioning(graph->colInd, vb, ve, graph->colInd, wb, we);
}

UINT_t tc_intersectPartition(GRAPH_TYPE *graph) {
  register UINT_t v, w;
  register UINT_t b, e;
  UINT_t num_triangles = 0;

  for (v = 0; v < graph->numVertices; v++) {
    b = graph->rowPtr[v  ];
    e = graph->rowPtr[v+1];
    for (UINT_t i=b ; i<e ; i++) {
      w  = graph->colInd[i];
      num_triangles += intersectSizePartition(graph, v, w);
    }
  }

  return (num_triangles/6);
}


 
void runTC(UINT_t (*f)(GRAPH_TYPE*), UINT_t scale, GRAPH_TYPE *originalGraph, GRAPH_TYPE *graph, char *name) {
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

  fprintf(stdout," scale: %2d \t tc: %12d \t %s: \t",
	  scale,numTriangles,name);
  if (strlen(name) <= 12) fprintf(stdout,"\t");
  fprintf(stdout, " %f\n",total_time);

}

int
main(int argc, char **argv) {
  GRAPH_TYPE 
    *originalGraph,
    *graph;
  UINT_t numTriangles;
  int scale;

  if (argc <= 1)
    scale = DEFAULT_SCALE;
  else
    scale = atoi(argv[1]);

  if (scale < SCALE_MIN) {
    fprintf(stderr, "Scale must be %d or greater.\n",SCALE_MIN);
    exit(-1);
  }
  
#if DEBUG
  fprintf(stdout,"Scale [%2d]\n",scale);
#endif
  
  
  originalGraph = (GRAPH_TYPE *)malloc(sizeof(GRAPH_TYPE));
  assert_malloc(originalGraph);
  allocate_graph_RMAT(scale, EDGE_FACTOR, originalGraph);
    
  create_graph_RMAT(originalGraph, scale);

#if DEBUG
  print_graph(originalGraph);
#endif

  graph = (GRAPH_TYPE *)malloc(sizeof(GRAPH_TYPE));
  assert_malloc(graph);
  allocate_graph_RMAT(scale, EDGE_FACTOR, graph);
  
  copy_graph(originalGraph, graph);
  numTriangles = tc_wedge(graph);
  correctTriangleCount = numTriangles;

  runTC(tc_wedge, scale, originalGraph, graph, "tc_wedge");
  runTC(tc_wedge_DO, scale, originalGraph, graph, "tc_wedge_DO");
  runTC(tc_intersectLin, scale, originalGraph, graph, "tc_intersectLin");
  runTC(tc_intersectLin_DO, scale, originalGraph, graph, "tc_intersectLin_DO");
  runTC(tc_intersectLog, scale, originalGraph, graph, "tc_intersectLog");
  runTC(tc_intersectLog_DO, scale, originalGraph, graph, "tc_intersectLog_DO");
  /*  runTC(tc_intersectPartition, scale, originalGraph, graph, "tc_intersectPartition"); */
  runTC(tc_triples, scale, originalGraph, graph, "tc_triples");
  runTC(tc_triples_DO, scale, originalGraph, graph, "tc_triples_DO");
  
  free_graph(originalGraph);
  free_graph(graph);
  return(0);
}

