#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <limits.h>
#include <strings.h>

#define DEFAULT_SCALE  10
#define EDGE_FACTOR    16
#define LOOP_CNT       25
#define SCALE_MIN       6
#define DEBUG           0

#ifdef GCC
#define INLINE inline
/* #define INLINE */
#else
#define INLINE
#endif

#define INT_t uint32_t

typedef struct {
    INT_t numVertices;
    INT_t numEdges;
    INT_t* rowPtr;
    INT_t* colInd;
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

static INT_t correctTriangleCount;
/* Check the correctness of the triangle count.
   Return 1 if worked, 0 if failed */
int check_triangleCount(GRAPH_TYPE *graph, INT_t numTriangles) {
  return (numTriangles==correctTriangleCount);
}


void copy_graph(GRAPH_TYPE *srcGraph, GRAPH_TYPE *dstGraph) {
  dstGraph->numVertices = srcGraph->numVertices;
  dstGraph->numEdges = srcGraph->numEdges;
  memcpy(dstGraph->rowPtr, srcGraph->rowPtr, (srcGraph->numVertices + 1) * sizeof(INT_t));
  memcpy(dstGraph->colInd, srcGraph->colInd, srcGraph->numEdges * sizeof(INT_t));
}

void allocate_graph(GRAPH_TYPE* graph) {
  graph->rowPtr = (INT_t*)calloc((graph->numVertices + 1), sizeof(INT_t));
    assert_malloc(graph->rowPtr);
    graph->colInd = (INT_t*)calloc(graph->numEdges, sizeof(INT_t));
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
  INT_t src;
  INT_t dst;
} edge_t;

int compareInt_t(const void *a, const void *b) {
    INT_t arg1 = *(const INT_t *)a;
    INT_t arg2 = *(const INT_t *)b;
    if (arg1 < arg2) return -1;
    if (arg1 > arg2) return 1;
    return 0;
}

void create_graph_RMAT(GRAPH_TYPE* graph, int scale) {

    register int good;
    register INT_t src, dst;
    INT_t connectedGraph;

    edge_t* edges = (edge_t*)calloc(graph->numEdges, sizeof(edge_t));
    assert_malloc(edges);

#if 1
    connectedGraph = 0;
#else
    connectedGraph = graph->numVertices-1;
#endif
    

    for (INT_t e = 0; e < connectedGraph ; e++) {
      edges[e].src = e;
      edges[e].dst = e+1;
    }

    for (INT_t e = connectedGraph ; e < graph->numEdges ; e+=2) {

      good = 0;

      while (!good) {
	src = 0;
	dst = 0;

	for (INT_t level = 0; level < scale; level++) {
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
	for (INT_t i = 0; i<e ; i++)
	  if ((edges[i].src == src) && (edges[i].dst == dst)) good = 0;
	/* Do not keep self-loops */
	if (src == dst) good = 0;
      }

      edges[e].src = src;
      edges[e].dst = dst;
      edges[e+1].src = dst;
      edges[e+1].dst = src;
#if DEBUG
      fprintf(stdout,"Edge[%5d]: (%5d, %5d)\n",e, edges[e].src, edges[e].dst);
#endif
    }

    // Count the number of edges incident to each vertex
    for (INT_t i = 0; i < graph->numEdges; i++) {
        INT_t vertex = edges[i].src;
        graph->rowPtr[vertex + 1]++;
    }

    // Compute the prefix sum of the rowPtr array
    for (INT_t i = 1; i <= graph->numVertices; i++) {
      graph->rowPtr[i] += graph->rowPtr[i - 1];
    }

    // Populate the col_idx array with the destination vertices
    INT_t *current_row = (INT_t *)calloc(graph->numVertices, sizeof(INT_t));
    assert_malloc(current_row);
    for (INT_t i = 0; i < graph->numEdges; i++) {
        INT_t src_vertex = edges[i].src;
        INT_t dst_vertex = edges[i].dst;
        INT_t index = graph->rowPtr[src_vertex] + current_row[src_vertex];
        graph->colInd[index] = dst_vertex;
        current_row[src_vertex]++;
    }

    // Sort the column indices within each row
    for (INT_t i = 0; i < graph->numVertices; i++) {
      INT_t start = graph->rowPtr[i];
      INT_t end = graph->rowPtr[i + 1];
      INT_t size = end - start;
      INT_t *row_indices = &graph->colInd[start];
      qsort(row_indices, size, sizeof(INT_t), compareInt_t);
    }

    free(current_row);

    free(edges);

}



void print_graph(const GRAPH_TYPE* graph) {
    printf("Number of Vertices: %u\n", graph->numVertices);
    printf("Number of Edges: %u\n", graph->numEdges);
    printf("RowPtr: ");
    for (INT_t i = 0; i <= graph->numVertices; i++)
        printf("%u ", graph->rowPtr[i]);
    printf("\n");
    printf("ColInd: ");
    for (INT_t i = 0; i < graph->numEdges; i++)
        printf("%u ", graph->colInd[i]);
    printf("\n");
}



int tc_base(GRAPH_TYPE *graph) {
  int num_triangles = 0;

  INT_t* row_ptr = graph->rowPtr;
  INT_t* col_ind = graph->colInd;
  INT_t num_vertices = graph->numVertices;

  for (INT_t i = 0; i < num_vertices; i++) {
    INT_t start = row_ptr[i];
    INT_t end = row_ptr[i + 1];

    for (INT_t j = start; j < end; j++) {
      INT_t neighbor1 = col_ind[j];

      for (INT_t k = start; k < end; k++) {
	INT_t neighbor2 = col_ind[k];

	if (neighbor1 != neighbor2) {
	  INT_t start_n1 = row_ptr[neighbor1];
	  INT_t end_n1 = row_ptr[neighbor1 + 1];
	  
	  for (INT_t l = start_n1; l < end_n1; l++) {
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


  
int check_edge(GRAPH_TYPE *graph, INT_t src, INT_t dst) {
  int exists = 0;

  INT_t start = graph->rowPtr[src];
  INT_t end = graph->rowPtr[src + 1];
  for (INT_t i = start; i < end; i++)
    if (graph->colInd[i] == dst) exists = 1;

  return(exists);
}

int tc_bruteforce(GRAPH_TYPE *graph) {
  register INT_t i, j, k;
  int num_triangles = 0;


  for (i = 0; i < graph->numVertices; i++)
    for (j = 0; j < graph->numVertices; j++)
      for (k = 0; k < graph->numVertices; k++) {

	if (check_edge(graph, i, j) && check_edge(graph, j, k) && check_edge(graph, k, i))
	  num_triangles++;

      }

  return (num_triangles/6);
}

int tc_oriented(GRAPH_TYPE *graph) {
  register INT_t i, j, k;
  int num_triangles = 0;


  for (i = 0; i < graph->numVertices; i++)
    for (j = i; j < graph->numVertices; j++)
      for (k = j; k < graph->numVertices; k++) {

	if (check_edge(graph, i, j) && check_edge(graph, j, k) && check_edge(graph, k, i))
	  num_triangles++;

      }

  return num_triangles;
}

int tc_orientIntersect(GRAPH_TYPE *graph) {
  register INT_t i, j, k, n1, n2;
  register INT_t start, end;
  int num_triangles = 0;

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

int tc_intersect(GRAPH_TYPE *graph) {
  register INT_t v1, v2, i;
  register INT_t start_v1, end_v1, start_v2, end_v2;
  register INT_t ptr_v1, ptr_v2;
  int num_triangles = 0;


  for (v1 = 0; v1 < graph->numVertices; v1++) {
    start_v1 = graph->rowPtr[v1];
    end_v1 = graph->rowPtr[v1+1];
    for (i = start_v1 ; i<end_v1 ; i++) {
      v2 = graph->colInd[i];
      /* Edge (v1, v2) */
      start_v2 = graph->rowPtr[v2];
      end_v2 = graph->rowPtr[v2+1];
      ptr_v1 = start_v1;
      ptr_v2 = start_v2;
      while ((ptr_v1 < end_v1) && (ptr_v2 < end_v2)) {
	if (graph->colInd[ptr_v1] == graph->colInd[ptr_v2]) {
	  num_triangles++;
	  ptr_v1++;
	  ptr_v2++;
	}
	else
	  if (graph->colInd[ptr_v1] < graph->colInd[ptr_v2])
	    ptr_v1++;
	  else
	    ptr_v2++;
      }
    }
  }

  return (num_triangles/6);
}


void runTC(int (*f)(GRAPH_TYPE*), INT_t scale, GRAPH_TYPE *originalGraph, GRAPH_TYPE *graph, char *name) {
  int loop, err;
  double 
    total_time,
    over_time;
  INT_t numTriangles;
  
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
  INT_t numTriangles;
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
  numTriangles = tc_base(graph);
  correctTriangleCount = numTriangles;

  runTC(tc_base, scale, originalGraph, graph, "tc_base");
  runTC(tc_oriented, scale, originalGraph, graph, "tc_oriented");
  runTC(tc_orientIntersect, scale, originalGraph, graph, "tc_orientIntersect");
  runTC(tc_intersect, scale, originalGraph, graph, "tc_intersect");
  runTC(tc_bruteforce, scale, originalGraph, graph, "tc_bruteforce");
  
  free_graph(originalGraph);
  free_graph(graph);
  return(0);
}

