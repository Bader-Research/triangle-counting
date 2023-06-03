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

/* Check the correctness of the triangle count.
   Return 1 if worked, 0 if failed */
int check_triangleCount(GRAPH_TYPE *graph, INT_t numTriangles) {
  int val = 1;

  /* if bad, val = 0 */

  return val;
}


void copy_graph(GRAPH_TYPE *srcGraph, GRAPH_TYPE *dstGraph) {
  dstGraph->numVertices = srcGraph->numVertices;
  dstGraph->numEdges = srcGraph->numEdges;
  memcpy(dstGraph->rowPtr, srcGraph->rowPtr, (srcGraph->numVertices + 1) * sizeof(INT_t));
  memcpy(dstGraph->colInd, srcGraph->colInd, srcGraph->numEdges * sizeof(INT_t));
}

void allocate_graph(GRAPH_TYPE* graph) {
    graph->rowPtr = (INT_t*)malloc((graph->numVertices + 1) * sizeof(INT_t));
    assert_malloc(graph->rowPtr);
    graph->colInd = (INT_t*)malloc(graph->numEdges * sizeof(INT_t));
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

    int found;

    edge_t* edges = (edge_t*)calloc(graph->numEdges, sizeof(edge_t));
    assert_malloc(edges);

    for (INT_t e = 0; e < graph->numEdges ; e+=2) {

      found = 1;

      while (found) {
	INT_t src = 0, dest = 0;

	for (INT_t level = 0; level < scale; level++) {
	  double randNum = (double)rand() / RAND_MAX;

	  double a = 0.57, b = 0.19, c = 0.19; /* d = 1 - a - b - c */

	  if (randNum < a)
	    continue;
	  else if (randNum < a + b)
	    dest |= 1 << level;
	  else if (randNum < a + b + c)
	    src |= 1 << level;
	  else {
	    src |= 1 << level;
	    dest |= 1 << level;
	  }
	}

	edges[e].src = src;
	edges[e].dst = dest;
	edges[e+1].src = dest;
	edges[e+1].dst = src;
	
	/* Only store unique edges */
	found = 0;
	for (INT_t i = 0; i<e ; i++)
	  if ((edges[i].src == src) && (edges[i].dst == dest)) found = 1;
	if (src == dest) found = 1; /* Do not add self-loops */
      }

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


int tc_chatGPT(GRAPH_TYPE *graph) {
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

  return (num_triangles/3);
}

int
main(int argc, char **argv) {
  GRAPH_TYPE 
    *originalGraph,
    *graph;
  INT_t numTriangles;
  int 
    scale,
    loop;
  double 
    total_time,
    over_time;

  int err;

  if (argc <= 1)
    scale = DEFAULT_SCALE;
  else
    scale = atoi(argv[1]);

  if (scale <= 5) {
    fprintf(stderr, "Scale must be 6 or greater.\n");
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
  

/******************************************************************/

  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    copy_graph(originalGraph, graph);
    numTriangles = tc_chatGPT(graph);
  }
  total_time = get_seconds() - total_time;
  err = check_triangleCount(graph,numTriangles);
  if (!err) fprintf(stderr,"ERROR with tc_chatGPT\n");

  over_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    copy_graph(originalGraph, graph);
  }
  over_time = get_seconds() - over_time;

  total_time -= over_time;
  total_time /= (double)LOOP_CNT;

  fprintf(stdout," scale: %2d \t tc: %12d \t tc_chatGPT: %f\n",
	  scale,numTriangles,total_time);

/******************************************************************/

  
  free_graph(originalGraph);
  free_graph(graph);
  return(0);
}

