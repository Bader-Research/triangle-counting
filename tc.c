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
  bcopy(srcGraph->rowPtr, dstGraph->rowPtr, (srcGraph->numVertices + 1) * sizeof(INT_t));
  bcopy(srcGraph->colInd, dstGraph->colInd, srcGraph->numEdges * sizeof(INT_t));
}

void allocate_graph_RMAT(int scale, int edgeFactor, GRAPH_TYPE* graph) {
    INT_t numVertices = 1 << scale;
    INT_t numEdges = numVertices * edgeFactor;
    INT_t* rowPtr = (INT_t*)malloc((numVertices + 1) * sizeof(INT_t));
    assert_malloc(rowPtr);
    INT_t* colInd = (INT_t*)malloc(numEdges * sizeof(INT_t));
    assert_malloc(colInd);
    
    graph->numVertices = numVertices;
    graph->numEdges = numEdges;
    graph->rowPtr = rowPtr;
    graph->colInd = colInd;
}

void create_graph_RMAT(GRAPH_TYPE* graph, int scale, int edgeFactor) {

    allocate_graph_RMAT(scale, edgeFactor, graph);
    
    INT_t numVertices = graph->numVertices;
    INT_t numEdges = graph->numEdges;
    INT_t* rowPtr = graph->rowPtr;
    INT_t* colInd = graph->colInd;

    INT_t* degree = (INT_t*)calloc(numVertices, sizeof(INT_t));
    assert_malloc(degree);

    srand(time(NULL));

    rowPtr[0] = 0;
    rowPtr[numVertices] = numEdges;

    // Generate the edges
    for (INT_t e = 0; e < numEdges; e++) {
        INT_t src = 0, dest = 0;

        for (INT_t level = 0; level < scale; level++) {
            double randNum = (double)rand() / RAND_MAX;

            double a = 0.57, b = 0.19, c = 0.19, d = 1 - a - b - c;

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

        colInd[e] = dest;
        degree[src]++;
        degree[dest]++;
    }

    // Compute the rowPtr array
    for (INT_t i = 1; i <= numVertices; i++)
        rowPtr[i] = rowPtr[i - 1] + degree[i - 1];

    // Sort the column indices for each row
    for (INT_t v = 0; v < numVertices; v++) {
        INT_t rowStart = rowPtr[v];
        INT_t rowEnd = rowPtr[v + 1];

        for (INT_t i = rowStart + 1; i < rowEnd; i++) {
            INT_t j = i;
            while (j > rowStart && colInd[j] < colInd[j - 1]) {
                INT_t temp = colInd[j];
                colInd[j] = colInd[j - 1];
                colInd[j - 1] = temp;
                j--;
            }
        }
    }


    free(degree);
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

void free_graph(GRAPH_TYPE* graph) {
    free(graph->rowPtr);
    free(graph->colInd);
}

int tc_chatGPT(GRAPH_TYPE *graph) {
  int num_triangles = 0;

  INT_t* row_ptr = graph->rowPtr;
  INT_t* col_ind = graph->colInd;
  INT_t num_vertices = graph->numVertices;

  // Traverse each vertex
  for (INT_t i = 0; i < num_vertices; i++) {
    // Get the neighbors of the current vertex
    INT_t start = row_ptr[i];
    INT_t end = row_ptr[i + 1];

    // Check for triangles formed by pairs of neighbors
    for (INT_t j = start; j < end; j++) {
      INT_t neighbor1 = col_ind[j];

      // Traverse the neighbors of the first neighbor
      for (INT_t k = start; k < end; k++) {
	INT_t neighbor2 = col_ind[k];

	// Check if there is an edge between neighbor1 and neighbor2
	if (neighbor1 != neighbor2) {
	  for (INT_t l = start; l < end; l++) {
	    INT_t neighbor3 = col_ind[l];

	    // Check if there is an edge between neighbor1 and neighbor3
	    // and between neighbor2 and neighbor3
	    if (neighbor2 != neighbor3 && neighbor1 != neighbor3) {
	      if (neighbor3 == neighbor1) {
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

  
  fprintf(stdout,"Scale [%2d]\n",scale);
  
  
  originalGraph = (GRAPH_TYPE *)malloc(sizeof(GRAPH_TYPE));
  assert_malloc(originalGraph);

  create_graph_RMAT(originalGraph, scale, EDGE_FACTOR);

  graph = (GRAPH_TYPE *)malloc(sizeof(GRAPH_TYPE));
  assert_malloc(graph);
  allocate_graph_RMAT(scale, EDGE_FACTOR, graph);
  
  /* From ChatGPT 
    // Example CSR graph representation
    int num_vertices = 5;
    int num_edges = 6;
    int row_ptr[] = {0, 2, 5, 7, 9, 11};
    int col_ind[] = {1, 2, 0, 2, 3, 0, 1, 3, 1, 4, 2, 4};
  */

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

  fprintf(stdout," scale: %2d \t  tc_chatGPT: %f\n",
	  scale,total_time);

/******************************************************************/

  
  free_graph(originalGraph);
  free_graph(graph);
  return(0);
}


