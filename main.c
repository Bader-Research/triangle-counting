#include "types.h"
#include "queue.h"
#include "graph.h"
#include "bfs.h"
#include "tc.h"
#ifdef PARALLEL
#include "tc_parallel.h"
#include <omp.h>
#endif

#define DEFAULT_SCALE  10
#define EDGE_FACTOR    16
#ifndef LOOP_CNT
#define LOOP_CNT       10
#endif
#define SCALE_MIN       6
#define DEBUG           0

bool QUIET  = false;
bool PRINT  = false;
bool NCUBED = true;

#ifdef PARALLEL
bool BENCHMARK_BFS = false;
bool PARALLEL_MAX = false;
int  PARALLEL_PROCS = 0;
#endif


static void benchmarkTC(UINT_t (*f)(const GRAPH_TYPE*), const GRAPH_TYPE *, GRAPH_TYPE *, const char *);
#ifdef PARALLEL
static void benchmarkTC_P(UINT_t (*f)(const GRAPH_TYPE*), const GRAPH_TYPE *, GRAPH_TYPE *, const char *);
#endif

static FILE *infile = NULL, *outfile = NULL;
static char *INFILENAME = NULL;
static int SCALE = 0;
static bool input_selected = 0;

static void usage(void) {

  printf("Triangle Counting\n\n");
  printf("Usage:\n\n");
  printf("Either one of these two must be selected:\n");
  printf(" -f <filename>   [Input Graph in Matrix Market format]\n");
  printf(" -r SCALE        [Use RMAT graph of size SCALE] (SCALE must be >= %d) \n", SCALE_MIN);
  printf("Optional arguments:\n");
  printf(" -o <filename>   [Output File]\n");
#ifdef PARALLEL
  printf(" -p #            [Parallel: Use # threads/cores]\n");
  printf(" -P              [Parallel: Use maximum number of cores]\n");
#endif
  printf(" -d              [Display/Print Input Graph]\n");
  printf(" -q              [Turn on Quiet mode]\n");
  printf(" -x              [Do not run N^3 algorithms]\n");
  printf(" -B              [Benchmark BFS algorithms]\n");
  exit (8);
}

static void parseFlags(int argc, char **argv) {

  if (argc < 1) usage();
  infile = NULL;

  while ((argc > 1) && (argv[1][0] == '-')) {

    switch (argv[1][1]) {

    case 'f':
      if (argc < 3) usage();
      if (!QUIET)
	printf("Input Graph: %s\n",argv[2]);
      infile = fopen(argv[2], "r");
      if (infile == NULL) usage();
      INFILENAME = argv[2];
      input_selected = true;
      argv+=2;
      argc-=2;
      break;

    case 'o':
      if (argc < 3) usage();
      if (!QUIET)
	printf("Output file: %s\n",argv[2]);
      outfile = fopen(argv[2], "a");
      if (outfile == NULL) usage();
      argv+=2;
      argc-=2;
      break;

    case 'r':
      if (argc < 3) usage();
      SCALE = atoi(argv[2]);
      if (!QUIET)
	printf("RMAT Scale: %d\n",SCALE);
      INFILENAME = "RMAT";
      if (SCALE >= SCALE_MIN) input_selected = true;
      argv+=2;
      argc-=2;
      break;

    case 'q':
      QUIET = true;
      argv++;
      argc--;
      break;
	
    case 'd':
      PRINT = true;
      argv++;
      argc--;
      break;
	
#ifdef PARALLEL
    case 'B':
      BENCHMARK_BFS = true;
      argv++;
      argc--;
      break;
      
    case 'P':
      PARALLEL_MAX = true;
      argv++;
      argc--;
      break;
      
    case 'p':
      if (argc < 3) usage();
      PARALLEL_PROCS = atoi(argv[2]);
      argv+=2;
      argc-=2;
      break;
#endif
	
    case 'x':
      NCUBED = false;
      argv++;
      argc--;
      break;
	
    default:
      fprintf(stderr,"Wrong Argument: %s\n", argv[1]);
      usage();
    }

  }

  if (!input_selected) usage();
  
  return;
}

static UINT_t correctTriangleCount;
/* Check the correctness of the triangle count.
   Return 1 if worked, 0 if failed */
bool check_triangleCount(const GRAPH_TYPE *graph, const UINT_t numTriangles) {
  return (numTriangles==correctTriangleCount);
}


static void benchmarkTC(UINT_t (*f)(const GRAPH_TYPE*), const GRAPH_TYPE *originalGraph, GRAPH_TYPE *graph, const char *name) {
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

  fprintf(outfile,"TC\t%s\t%12d\t%12d\t%-30s\t%9.6f\t%12d\n",
	  INFILENAME,
	  graph->numVertices, (graph->numEdges)/2,
	  name, total_time, numTriangles);
  fflush(outfile);

}

#ifdef PARALLEL
static void benchmarkBFS(void (*f)(const GRAPH_TYPE*, const UINT_t, UINT_t*, bool*), const GRAPH_TYPE *originalGraph, const char *name) {
  int loop, err;
  double 
    total_time,
    over_time;
  bool *visited;
  UINT_t *level;
  UINT_t i, n;

  n = originalGraph->numVertices;

  visited = (bool *)malloc(n * sizeof(bool));
  assert_malloc(visited);
  level = (UINT_t *)malloc(n * sizeof(UINT_t));
  assert_malloc(level);
  
  total_time = get_seconds();
  for (loop=0 ; loop<LOOP_CNT ; loop++) {
    for(i=0; i<n ; i++) {
      visited[i] = false;
      level[i]   = n+1;
    }
    for(i=0; i<n ; i++) {
      if (!visited[i])
	(*f)(originalGraph, i, level, visited);
    }
  }
  total_time = get_seconds() - total_time;

  total_time /= (double)LOOP_CNT;

  if (name[strlen(name)-1] != 'P') {
    fprintf(outfile,"BFS\t%s\t%12d\t%12d\t%-30s\t%9.6f\n",
	    INFILENAME,
	    originalGraph->numVertices, (originalGraph->numEdges)/2,
	    name, total_time);
  }
  else {
#pragma omp parallel
#pragma omp master
    fprintf(outfile,"BFS\t%s\t%12d\t%12d\t%-30s\t%9.6f\t%12d\n",
	    INFILENAME,
	    originalGraph->numVertices, (originalGraph->numEdges)/2,
	    name, total_time,
	    omp_get_num_threads());
  }
  fflush(outfile);

  free(level);
  free(visited);
}


static void benchmarkTC_P(UINT_t (*f)(const GRAPH_TYPE*), const GRAPH_TYPE *originalGraph, GRAPH_TYPE *graph, const char *name) {
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

#pragma omp parallel
#pragma omp master
  fprintf(outfile,"TC_P\t%s\t%12d\t%12d\t%-30s\t%9.6f\t%12d\t%12d\n",
	  INFILENAME,
	  graph->numVertices, (graph->numEdges)/2,
	  name, total_time, numTriangles,
	  omp_get_num_threads());
  fflush(outfile);

}
#endif

static int compareEdge_t(const void *a, const void *b) {
    edge_t arg1 = *(const edge_t *)a;
    edge_t arg2 = *(const edge_t *)b;
    if (arg1.src < arg2.src) return -1;
    if (arg1.src > arg2.src) return 1;
    if ((arg1.src == arg2.src) && (arg1.dst < arg2.dst)) return -1;
    if ((arg1.src == arg2.src) && (arg1.dst > arg2.dst)) return 1;
    return 0;
}


static void readMatrixMarketFile(const char *filename, GRAPH_TYPE* graph) {
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
  graph->numEdges = edgeCountNoDup;
  allocate_graph(graph);

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


  if (SCALE) {
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
    fprintf(outfile,"Graph has %d vertices and %d undirected edges. Timing loop count %d.\n", originalGraph->numVertices, originalGraph->numEdges/2, LOOP_CNT);

  if (PRINT)
    print_graph(originalGraph, outfile);

  if (!QUIET)
    fprintf(outfile,"%% of horizontal edges from bfs (k): %9.6f\n",tc_bader_compute_k(originalGraph));

  copy_graph(originalGraph, graph);
  numTriangles = tc_wedge(graph);
  correctTriangleCount = numTriangles;

#ifdef PARALLEL
  if (BENCHMARK_BFS) {
    benchmarkBFS(bfs_visited, originalGraph, "bfs_visited");
    benchmarkBFS(bfs_visited_P, originalGraph, "bfs_visited_P");
    benchmarkBFS(bfs_hybrid_visited, originalGraph, "bfs_hybrid_visited");
    benchmarkBFS(bfs_hybrid_visited_P, originalGraph, "bfs_hybrid_visited_P");
    benchmarkBFS(bfs_chatgpt_P, originalGraph, "bfs_chatgpt_P");
    benchmarkBFS(bfs_locks_P, originalGraph, "bfs_locks_P");
    //    benchmarkBFS(bfs_beamerGAP_P, originalGraph, "bfs_beamerGAP_P");
    goto done;
  }
#endif

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
  benchmarkTC(tc_forward, originalGraph, graph, "tc_forward");
  benchmarkTC(tc_forward_hash, originalGraph, graph, "tc_forward_hash");
  benchmarkTC(tc_forward_hash_skip, originalGraph, graph, "tc_forward_hash_skip");
  benchmarkTC(tc_forward_hash_degreeOrder, originalGraph, graph, "tc_forward_hash_degreeOrder");
  benchmarkTC(tc_forward_hash_degreeOrderReverse, originalGraph, graph, "tc_forward_hash_degreeOrderRev");
  benchmarkTC(tc_compact_forward, originalGraph, graph, "tc_compact_forward");
  benchmarkTC(tc_davis, originalGraph, graph, "tc_davis");
  benchmarkTC(tc_low, originalGraph, graph, "tc_low");
  benchmarkTC(tc_bader, originalGraph, graph, "tc_bader");
  benchmarkTC(tc_bader2, originalGraph, graph, "tc_bader2");
  benchmarkTC(tc_bader3, originalGraph, graph, "tc_bader3");
  benchmarkTC(tc_bader4, originalGraph, graph, "tc_bader4");
  benchmarkTC(tc_bader5, originalGraph, graph, "tc_bader5");
  benchmarkTC(tc_bader4_degreeOrder, originalGraph, graph, "tc_bader4_degreeOrder");
  benchmarkTC(tc_bader_forward_hash, originalGraph, graph, "tc_bader_forward_hash");
  benchmarkTC(tc_bader_forward_hash_degreeOrder, originalGraph, graph, "tc_bader_forward_hash_degOrd");
  benchmarkTC(tc_bader_recursive, originalGraph, graph, "tc_bader_recursive");
  benchmarkTC(tc_bader_hybrid, originalGraph, graph, "tc_bader_hybrid");
  benchmarkTC(tc_bader_new_bfs, originalGraph, graph, "tc_bader_new_bfs");
  benchmarkTC(tc_treelist, originalGraph, graph, "tc_treelist");
  benchmarkTC(tc_treelist2, originalGraph, graph, "tc_treelist2");
  if (NCUBED)
    benchmarkTC(tc_triples, originalGraph, graph, "tc_triples");
  if (NCUBED)
    benchmarkTC(tc_triples_DO, originalGraph, graph, "tc_triples_DO");
  
#ifdef PARALLEL

  omp_set_num_threads(PARALLEL_MAX ? omp_get_max_threads() : PARALLEL_PROCS);
#pragma omp parallel
  if (!QUIET) {
#pragma omp master
    fprintf(outfile,"OpenMP threads: %12d\n",omp_get_num_threads());
  }

  
  benchmarkTC_P(tc_wedge_P, originalGraph, graph, "tc_wedge_P");
  benchmarkTC_P(tc_wedge_DO_P, originalGraph, graph, "tc_wedge_DO_P");
  benchmarkTC_P(tc_intersectMergePath_P, originalGraph, graph, "tc_intersect_MergePath_P");
  benchmarkTC_P(tc_intersectMergePath_DO_P, originalGraph, graph, "tc_intersect_MergePath_DO_P");
  benchmarkTC_P(tc_intersectBinarySearch_P, originalGraph, graph, "tc_intersect_BinarySearch_P");
  benchmarkTC_P(tc_intersectBinarySearch_DO_P, originalGraph, graph, "tc_intersect_BinarySearch_DO_P");
  benchmarkTC_P(tc_intersectPartition_P, originalGraph, graph, "tc_intersect_Partition_P");
  benchmarkTC_P(tc_intersectPartition_DO_P, originalGraph, graph, "tc_intersect_Partition_DO_P");
  benchmarkTC_P(tc_intersectHash_P, originalGraph, graph, "tc_intersect_Hash_P");
  benchmarkTC_P(tc_intersectHash_DO_P, originalGraph, graph, "tc_intersect_Hash_DO_P");
  benchmarkTC_P(tc_bader_bfs1_P, originalGraph, graph, "tc_bader_bfs1_P");
  benchmarkTC_P(tc_bader_bfs3_P, originalGraph, graph, "tc_bader_bfs3_P");
  benchmarkTC_P(tc_bader_bfs_visited_P, originalGraph, graph, "tc_bader_bfs_visited_P");
  benchmarkTC_P(tc_bader_bfs_hybrid_P, originalGraph, graph, "tc_bader_bfs_hybrid_P");
  benchmarkTC_P(tc_bader_bfs_hybrid2_P, originalGraph, graph, "tc_bader_bfs_hybrid2_P");
  benchmarkTC_P(tc_bader_bfs_chatgpt_P, originalGraph, graph, "tc_bader_bfs_chatgpt_P");
  benchmarkTC_P(tc_bader_bfs_locks_P, originalGraph, graph, "tc_bader_bfs_locks_P");
  benchmarkTC_P(tc_MapJIK_P, originalGraph, graph, "tc_MapJIK_P");
  benchmarkTC_P(tc_forward_hash_P, originalGraph, graph, "tc_forward_hash_P");
  if (NCUBED)
    benchmarkTC_P(tc_triples_P, originalGraph, graph, "tc_triples_P");
  if (NCUBED)
    benchmarkTC_P(tc_triples_DO_P, originalGraph, graph, "tc_triples_DO_P");
#endif

 done:
  
  free_graph(originalGraph);
  free_graph(graph);

#if 0
  if (!QUIET)
    fprintf(outfile,"Number of Triangles: %12d\n",numTriangles);
#endif

  fclose(outfile);
  
  return(0);
}

