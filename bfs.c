#include "types.h"
#include "bfs.h"
#ifdef PARALLEL
#include "omp.h"
#endif

void bfs(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level) {
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

void bfs_visited(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level, bool *visited) {

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

  free_queue(queue);
}


void bfs_visited_P(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level, bool *visited) {

  omp_lock_t qlock;

  Queue *queue = createQueue(graph->numVertices);

  visited[startVertex] = true;
  enqueue(queue, startVertex);
  level[startVertex] = 0;

  omp_init_lock(&qlock);

#pragma omp parallel
  {
    UINT_t v;

    while (!isEmpty(queue)) {
      omp_set_lock(&qlock);
      if (!isEmpty(queue)) {
	v = dequeue(queue);
	omp_unset_lock(&qlock);
      }
      else {
	omp_unset_lock(&qlock);
	continue;
      }
      omp_set_lock(&qlock);
      for (UINT_t i = graph->rowPtr[v]; i < graph->rowPtr[v + 1]; i++) {
	UINT_t w = graph->colInd[i];
	if (!visited[w])  {
	  visited[w] = true;
	  enqueue(queue, w);
	  level[w] = level[v] + 1;
	}
      }
      omp_unset_lock(&qlock);
    }
  }

  free_queue(queue);
  omp_destroy_lock(&qlock);
}


void bfs_visited_P_DEBUG(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level, bool *visited) {

  omp_lock_t qlock;
  UINT_t qsize;
  UINT_t totalenq;
  UINT_t totaldeq;

  printf("BFS_Visited_P (%d):\n",startVertex);

  int *vstate = (int *)calloc(graph->numVertices, sizeof(int));
  assert_malloc(vstate);
  
  Queue *queue = createQueue(graph->numVertices);

  visited[startVertex] = true;
  enqueue(queue, startVertex);
  if (vstate[startVertex] != 0) printf("T%2d: ERROR: ENQ %3d\n",omp_get_thread_num(), startVertex);
  vstate[startVertex] = 1;
  qsize = 1;
  totalenq = 1;
  level[startVertex] = 0;

  omp_init_lock(&qlock);

#pragma omp parallel
  {
    UINT_t v;

    while ( (!isEmpty(queue)) || (qsize>0)) {
      omp_set_lock(&qlock);
      if (!isEmpty(queue)) {
	v = dequeue(queue);
	printf("T%2d: DEQ %d\n",omp_get_thread_num(), v);
	if (vstate[v] != 1) printf("T%2d: ERROR: DEQ %3d\n",omp_get_thread_num(), v);
	vstate[v] = 2;
	qsize--;
	totaldeq++;
	omp_unset_lock(&qlock);
      }
      else {
	omp_unset_lock(&qlock);
	continue;
      }
      omp_set_lock(&qlock);
      for (UINT_t i = graph->rowPtr[v]; i < graph->rowPtr[v + 1]; i++) {
	UINT_t w = graph->colInd[i];
	if (!visited[w])  {
	  /* omp_set_lock(&qlock); */
	  if (!visited[w]) {
	    visited[w] = true;
	    printf("T%2d: ENQ %d\n",omp_get_thread_num(), w);
	    if (vstate[w] != 0) printf("T%2d: ERROR: ENQ %3d\n",omp_get_thread_num(), w);
	    vstate[w] = 1;
	    enqueue(queue, w);
	    qsize++;
	    totalenq++;
	    
	    level[w] = level[v] + 1;
	    printf("T%2d: Level[%3d] <-- %3d  treeedge: (%3d, %3d) \n",omp_get_thread_num(), w, level[w],v, w);
	  }
	  /* omp_unset_lock(&qlock); */
	}
      }
      omp_unset_lock(&qlock);

    }
    fflush(stdout);
    #pragma omp barrier
  }

  printf("q empty: %d  qsize: %d  total enq: %d total deq: %d  n: %d\n",isEmpty(queue)?1:0, qsize, totalenq, totaldeq, graph->numVertices);


  /* CHECK */
  UINT_t *checkLevel = (UINT_t *)malloc(graph->numVertices * sizeof(UINT_t));
  assert_malloc(checkLevel);
  bool *checkVisited = (bool *)malloc(graph->numVertices * sizeof(bool));
  assert_malloc(checkVisited);

  for (int i=0 ; i<graph->numVertices; i++) checkLevel[i] = 0;
  for (int i=0 ; i<graph->numVertices; i++) checkVisited[i] = false;
  bfs_visited(graph, startVertex, checkLevel, checkVisited);

  for (int i=0 ; i<graph->numVertices; i++) {
    if (checkVisited[i]) {
      if (level[i] != checkLevel[i])
	printf("ERROR: Level[%3d]: %3d (not %3d)\n",i, level[i], checkLevel[i]);
    }
  }


  free(checkVisited);
  free(checkLevel);
  /********/

  free(vstate);
  free_queue(queue);
  omp_destroy_lock(&qlock);
}


void bfs_mark_horizontal_edges(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* restrict level, Queue* queue, bool* visited, bool* horiz) {
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


// Beamer

//    function top-down-step(frontier, next, parents)
//        for v ∈ frontier do
//            for n ∈ neighbors[v] do
//                if parents[n] = -1 then
//                        parents[n] ← v
//                        next ← next ∪ {n}
//                end if
//            end for
//        end for

void top_down_step(UINT_t* frontier, UINT_t* next, bool* visited, const GRAPH_TYPE* graph, UINT_t frontier_size, UINT_t *level) {

  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  UINT_t next_size = 0; // track number of elements in the next array
  for (UINT_t i = 0; i < frontier_size; i++) {
    UINT_t v = frontier[i];
    for (UINT_t j = Ap[v]; j < Ap[v + 1]; j++) {
      UINT_t w = Ai[j];
      if (!visited[w]) {
	visited[w] = true;
	level[w] = level[v] + 1;
	next[next_size++] = w;
      }
    }
  }
}

//
//    function bottom-up-step(frontier, next, parents)
//        for v ∈ vertices do
//            if parents[v] = -1 then
//                for n ∈ neighbors[v] do
//                    if n ∈ frontier then
//                        parents[v] ← n
//                                next ← next ∪ {v}
//                        break
//                    end if
//                end for
//            end if
//        end for

void bottom_up_step(UINT_t* frontier, UINT_t* next, bool *visited, const GRAPH_TYPE* graph, UINT_t frontier_size, UINT_t *level) {

  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  UINT_t next_size = 0;
  for (UINT_t i = 0; i < frontier_size; i++) {
    UINT_t v = frontier[i];
    if (!visited[v]) {
      for (UINT_t j = Ap[v]; j < Ap[v + 1]; j++) {
	UINT_t w = Ai[j];
	if (w == frontier[j]) {
	  visited[v] = true;
	  level[v] = level[w] + 1;
	  next[next_size++] = v;
	  break;
	}
      }
    }
  }
}

#define ALPHA 14.0
#define BETA 24.0

//function breadth-first-search(graph, source)
//    frontier ← {source}
//    next ← {}
//    parents ← [-1,-1,...-1]
//    while frontier ̸= {} do
//        top-down-step(frontier, next, parents)
//        frontier ← next
//        next ← {}
//    end while
//    return tree

// frontier - vertices that considered for exploration
void bfs_hybrid_visited(const GRAPH_TYPE* graph, const UINT_t startVertex, UINT_t* level, bool* visited) {

  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  UINT_t* frontier = (UINT_t*)malloc(n * sizeof(UINT_t));
  assert_malloc(frontier);
  UINT_t* next = (UINT_t*)malloc(n * sizeof(UINT_t));
  assert_malloc(next);

  UINT_t frontierSize = 0, nextSize = 0;
    
  // Set the initial frontier to vertex startVertex
  frontier[frontierSize++] = startVertex;
  level[startVertex] = 0;

  while (frontierSize > 0) {
    UINT_t numEdgesFrontier = 0; // Number of edges in the frontier
    for (UINT_t i = 0; i < frontierSize; i++) {
      UINT_t v = frontier[i];
      numEdgesFrontier += Ap[v + 1] - Ap[v];
    }

    UINT_t numEdgesUnexplored = 0; // Number of edges to check from unexplored vertices
    for (UINT_t v = 0; v < n; v++) {
      if (!visited[v]) {
	numEdgesUnexplored += Ap[v + 1] - Ap[v];
      }
    }

    if (numEdgesFrontier > numEdgesUnexplored / ALPHA) {
      // Use bottom-up approach
      bottom_up_step(frontier, next, visited, graph, frontierSize, level);
#if 0
      printf("USING: bottom_up_step\n");
#endif
    } else {
      // Use top-down approach
      top_down_step(frontier, next, visited, graph, frontierSize, level);
#if 0
	printf("USING: top_down_step\n");
#endif
    }

    // Swap frontier and next arrays for the next iteration
    UINT_t* temp = frontier;
    frontier = next;
    next = temp;
    frontierSize = nextSize;
    nextSize = 0; // Reset nextSize for the next iteration
    
    if (frontierSize <= n / BETA) {
      top_down_step(frontier, next, visited, graph, frontierSize, level);
#if 0
	printf("USING: top_down_step\n");
#endif
    }
  }

  free(frontier);
  free(next);
}

#ifdef PARALLEL

void top_down_step_P(UINT_t* frontier, UINT_t* next, bool* visited, const GRAPH_TYPE* graph, UINT_t frontier_size, UINT_t *level) {

  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  UINT_t next_size = 0; // track number of elements in the next array
#pragma omp for schedule(dynamic)
  for (UINT_t i = 0; i < frontier_size; i++) {
    UINT_t v = frontier[i];
    for (UINT_t j = Ap[v]; j < Ap[v + 1]; j++) {
      UINT_t w = Ai[j];
      if (!visited[w]) {
	visited[w] = true;
	level[w] = level[v] + 1;
	next[next_size++] = w;
      }
    }
  }
}

void bottom_up_step_P(UINT_t* frontier, UINT_t* next, bool *visited, const GRAPH_TYPE* graph, UINT_t frontier_size, UINT_t *level) {

  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  UINT_t next_size = 0;
#pragma omp for schedule(dynamic)
  for (UINT_t i = 0; i < frontier_size; i++) {
    UINT_t v = frontier[i];
    if (!visited[v]) {
      for (UINT_t j = Ap[v]; j < Ap[v + 1]; j++) {
	UINT_t w = Ai[j];
	if (w == frontier[j]) {
	  visited[v] = true;
	  level[v] = level[w] + 1;
	  next[next_size++] = v;
	  break;
	}
      }
    }
  }
}

void bfs_hybrid_visited_P(const GRAPH_TYPE* graph, const UINT_t startVertex, UINT_t* level, bool* visited) {

  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  UINT_t* frontier = (UINT_t*)malloc(n * sizeof(UINT_t));
  assert_malloc(frontier);
  UINT_t* next = (UINT_t*)malloc(n * sizeof(UINT_t));
  assert_malloc(next);

  UINT_t frontierSize = 0, nextSize = 0;
    
  // Set the initial frontier to vertex startVertex
  frontier[frontierSize++] = startVertex;
  level[startVertex] = 0;

  UINT_t numEdgesFrontier; // Number of edges in the frontier
  UINT_t numEdgesUnexplored; // Number of edges to check from unexplored vertices

#pragma omp parallel
    {
      while (frontierSize > 0) {
	numEdgesFrontier = 0; // Number of edges in the frontier
	numEdgesUnexplored = 0; // Number of edges to check from unexplored vertices
#pragma omp for reduction (+:numEdgesFrontier)
	for (UINT_t i = 0; i < frontierSize; i++) {
	  UINT_t v = frontier[i];
	  numEdgesFrontier += Ap[v + 1] - Ap[v];
	}

#pragma omp for reduction (+:numEdgesUnexplored)
	for (UINT_t v = 0; v < n; v++) {
	  if (!visited[v]) {
	    numEdgesUnexplored += Ap[v + 1] - Ap[v];
	  }
	}

	if (numEdgesFrontier > numEdgesUnexplored / ALPHA) {
	  bottom_up_step_P(frontier, next, visited, graph, frontierSize, level);
	} else {
	  top_down_step_P(frontier, next, visited, graph, frontierSize, level);
	}

      // Swap frontier and next arrays for the next iteration
#pragma omp single
	{
	  UINT_t* temp = frontier;
	  frontier = next;
	  next = temp;
	  frontierSize = nextSize;
	  nextSize = 0; // Reset nextSize for the next iteration
	}
    
	if (frontierSize <= n / BETA) {
	  top_down_step_P(frontier, next, visited, graph, frontierSize, level);
	}
      }
    }

  free(frontier);
  free(next);
}



void bfs_chatgpt_P(const GRAPH_TYPE* graph, const UINT_t startVertex, UINT_t* level, bool* visited) {

  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  visited[startVertex] = true;

  UINT_t *current_level_vertices = (UINT_t *)malloc(n * sizeof(UINT_t));
  assert_malloc(current_level_vertices);

  UINT_t current_level = 0;
  UINT_t current_level_size = 1;
  UINT_t curStart = 0;
  UINT_t curEnd = 1;

#pragma omp parallel
  {

#pragma omp single
    {
      current_level_vertices[0] = startVertex;
    }

    while (curEnd > curStart) {
#pragma omp for
      for (UINT_t i = curStart; i < curEnd; i++) {
	UINT_t v = current_level_vertices[i];
	UINT_t s = Ap[v];
	UINT_t e = Ap[v + 1];

	for (UINT_t j = s; j < e; j++) {
	  UINT_t w = Ai[j];
	  if (!visited[w]) {
#pragma omp critical
	    {
	      visited[w] = true;
	    }
#pragma omp critical
	    {
	      level[w] = current_level + 1;
	    }
#pragma omp critical
	    {
	      if (current_level_size >= n)
		printf("ERROR: current_level_size: %d  n: %d\n",current_level_size,n);
	      current_level_vertices[current_level_size++] = w;
	    }
	  }
	}
      }

#pragma omp single
      {
	current_level++;
	curStart = curEnd;
	curEnd = current_level_size;
      }
    }

  }

  free(current_level_vertices);
}


void bfs_locks_P(const GRAPH_TYPE* graph, const UINT_t startVertex, UINT_t* level, bool* visited) {

  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  visited[startVertex] = true;

  UINT_t *current_level_vertices = (UINT_t *)malloc(n * sizeof(UINT_t));
  assert_malloc(current_level_vertices);

  UINT_t current_level = 0;
  UINT_t current_level_size = 1;
  UINT_t curStart = 0;
  UINT_t curEnd = 1;

  omp_lock_t *vlocks = (omp_lock_t *)malloc(n * sizeof(omp_lock_t));
  assert_malloc(vlocks);
  for (int i=0 ; i<n ; i++)
    omp_init_lock(&vlocks[i]);

#pragma omp parallel
  {

#pragma omp single
    {
      current_level_vertices[0] = startVertex;
    }

    while (curEnd > curStart) {
#pragma omp for
      for (UINT_t i = curStart; i < curEnd; i++) {
	UINT_t v = current_level_vertices[i];
	UINT_t s = Ap[v];
	UINT_t e = Ap[v + 1];

	for (UINT_t j = s; j < e; j++) {
	  UINT_t w = Ai[j];
	  omp_set_lock(&vlocks[w]);
	  if (!visited[w]) {
	    visited[w] = true;
	    level[w] = current_level + 1;
#pragma omp critical
	    {
	      if (current_level_size >= n)
		printf("ERROR: current_level_size: %d  n: %d\n",current_level_size,n);
	      current_level_vertices[current_level_size++] = w;
	    }
	  }
	  omp_unset_lock(&vlocks[w]);
	}
      }

#pragma omp single
      {
	current_level++;
	curStart = curEnd;
	curEnd = current_level_size;
      }
    }

  }

  for (int i=0 ; i<n ; i++)
    omp_destroy_lock(&vlocks[i]);
  free(vlocks);
  free(current_level_vertices);
}

#if 1



/* BEAMER GAP BENCHMARK */

typedef struct {
  UINT_t *shared;
  UINT_t shared_in;
  UINT_t shared_out_start;
  UINT_t shared_out_end;
} SlidingQueue;

SlidingQueue *SQ_init(UINT_t shared_size);
void SQ_destroy(SlidingQueue *queue);
void SQ_push_back(SlidingQueue *queue, UINT_t to_add);
bool SQ_empty(const SlidingQueue *queue);
void SQ_reset(SlidingQueue *queue);
void SQ_slide_window(SlidingQueue *queue);
UINT_t SQ_begin(const SlidingQueue *queue);
UINT_t SQ_end(const SlidingQueue *queue);
UINT_t SQ_size(const SlidingQueue *queue);

SlidingQueue *SQ_init(UINT_t shared_size) {
  SlidingQueue *queue = (SlidingQueue *)malloc(sizeof(SlidingQueue));
  assert_malloc(queue);
  queue->shared = (UINT_t*)malloc(shared_size*sizeof(UINT_t));
  assert_malloc(queue->shared);
  SQ_reset(queue);
  return queue;
}

void SQ_destroy(SlidingQueue *queue) {
  free(queue->shared);
  free(queue);
  return;
}

void SQ_push_back(SlidingQueue *queue, UINT_t to_add) {
  queue->shared[queue->shared_in++] = to_add;
  return;
}

bool SQ_empty(const SlidingQueue *queue) {
  return queue->shared_out_start == queue->shared_out_end;
}

void SQ_reset(SlidingQueue *queue) {
  queue->shared_out_start = 0;
  queue->shared_out_end = 0;
  queue->shared_in = 0;
  return;
}

void SQ_slide_window(SlidingQueue *queue) {
  queue->shared_out_start = queue->shared_out_end;
  queue->shared_out_end = queue->shared_in;
  return;
}


UINT_t SQ_begin(const SlidingQueue *queue) {
  return queue->shared[queue->shared_out_start];
}

UINT_t SQ_end(const SlidingQueue *queue) {
  return queue->shared[queue->shared_out_end];
}

UINT_t SQ_size(const SlidingQueue *queue) {
  return SQ_end(queue) - SQ_begin(queue);
}

UINT_t BUStep(const GRAPH_TYPE *graph, INT_t *parent, bool *visited, bool *front, bool *next) {
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  UINT_t awake_count = 0;
  for (UINT_t i=0 ; i<n ; i++)
    next[i] = false;
#pragma omp parallel for reduction(+ : awake_count) schedule(dynamic, 1024)
  for (UINT_t u=0; u < n; u++) {
    if (parent[u] < 0) {
      UINT_t s = Ap[u];
      UINT_t e = Ap[u+1];
      for (UINT_t j=s ; j<e ; j++) {
	UINT_t v = Ai[j];
        if (front[v]) {
	  parent[u] = v;
	  visited[u] = true;
          awake_count++;
          next[u] = true;
          break;
        }
      }
    }
  }
  return awake_count;
}


UINT_t TDStep(const GRAPH_TYPE *graph, INT_t *parent, bool *visited, SlidingQueue *queue) {
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  UINT_t scout_count = 0;
#pragma omp parallel
  {
    /* QueueBuffer<UINT_t> lqueue(queue); */
    queue = SQ_init(n);
#pragma omp for reduction(+ : scout_count) nowait
    for (UINT_t i = queue->shared_out_start; i < queue->shared_out_end; i++) {
      UINT_t u = queue->shared[i];
      UINT_t s = Ap[u];
      UINT_t e = Ap[u+1];
      for (UINT_t j=s ; j<e ; j++) {
	UINT_t v = Ai[j];
        UINT_t curr_val = parent[v];
        if (curr_val < 0) {

#if 1
	  //#pragma omp atomic compare capture
#pragma omp critical
	  {
	    if (parent[v] == curr_val) {
	      parent[v] = u;
	      visited[v] = true;
	      SQ_push_back(queue, v);
	      scout_count += -curr_val;
	    }
	  }
#else	  
          if (compare_and_swap(parent[v], curr_val, u)) {
            SQ_push_back(queue, v);
            scout_count += -curr_val;
          }
#endif
        }
      }
    }
    SQ_destroy(queue);
  }
  return scout_count;
}


void QueueToBitmap(const SlidingQueue *queue, bool *bm) {
#pragma omp parallel for
  for (UINT_t i = queue->shared_out_start; i < queue->shared_out_end; i++) {
    UINT_t u = queue->shared[i];
    bm[u] = true;
  }
}

void BitmapToQueue(const GRAPH_TYPE *graph, const bool *bm, SlidingQueue *queue) {
  const UINT_t n = graph->numVertices;
#pragma omp parallel
  {
#pragma omp for nowait
    for (UINT_t i=0; i < n; i++)
      if (bm[i])
        SQ_push_back(queue, i);
    SQ_destroy(queue);
  }
  SQ_slide_window(queue);
}




#define BEAMERGAP_ALPHA 15
#define BEAMERGAP_BETA  18


void bfs_beamerGAP_P(const GRAPH_TYPE* graph, const UINT_t startVertex, UINT_t* level, bool* visited) {
  const UINT_t n = graph->numVertices;
  const UINT_t m = graph->numEdges;
  UINT_t* Ap = graph->rowPtr;
  UINT_t* Ai = graph->colInd;

  UINT_t alpha = BEAMERGAP_ALPHA;
  UINT_t beta = BEAMERGAP_BETA;
  
  UINT_t source = startVertex;

  SlidingQueue *queue = SQ_init(n);
  
  SQ_push_back(queue, source);
  SQ_slide_window(queue);

  INT_t *parent = (INT_t *)malloc(n * sizeof(INT_t));
  assert_malloc(parent);

  printf("here 1 startVertex: %d\n",startVertex);

#pragma omp parallel for
  for (UINT_t i=0 ; i<n ; i++) {
    UINT_t d = Ap[i+1] = Ap[i];
    parent[i] = (d != 0 ? -d : -1);
  }
  parent[startVertex] = startVertex;
    
  printf("here 2\n");
  bool *curr = (bool *)malloc(n * sizeof(bool));
  assert_malloc(curr);
  for (UINT_t i=0 ; i<n ; i++)
    curr[i] = false;
  bool *front = (bool *)malloc(n * sizeof(bool));
  assert_malloc(front);
  for (UINT_t i=0 ; i<n ; i++)
    front[i] = false;

  UINT_t edges_to_check = m;
  UINT_t scout_count = Ap[source+1] - Ap[source];
  while (!SQ_empty(queue)) {
    printf("here 3\n");
    if (scout_count > edges_to_check / alpha) {
      printf("here 3a\n");
      UINT_t awake_count, old_awake_count;
      QueueToBitmap(queue, front);
      awake_count = SQ_size(queue);
      SQ_slide_window(queue);
      do {
        old_awake_count = awake_count;
        awake_count = BUStep(graph, parent, visited, front, curr);
	bool *temp = front;
	front = curr;
	curr = temp;
      } while ((awake_count >= old_awake_count) ||
               (awake_count > n / beta));
      BitmapToQueue(graph, front, queue);
      scout_count = 1;
    } else {
      printf("here 3b\n");
      edges_to_check -= scout_count;
      scout_count = TDStep(graph, parent, visited, queue);
      SQ_slide_window(queue);

    }
  }

  printf("here x\n");

  SQ_destroy(queue);
  free(parent);
  free(front);
  free(curr);

  printf("here y\n");
  return;
}

#endif

#ifdef PARALLEL

// Thread-safe concurrent queue implementation
typedef struct ConcurrentQueue {
    UINT_t* items;
    UINT_t capacity;
    volatile UINT_t head;
    volatile UINT_t tail;
    omp_lock_t lock;
} ConcurrentQueue;

ConcurrentQueue* create_concurrent_queue(UINT_t capacity) {
    ConcurrentQueue* q = (ConcurrentQueue*)malloc(sizeof(ConcurrentQueue));
    assert_malloc(q);
    q->items = (UINT_t*)malloc(capacity * sizeof(UINT_t));
    assert_malloc(q->items);
    q->capacity = capacity;
    q->head = 0;
    q->tail = 0;
    omp_init_lock(&q->lock);
    return q;
}

void destroy_concurrent_queue(ConcurrentQueue* q) {
    omp_destroy_lock(&q->lock);
    free(q->items);
    free(q);
}

bool concurrent_enqueue(ConcurrentQueue* q, UINT_t item) {
    omp_set_lock(&q->lock);
    UINT_t next_tail = (q->tail + 1) % q->capacity;
    if (next_tail == q->head) {
        omp_unset_lock(&q->lock);
        return false; // Queue full
    }
    q->items[q->tail] = item;
    q->tail = next_tail;
    omp_unset_lock(&q->lock);
    return true;
}

bool concurrent_try_dequeue(ConcurrentQueue* q, UINT_t* item) {
    omp_set_lock(&q->lock);
    if (q->head == q->tail) {
        omp_unset_lock(&q->lock);
        return false; // Queue empty
    }
    *item = q->items[q->head];
    q->head = (q->head + 1) % q->capacity;
    omp_unset_lock(&q->lock);
    return true;
}

bool concurrent_is_empty(ConcurrentQueue* q) {
    omp_set_lock(&q->lock);
    bool empty = (q->head == q->tail);
    omp_unset_lock(&q->lock);
    return empty;
}

void bfs_claude_chaotic_P(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level, bool* visited) {
    const UINT_t n = graph->numVertices;
    UINT_t* Ap = graph->rowPtr;
    UINT_t* Ai = graph->colInd;
    
    INT_t* parent = (INT_t*)malloc(n * sizeof(INT_t));
    assert_malloc(parent);
    ConcurrentQueue* queue = create_concurrent_queue(n);
    
    for (UINT_t i = 0; i < n; i++) {
        parent[i] = -1;
        level[i] = UINT_MAX;
    }
    
    concurrent_enqueue(queue, startVertex);
    visited[startVertex] = true;
    parent[startVertex] = startVertex;
    level[startVertex] = 0;
    
    volatile bool done = false;
    
    #pragma omp parallel
    {
        while (!done) {
            UINT_t u;
            if (concurrent_try_dequeue(queue, &u)) {
                for (UINT_t i = Ap[u]; i < Ap[u + 1]; i++) {
                    UINT_t v = Ai[i];
                    
                    // Try to claim this vertex using compare-and-swap
                    bool expected = false;
                    if (__sync_bool_compare_and_swap(&visited[v], expected, true)) {
                        parent[v] = u;
                        level[v] = level[u] + 1;
                        concurrent_enqueue(queue, v);
                    }
                }
            }
            
            #pragma omp barrier
            #pragma omp single
            {
                done = concurrent_is_empty(queue);
            }
        }
    }
    
    destroy_concurrent_queue(queue);
    free(parent);
}
#endif


#ifdef PARALLEL

void bfs_claude_fine_lock_P(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level, bool* visited) {
    const UINT_t n = graph->numVertices;
    UINT_t* Ap = graph->rowPtr;
    UINT_t* Ai = graph->colInd;
    
    // Initialize arrays
    for (UINT_t i = 0; i < n; i++) {
        visited[i] = false;
        level[i] = UINT_MAX;
    }
    
    visited[startVertex] = true;
    level[startVertex] = 0;
    
    // Multiple queues with locks, one per thread
    int max_threads = omp_get_max_threads();
    Queue** queues = (Queue**)malloc(max_threads * sizeof(Queue*));
    assert_malloc(queues);
    omp_lock_t* queue_locks = (omp_lock_t*)malloc(max_threads * sizeof(omp_lock_t));
    assert_malloc(queue_locks);
    
    for (int i = 0; i < max_threads; i++) {
        queues[i] = createQueue(n);
        omp_init_lock(&queue_locks[i]);
    }
    
    // Add start vertex to first queue
    omp_set_lock(&queue_locks[0]);
    enqueue(queues[0], startVertex);
    omp_unset_lock(&queue_locks[0]);
    
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        int rounds_without_work = 0;
        
        while (rounds_without_work < num_threads * 2) {
            bool did_work = false;
            
            // Check all queues starting from own
            for (int q = 0; q < num_threads; q++) {
                int qid = (tid + q) % num_threads;
                
                // Try to get work from this queue
                UINT_t u;
                bool got_vertex = false;
                
                omp_set_lock(&queue_locks[qid]);
                if (!isEmpty(queues[qid])) {
                    u = dequeue(queues[qid]);
                    got_vertex = true;
                }
                omp_unset_lock(&queue_locks[qid]);
                
                if (got_vertex) {
                    did_work = true;
                    rounds_without_work = 0;
                    
                    // Process neighbors
                    for (UINT_t i = Ap[u]; i < Ap[u + 1]; i++) {
                        UINT_t v = Ai[i];
                        
                        // Use compare-and-swap for thread safety on visited array
                        bool expected = false;
                        if (__sync_bool_compare_and_swap(&visited[v], expected, true)) {
                            level[v] = level[u] + 1;
                            
                            // Add to a queue (distribute by vertex number)
                            int target_queue = v % num_threads;
                            omp_set_lock(&queue_locks[target_queue]);
                            enqueue(queues[target_queue], v);
                            omp_unset_lock(&queue_locks[target_queue]);
                        }
                    }
                    break;  // Process one vertex at a time before checking queues again
                }
            }
            
            if (!did_work) {
                rounds_without_work++;
            }
            
            // Synchronization point to check if all threads are done
            #pragma omp barrier
        }
    }
    
    // Cleanup
    for (int i = 0; i < max_threads; i++) {
        omp_destroy_lock(&queue_locks[i]);
        free_queue(queues[i]);
    }
    free(queue_locks);
    free(queues);
}

#endif


#ifdef PARALLEL

typedef struct Bag {
    UINT_t* vertices;
    UINT_t size;
    UINT_t capacity;
} Bag;

Bag* create_bag(UINT_t capacity) {
    Bag* bag = (Bag*)malloc(sizeof(Bag));
    assert_malloc(bag);
    bag->vertices = (UINT_t*)malloc(capacity * sizeof(UINT_t));
    assert_malloc(bag->vertices);
    bag->size = 0;
    bag->capacity = capacity;
    return bag;
}

void destroy_bag(Bag* bag) {
    free(bag->vertices);
    free(bag);
}

void add_to_bag(Bag* bag, UINT_t v) {
    if (bag->size >= bag->capacity) {
        bag->capacity *= 2;
        bag->vertices = (UINT_t*)realloc(bag->vertices, bag->capacity * sizeof(UINT_t));
        assert_malloc(bag->vertices);
    }
    bag->vertices[bag->size++] = v;
}

void merge_bags(Bag* dest, Bag* src) {
    for (UINT_t i = 0; i < src->size; i++) {
        add_to_bag(dest, src->vertices[i]);
    }
    src->size = 0;
}

void clear_bag(Bag* bag) {
    bag->size = 0;
}

void bfs_claude_bags_P(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level, bool* visited) {
    const UINT_t n = graph->numVertices;
    UINT_t* Ap = graph->rowPtr;
    UINT_t* Ai = graph->colInd;
    
    Bag* current_bag = create_bag(n);
    Bag* next_bag = create_bag(n);
    
    add_to_bag(current_bag, startVertex);
    visited[startVertex] = true;
    level[startVertex] = 0;
    
    UINT_t current_level = 0;
    
    while (current_bag->size > 0) {
        #pragma omp parallel
        {
            // Thread-local bag
            Bag* thread_local_bag = create_bag(n / omp_get_num_threads() + 1);
            
            #pragma omp for schedule(dynamic)
            for (UINT_t i = 0; i < current_bag->size; i++) {
                UINT_t u = current_bag->vertices[i];
                
                for (UINT_t j = Ap[u]; j < Ap[u + 1]; j++) {
                    UINT_t v = Ai[j];
                    
                    bool expected = false;
                    if (__sync_bool_compare_and_swap(&visited[v], expected, true)) {
                        level[v] = current_level + 1;
                        add_to_bag(thread_local_bag, v);
                    }
                }
            }
            
            // Merge thread-local bags into next_bag
            #pragma omp critical
            {
                merge_bags(next_bag, thread_local_bag);
            }
            
            destroy_bag(thread_local_bag);
        }
        
        // Swap bags
        Bag* temp = current_bag;
        current_bag = next_bag;
        next_bag = temp;
        clear_bag(next_bag);
        current_level++;
    }
    
    destroy_bag(current_bag);
    destroy_bag(next_bag);
}
#endif

#ifdef PARALLEL

typedef struct Deque {
    UINT_t* items;
    volatile UINT_t head;
    volatile UINT_t tail;
    UINT_t capacity;
    omp_lock_t lock;
} Deque;

Deque* create_deque(UINT_t capacity) {
    Deque* d = (Deque*)malloc(sizeof(Deque));
    assert_malloc(d);
    d->items = (UINT_t*)malloc(capacity * sizeof(UINT_t));
    assert_malloc(d->items);
    d->head = 0;
    d->tail = 0;
    d->capacity = capacity;
    omp_init_lock(&d->lock);
    return d;
}

void destroy_deque(Deque* d) {
    omp_destroy_lock(&d->lock);
    free(d->items);
    free(d);
}

bool push_back_deque(Deque* d, UINT_t item) {
    omp_set_lock(&d->lock);
    UINT_t next_tail = (d->tail + 1) % d->capacity;
    if (next_tail == d->head) {
        // Deque is full - resize it
        UINT_t old_capacity = d->capacity;
        UINT_t new_capacity = old_capacity * 2;
        UINT_t* new_items = (UINT_t*)malloc(new_capacity * sizeof(UINT_t));
        assert_malloc(new_items);
        
        // Copy elements to new array
        UINT_t count = 0;
        UINT_t curr = d->head;
        while (curr != d->tail) {
            new_items[count++] = d->items[curr];
            curr = (curr + 1) % old_capacity;
        }
        
        free(d->items);
        d->items = new_items;
        d->capacity = new_capacity;
        d->head = 0;
        d->tail = count;
        next_tail = count + 1;
    }
    d->items[d->tail] = item;
    d->tail = next_tail;
    omp_unset_lock(&d->lock);
    return true;
}

bool pop_front_deque(Deque* d, UINT_t* item) {
    omp_set_lock(&d->lock);
    if (d->head == d->tail) {
        omp_unset_lock(&d->lock);
        return false;
    }
    *item = d->items[d->head];
    d->head = (d->head + 1) % d->capacity;
    omp_unset_lock(&d->lock);
    return true;
}

bool pop_back_deque(Deque* d, UINT_t* item) {
    omp_set_lock(&d->lock);
    if (d->head == d->tail) {
        omp_unset_lock(&d->lock);
        return false;
    }
    d->tail = (d->tail - 1 + d->capacity) % d->capacity;
    *item = d->items[d->tail];
    omp_unset_lock(&d->lock);
    return true;
}

bool is_deque_empty(Deque* d) {
    omp_set_lock(&d->lock);
    bool empty = (d->head == d->tail);
    omp_unset_lock(&d->lock);
    return empty;
}

// Simple thread-safe linear congruential generator
static inline unsigned int thread_safe_rand(unsigned int* seed) {
    *seed = (*seed * 1103515245u + 12345u) & 0x7fffffff;
    return *seed;
}

void bfs_claude_work_stealing_P(const GRAPH_TYPE *graph, const UINT_t startVertex, UINT_t* level, bool* visited) {
    const UINT_t n = graph->numVertices;
    UINT_t* Ap = graph->rowPtr;
    UINT_t* Ai = graph->colInd;
    
    // Initialize arrays
    for (UINT_t i = 0; i < n; i++) {
        visited[i] = false;
        level[i] = UINT_MAX;
    }
    
    int num_threads = omp_get_max_threads();
    Deque** local_deques = (Deque**)malloc(num_threads * sizeof(Deque*));
    assert_malloc(local_deques);
    
    // Create deques with reasonable initial capacity
    for (int i = 0; i < num_threads; i++) {
        local_deques[i] = create_deque(1024);  // Start small, will resize as needed
    }
    
    visited[startVertex] = true;
    level[startVertex] = 0;
    push_back_deque(local_deques[0], startVertex);
    
    // Shared termination flag
    volatile int active_threads = num_threads;
    
    #pragma omp parallel
    {
        int tid = omp_get_thread_num();
        unsigned int seed = (unsigned int)(tid * 1000 + 1);  // Thread-specific seed
        int consecutive_steal_fails = 0;
        bool thread_active = true;
        
        while (thread_active) {
            UINT_t u;
            bool got_work = false;
            
            // Try local work first
            if (pop_front_deque(local_deques[tid], &u)) {
                got_work = true;
                consecutive_steal_fails = 0;
            } else if (consecutive_steal_fails < num_threads * 2) {
                // Try to steal from a random victim
                int attempts = 0;
                while (!got_work && attempts < num_threads) {
                    int victim = thread_safe_rand(&seed) % num_threads;
                    if (victim != tid && pop_back_deque(local_deques[victim], &u)) {
                        got_work = true;
                        consecutive_steal_fails = 0;
                    }
                    attempts++;
                }
                
                if (!got_work) {
                    consecutive_steal_fails++;
                }
            }
            
            if (got_work) {
                // Process the vertex
                for (UINT_t i = Ap[u]; i < Ap[u + 1]; i++) {
                    UINT_t v = Ai[i];
                    
                    bool expected = false;
                    if (__sync_bool_compare_and_swap(&visited[v], expected, true)) {
                        level[v] = level[u] + 1;
                        push_back_deque(local_deques[tid], v);
                    }
                }
            }
            
            // Check for global termination
            if (!got_work && consecutive_steal_fails >= num_threads * 2) {
                #pragma omp barrier
                
                // Check if any thread has work
                #pragma omp single
                {
                    bool any_work = false;
                    for (int t = 0; t < num_threads; t++) {
                        if (!is_deque_empty(local_deques[t])) {
                            any_work = true;
                            break;
                        }
                    }
                    if (!any_work) {
                        active_threads = 0;  // Signal all threads to stop
                    }
                }
                
                #pragma omp barrier
                
                if (active_threads == 0) {
                    thread_active = false;
                } else {
                    consecutive_steal_fails = 0;  // Reset and try again
                }
            }
        }
    }
    
    // Cleanup
    for (int i = 0; i < num_threads; i++) {
        destroy_deque(local_deques[i]);
    }
    free(local_deques);
}

#endif

#endif

