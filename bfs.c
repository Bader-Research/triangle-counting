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

#endif

