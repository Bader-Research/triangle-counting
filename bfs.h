#ifndef _BFS_H
#define _BFS_H

#include "queue.h"

void bfs(const GRAPH_TYPE *, const UINT_t, UINT_t*);
void bfs_visited(const GRAPH_TYPE *, const UINT_t, UINT_t*, bool *);
void bfs_hybrid_visited(const GRAPH_TYPE *, const UINT_t, UINT_t*, bool *);
#ifdef PARALLEL
void bfs_visited_P(const GRAPH_TYPE *, const UINT_t, UINT_t*, bool *);
void bfs_hybrid_visited_P(const GRAPH_TYPE *, const UINT_t, UINT_t*, bool *);
void bfs_chatgpt_P(const GRAPH_TYPE*, const UINT_t, UINT_t *, bool *);
void bfs_locks_P(const GRAPH_TYPE*, const UINT_t, UINT_t *, bool *);
#endif
void bfs_mark_horizontal_edges(const GRAPH_TYPE *, const UINT_t, UINT_t* restrict, Queue*, bool*, bool*);
void bfs_beamerGAP_P(const GRAPH_TYPE *, const UINT_t, UINT_t*, bool*);

#endif

