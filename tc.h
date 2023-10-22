#ifndef _TC_H
#define _TC_H

#include "types.h"
#include "graph.h"

double tc_bader_compute_k(const GRAPH_TYPE *);

UINT_t tc_davis(const GRAPH_TYPE *);
UINT_t tc_wedge(const GRAPH_TYPE *);
UINT_t tc_wedge_DO(const GRAPH_TYPE *);
UINT_t tc_triples(const GRAPH_TYPE *);
UINT_t tc_triples_DO(const GRAPH_TYPE *);
UINT_t tc_intersectMergePath(const GRAPH_TYPE *);
UINT_t tc_intersectMergePath_DO(const GRAPH_TYPE *);
UINT_t tc_intersectBinarySearch(const GRAPH_TYPE *);
UINT_t tc_intersectBinarySearch_DO(const GRAPH_TYPE *);
UINT_t tc_intersectPartition(const GRAPH_TYPE *);
UINT_t tc_intersectPartition_DO(const GRAPH_TYPE *);
UINT_t tc_intersectHash(const GRAPH_TYPE *);
UINT_t tc_intersectHash_DO(const GRAPH_TYPE *);
UINT_t tc_low(const GRAPH_TYPE *);
UINT_t tc_treelist(const GRAPH_TYPE *);
UINT_t tc_treelist2(const GRAPH_TYPE *);
UINT_t tc_forward(const GRAPH_TYPE *);
UINT_t tc_forward_hash(const GRAPH_TYPE *);
UINT_t tc_forward_hash_skip(const GRAPH_TYPE *);
UINT_t tc_forward_hash_degreeOrder(const GRAPH_TYPE *);
UINT_t tc_forward_hash_degreeOrderReverse(const GRAPH_TYPE *);
UINT_t tc_compact_forward(const GRAPH_TYPE *);
UINT_t tc_bader(const GRAPH_TYPE *);
UINT_t tc_bader2(const GRAPH_TYPE *);
UINT_t tc_bader3(const GRAPH_TYPE *);
UINT_t tc_bader4(const GRAPH_TYPE *);
UINT_t tc_bader4_degreeOrder(const GRAPH_TYPE *);
UINT_t tc_bader5(const GRAPH_TYPE *);
UINT_t tc_bader_forward_hash(const GRAPH_TYPE *);
UINT_t tc_bader_forward_hash_degreeOrder(const GRAPH_TYPE *);
UINT_t tc_bader_recursive(const GRAPH_TYPE *);
UINT_t tc_bader_hybrid(const GRAPH_TYPE *);
UINT_t tc_bader_new_bfs(const GRAPH_TYPE *);

#define _BADER_RECURSIVE_BASE 100000
#endif


