#ifdef PARALLEL

#ifndef _TC_PARALLEL_H
#define _TC_PARALLEL_H

UINT_t tc_triples_P(const GRAPH_TYPE *);
UINT_t tc_triples_DO_P(const GRAPH_TYPE *);
UINT_t tc_wedge_P(const GRAPH_TYPE *);
UINT_t tc_wedge_DO_P(const GRAPH_TYPE *);
UINT_t tc_intersectMergePath_P(const GRAPH_TYPE *);
UINT_t tc_intersectMergePath_DO_P(const GRAPH_TYPE *);
UINT_t tc_intersectBinarySearch_P(const GRAPH_TYPE *);
UINT_t tc_intersectBinarySearch_DO_P(const GRAPH_TYPE *);
UINT_t tc_intersectPartition_P(const GRAPH_TYPE *);
UINT_t tc_intersectPartition_DO_P(const GRAPH_TYPE *);
UINT_t tc_intersectHash_P(const GRAPH_TYPE *);
UINT_t tc_intersectHash_DO_P(const GRAPH_TYPE *);
UINT_t tc_forward_hash_P(const GRAPH_TYPE *);
UINT_t tc_bader_bfs1_P(const GRAPH_TYPE *);
UINT_t tc_bader_bfs3_P(const GRAPH_TYPE *);
UINT_t tc_bader_bfs_visited_P(const GRAPH_TYPE *);
UINT_t tc_bader_bfs_hybrid_P(const GRAPH_TYPE *);
UINT_t tc_bader_bfs_hybrid2_P(const GRAPH_TYPE *);
UINT_t tc_bader_bfs_chatgpt_P(const GRAPH_TYPE *);
UINT_t tc_bader_bfs_locks_P(const GRAPH_TYPE *);
UINT_t tc_MapJIK_P(const GRAPH_TYPE *);

  
#endif

#endif
