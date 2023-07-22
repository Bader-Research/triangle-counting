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

#endif

#endif
