#ifndef _GRAPH_H
#define _GRAPH_H

void print_graph(const GRAPH_TYPE*, FILE*);
void convert_edges_to_graph(const edge_t*, GRAPH_TYPE*);
void copy_graph(const GRAPH_TYPE *, GRAPH_TYPE *);
bool check_triangleCount(const GRAPH_TYPE *, const UINT_t);
void allocate_graph(GRAPH_TYPE*);
void free_graph(GRAPH_TYPE*);
void allocate_graph_RMAT(const int, const int, GRAPH_TYPE*);
void create_graph_RMAT(GRAPH_TYPE*, const UINT_t);
bool check_edge(const GRAPH_TYPE *, const UINT_t, const UINT_t);
enum reorderDegree_t { REORDER_HIGHEST_DEGREE_FIRST = 0, REORDER_LOWEST_DEGREE_FIRST };
GRAPH_TYPE *reorder_graph_by_degree(const GRAPH_TYPE *, enum reorderDegree_t reorderDegree);

UINT_t intersectSizeMergePath(const GRAPH_TYPE*, const UINT_t, const UINT_t);
UINT_t intersectSizeBinarySearch(const GRAPH_TYPE*, const UINT_t, const UINT_t);
UINT_t searchLists_with_partitioning(const UINT_t*, const INT_t, const INT_t, const UINT_t*, const INT_t, const INT_t);
UINT_t intersectSizeHash(const GRAPH_TYPE *, bool *, const UINT_t, const UINT_t);

UINT_t intersectSizeMergePath_forward(const GRAPH_TYPE*, const UINT_t, const UINT_t, const UINT_t*, const UINT_t*);
UINT_t intersectSizeHash_forward(const GRAPH_TYPE *, bool *, const UINT_t, const UINT_t, const UINT_t*, const UINT_t*);
#ifdef PARALLEL
UINT_t intersectSizeHash_forward_P(const GRAPH_TYPE *, bool *, const UINT_t, const UINT_t, const UINT_t*, const UINT_t*);
#endif
UINT_t intersectSizeHashSkip_forward(const GRAPH_TYPE *, bool *, const UINT_t, const UINT_t, const UINT_t*, const UINT_t*);


#endif

