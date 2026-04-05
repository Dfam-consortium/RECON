/*
 * redef_edges.h  --  Edge-building and image-consistency interface
 *
 * Declarations for the edge-graph construction and partial-primary-image
 * consistency-tree functions used during element redefinition (Stage 3).
 * Implementations are in redef_edges.c.
 *
 * Callers must include ele.h and eleredef.h before this header so that
 * ELE_INFO_t, IMAGE_t, EDGE_t, EDGE_TREE_t, and IMG_NODE_t are defined.
 */
#ifndef REDEF_EDGES_H
#define REDEF_EDGES_H

/* ---- Sort comparator (used by edges_and_cps and external callers) --------- */
int partner_cmp(const void *, const void *);

/* ---- Main edge / CP driver ----------------------------------------------- */
void edges_and_cps(ELE_INFO_t *, IMAGE_t **);

/* ---- Edge helpers --------------------------------------------------------- */
int  full_length(IMAGE_t *, float);
void add_edge(ELE_INFO_t *, ELE_INFO_t *, char, int32_t, short);
void adjust_edge_tree(ELE_INFO_t *);
int  charge_edge_array(EDGE_t **, EDGE_TREE_t *, int);

/* ---- Consistency tree ----------------------------------------------------- */
int          consis_tree_build(IMG_NODE_t *, IMAGE_t *, int);
int          print_consis_tree(IMG_NODE_t *);
int          consis(IMAGE_t *, IMAGE_t *, float);
IMG_NODE_t **node_entry(IMG_NODE_t **);
void         consis_tree_free(IMG_NODE_t *);
int          find_prim(IMG_NODE_t *, float,
                       int32_t, int32_t, int32_t, int32_t,
                       int32_t, int32_t, int32_t,
                       int *, int32_t *, short *);

#endif /* REDEF_EDGES_H */
