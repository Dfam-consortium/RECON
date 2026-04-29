/*
 * redef_boundary.h  --  Potential-cut-point and boundary-detection interface
 *
 * Declarations for the CP/BD clustering functions used during element
 * redefinition (Stage 3).  Implementations are in redef_boundary.c.
 *
 * Callers must include ele.h before this header so that CP_t, BD_t,
 * ELEMENT_t, and ELE_INFO_t are already defined.
 */
#ifndef REDEF_BOUNDARY_H
#define REDEF_BOUNDARY_H

/* ---- CP/BD sort comparators (also used by qsort callers in other units) --- */
int CP_cmp(const void *, const void *);
int BD_cmp(const void *, const void *);

/* ---- Core clustering pipeline -------------------------------------------- */
void CP_sort(CP_t **);
void BD_sort(BD_t **);
BD_t *CP_cluster(CP_t *);
int  span(ELEMENT_t *, int32_t);
void TBD_merge(ELEMENT_t *);
void PCP_to_TBDs(ELEMENT_t *);

/* ---- CP list management --------------------------------------------------- */
void add_CP(CP_t **, int32_t, ELE_INFO_t *);
void CP_clean(CP_t **, ELE_INFO_t *);

#endif /* REDEF_BOUNDARY_H */
