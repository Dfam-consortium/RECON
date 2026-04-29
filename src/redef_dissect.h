/*
 * redef_dissect.h  --  Element splitting and dissection interface
 *
 * Declarations for the element re-definition, image dissection, combo/obs
 * bookkeeping, and comparison functions used during Stage 3 (eleredef).
 * Implementations are in redef_dissect.c.
 *
 * Callers must include ele.h before this header so that ELE_INFO_t,
 * ELE_DATA_t, ELEMENT_t, IMAGE_t, IMG_DATA_t, FRAG_t, MSP_t, EDGE_TREE_t,
 * and FILE are all defined.
 */
#ifndef REDEF_DISSECT_H
#define REDEF_DISSECT_H

/* ---- Element construction ------------------------------------------------- */
ELE_DATA_t  *ele_def(IMG_DATA_t **, float);
void         generate_img_tree(ELEMENT_t *);
ELE_INFO_t  *new_element(void);
void         add_ele_info(ELE_INFO_t *);
IMG_DATA_t  *img_data_sort(IMG_DATA_t *, int);

/* ---- Element redefinition ------------------------------------------------- */
void ele_redef(ELE_INFO_t *, IMAGE_t **);

/* ---- Image dissection ----------------------------------------------------- */
void  dissect(ELE_INFO_t *);
int   too_short(FRAG_t *);
MSP_t *add_msp(MSP_t *);
void  register_image(IMAGE_t *, ELEMENT_t *);
void  put_image(IMG_DATA_t **, IMAGE_t *);
void  dump_image(IMAGE_t *);
void  remove_image(IMAGE_t *);

/* ---- Combo / obs bookkeeping ---------------------------------------------- */
void combo_update(ELE_INFO_t *);
void combo_edge_update(ELE_INFO_t *, EDGE_TREE_t **);
void combo_output(ELE_INFO_t *);
void obs_output(ELE_INFO_t *);
void fprint_ele_obs(FILE *, ELE_INFO_t *);

/* ---- Comparison functions for qsort --------------------------------------- */
int frag_cmp(const void *, const void *);
int fam_cmp(const void *, const void *);

#endif /* REDEF_DISSECT_H */
