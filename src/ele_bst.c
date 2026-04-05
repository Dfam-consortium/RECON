/*
 * ele_bst.c  --  Binary search tree operations for image and edge trees
 *
 * Implements the BST insert/delete/find/build operations for IMG_TREE_t and
 * EDGE_TREE_t.  All prototypes are declared in ele.h.
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#include "ele.h"
#include "ele_db.h"


/* ============================================================
 * Image tree operations
 * ============================================================ */

void build_img_tree(IMG_TREE_t **rt_ptr, IMAGE_t **img_array, int ori, int end) {
  int mid = (ori+end)/2;
  IMG_TREE_t *new = (IMG_TREE_t *) malloc(sizeof(IMG_TREE_t));
  new->to_image = *(img_array+mid);
  new->l = NULL;
  new->r = NULL;
  new->p = NULL;
  *rt_ptr = new;
  if (mid > ori) {
    build_img_tree(&new->l, img_array, ori, mid-1);
    new->l->p = new;
  }
  if (mid < end) {
    build_img_tree(&new->r, img_array, mid+1, end);
    new->r->p = new;
  }
}


int count_img_nodes(IMG_TREE_t *rt) {
  int ct=0;

  if (!rt) return 0;
  if (rt->l) ct += count_img_nodes(rt->l);
  if (rt->r) ct += count_img_nodes(rt->r);

  ct ++;
  return ct;
}


void insert_image(IMG_TREE_t **rt_ptr, IMAGE_t *i) {
  IMG_TREE_t *new, *x, *y;

  new = (IMG_TREE_t *) malloc(sizeof(IMG_TREE_t));
  new->to_image = i;
  new->p = NULL;
  new->l = NULL;
  new->r = NULL;

  y = NULL;
  x = *rt_ptr;
  while (x) {
    y = x;
    if (new->to_image->index < x->to_image->index) x = x->l;
    else x = x->r;
  }
  new->p = y;
  if (!y) *rt_ptr = new;
  else {
    if (new->to_image->index < y->to_image->index) y->l = new;
    else y->r = new;
  }
}


void delete_image(IMG_TREE_t **rt_ptr, IMAGE_t *i) {
  IMG_TREE_t *x, *y, *z;

  if (!(*rt_ptr)) {
    err_no ++;
    fprintf(log_file, "error:  Image tree does not exist: %d\n", i->ele_info->index);
    fflush(log_file);
    exit(2);
    return;
  }
  z= find_image(*rt_ptr, i->index);
  if (z) {
    if (!z->l || !z->r) y = z;
    else y = minimal_image(z->r);
    if (y->l) x = y->l;
    else x = y->r;
    if (x) x->p = y->p;
    if (!y->p) *rt_ptr = x;
    else {
      if (y == y->p->l) y->p->l = x;
      else y->p->r = x;
    }
    if (y != z) z->to_image = y->to_image;
    free(y);
  } else {
    err_no ++;
    fprintf(log_file, "error:  Can not find the image in the tree: %d %d\n", i->ele_info->index, i->index);
    fflush(log_file);
    exit(2);
    return;
  }
}


IMG_TREE_t *minimal_image(IMG_TREE_t *rt) {
  while (rt->l) {
    rt = rt->l;
  }

  return rt;
}


IMG_TREE_t *find_image(IMG_TREE_t *rt, int index) {
  while (rt && index != rt->to_image->index) {
    if (index < rt->to_image->index) rt = rt->l;
    else rt = rt->r;
  }

  return rt;
}


IMG_DATA_t **listify(IMG_TREE_t *rt, IMG_DATA_t **tail_ptr) {
  IMG_DATA_t *new, **new_tail_ptr;

  if (!rt) return NULL;
  if (rt->l) tail_ptr = listify(rt->l, tail_ptr);
  new = (IMG_DATA_t *) malloc(sizeof(IMG_DATA_t));
  new->to_image = rt->to_image;
  new->next = NULL;
  *tail_ptr = new;
  tail_ptr = &new->next;
  if (rt->r) tail_ptr = listify(rt->r, tail_ptr);
  return tail_ptr;
}


/* ============================================================
 * Edge tree operations
 * ============================================================ */

void build_edge_tree(EDGE_TREE_t **rt_ptr, EDGE_t **edge_array, int ori, int end) {
  int mid = (ori+end)/2;
  EDGE_TREE_t *new = (EDGE_TREE_t *) malloc(sizeof(EDGE_TREE_t));
  new->to_edge = *(edge_array+mid);
  new->l = NULL;
  new->r = NULL;
  new->p = NULL;
  *rt_ptr = new;
  if (mid > ori) {
    build_edge_tree(&new->l, edge_array, ori, mid-1);
    new->l->p = new;
  }
  if (mid < end) {
    build_edge_tree(&new->r, edge_array, mid+1, end);
    new->r->p = new;
  }
}


int count_edge_nodes(EDGE_TREE_t *rt) {
  int ct=0;

  if (!rt) return 0;
  if (rt->l) ct += count_edge_nodes(rt->l);
  if (rt->r) ct += count_edge_nodes(rt->r);

  if (rt->to_edge->type != 'c' && rt->to_edge->ele1_info->index != rt->to_edge->ele2_info->index) ct ++;
  return ct;
}


int count_total_edges(EDGE_TREE_t *rt) {
  int ct=0;

  if (!rt) return 0;
  if (rt->l) ct += count_total_edges(rt->l);
  if (rt->r) ct += count_total_edges(rt->r);

  ct ++;
  return ct;
}


void insert_edge(EDGE_TREE_t **rt_ptr, EDGE_t *ed) {
  EDGE_TREE_t *new, *x, *y;

  new = (EDGE_TREE_t *) malloc(sizeof(EDGE_TREE_t));
  new->to_edge = ed;
  new->p = NULL;
  new->l = NULL;
  new->r = NULL;

  y = NULL;
  x = *rt_ptr;
  while (x) {
    y = x;
    if (new->to_edge->index < x->to_edge->index) x = x->l;
    else x = x->r;
  }
  new->p = y;
  if (!y) *rt_ptr = new;
  else {
    if (new->to_edge->index < y->to_edge->index) y->l = new;
    else y->r = new;
  }
}


void delete_edge(EDGE_TREE_t **rt_ptr, EDGE_t *ed) {
  EDGE_TREE_t *x, *y, *z;

  if (!(*rt_ptr)) {
    err_no ++;
    fprintf(log_file, "error:  edge tree empty: %d %d %d\n",  ed->index, ed->ele1_info->index, ed->ele2_info->index);
    fflush(log_file);
    exit(4);
  }
  z= find_edge(*rt_ptr, ed->index);
  if (z) {
    if (!z->l || !z->r) y = z;
    else y = minimal_edge(z->r);
    if (y->l) x = y->l;
    else x = y->r;
    if (x) x->p = y->p;
    if (!y->p) *rt_ptr = x;
    else {
      if (y == y->p->l) y->p->l = x;
      else y->p->r = x;
    }
    if (y != z) z->to_edge = y->to_edge;
    free(y);
  } else {
    err_no ++;
    fprintf(log_file, "error:  Can not find edge %d between ele %d and %d\n", ed->index, ed->ele1_info->index, ed->ele2_info->index);
    fflush(log_file);
    exit(4);
  }
}


EDGE_TREE_t *minimal_edge(EDGE_TREE_t *rt) {
  while (rt->l) {
    rt = rt->l;
  }

  return rt;
}


EDGE_TREE_t *find_edge(EDGE_TREE_t *rt, int index) {
  while (rt && index != rt->to_edge->index) {
    if (index < rt->to_edge->index) rt = rt->l;
    else rt = rt->r;
  }

  return rt;
}
