/*
 * redef_boundary.c  --  Potential-cut-point and boundary-detection functions
 *
 * This module implements the CP/BD clustering pipeline used during Stage 3
 * (eleredef) to identify positions where an element should be split into two
 * or more child elements.
 *
 * Terminology:
 *   PCP  (Potential Cut Point)  -- an element-boundary candidate derived from
 *        the endpoints of full-length MSP images.  Stored as a CP_t list on
 *        ELEMENT_t.PCP during edges_and_cps().
 *
 *   BD   (Boundary / cluster)   -- a PCP cluster produced by CP_cluster().
 *        Represented as a BD_t list; each node carries the cluster centroid
 *        and a support count.
 *
 *   TBD  (To-Be-Determined)     -- a BD that has passed the span() threshold
 *        and is promoted to ELEMENT_t.TBD for use by dissect().
 *
 * Public interface: see redef_boundary.h
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#include "ele.h"
#include "redef_boundary.h"


/* =========================================================================
 * Comparison functions (used by qsort callers in this and other units)
 * ========================================================================= */

int CP_cmp(const void *cp1, const void *cp2) {
  return (*((CP_t **) cp1))->cp - (*((CP_t **) cp2))->cp;
}

int BD_cmp(const void *bd1, const void *bd2) {
  return (*((BD_t **) bd1))->bd - (*((BD_t **) bd2))->bd;
}


/* =========================================================================
 * CP / BD sort
 * ========================================================================= */

void CP_sort(CP_t **CP_ptr) {
  int i=0, j;
  CP_t *cur, **CPs;

  cur = *CP_ptr;
  while (cur != NULL) {
    i ++;
    cur = cur->next;
  }
  CPs = (CP_t **) malloc(i*sizeof(CP_t *));
  i = 0;
  cur = *CP_ptr;
  while (cur != NULL) {
    *(CPs+i) = cur;
    i ++;
    cur = cur->next;
  }
  qsort(CPs, i, sizeof(CP_t *), CP_cmp);
  i --;
  for (j=0; j<i; j++) {
    (*(CPs+j))->next = *(CPs+j+1);
  }
  (*(CPs+i))->next = NULL;
  *CP_ptr = *CPs;

  free(CPs);
}


void BD_sort(BD_t **BD_ptr) {
  int i=0, j;
  BD_t *cur, **BDs;

  cur = *BD_ptr;
  while (cur != NULL) {
    i ++;
    cur = cur->next;
  }

  BDs = (BD_t **) malloc(i*sizeof(BD_t*));
  i = 0;
  cur = *BD_ptr;;
  while (cur != NULL) {
    *(BDs+i) = cur;
    i ++;
    cur = cur->next;
  }

  qsort(BDs, i, sizeof(BD_t *), BD_cmp);

  i --;
  for (j=0; j<i; j++) {
    (*(BDs+j))->next = *(BDs+j+1);
  }
  (*(BDs+i))->next = NULL;;
  *BD_ptr = *BDs;

  free(BDs);
}


/* =========================================================================
 * CP_cluster  --  group sorted CPs into BD boundary candidates
 *
 * CP_t is a linked list holding an integer coordinate (cp) and a pointer
 * to the contributing element.  Consecutive CPs within 20 bp of the cluster
 * start and 10 bp of the running last are merged into one BD_t node.
 *
 * Example CP values (endpoints of full-length images):
 *   cp=1,    contributor ele 2
 *   cp=1597, contributor ele 2
 *   cp=2,    contributor ele 1
 *   cp=1600, contributor ele 1
 *   ...
 * ========================================================================= */
BD_t *CP_cluster(CP_t *cps) {
  int32_t first = cps->cp, last = cps->cp, sum = 0;
  CP_t *begin = cps;
  int cpct = 0;
  BD_t *bds = NULL, *bd_tmp;

  while (cps != NULL) {
    if (cps->cp - first <= 20 && cps->cp - last <= 10) {
      sum += cps->cp - first;
      cpct ++;
      last = cps->cp;
      cps = cps->next;
    } else {
      bd_tmp = (BD_t *) malloc(sizeof(BD_t));
      bd_tmp->bd = sum/cpct + first;
      bd_tmp->support = cpct;
      bd_tmp->next = bds;
      bds = bd_tmp;
      if (cps->cp - last <= 10)	begin = cps;
      else begin = begin->next;
      if (begin) {
	first = begin->cp;
	last = first;
      }
      sum = 0;
      cpct = 0;
      cps = begin;
    }
  }

  bd_tmp = (BD_t *) malloc(sizeof(BD_t));
  bd_tmp->bd = sum/cpct + first;
  bd_tmp->support = cpct;
  bd_tmp->next = bds;
  bds = bd_tmp;

  return bds;
}


/* =========================================================================
 * span  --  count images that span a candidate cut point
 *
 * Returns the number of images that start before (cut-10) but end before
 * (cut+10), multiplied by FUDGE.  Used by PCP_to_TBDs() as the minimum
 * support threshold for promoting a BD to a TBD.
 *
 * Requires:  PCP sorted by CP_cmp; images sorted by frag_cmp.
 * ========================================================================= */
int span(ELEMENT_t *ele, int32_t cut) {
  int left=0, rite=0;
  IMG_DATA_t *id;

  id = ele->to_img_data;
  while (id && id->to_image->frag.lb <= cut-10) {
    left ++;
    if (id->to_image->frag.rb <= cut+10) rite ++;
    id = id->next;
  }
  return (left-rite)*FUDGE;
}


/* =========================================================================
 * TBD_merge  --  collapse TBD nodes that are within 10 bp of each other
 *
 * After BD_sort(), adjacent boundaries closer than 10 bp are merged by
 * keeping the one with greater support.
 * ========================================================================= */
void TBD_merge(ELEMENT_t *ele) {
  BD_t *prev=NULL, *cur, *next;

  BD_sort(&ele->TBD);
  cur = ele->TBD;
  next = cur->next;
  while (next != NULL) {
    if (next->bd - cur->bd <= 10) {
      if (cur->support < next->support) {
	if (prev == NULL) {
	  ele->TBD = next;
	} else {
	  prev->next = next;
	}
	free(cur);
	cur = next;
      } else {
	cur->next = next->next;
	free(next);
      }
      if (cur) next = cur->next;
      else next = NULL;
    } else {
      prev = cur;
      cur = next;
      next = cur->next;
    }
  }
}


/* =========================================================================
 * PCP_to_TBDs  --  promote well-supported boundary candidates to TBDs
 *
 * 1. Sorts the element's PCP list by coordinate.
 * 2. Clusters them into BD nodes via CP_cluster().
 * 3. Promotes each BD whose support >= span() to ELEMENT_t.TBD.
 * ========================================================================= */
void PCP_to_TBDs(ELEMENT_t *ele) {
  RLOG_DBG("PCP_to_TBDs: ele %d, first PCP contributor ele %d\n",
           ele->index, ele->PCP->contributor->index);
  int s = 0, left;
  BD_t *pbd_tmp, *pbd_prev, *pbd, *pbds;
  CP_t *cp;

  /* sort the PCP list according to cp */
  CP_sort(&ele->PCP);
  /* cluster the PCPs into PBDs */
  pbds = CP_cluster(ele->PCP);
  /* identify TBDs from PBDs:
   * TBDs are removed from pbds; what remains are unsuccessful ones */
  pbd_tmp = pbds;
  pbd_prev = NULL;
  while (pbd_tmp != NULL) {
    /* s is the KEY: minimum support required */
    s = span(ele, pbd_tmp->bd);
    if (pbd_tmp->support >= s) {
	if (pbd_prev == NULL) {
	  pbds = pbd_tmp->next;
	  pbd_tmp->next = ele->TBD;
	  ele->TBD = pbd_tmp;
	  pbd_tmp = pbds;
	} else {
	  pbd_prev->next = pbd_tmp->next;
	  pbd_tmp->next = ele->TBD;
	  ele->TBD = pbd_tmp;
	  pbd_tmp = pbd_prev->next;
	}
    } else {
      pbd_prev = pbd_tmp;
      pbd_tmp = pbd_tmp->next;
    }
  }

  BD_free(&pbds);
}


/* =========================================================================
 * add_CP  --  prepend a new cut-point to a CP list
 * ========================================================================= */
void add_CP(CP_t **CP_ptr, int32_t cp, ELE_INFO_t *cont) {
  CP_t *new = (CP_t *) malloc(sizeof(CP_t));
  new->cp = cp;
  new->contributor = cont;
  new->next = *CP_ptr;
  *CP_ptr = new;
}


/* =========================================================================
 * CP_clean  --  remove all CPs contributed by a specific element
 *
 * Called from combo_edge_update() when an element is dismissed, to purge
 * its contribution from its partners' PCP lists.
 * ========================================================================= */
void CP_clean(CP_t **cps_p, ELE_INFO_t *cont) {
  CP_t *cp_cur, *cp_prev;

  cp_cur = *cps_p;
  cp_prev = NULL;
  while (cp_cur != NULL) {
    if (cp_cur->contributor == cont) {
      if (cp_prev == NULL) {
	*cps_p = cp_cur->next;
	free(cp_cur);
	cp_cur = *cps_p;
      } else {
	cp_prev->next = cp_cur->next;
	free(cp_cur);
	cp_cur = cp_prev->next;
      }
    } else {
      cp_prev = cp_cur;
      cp_cur = cp_cur->next;
    }
  }
}
