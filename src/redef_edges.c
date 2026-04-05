/*
 * redef_edges.c  --  Edge-building and image-consistency functions
 *
 * This module implements two related subsystems used during Stage 3
 * (eleredef):
 *
 * 1. Edge graph construction (edges_and_cps)
 *    For each element, examines its non-self images, determines which
 *    are full-length (primary) or partial, and builds:
 *      - edges between element pairs (primary 'p' or secondary 's')
 *      - PCP lists on each element from the endpoints of primary images
 *
 *    The element graph:
 *      E1 ---primary---> E2   when a nearly full-length alignment exists
 *                              (within 10 bp of both endpoints) for E1 or E2.
 *      E1 --secondary--> E2   when only partial alignments are found.
 *
 * 2. Consistency tree (consis_tree_build / find_prim)
 *    Partial images that could collectively tile a full-length alignment
 *    are detected via a tree structure where each level holds non-overlapping
 *    images.  find_prim() traverses this tree to identify paths whose
 *    combined coverage qualifies as a primary connection.
 *
 * Public interface: see redef_edges.h
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#include "ele.h"
#include "eleredef.h"
#include "redef_boundary.h"
#include "redef_edges.h"


/* =========================================================================
 * partner_cmp  --  sort IMAGE_t* by partner element index, direction, lb, rb
 *
 * Used by edges_and_cps() to group images by their partner element before
 * the edge/CP classification pass.
 * ========================================================================= */
int partner_cmp(const void *i1, const void *i2) {
  int res = partner(*((IMAGE_t **)i1))->ele_info->index -
            partner(*((IMAGE_t **)i2))->ele_info->index;
  if (res == 0) {
    res = (*((IMAGE_t **)i1))->to_msp->direction -
          (*((IMAGE_t **)i2))->to_msp->direction;
  }
  if (res == 0) {
    res = (*((IMAGE_t **)i1))->frag.lb - (*((IMAGE_t **)i2))->frag.lb;
  }
  if (res == 0) {
    res = (*((IMAGE_t **)i1))->frag.rb - (*((IMAGE_t **)i2))->frag.rb;
  }
  return res;
}


/* =========================================================================
 * full_length  --  test whether an image spans nearly all of its element
 *
 * Returns 1 if the image starts within 10 bp of the element's left bound,
 * ends within 10 bp of the element's right bound, and covers at least
 * 'cutoff' fraction of the element's length.
 * ========================================================================= */
int full_length(IMAGE_t *i, float cutoff) {
  if (!i->ele_info->ele) {
    err_no ++;
    fprintf(log_file, "error:  element %d not in memory\n", i->ele_info->index);
    fflush(log_file);
    /* RMH: Looks like this could have recovered here, but exited instead. */
    exit(3);
    ele_read_in(i->ele_info, 1);
  }
  /* For cutoff=0.9 this is only more stringent for sequences < 200 bp. */
  if ( i->frag.lb - i->ele_info->ele->frag.lb < 10 &&
       i->frag.rb - i->ele_info->ele->frag.rb > -10 &&
       ((float)i->frag.rb - i->frag.lb) /
           (i->ele_info->ele->frag.rb - i->ele_info->ele->frag.lb) > cutoff) {
    return 1;
  }
  return 0;
}


/* =========================================================================
 * add_edge  --  create an edge between two elements and insert into both BSTs
 * ========================================================================= */
void add_edge(ELE_INFO_t *ele1_info, ELE_INFO_t *ele2_info,
              char type, int32_t score, short dir) {
  EDGE_t *new = EDGE_malloc();

  edge_index ++;
  new->index     = edge_index;
  new->ele1_info = ele1_info;
  new->ele2_info = ele2_info;
  new->type      = type;
  new->score     = score;
  new->direction = dir;

  ele1_info->ele->edge_no ++;
  insert_edge(&ele1_info->ele->edges, new);

  if (ele1_info->index != ele2_info->index) {
    ele2_info->ele->edge_no ++;
    insert_edge(&ele2_info->ele->edges, new);
  }
}


/* =========================================================================
 * adjust_edge_tree  --  rebuild the edge BST from a flat array
 *
 * The edge BST can become unbalanced after many insertions.  This function
 * serialises it to an array, frees the old tree, and rebuilds it balanced.
 * ========================================================================= */
void adjust_edge_tree(ELE_INFO_t *ele_info) {
  EDGE_t **edge_array;
  int ct, ct1;

  ct = ele_info->ele->edge_no;
  edge_array = (EDGE_t **) malloc(ct*sizeof(EDGE_t *));
  ct1 = charge_edge_array(edge_array, ele_info->ele->edges, 0);
  if (ct1 != ct) {
    err_no ++;
    fprintf(log_file,
            "error:  trouble charging the edge array in ele %d: "
            "%d charged, %d expected\n",
            ele_info->index, ct1, ct);
    fflush(log_file);
    exit(4);
  }
  edge_tree_free(&ele_info->ele->edges);
  build_edge_tree(&ele_info->ele->edges, edge_array, 0, ct-1);
  free(edge_array);
}


/* =========================================================================
 * charge_edge_array  --  in-order serialisation of the edge BST
 * ========================================================================= */
int charge_edge_array(EDGE_t **edge_array, EDGE_TREE_t *rt, int pos) {
  if (rt->l) pos = charge_edge_array(edge_array, rt->l, pos);

  *(edge_array+pos) = rt->to_edge;
  pos ++;

  if (rt->r) pos = charge_edge_array(edge_array, rt->r, pos);

  return pos;
}


/* =========================================================================
 * edges_and_cps  --  core function: build edges and CP lists for one element
 *
 * For each partner element connected to ele_info by at least one image,
 * determines whether the connection is primary or secondary and records:
 *   - an edge in both elements' edge BSTs
 *   - CP entries on both elements for the endpoints of primary images
 *
 * When no single full-length image exists, the consistency tree is used to
 * check whether a chain of partial images could collectively serve as a
 * full-length connection (find_prim).
 *
 * On exit, ele_info->stat is set to 't'.
 * ========================================================================= */
void edges_and_cps(ELE_INFO_t *ele_info, IMAGE_t **img_ptr) {
    int eff_img_ct, i, prim;
    int j, riteplace;
    short ritetime;
    MSP_t *prim_p;

    IMG_DATA_t *cur_img_data;
    ELEMENT_t *ele_partner, *cur_ele = ele_info->ele;
    ELE_INFO_t *epi;
    IMAGE_t *img_partner, *cur_img;

    IMAGE_t *token_image;
    int token_mark=0;
    IMG_NODE_t *consis_rt, *cur_consis_nd;
    ELE_DATA_t *ele_data_tmp;
    IMG_DATA_t *scovs;

    int32_t max_score=0;
    short dir;

    RLOG_DBG("edges_and_cps: ele_info->index = %d, ele_info->ele->index = %d\n",
             ele_info->index, ele_info->ele->index);

    /* Sort unprocessed images by partner element, then left bound.
     * When update=1 all images need processing; when update=0, images whose
     * partner is already 'v' or 'w' are skipped.  Self-images are always
     * excluded from CP/edge generation. */
    if (!cur_ele->to_img_data) listify(cur_ele->to_img_tree, &cur_ele->to_img_data);
    cur_img_data = cur_ele->to_img_data;

    /* First pass: count effective images (non-self, eligible) */
    /* TODO: This is redundant with the loop below.  Could be combined. */
    eff_img_ct = 0;
    while(cur_img_data != NULL) {
      epi = partner(cur_img_data->to_image)->ele_info;
      if (epi->index != ele_info->index) {
        if (cur_ele->update || epi->stat == 'z') eff_img_ct ++;
      }
      cur_img_data = cur_img_data->next;
    }
    RLOG_DBG("edges_and_cps: img_no = %d, eff_img_ct = %d\n",
             cur_ele->img_no, eff_img_ct);

  if (eff_img_ct) {
    token_image = (IMAGE_t *) malloc(sizeof(IMAGE_t));
    token_image->frag.lb = 0;
    token_image->frag.rb = 0;
    token_image->to_msp = NULL;
    token_image->ele_info = NULL;
    consis_rt = (IMG_NODE_t *) malloc(sizeof(IMG_NODE_t));
    consis_rt->to_image = token_image;
    consis_rt->children = NULL;
    consis_rt->sib = NULL;

    if (eff_img_ct > MAX_IMG) {
      img_ptr = (IMAGE_t **) malloc(eff_img_ct*sizeof(IMAGE_t *));
    }

    /* Second pass: populate img_ptr with eligible images */
    cur_img_data = cur_ele->to_img_data;
    eff_img_ct = 0;
    while(cur_img_data != NULL) {
      epi = partner(cur_img_data->to_image)->ele_info;
      if (epi->index != ele_info->index) {
        if (cur_ele->update || epi->stat == 'z') {
  	  *(img_ptr+eff_img_ct) = cur_img_data->to_image;
	  eff_img_ct ++;
        }
      }
      cur_img_data = cur_img_data->next;
    }
    /* Sort by: partner element index, MSP direction, left bound, right bound */
    qsort(img_ptr, eff_img_ct, sizeof(IMAGE_t *), partner_cmp);

    /* Recognise full-length images and build the consistency tree for
     * partial images, grouped by partner element. */
    ritetime = 1;
    for (i=0; i<eff_img_ct; i++) {
      cur_img = *(img_ptr+i);
      img_partner = partner(cur_img);

      RLOG_DBG("edges_and_cps: image cur_ele_info->index=%d, par_ele_info->index=%d\n",
               cur_img->ele_info->index, img_partner->ele_info->index);

      if (ritetime) { /* new partner element begins */
	epi = img_partner->ele_info;
	ritetime = 0;
	if (!epi->ele) ele_read_in(epi, 1);
	ele_partner = epi->ele;
	riteplace = i;
	prim = 0;
	prim_p = NULL;
	max_score = 0;
      }
      if (img_partner->ele_info->index == epi->index) {
	/* full length check for current element */
	if (full_length(cur_img, CUTOFF2)) {
          RLOG_DBG("edges_and_cps:     image is primary (full_length with current element)\n");
	  prim = 1;
	  cur_ele->flimg_no ++;
	}
	/* full length check for partner element */
	if (full_length(img_partner, CUTOFF2)) {
          RLOG_DBG("edges_and_cps:     image is primary (full_length with partner element)\n");
	  prim = 1;
	  ele_partner->flimg_no ++;
	}
	/* track the highest-scoring primary image */
	if (prim == 1) {
	  prim = 0;
	  if (cur_img->to_msp->iden > max_score) {
            RLOG_DBG("edges_and_cps:     **** new high-scoring primary image\n");
	    max_score = cur_img->to_msp->iden;
	    dir = cur_img->to_msp->direction;
	    prim_p = cur_img->to_msp;
	  }
	}
	/* accumulate partial images into the consistency tree */
	if (!prim_p) {
          RLOG_DBG("edges_and_cps:      Adding to consis tree: "
                   "e%d im=%s:%d-%d  e%d pt=%s:%d-%d\n",
                   epi->index,
                   cur_img->frag.seq_name, cur_img->frag.lb, cur_img->frag.rb,
                   img_partner->ele_info->index,
                   img_partner->frag.seq_name, img_partner->frag.lb,
                   img_partner->frag.rb);
	  consis_tree_build(consis_rt, cur_img, 1);
	}
      }
      if (img_partner->ele_info->index != epi->index || i == eff_img_ct-1) {
        RLOG_DBG("edges_and_cps:     End of current partner element\n");
	ritetime = 1;
	/* back up loop so this image is reprocessed with the next partner */
	if (img_partner->ele_info->index != epi->index) i --;
	/* finish the current partner */
	if (prim_p) {
	  prim_p->stat = 'p';
	  prim = 1;
	} else {
	  /* no full-length image; search for partial primaries in consis tree */
	  prim = find_prim(consis_rt->children, CUTOFF2,
                           ele_info->ele->frag.lb, -1, 0, 0, 0, 0, 0,
                           &token_mark, &max_score, &dir);
          if (prim)
            RLOG_DBG("edges_and_cps:     Identified partial primary image\n");
	}
	/* build edge */
	if (prim) {
          RLOG_DBG("edges_and_cps:     Adding primary edge\n");
	  if (ele_info->index != epi->index) {
	    add_edge(ele_info, epi, 'p', max_score, dir);
	  } else {
	    err_no ++;
	    fprintf(log_file, "error:  self edge seen: ele %d\n", ele_info->index);
	    fflush(log_file);
	  }
	  /* record CP endpoints from all primary images in this group */
	  for (j=riteplace; j<=i; j++) {
	    if ((*(img_ptr+j))->to_msp->stat == 'p') {
	      cur_img = *(img_ptr+j);
	      img_partner = partner(cur_img);
	      add_CP(&cur_ele->PCP, cur_img->frag.lb, epi);
	      add_CP(&cur_ele->PCP, cur_img->frag.rb, epi);
	      add_CP(&ele_partner->PCP, img_partner->frag.lb, ele_info);
	      add_CP(&ele_partner->PCP, img_partner->frag.rb, ele_info);
	    }
	  }
	} else {
          RLOG_DBG("edges_and_cps:     Adding secondary edge\n");
	  if (ele_info->index != epi->index) {
	    add_edge(ele_info, epi, 's', 0, 0);
	  } else {
	    err_no ++;
	    fprintf(log_file, "error:  self edge seen: ele %d\n", ele_info->index);
	    fflush(log_file);
	  }
	}
	if (consis_rt->children != NULL) {
          RLOG_DBG("edges_and_cps:     Freeing consis tree\n");
	  consis_tree_free(consis_rt->children);
	  consis_rt->children = NULL;
	}
      } else {
        RLOG_DBG("edges_and_cps:     image fell through all cases\n");
      }
    }

    if (eff_img_ct > MAX_IMG) {
      free(img_ptr);
    }
    free(consis_rt);
    free(token_image);
  }

  if (ele_info->ele->edges) adjust_edge_tree(ele_info);

  cur_ele->update = 0;
  ele_info->stat = 't';
}


/* =========================================================================
 * Consistency tree
 *
 * The tree organises partial images by their non-overlapping relationship:
 * each level holds images that are consistent (non-overlapping, same pair,
 * same direction) with the node at the level above.  find_prim() then
 * searches for a root-to-leaf path whose combined length qualifies as a
 * primary connection.
 * ========================================================================= */

/* print_consis_tree  --  debug dump of the consistency tree */
int print_consis_tree(IMG_NODE_t *rt) {
  int level = 0;
  IMG_NODE_t *child_rt;
  IMG_NODE_t *sib_rt;
  if ( rt->sib != NULL )
    printf("c_tree: level=%d [(null):0:0-0] ( has sibs )\n",level);
  else
    printf("c_tree: level=%d [(null):0:0-0]\n",level);
  level++;
  child_rt = rt->children;
  while ( child_rt != NULL ) {
    printf("c_tree: level=%d [%s:%d:%d-%d]",level,
           child_rt->to_image->frag.seq_name,
           child_rt->to_image->to_msp->direction,
           child_rt->to_image->frag.lb,
           child_rt->to_image->frag.rb);
    sib_rt = child_rt->sib;
    while ( sib_rt != NULL ) {
      if ( sib_rt->children == NULL )
        printf(", image %s:%d:%d-%d",
               sib_rt->to_image->frag.seq_name,
               sib_rt->to_image->to_msp->direction,
               sib_rt->to_image->frag.lb,
               sib_rt->to_image->frag.rb);
      else
        printf(", image %s:%d:%d-%d(*)",
               sib_rt->to_image->frag.seq_name,
               sib_rt->to_image->to_msp->direction,
               sib_rt->to_image->frag.lb,
               sib_rt->to_image->frag.rb);
      sib_rt = sib_rt->sib;
    }
    printf("\n");
    level++;
    child_rt = child_rt->children;
  }
  return 0;
}


/* consis_tree_build  --  insert an image into the consistency tree
 *
 * Non-overlapping images (consistent with the parent node) are placed on
 * a new child level; overlapping ones are placed as siblings on the same
 * level.  prequal=1 on the first call (root) to bypass the overlap check.
 */
int consis_tree_build(IMG_NODE_t *rt, IMAGE_t *im, int prequal) {
  int sum=0;
  IMG_NODE_t *nex_rt, *node;
  if (prequal || consis(rt->to_image, im, CUTOFF2)) {
    nex_rt = rt->children;
    while (nex_rt != NULL) {
      sum += consis_tree_build(nex_rt, im, 0);
      nex_rt = nex_rt->sib;
    }
    if (!sum) {
      node = (IMG_NODE_t *) malloc(sizeof(IMG_NODE_t));
      node->recorded = 0;
      node->to_image = im;
      node->sib = NULL;
      node->children = NULL;
      *node_entry(&rt->children) = node;
      sum = 1;
    }
  }
  return sum;
}


/* consis  --  test whether two images are consistent (non-overlapping, co-linear)
 *
 * Returns 1 if:
 *   - Both images belong to the same element pair.
 *   - Both images have the same MSP direction.
 *   - Neither image nor its partner overlaps the other by >= (1-cutoff) of
 *     either sequence's length.
 */
int consis(IMAGE_t *i1, IMAGE_t *i2, float cutoff) {
  int res = 0;
  IMAGE_t *ip1 = partner(i1), *ip2 = partner(i2);
  if (i1->ele_info->index == i2->ele_info->index &&
      ip1->ele_info->index == ip2->ele_info->index &&
      i1->to_msp->direction == i2->to_msp->direction) {
    if (i1->to_msp->direction == 1) {
      /* RMH: If either image or its partner starts at the same position
       * they probably overlap too much.  This is an optimisation that
       * avoids calling sing_cov.  Assumes sequences are the same when
       * left bounds collide -- acceptable in practice but worth noting. */
      if ((i1->frag.lb - i2->frag.lb)*(ip1->frag.lb - ip2->frag.lb) > 0) {
	if (!sing_cov(&i1->frag, &i2->frag, 1.0-cutoff) &&
            !sing_cov(&ip1->frag, &ip2->frag, 1.0-cutoff)) {
	  res = 1;
	}
      }
    } else {
      if ((i1->frag.lb - i2->frag.lb)*(ip1->frag.rb - ip2->frag.rb) < 0) {
	if (!sing_cov(&i1->frag, &i2->frag, 1.0-cutoff) &&
            !sing_cov(&ip1->frag, &ip2->frag, 1.0-cutoff)) {
	  res = 1;
	}
      }
    }
  }
  return res;
}


/* node_entry  --  find the last NULL sib pointer in a sibling chain */
IMG_NODE_t **node_entry(IMG_NODE_t **node_pp) {
  if (*node_pp != NULL) {
    return node_entry(&(*node_pp)->sib);
  }
  return node_pp;
}


/* consis_tree_free  --  recursively free a consistency (sub)tree */
void consis_tree_free(IMG_NODE_t *rt) {
  if (rt->sib != NULL) {
    consis_tree_free(rt->sib);
  }
  if (rt->children != NULL) {
    consis_tree_free(rt->children);
  }
  free(rt);
}


/* find_prim  --  search the consistency tree for a primary-image path
 *
 * Traverses the tree recursively.  At each leaf, checks whether the
 * accumulated alignment (al1/al2) covers enough of the inferred full-length
 * extent (efl1/efl2) to qualify as a primary connection.
 *
 * Initially called with:
 *   end1 = ele->frag.lb, end2 = -1
 *   efl1 = efl2 = al1 = al2 = score = 0
 *
 * Note: int32_t may overflow for al1/al2 on very long elements; this is a
 * known limitation of the original algorithm.
 */
int find_prim(IMG_NODE_t *nd, float cutoff,
              int32_t end1, int32_t end2,
              int32_t efl1, int32_t efl2,
              int32_t al1,  int32_t al2,
              int32_t score,
              int *pmarkp, int32_t *sc, short *d) {
  int sum = 0, mark=0;
  int32_t skip1, skip2, len1, len2;
  IMAGE_t *ipt;

  if (nd->sib)
    sum += find_prim(nd->sib, cutoff, end1, end2, efl1, efl2,
                     al1, al2, score, pmarkp, sc, d);

  ipt = partner(nd->to_image);
  /* end2 is the partner element's boundary for the first image in the group */
  if (end2 < 0) {
    if (nd->to_image->to_msp->direction == 1) end2 = ipt->ele_info->ele->frag.lb;
    else end2 = ipt->ele_info->ele->frag.rb;
  }

  /* skip1: gap between element left bound and image left bound (>0 = image starts late)
   * skip2: gap between partner's boundary and partner image's boundary */
  skip1 = nd->to_image->frag.lb - end1;
  if (nd->to_image->to_msp->direction == 1) skip2 = ipt->frag.lb - end2;
  else skip2 = end2 - ipt->frag.rb;

  /* If both image and partner start with a gap > 10 bp, accumulate them */
  if (skip1>10 && skip2>10) {
    efl1 += skip1;
    efl2 += skip2;
  }
  len1 = nd->to_image->frag.rb - nd->to_image->frag.lb;
  len2 = ipt->frag.rb - ipt->frag.lb;
  efl1 += len1;
  efl2 += len2;
  al1  += len1;
  al2  += len2;
  score += ((int32_t) nd->to_image->to_msp->iden)*(len1+len2)/2;

  if (nd->children) {
    end1 = nd->to_image->frag.rb;
    if (nd->to_image->to_msp->direction == 1) end2 = ipt->frag.rb;
    else end2 = ipt->frag.lb;
    sum += find_prim(nd->children, cutoff, end1, end2, efl1, efl2,
                     al1, al2, score, &mark, sc, d);
  } else {
    /* last image in group: check trailing gaps */
    skip1 = nd->to_image->ele_info->ele->frag.rb - nd->to_image->frag.rb;
    if (nd->to_image->to_msp->direction == 1)
      skip2 = ipt->ele_info->ele->frag.rb - ipt->frag.rb;
    else
      skip2 = ipt->frag.lb - ipt->ele_info->ele->frag.lb;
    if (skip1>10 && skip2>10) {
      efl1 += skip1;
      efl2 += skip2;
    }

    RLOG_DBG("find_prim: al1=%d al2=%d efl1=%d efl2=%d:  "
             "al/ef %f, %f   ef-al %d, %d  ele:%d-%d, ptn:%d-%d\n",
             al1, al2, efl1, efl2,
             (1.0*al1/efl1), (1.0*al2/efl2),
             (efl1-al1), (efl2-al2),
             nd->to_image->ele_info->ele->frag.lb,
             nd->to_image->ele_info->ele->frag.rb,
             ipt->ele_info->ele->frag.lb,
             ipt->ele_info->ele->frag.rb);

    if ( (1.0*al1/efl1 > cutoff || 1.0*al2/efl2 > cutoff) &&
         (efl1-al1 < 30 || efl2-al2 < 30) ) {
      sum = 1;
      mark = 1;
      if ( (al1+al2) == 0 ) {
        /* RMH: Divide-by-zero guard.  Root cause not yet identified;
         * setting al1=1 avoids the crash without affecting the primary
         * classification result. */
        printf("eleredef warning: Divide by zero averted -- setting al1 to 1.\n");
        al1 = 1;
      }
      score = score*2/(al1+al2);
      if (score > *sc) {
	*sc = score;
	*d = nd->to_image->to_msp->direction;
      }
    }
  }
  if (mark) {
      nd->to_image->to_msp->stat = 'p';
      *pmarkp = 1;
  }
  return sum;
}
