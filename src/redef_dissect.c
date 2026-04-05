/*
 * redef_dissect.c  --  Element construction, redefinition, and dissection
 *
 * This module implements four related areas used during Stage 3 (eleredef):
 *
 * 1. Element construction (ele_def, new_element, generate_img_tree)
 *    Partitions a sorted list of images into elements using single-coverage
 *    clustering (sing_cov).  Creates ELE_INFO_t / ELEMENT_t structures and
 *    builds image BSTs for each resulting element.
 *
 * 2. Element redefinition (ele_redef)
 *    Drives the split of an existing element at its TBD positions by calling
 *    PCP_to_TBDs(), TBD_merge(), dissect(), and ele_def() on the resulting
 *    image subsets.  Recursively redefines child elements.
 *
 * 3. Image dissection (dissect, register_image, dump_image)
 *    Splits MSPs that span a TBD boundary into left and right halves,
 *    creating new MSP records and updating both elements' image trees.
 *
 * 4. Combo/obs bookkeeping (combo_update, combo_output, obs_output)
 *    Handles the cleanup of elements that are subsumed (combo) or
 *    eliminated (obsolete) as a result of dissection.
 *
 * Public interface: see redef_dissect.h
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#include <time.h>

#include "ele.h"
#include "eleredef.h"      /* IMG_NODE_t -- needed by redef_edges.h */
#include "redef_boundary.h"
#include "redef_edges.h"
#include "redef_dissect.h"


/* =========================================================================
 * Comparison functions for qsort
 * ========================================================================= */

/* frag_cmp  --  sort IMG_DATA_t by sequence name pointer, then lb, then rb */
int frag_cmp(const void *i1, const void *i2) {
  int res = (*(IMG_DATA_t **)i1)->to_image->frag.seq_name -
            (*(IMG_DATA_t **)i2)->to_image->frag.seq_name;
  if (res == 0) {
    res = (*(IMG_DATA_t **)i1)->to_image->frag.lb -
          (*(IMG_DATA_t **)i2)->to_image->frag.lb;
    if (res == 0) {
      res = (*(IMG_DATA_t **)i1)->to_image->frag.rb -
            (*(IMG_DATA_t **)i2)->to_image->frag.rb;
    }
  }
  return res;
}


int fam_cmp(const void *fd1, const void *fd2) {
  return (*((FAM_DATA_t **) fd1))->to_family->index -
         (*((FAM_DATA_t **) fd2))->to_family->index;
}


/* =========================================================================
 * img_data_sort  --  sort a IMG_DATA_t linked list via a temporary array
 * ========================================================================= */
IMG_DATA_t *img_data_sort(IMG_DATA_t *img_data, int ct) {
  IMG_DATA_t **img_data_ptr;
  int i;

  if (!ct) return NULL;
  img_data_ptr = (IMG_DATA_t **) malloc(ct*sizeof(IMG_DATA_t *));
  for (i=0; i<ct; i++) {
    *(img_data_ptr+i) = img_data;
    img_data = img_data->next;
  }
  qsort(img_data_ptr, ct, sizeof(IMAGE_t *), frag_cmp);
  img_data = *img_data_ptr;
  ct --;
  for (i=0; i<ct; i++) {
    (*(img_data_ptr+i))->next = *(img_data_ptr+i+1);
  }
  (*(img_data_ptr+ct))->next = NULL;
  free(img_data_ptr);
  return img_data;
}


/* =========================================================================
 * Element construction
 * ========================================================================= */

/* new_element  --  allocate an ELE_INFO_t / ELEMENT_t pair for a new element
 *
 * The initial storage for elements is the pre-allocated all_ele[] array.
 * Once that is full (ele_ct > ele_array_size) new elements are allocated
 * individually and appended to the ele_info_data overflow linked list.
 */
ELE_INFO_t *new_element() {
  ELE_INFO_t *ele_info_tmp;

  ele_ct ++;

  if (ele_ct <= ele_array_size) {
    ele_info_tmp = *(all_ele+ele_ct-1);
  } else {
    ele_info_tmp = ele_info_init(ele_ct);
    add_ele_info(ele_info_tmp);
  }

  ele_info_tmp->ele = ele_init(ele_ct);
  ele_info_tmp->file_updated = 1;

  return ele_info_tmp;
}


/* add_ele_info  --  append an ELE_INFO_t to the overflow linked list */
void add_ele_info(ELE_INFO_t *ele_info_tmp) {
  if (!ele_info_data) ele_info_data = ele_info_tmp;
  else ele_info_tail->next = ele_info_tmp;
  ele_info_tail = ele_info_tmp;
}


/* generate_img_tree  --  build the image BST for an element from its img_data list */
void generate_img_tree(ELEMENT_t *ele) {
  RLOG_DBG("generate_img_tree: ele %d\n", ele->index);
  IMAGE_t **img_array = (IMAGE_t **) malloc(ele->img_no*sizeof(IMAGE_t *));
  IMG_DATA_t *cur;
  int ct=0;

  cur = ele->to_img_data;
  while (cur) {
    *(img_array+ct) = cur->to_image;
    ct ++;
    cur = cur->next;
  }

  qsort(img_array, ct, sizeof(IMAGE_t *), img_index_cmp);
  build_img_tree(&ele->to_img_tree, img_array, 0, ct-1);
  free(img_array);
}


/* ele_def  --  partition a sorted image list into elements (single-coverage)
 *
 * Assumes the images are sorted by frag_cmp.  Iterates through the remaining
 * list and greedily assigns each image to the current element if sing_cov()
 * confirms it falls within the element's current extent.  When a gap is
 * detected, a new element is started.
 *
 * Returns a linked list of ELE_DATA_t nodes, one per new element.
 */
ELE_DATA_t *ele_def(IMG_DATA_t **img_data_p, float cutoff) {
  clock_t t;
  t = clock();

  /* cur_img_data  -- image under inspection
   * *img_data_p   -- first image in the remaining (un-partitioned) list
   * prev_img_data -- last unpartitioned image before cur_img_data */
  IMG_DATA_t *cur_img_data, *prev_img_data;
  ELEMENT_t *ele_tmp;
  ELE_INFO_t *ele_info_tmp;
  ELE_DATA_t *ele_data=NULL, *ele_data_tmp;
  short ritetime;
  int ct;

 while (*img_data_p) {

  cur_img_data = *img_data_p;
  prev_img_data = NULL;
  ritetime = 1; /* ritetime marks the right time to start a new element */

  while(cur_img_data != NULL) {
    if (ritetime) { /* starting a new element */
      ele_info_tmp = new_element();
      ele_tmp = ele_info_tmp->ele;
      ele_tmp->img_no = 1;
      ele_tmp->frag = cur_img_data->to_image->frag;
      /* kick cur_img_data out of the original list */
      *img_data_p = cur_img_data->next;
      cur_img_data->to_image->ele_info = ele_info_tmp;
      cur_img_data->next = NULL;
      ele_tmp->to_img_data = cur_img_data;
      cur_img_data = *img_data_p;
      ritetime = 0;
      continue;
    }
    /* check images on the same sequence that overlap the current element */
    if (ele_tmp->frag.seq_name == cur_img_data->to_image->frag.seq_name &&
        ele_tmp->frag.rb - cur_img_data->to_image->frag.lb > 10) {
      if (sing_cov(&ele_tmp->frag, &cur_img_data->to_image->frag, cutoff)) {
	ele_tmp->img_no ++;
	if (prev_img_data != NULL) {
	  prev_img_data->next = cur_img_data->next;
	} else {
	  *img_data_p = cur_img_data->next;
	}
	cur_img_data->to_image->ele_info = ele_info_tmp;
	cur_img_data->next = ele_tmp->to_img_data;
	ele_tmp->to_img_data = cur_img_data;

	if (ele_tmp->frag.rb < cur_img_data->to_image->frag.rb) {
	  ele_tmp->frag.rb = cur_img_data->to_image->frag.rb;
	  /* rescan: some previously rejected images may now qualify */
	  cur_img_data = *img_data_p;
	  prev_img_data = NULL;
	} else {
	  if (prev_img_data != NULL) cur_img_data = prev_img_data->next;
	  else cur_img_data = *img_data_p;
	}
      } else {
	prev_img_data = cur_img_data;
	cur_img_data = cur_img_data->next;
      }
    } else {
      /* gap detected: finalise current element */
      ele_data_tmp = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
      ele_data_tmp->ele_info = ele_info_tmp;
      ele_data_tmp->next = ele_data;
      ele_data = ele_data_tmp;
      ritetime = 1;
      cur_img_data = *img_data_p;
      prev_img_data = NULL;
    }
  }
  /* add the last element */
  ele_data_tmp = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
  ele_data_tmp->ele_info = ele_info_tmp;
  ele_data_tmp->next = ele_data;
  ele_data = ele_data_tmp;

 }

  ele_data_tmp = ele_data;
  while (ele_data_tmp) {
    if (ele_data_tmp->ele_info->ele->frag.lb > ele_data_tmp->ele_info->ele->frag.rb) {
      fprintf(log_file, "error:  ele %d reversed after ele_def\n",
              ele_data_tmp->ele_info->index);
      fflush(log_file);
      exit(3);
    }
    generate_img_tree(ele_data_tmp->ele_info->ele);
    {
      clock_t u = clock() - t;
      double ele_defTIME = ((double)u)/CLOCKS_PER_SEC;
      (void)ele_defTIME; /* timing available; suppress unused-variable warning */
    }
    ele_data_tmp = ele_data_tmp->next;
  }

  return ele_data;
}


/* =========================================================================
 * Element redefinition
 * ========================================================================= */

/* ele_redef  --  redefine an element at its TBD cut points
 *
 * 1. Sorts images by position.
 * 2. Converts PCPs to TBDs and merges adjacent ones.
 * 3. If dissection is needed, splits images at TBD positions and
 *    re-clusters into child elements via ele_def().
 * 4. Runs edges_and_cps() on each child and recurses.
 * 5. Calls combo_update() to release the parent element.
 */
void ele_redef(ELE_INFO_t *ele_info, IMAGE_t **img_ptr) {
    ELE_DATA_t *new_ele_data;
    ELEMENT_t *cur_ele = ele_info->ele;
    BD_t *pbd;
    CP_t *cur;
    short to_dissect=0;
    IMG_DATA_t *cur_img_data;
    if (!cur_ele->to_img_data) listify(cur_ele->to_img_tree, &cur_ele->to_img_data);
    cur_ele->to_img_data = img_data_sort(cur_ele->to_img_data, cur_ele->img_no);

    /* find TBDs for the element */
    if (cur_ele->PCP) {
      PCP_to_TBDs(cur_ele);
    } else {
      RLOG_DBG("ele_redef: ele %d has no PCP\n", cur_ele->index);
    }
    if (cur_ele->TBD) {
      TBD_merge(cur_ele);
    }
    if (cur_ele->TBD) {
      /* count TBDs that are interior (not within FLURRY bp of either boundary) */
      pbd = cur_ele->TBD;
      while (pbd) {
	if (pbd->bd-cur_ele->frag.lb>FLURRY && pbd->bd-cur_ele->frag.rb<-FLURRY) {
	  to_dissect ++;
	}
	pbd = pbd->next;
      }

      if (!to_dissect) {
	/* TBDs near the boundary: adjust the element extent rather than split */
	pbd = cur_ele->TBD;
	while (pbd) {
	  if (pbd->bd-cur_ele->frag.lb<=FLURRY) {
	    cur_ele->frag.lb = pbd->bd;
	  } else if (pbd->bd-cur_ele->frag.rb>=-FLURRY) {
	    cur_ele->frag.rb = pbd->bd;
	  }
	  pbd = pbd->next;
	}
	BD_free(&cur_ele->TBD);
      } else {
	/* dissect all images at TBD positions */
	dissect(ele_info);
	/* re-cluster remaining images into child elements */
	if (cur_ele->img_no) {
	  cur_ele->to_img_data = img_data_sort(cur_ele->to_img_data, cur_ele->img_no);
	  cur_ele->redef = ele_def(&cur_ele->to_img_data, CUTOFF1);
	}
	/* mark parent as dismissed, update CPs for affected partners */
	combo_update(ele_info);
	/* redefine child elements before pulling in new partners */
	new_ele_data = cur_ele->redef;
	while (new_ele_data != NULL) {
	  new_ele_data->ele_info->ele->update = 1;
	  new_ele_data->ele_info->ele->l_hold = cur_ele->l_hold;
	  if (new_ele_data->ele_info->ele->to_img_tree) {
	    edges_and_cps(new_ele_data->ele_info, img_ptr);
	    if (new_ele_data->ele_info->ele->PCP)
	      ele_redef(new_ele_data->ele_info, img_ptr);
	    else
	      new_ele_data->ele_info->stat = 'v';
	  } else {
	   if (new_ele_data->ele_info->ele->img_no) {
	    err_no ++;
	    fprintf(log_file,
	            "error:  image tree missing in newly redefed offspring %d\n",
	            new_ele_data->ele_info->index);
	    fflush(log_file);
	    exit(2);
	   } else {
	    combo_output(new_ele_data->ele_info);
	   }
	  }
	  new_ele_data = new_ele_data->next;
	}
      }
    }
    if (cur_ele->img_no > 0) {
      ele_info->stat = 'v';
    } else if (!to_dissect) combo_update(ele_info);
}


/* =========================================================================
 * Image dissection
 * ========================================================================= */

/* too_short  --  return 1 if the fragment is shorter than MIN_ELEMENT_LEN_BP */
int too_short(FRAG_t *f) {
  if (f->rb - f->lb <= TOO_SHORT) return 1;
  return 0;
}


/* add_msp  --  allocate a new MSP as a copy of an existing one */
MSP_t *add_msp(MSP_t *m) {
    MSP_t *msp_tmp;

    msp_tmp = MSP_malloc();
    *msp_tmp = *m;
    msp_index ++;
    msp_tmp->query.to_msp = msp_tmp;
    msp_tmp->query.index = 2*msp_index;
    msp_tmp->sbjct.to_msp = msp_tmp;
    msp_tmp->sbjct.index = 2*msp_index+1;
    msp_tmp->stat = 's';

    return msp_tmp;
}


/* register_image  --  insert a new image into both elements' image trees
 *
 * If the image (or its partner) is shorter than TOO_SHORT after dissection,
 * the MSP is freed instead.
 */
void register_image(IMAGE_t *i, ELEMENT_t *ele) {
  IMAGE_t *ip = partner(i);

  if (too_short(&i->frag)) {
    MSP_free(i->to_msp);
  } else {
      insert_image(&i->ele_info->ele->to_img_tree, i);
      i->ele_info->ele->img_no ++;
      if (i->ele_info->ele->to_img_data) put_image(&i->ele_info->ele->to_img_data, i);
      insert_image(&ip->ele_info->ele->to_img_tree, ip);
      ip->ele_info->ele->img_no ++;
      if (ip->ele_info->ele->to_img_data) put_image(&ip->ele_info->ele->to_img_data, ip);
  }
}


/* put_image  --  prepend an IMAGE_t to an img_data linked list */
void put_image(IMG_DATA_t **idp, IMAGE_t *i) {
  IMG_DATA_t *img_data_tmp;

  img_data_tmp = (IMG_DATA_t *)malloc(sizeof(IMG_DATA_t));
  img_data_tmp->to_image = i;
  img_data_tmp->next = *idp;
  *idp = img_data_tmp;
}


/* dump_image  --  remove an image from both elements' trees and free the MSP */
void dump_image(IMAGE_t *i) {
    IMAGE_t *ip = partner(i);

    if (i->ele_info->ele->to_img_data) remove_image(i);
    if (ip->ele_info->ele->to_img_data) remove_image(ip);

    delete_image(&i->ele_info->ele->to_img_tree, i);
    i->ele_info->ele->img_no --;
    delete_image(&ip->ele_info->ele->to_img_tree, ip);
    ip->ele_info->ele->img_no --;

    MSP_free(i->to_msp);
}


/* remove_image  --  unlink an IMAGE_t from the to_img_data list */
void remove_image(IMAGE_t *i) {
    IMG_DATA_t *prev_img_data=NULL, *cur_img_data;

    cur_img_data = i->ele_info->ele->to_img_data;
    while (cur_img_data != NULL) {
      if (cur_img_data->to_image == i) {
	if (prev_img_data != NULL) {
	  prev_img_data->next = cur_img_data->next;
	} else {
	  i->ele_info->ele->to_img_data = cur_img_data->next;
	}
	free(cur_img_data);
	break;
      }
      prev_img_data = cur_img_data;
      cur_img_data = cur_img_data->next;
    }
}


/* dissect  --  split images in an element at each TBD boundary
 *
 * For each image that spans a TBD position, creates a new MSP for the left
 * half, adjusts the original MSP to represent the right half, and calls
 * register_image() to insert both into the appropriate element trees.
 *
 * Images that become too short after splitting are discarded via dump_image().
 */
void dissect(ELE_INFO_t *ele_info) {
    clock_t r;
    r = clock();
  IMG_DATA_t *cur_img_data, *img_data_tmp, *next;
  MSP_t *msp_tmp, *msp_ori;
  IMAGE_t *img_partner, *target_img, *target_partner;
  ELEMENT_t *ele_partner, *cur_ele = ele_info->ele;
  BD_t *tbd_tmp, *new_tbd;
  short dissected;

    cur_img_data = cur_ele->to_img_data;
    while (cur_img_data != NULL) {
      next = cur_img_data->next;
      dissected = 0;
      img_partner = partner(cur_img_data->to_image);
      ele_partner = img_partner->ele_info->ele;
      if (full_length(img_partner, CUTOFF2)){
        ele_partner->flimg_no --;
      }
      tbd_tmp = cur_ele->TBD;
      while (tbd_tmp != NULL) {
	if (tbd_tmp->bd > cur_img_data->to_image->frag.lb || !tbd_tmp->next) {
	  if (tbd_tmp->bd > cur_img_data->to_image->frag.lb &&
	      tbd_tmp->bd < cur_img_data->to_image->frag.rb) {
	    if (tbd_tmp->bd - cur_img_data->to_image->frag.lb <= TOO_SHORT) {
	      /* TBD near left edge: trim rather than split */
	      cur_img_data->to_image->to_msp->score =
	        (int32_t) (cur_img_data->to_image->frag.rb-tbd_tmp->bd+1.) /
	                  (cur_img_data->to_image->frag.rb-cur_img_data->to_image->frag.lb+1.) *
	                  cur_img_data->to_image->to_msp->score;
	      if (cur_img_data->to_image->to_msp->direction == 1)
	        img_partner->frag.lb += tbd_tmp->bd - cur_img_data->to_image->frag.lb;
	      else
	        img_partner->frag.rb -= tbd_tmp->bd - cur_img_data->to_image->frag.lb;
	      cur_img_data->to_image->frag.lb = tbd_tmp->bd;
	    } else if (tbd_tmp->bd - cur_img_data->to_image->frag.rb >= -TOO_SHORT) {
	      /* TBD near right edge: trim rather than split */
	      cur_img_data->to_image->to_msp->score =
	        (int32_t) (tbd_tmp->bd-cur_img_data->to_image->frag.lb+1.) /
	                  (cur_img_data->to_image->frag.rb-cur_img_data->to_image->frag.lb+1.) *
	                  cur_img_data->to_image->to_msp->score;
	      if (cur_img_data->to_image->to_msp->direction == 1)
	        img_partner->frag.rb -= cur_img_data->to_image->frag.rb - tbd_tmp->bd;
	      else
	        img_partner->frag.lb += cur_img_data->to_image->frag.rb - tbd_tmp->bd;
	      cur_img_data->to_image->frag.rb = tbd_tmp->bd;
	    } else {
	      /* interior split: create left MSP copy, adjust original to right half */
	      dissected = 1;
	      msp_tmp = add_msp(cur_img_data->to_image->to_msp);
	      if (cur_img_data->to_image == &cur_img_data->to_image->to_msp->query) {
		target_img     = &msp_tmp->query;
		target_partner = &msp_tmp->sbjct;
	      } else {
		target_img     = &msp_tmp->sbjct;
		target_partner = &msp_tmp->query;
	      }
	      msp_tmp->score =
	        (int32_t) (tbd_tmp->bd-target_img->frag.lb+1.) /
	                  (target_img->frag.rb-target_img->frag.lb+1.) *
	                  msp_tmp->score;
	      if (msp_tmp->direction == 1) {
		target_partner->frag.rb -= target_img->frag.rb - tbd_tmp->bd;
	      } else {
		target_partner->frag.lb += target_img->frag.rb - tbd_tmp->bd;
	      }
	      target_img->frag.rb = tbd_tmp->bd;
	      fprint_msp(new_msps, msp_tmp);
	      register_image(target_img, cur_ele);

	      /* adjust the original MSP to represent the right half */
	      cur_img_data->to_image->to_msp->score -= msp_tmp->score;
	      if (cur_img_data->to_image->to_msp->direction == 1)
	        img_partner->frag.lb += tbd_tmp->bd - cur_img_data->to_image->frag.lb + 1;
	      else
	        img_partner->frag.rb -= tbd_tmp->bd - cur_img_data->to_image->frag.lb + 1;
	      cur_img_data->to_image->frag.lb = tbd_tmp->bd + 1;
	    }
	  }
	  if (tbd_tmp->bd >= cur_img_data->to_image->frag.rb || !tbd_tmp->next) {
	    if (dissected) {
		fprint_msp(new_msps, cur_img_data->to_image->to_msp);
	    }
	    if (too_short(&cur_img_data->to_image->frag) ||
	        too_short(&img_partner->frag)) {
	      if (next && next->to_image->to_msp == cur_img_data->to_image->to_msp)
	        next = next->next;
	      dump_image(cur_img_data->to_image);
	    }
	    break;
	  }
	}
	tbd_tmp = tbd_tmp->next;
      }
      {
        clock_t s = clock() - r;
        double dissectTIME = ((double)s)/CLOCKS_PER_SEC;
        (void)dissectTIME; /* timing available; suppress unused-variable warning */
      }
      cur_img_data = next;
    }
}


/* =========================================================================
 * Combo / obs bookkeeping
 * ========================================================================= */

/* combo_update  --  dismiss a parent element after its images are redistributed
 *
 * Marks the element 'X', removes all its edges (updating partner stats),
 * cleans up CP lists and image trees, then records it in the combo or
 * obsolete output depending on whether it has child elements.
 */
void combo_update(ELE_INFO_t *ele_info) {
  RLOG_DBG("combo_update before stat: %d\n", ele_info->stat);
  if (ele_info->ele->img_no < 0) {
    err_no ++;
    fprintf(log_file, "error:  combo ele %d has %d images\n",
            ele_info->index, ele_info->ele->img_no);
    fflush(log_file);
    exit(2);
  }
    ele_info->stat = 'X';
    if (ele_info->ele->edges) combo_edge_update(ele_info, &ele_info->ele->edges);
    if (ele_info->ele->edge_no) {
	err_no ++;
	fprintf(log_file, "error:  combo_ele %d, %d edge_node left\n",
	        ele_info->index, ele_info->ele->edge_no);
	fflush(log_file);
	exit(4);
    }
    ele_info->ele->flimg_no = 0;
    if (ele_info->ele->PCP) CP_free(&ele_info->ele->PCP);
    if (ele_info->ele->redef) {
      if (ele_info->ele->to_img_data) {
	err_no ++;
	fprintf(log_file, "error re-defining ele %d, still images left\n",
	        ele_info->index);
	fflush(log_file);
	exit(5);
      }
      if (ele_info->ele->to_img_tree)
        img_tree_free(&ele_info->ele->to_img_tree, ele_info);
      if (ele_info->ele->img_no) {
	err_no ++;
	fprintf(log_file, "error:  combo_ele %d, %d img_node left\n",
	        ele_info->index, ele_info->ele->img_no);
	fflush(log_file);
	ele_info->ele->img_no = 0;
	exit(2);
      }
      combo_output(ele_info);
    } else {
      if (ele_info->ele->img_no || ele_info->ele->to_img_data ||
          ele_info->ele->to_img_tree) {
	err_no ++;
	fprintf(log_file, "error:  images not cleaned in obs ele %d\n",
	        ele_info->index);
	fflush(log_file);
	exit(2);
      }
      else obs_output(ele_info);
    }
  RLOG_DBG("combo_update after stat: %d\n", ele_info->stat);
}


/* combo_edge_update  --  recursively remove all edges from a dismissed element
 *
 * For each edge, marks the partner as 'w' (needs combo update), cleans
 * the partner's PCP list, and removes the edge from the partner's BST.
 */
void combo_edge_update(ELE_INFO_t *ele_info, EDGE_TREE_t **edge_node_p) {
  ELE_INFO_t *epi;
  EDGE_TREE_t *edge_node = *edge_node_p;
  ELEMENT_t *ele = ele_info->ele;

  if (edge_node->l) combo_edge_update(ele_info, &(*edge_node_p)->l);
  if (edge_node->r) combo_edge_update(ele_info, &(*edge_node_p)->r);

  epi = linked_ele(ele_info, edge_node->to_edge);
  if (epi->stat == 'v') epi->stat = 'w';
  if (epi->ele->PCP) CP_clean(&epi->ele->PCP, ele_info);
  delete_edge(&epi->ele->edges, edge_node->to_edge);
  epi->ele->edge_no --;

  EDGE_free(edge_node->to_edge);

  ele_info->ele->edge_no --;
  free(edge_node);
  *edge_node_p = NULL;
}


/* =========================================================================
 * Output helpers
 * ========================================================================= */

/* combo_output  --  record a combo (parent of a dissection) in summary/combo */
void combo_output(ELE_INFO_t *ele_info) {
  fprintf(combo, "%d %s %d %d \n",
          ele_info->index, ele_info->ele->frag.seq_name,
          ele_info->ele->frag.lb, ele_info->ele->frag.rb);
  fflush(combo);
}


/* obs_output  --  record an obsolete (eliminated) element in summary/obsolete */
void obs_output(ELE_INFO_t *ele_info) {
  fprintf(obs, "%d %s %d %d \n",
          ele_info->index, ele_info->ele->frag.seq_name,
          ele_info->ele->frag.lb, ele_info->ele->frag.rb);
  fflush(obs);
}


/* fprint_ele_obs  --  write a detailed obsolete-element record to a FILE */
void fprint_ele_obs(FILE *fp, ELE_INFO_t *ele_info) {
  int i;
  BD_t *cur_bd;
  ELEMENT_t *ele=ele_info->ele;

  fprintf(fp, "ele %d\n", ele_info->index);
  fprintf(fp, "%s %d %d \n", ele->frag.seq_name, ele->frag.lb, ele->frag.rb);
  if (ele->TBD != NULL) {
    fprintf(fp, "cutting points\n");
    i = 0;
    cur_bd = ele->TBD;
    while (cur_bd != NULL) {
      fprintf(fp, "%d ", cur_bd->bd);
      i ++;
      if (i%10 == 0) fprintf(fp, "\n");
      cur_bd = cur_bd->next;
    }
  }
}
