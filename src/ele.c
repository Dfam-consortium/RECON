/*
 * ele.c  --  Element, edge, and family operations: globals, initializers, I/O
 *
 * Contains the definitions of all pipeline-wide globals declared extern in
 * ele.h, and the implementations of the core initializer, I/O, and memory-
 * cleanup functions.  BST operations are in ele_bst.c; debug print helpers
 * are in ele_print.c.  Data structure definitions (ELEMENT_t, ELE_INFO_t,
 * EDGE_t, FAMILY_t, etc.) remain in ele.h.
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#include "ele.h"
#include "ele_db.h"

/*
 * Pipeline-wide global definitions.
 * All of these are declared extern in ele.h; the linker resolves them
 * to the single instance defined here.  Each program that links ele.o
 * (eleredef, edgeredef, famdef) gets its own copy at run time (separate
 * processes, separate address spaces).
 */

/* all_ele        -- pre-allocated array of ELE_INFO_t pointers */
ELE_INFO_t **all_ele;

/* ele_ct, ele_array_size, fam_ct */
int ele_ct, ele_array_size, fam_ct;

/* clan_size, clan_core_size */
int clan_size, clan_core_size;

/* MSP memory-tracking counters */
int32_t msp_in_mem;
int32_t msp_left;
int32_t msp_ct;
int32_t msp_index;

/* Edge memory-tracking counters */
int32_t edge_index, edge_in_mem, edge_left, edge_ct;

/* Miscellaneous run-time counters */
int32_t files_read;
int32_t clan_ct;
int32_t err_no;

/* Overflow linked list for ele_info beyond ele_array_size */
ELE_INFO_t *ele_info_data, *ele_info_tail;

/* Pipeline output files -- opened by each program's main() */
FILE *new_msps;
FILE *eles;
FILE *unproc;
FILE *combo;
FILE *obs;
FILE *fams;
FILE *log_file;


/* ============================================================
 * Initializers
 * ============================================================ */

ELEMENT_t *ele_init(int index) {
  ELEMENT_t *ele_tmp = (ELEMENT_t *) malloc(sizeof(ELEMENT_t));

  ele_tmp->index = index;
  ele_tmp->direction = 1;
  ele_tmp->update = 0;
  ele_tmp->l_hold = 0;
  ele_tmp->img_no = 0;
  ele_tmp->flimg_no = 0;
  ele_tmp->edge_no = 0;
  ele_tmp->to_img_data = NULL;
  ele_tmp->to_img_tree = NULL;
  ele_tmp->edges = NULL;
  ele_tmp->PCP = NULL;
  ele_tmp->TBD = NULL;
  ele_tmp->redef =  NULL;

  return ele_tmp;
}


ELE_INFO_t *ele_info_init(int index) {
  ELE_INFO_t *ele_info_tmp;

  ele_info_tmp = (ELE_INFO_t *) malloc(sizeof(ELE_INFO_t));

  ele_info_tmp->index = index;
  ele_info_tmp->file_updated = 0;
  ele_info_tmp->stat = 'z';
  ele_info_tmp->ele = NULL;
  ele_info_tmp->to_family = NULL;
  ele_info_tmp->next = NULL;

  return ele_info_tmp;
}


MSP_t *MSP_malloc() {
  MSP_t *msp_tmp = (MSP_t *) malloc(sizeof(MSP_t));

  msp_in_mem ++;
  msp_ct ++;

  return msp_tmp;
}


void MSP_free(MSP_t *m) {
  msp_in_mem --;
  free(m);
}


EDGE_t *EDGE_malloc() {
  EDGE_t *edge_tmp = (EDGE_t *) malloc(sizeof(EDGE_t));

  edge_in_mem ++;
  edge_ct ++;

  return edge_tmp;
}


void EDGE_free(EDGE_t *ed) {
  edge_in_mem --;
  free(ed);
}


/* ============================================================
 * Small tools
 * ============================================================ */

ELE_INFO_t *linked_ele(ELE_INFO_t *ele_info, EDGE_t *edge) {
  if (ele_info->index == edge->ele1_info->index) return edge->ele2_info;
  return edge->ele1_info;
}


/*
 * outthrow_big_tandems  --  filter elements that are likely tandem repeats
 *
 * An element is classified as a tandem repeat if:
 *   ino  > TANDEM_SIZE_LIMIT (was: SIZE_LIMIT)  AND
 *   ino / p_ct > TANDEM_IMG_RATIO (was: TANDEM)
 *
 * The logic is that tandem repeats produce many images (ino) that all map
 * back to the same small set of nearby partner elements (p_ct), giving a
 * very high images-per-partner ratio.  Dispersed repeats, by contrast,
 * have many distinct partner elements across the genome.
 *
 * Filtered elements are marked 'O' (orphan) in all_ele and written to the
 * unproc summary file.
 *
 * Parameters
 *   size_list  FILE* positioned at the beginning of ele_def_res/size_list,
 *              which contains one "<ei> <ino>" (element_index image_count)
 *              line per element.
 */
void outthrow_big_tandems(FILE *size_list) {
  int ei, ino;
  FILE *ele_file;
  char line[150], head[10], *msp = "msp";
  int id, ele1, ele2;
  int *partners, i_ct, i, p_id, p_ct;
  char *db_buf;
  int   db_len;

  while (fgets(line, 25, size_list)) {
    sscanf(line, "%d %d\n", &ei, &ino);
    if (ino > SIZE_LIMIT) {
      db_buf = ele_db_read(ei, &db_len);
      if (!db_buf) {
	fprintf(log_file, "Can not open ele record %d.  Exit.\n", ei);
	exit(2);
      }
      ele_file = fmemopen(db_buf, db_len, "r");
      if (!ele_file) {
	fprintf(log_file, "fmemopen failed for ele %d.  Exit.\n", ei);
	free(db_buf);
	exit(2);
      }
      partners = (int *) malloc(ino*sizeof(int));
      p_ct = 0;
      i_ct = -1;
      p_id = -1;
      while (fgets(line, 150, ele_file)) {
	sscanf(line, "%s %*s", head);
	if (!strncmp(head, msp, 3)) {
	  sscanf(line, "msp %d %*c %*d %*f %*d %d %*s %*d %*d %d %*s %*d %*d\n", &id, &ele1, &ele2);
	  i_ct ++;
	  if (id%2) *(partners+i_ct) = ele1;
	  else *(partners+i_ct) = ele2;
	}
      }
      qsort(partners, i_ct+1, sizeof(int), int_cmp);
      for (i=0; i<=i_ct; i++) {
	if (*(partners+i) != p_id) {
	  p_id = *(partners+i);
	  p_ct ++;
	}
      }
      if (ino/p_ct > TANDEM) {
	(*(all_ele+ei-1))->stat = 'O';
	spit_out_ele(*(all_ele+ei-1));
      }
      free(partners);
      fclose(ele_file);
      free(db_buf);
    }
  }
}


int int_cmp(const void *i1, const void *i2) {
  return *((int *) i1) - *((int *) i2);
}


void spit_out_ele(ELE_INFO_t *ele_info) {
  ele_info->file_updated = 1;
  fprintf(unproc, "%d\n", ele_info->index);
}


/* ============================================================
 * Read in an element
 * ============================================================ */

/*
 * ele_read_in  --  deserialize an element from its file on disk
 *
 * Reads the element file (tmp/e<N> or tmp2/e<N>) into a newly allocated
 * ELEMENT_t and attaches it to ele_info->ele.  The file is selected by
 * ele_info->file_updated: 0 => tmp/, 1 => tmp2/.
 *
 * The 'stage' parameter controls how much of the file is parsed:
 *   stage 1  -- full load: frag, images (msp lines), edges, PCPs, TBDs
 *   stage 2  -- skip PCPs (eleredef second pass)
 *   stage 3  -- skip images and PCPs (edgeredef / famdef)
 *   stage 4  -- full load without PCPs (famdef variant)
 */
ELEMENT_t *ele_read_in(ELE_INFO_t *ele_info, int stage) {
  char line[200], head[10], rest[150];
  char fragname[SEQ_NAME_MAX_LEN];
  int pos;
  FILE *fp;
  IMAGE_t **img_array;
  EDGE_t **edge_array;
  int holder, img_encountered=0, edge_encountered=0, ct;
  CP_t *new_pcp;
  BD_t *new_tbd;
  ELE_DATA_t *cur_ele_data;
  char *db_buf;
  int   db_len;

  char *index="index", *stat = "stat", *file_updated="file_updated", *family="family", *direc="direc", *update="update", *l_hold="l_hold", *img_no="img_no", *flimg_no="flimg_no", *edge_no="edge_no", *frag="frag", *pcp="pcp", *tbd="tbd", *redef="redef", *msp="msp", *edge="edge";

  db_buf = ele_db_read(ele_info->index, &db_len);
  if (!db_buf) {
    fprintf(log_file, "error: ele_read_in: no record for element %d\n", ele_info->index);
    fflush(log_file);
    exit(2);
  }
  fp = fmemopen(db_buf, db_len, "r");
  if (!fp) {
    fprintf(log_file, "error: ele_read_in: fmemopen failed for element %d\n", ele_info->index);
    fflush(log_file);
    exit(2);
  }
  files_read++;

  ele_info->ele = ele_init(ele_info->index);

  while (fgets(line, 200, fp)) {
    sscanf(line, "%s %*s", head);

    if (!strncmp(head, update, 10)) {sscanf(line, "%*s %hd\n", &ele_info->ele->update); continue;}
    if (!strncmp(head, direc, 10)) {sscanf(line, "%*s %d\n", &ele_info->ele->direction); continue;}
    if (!strncmp(head, flimg_no, 10)) {sscanf(line, "%*s %d\n", &ele_info->ele->flimg_no); continue;}
    if (!strncmp(head, frag, 10)) {
      sscanf(line, "%*s %s %d %d\n", fragname, &ele_info->ele->frag.lb, &ele_info->ele->frag.rb);
      pos = GetSeqIndex(0, seq_count-1, fragname);
      ele_info->ele->frag.seq_name = seq_name_table[pos];
      if (ele_info->ele->frag.lb > ele_info->ele->frag.rb) {
	fprintf(log_file, "error:  ele %d reversed from read_in\n", ele_info->index);
	fflush(log_file);
	exit(3);
      }
      continue;
    }
    if (!strncmp(head, edge_no, 10)) {
      sscanf(line, "%*s %d\n", &ele_info->ele->edge_no);
      if (ele_info->ele->edge_no) edge_array = (EDGE_t **) malloc(ele_info->ele->edge_no*sizeof(EDGE_t *));
      continue;
    }
    if (!strncmp(head, edge, 10)) {
      if (edge_encountered < ele_info->ele->edge_no) {
	edge_scan(line, edge_array, ele_info, &edge_encountered);
      } else {
	err_no ++;
	fprintf(log_file, "error:  more edges in file for ele %d\n", ele_info->index);
	fflush(log_file);
	exit(4);
      }
      continue;
    }
    if (!strncmp(head, img_no, 10)) {
      sscanf(line, "%*s %d\n", &ele_info->ele->img_no);
      if (stage != 3 && ele_info->ele->img_no) img_array = (IMAGE_t **) malloc(ele_info->ele->img_no*sizeof(IMAGE_t *));
      continue;
    }
    if (stage != 3 && !strncmp(head, msp, 10)) {
      if (img_encountered < ele_info->ele->img_no) {
        img_scan(ele_info, line, img_array, &img_encountered);
      } else {
	err_no ++;
	fprintf(log_file, "error:  more image in file for ele %d\n", ele_info->index);
	fflush(log_file);
	exit(2);
      }
      continue;
    }
    if (stage < 2 && !strncmp(head, pcp, 10)) {
      new_pcp = (CP_t *) malloc(sizeof(CP_t));
      sscanf(line, "%*s %d %d\n", &new_pcp->cp, &holder);
      new_pcp->contributor = get_ele_info(holder);
      new_pcp->next = ele_info->ele->PCP;
      ele_info->ele->PCP = new_pcp;
      continue;
    }
    if (!strncmp(head, tbd, 10)) {
      new_tbd = (BD_t *) malloc(sizeof(BD_t));
      sscanf(line, "%*s %d %d\n", &new_tbd->bd, &new_tbd->support);
      new_tbd->next = ele_info->ele->TBD;
      ele_info->ele->TBD = new_tbd;
      continue;
    }
    if (!strncmp(head, redef, 10)) {
      sscanf(line, "%*s %d\n", &holder);
      cur_ele_data = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
      cur_ele_data->ele_info = get_ele_info(holder);
      cur_ele_data->next = ele_info->ele->redef;
#ifdef ORIGINAL_BUGS
      ele_info->ele->redef = cur_ele_data->next;  /* original: assigns old head, cur_ele_data lost */
#else
      ele_info->ele->redef = cur_ele_data;         /* fix: prepend cur_ele_data to redef list */
#endif
      continue;
    }
  }

  fclose(fp);
  free(db_buf);
  if (stage != 3 && img_encountered != ele_info->ele->img_no) {
    if (ele_info->file_updated) {
      err_no ++;
      fprintf(log_file, "error:  %d: image number not matched.\n", ele_info->index);
      fflush(log_file);
      exit(2);
    }
    ele_info->ele->img_no = img_encountered;
  }

  if (edge_encountered != ele_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  %d: edge number not matched.\n", ele_info->index);
    fflush(log_file);
    exit(4);
    ele_info->ele->edge_no = edge_encountered;
  }

  if (edge_encountered) {
    qsort(edge_array, edge_encountered, sizeof(EDGE_t *), edge_index_cmp);
    build_edge_tree(&ele_info->ele->edges, edge_array, 0, edge_encountered-1);
    free(edge_array);
  }
  if (stage != 3 && ele_info->ele->img_no) {
    qsort(img_array, img_encountered, sizeof(IMAGE_t *), img_index_cmp);
    build_img_tree(&ele_info->ele->to_img_tree, img_array, 0, img_encountered-1);
    free(img_array);
  }
  ele_info->file_updated = 1;

  return ele_info->ele;
}


void img_scan(ELE_INFO_t *ele_info, char *line, IMAGE_t **img_array, int *img_encountered_ptr) {
  MSP_t *msp_tmp = MSP_malloc();
  int id, ele1, ele2, id2;
  ELE_INFO_t *epi;
  IMAGE_t *ip;
  char qname[NAME_LEN], sname[NAME_LEN];
  int pos;

  sscanf(line, "msp %d %c %d %f %d %d %s %d %d %d %s %d %d\n", &id, &msp_tmp->stat, &msp_tmp->score, &msp_tmp->iden, &msp_tmp->direction, &ele1, qname, &msp_tmp->query.frag.lb, &msp_tmp->query.frag.rb, &ele2, sname, &msp_tmp->sbjct.frag.lb, &msp_tmp->sbjct.frag.rb);
  pos = GetSeqIndex(0, seq_count-1, qname);
  msp_tmp->query.frag.seq_name = seq_name_table[pos];
  pos = GetSeqIndex(0, seq_count-1, sname);
  msp_tmp->sbjct.frag.seq_name = seq_name_table[pos];

  if (id%2) {
    epi = get_ele_info(ele1);
    id2 = id - 1;
  } else {
    epi = get_ele_info(ele2);
    id2 = id + 1;
  }
  if (epi->stat == 'O') {
    MSP_free(msp_tmp);
  } else {
    if (ele1 != ele2) {
      if (epi->ele) {
	ip = find_image(epi->ele->to_img_tree, id2)->to_image;
	if (id%2) *(img_array+*img_encountered_ptr) = &ip->to_msp->sbjct;
	else *(img_array+*img_encountered_ptr) = &ip->to_msp->query;
	MSP_free(msp_tmp);
      } else {
	if (id%2) {
	  msp_tmp->query.index = id2;
	  msp_tmp->sbjct.index = id;
	  *(img_array+*img_encountered_ptr) = &msp_tmp->sbjct;
	  msp_tmp->query.ele_info = epi;
	  msp_tmp->sbjct.ele_info = ele_info;
	} else {
	  msp_tmp->query.index = id;
	  msp_tmp->sbjct.index = id2;
	  *(img_array+*img_encountered_ptr) = &msp_tmp->query;
	  msp_tmp->query.ele_info = ele_info;
	  msp_tmp->sbjct.ele_info = epi;
	}
	msp_tmp->query.to_msp = msp_tmp;
	msp_tmp->sbjct.to_msp = msp_tmp;
      }
      (*img_encountered_ptr) ++;
    } else {
      if (id%2) {
	msp_tmp->query.index = id2;
	msp_tmp->sbjct.index = id;
	msp_tmp->query.ele_info = epi;
	msp_tmp->sbjct.ele_info = ele_info;
      } else {
	msp_tmp->query.index = id;
	msp_tmp->sbjct.index = id2;
	msp_tmp->query.ele_info = ele_info;
	msp_tmp->sbjct.ele_info = epi;
      }
      msp_tmp->query.to_msp = msp_tmp;
      msp_tmp->sbjct.to_msp = msp_tmp;
      *(img_array+*img_encountered_ptr) = &msp_tmp->sbjct;
      (*img_encountered_ptr) ++;
      *(img_array+*img_encountered_ptr) = &msp_tmp->query;
      (*img_encountered_ptr) ++;
    }
  }
}


void edge_scan(char *line, EDGE_t **edge_array, ELE_INFO_t *ele_info, int *edge_encountered_ptr) {
  EDGE_t *new = EDGE_malloc();
  int ele1, ele2;
  ELE_INFO_t *epi;

  sscanf(line, "edge %d %c %d %d %d %d\n", &new->index, &new->type, &new->direction, &new->score, &ele1, &ele2);

  if (ele1 != ele2 && new->type != 'c') {
    if (ele_info->index == ele1) {
      epi = get_ele_info(ele2);
      new->ele1_info = ele_info;
      new->ele2_info = epi;
    } else {
      epi = get_ele_info(ele1);
      new->ele1_info = epi;
      new->ele2_info = ele_info;
    }
    if (epi->ele) {
      *(edge_array+*edge_encountered_ptr) = find_edge(epi->ele->edges, new->index)->to_edge;
      EDGE_free(new);
    } else {
      *(edge_array+*edge_encountered_ptr) = new;
    }
    (*edge_encountered_ptr) ++;
  } else EDGE_free(new);
}


ELE_INFO_t *get_ele_info(int index) {
  ELE_INFO_t *cur;
  if (index <= ele_array_size) return *(all_ele+index-1);
  cur = ele_info_data;
  while (cur && cur->index < index) {
    cur = cur->next;
  }
  if (cur && cur->index == index) return cur;
  return NULL;
}


int img_index_cmp(const void *i1, const void *i2) {
  return (*((IMAGE_t **) i1))->index - (*((IMAGE_t **) i2))->index;
}


int edge_index_cmp(const void *ed1, const void *ed2) {
  return (*((EDGE_t **) ed1))->index - (*((EDGE_t **) ed2))->index;
}


/* ============================================================
 * Write out an element
 * ============================================================ */

void ele_write_out(ELE_INFO_t *ele_info, int stage) {
  char   *buf  = NULL;
  size_t  bsz  = 0;
  FILE   *fp;
  CP_t   *cur_pcp;
  BD_t   *cur_bd;
  ELE_DATA_t *cur_ele_data;

  fp = open_memstream(&buf, &bsz);
  if (!fp) { perror("ele_write_out: open_memstream"); exit(1); }

  fprintf(fp, "frag %s %d %d\n", ele_info->ele->frag.seq_name, ele_info->ele->frag.lb, ele_info->ele->frag.rb);
  if (ele_info->ele->direction != 1) fprintf(fp, "direc %d\n", ele_info->ele->direction);
  if (ele_info->ele->update) fprintf(fp, "update %d\n", ele_info->ele->update);
  fprintf(fp, "img_no %d\n", ele_info->ele->img_no);
  if (ele_info->ele->flimg_no) fprintf(fp, "flimg_no %d\n", ele_info->ele->flimg_no);
  fprintf(fp, "edge_no %d\n", ele_info->ele->edge_no);
  if (ele_info->ele->img_no) write_out_msps(fp, ele_info->ele->to_img_tree);
  if (ele_info->ele->edge_no) write_out_edges(fp, ele_info->ele->edges, ele_info);
  if (stage < 2 && ele_info->ele->PCP) {
    cur_pcp = ele_info->ele->PCP;
    while (cur_pcp) {
      fprintf(fp, "pcp %d %d\n", cur_pcp->cp, cur_pcp->contributor->index);
      cur_pcp = cur_pcp->next;
    }
  }
  if (ele_info->ele->TBD) {
    cur_bd = ele_info->ele->TBD;
    while (cur_bd) {
      fprintf(fp, "tbd %d %d\n", cur_bd->bd, cur_bd->support);
      cur_bd = cur_bd->next;
    }
  }
  if (ele_info->ele->redef) {
    cur_ele_data = ele_info->ele->redef;
    while (cur_ele_data) {
      fprintf(fp, "redef %d\n", cur_ele_data->ele_info->index);
      cur_ele_data = cur_ele_data->next;
    }
  }
  fclose(fp);

  ele_db_write(ele_info->index, buf, (int)bsz);
  free(buf);
}


void write_out_msps(FILE *fp, IMG_TREE_t *rt) {
  if (rt->l) write_out_msps(fp, rt->l);

  if (rt->to_image->to_msp->query.ele_info->index != rt->to_image->to_msp->sbjct.ele_info->index || rt->to_image == &rt->to_image->to_msp->sbjct) {
    fprintf(fp, "msp %d ", rt->to_image->index);
    fprintf(fp, "%c ", rt->to_image->to_msp->stat);
    fprintf(fp, "%d %3.1f %d ", rt->to_image->to_msp->score, rt->to_image->to_msp->iden, rt->to_image->to_msp->direction);
    fprintf(fp, "%d %s %d %d ", rt->to_image->to_msp->query.ele_info->index, rt->to_image->to_msp->query.frag.seq_name, rt->to_image->to_msp->query.frag.lb, rt->to_image->to_msp->query.frag.rb);
    fprintf(fp, "%d %s %d %d\n", rt->to_image->to_msp->sbjct.ele_info->index, rt->to_image->to_msp->sbjct.frag.seq_name, rt->to_image->to_msp->sbjct.frag.lb, rt->to_image->to_msp->sbjct.frag.rb);
  }
  if (rt->r) write_out_msps(fp, rt->r);
}


void write_out_edges(FILE *fp, EDGE_TREE_t *rt, ELE_INFO_t *ele_info) {
  if (rt->l) write_out_edges(fp, rt->l, ele_info);

  if (rt->to_edge->type != 'c' && rt->to_edge->ele1_info->index != rt->to_edge->ele2_info->index) {
    fprintf(fp, "edge %d %c %d %d %d %d\n", rt->to_edge->index, rt->to_edge->type, rt->to_edge->direction, rt->to_edge->score, rt->to_edge->ele1_info->index, rt->to_edge->ele2_info->index);
  } else {
    fprintf(log_file, "error:  illegitimate edges in ele %d\n", ele_info->index);
    fflush(log_file);
    exit(4);
  }

  if (rt->r) write_out_edges(fp, rt->r, ele_info);
}


/* ============================================================
 * Memory cleanup
 * ============================================================ */

void ele_info_free(ELE_INFO_t **ei) {
  ELE_INFO_t *cur, *next;

  cur = *ei;
  while (cur) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *ei = NULL;
}


void ele_data_free(ELE_DATA_t **e) {
  ELE_DATA_t *cur, *next;

  cur = *e;
  while (cur != NULL) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *e = NULL;
}


void ele_cleanup(ELEMENT_t **e) {
  ELE_INFO_t *ele_info;

  if (!(*e)) return;

  if ((*e)->frag.lb > (*e)->frag.rb) {
    fprintf(log_file, "error:  ele %d reversed before clean up\n", (*e)->index);
    fflush(log_file);
    exit(3);
  }
  ele_info = get_ele_info((*e)->index);
  if ((*e)->to_img_data) img_data_free(&(*e)->to_img_data);
  if ((*e)->PCP) CP_free(&(*e)->PCP);
  if ((*e)->TBD) BD_free(&(*e)->TBD);
  if ((*e)->edges) edge_tree_cleanup(ele_info, &(*e)->edges);
  if ((*e)->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  element %d, %d edges left\n", (*e)->index, (*e)->edge_no);
    fflush(log_file);
    exit(4);
  }
  if ((*e)->to_img_tree) {
    img_tree_cleanup(&(*e)->to_img_tree, ele_info);
    if ((*e)->img_no) {
      err_no ++;
      fprintf(log_file, "error:  element %d, %d images left\n", (*e)->index, (*e)->img_no);
      fflush(log_file);
      exit(2);
    }
  }
  if ((*e)->redef) ele_data_free(&(*e)->redef);
  free(*e);
  *e = NULL;
}


void BD_free(BD_t **bds) {
  BD_t *cur, *next;

  cur = *bds;
  while (cur != NULL) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *bds = NULL;
}


void CP_free(CP_t **cps) {
  CP_t *cur, *next;

  cur = *cps;
  while (cur != NULL) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *cps = NULL;
}


void img_data_free(IMG_DATA_t **i) {
  IMG_DATA_t *cur, *next;

  cur = *i;
  while (cur != NULL) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *i = NULL;
}


void img_tree_cleanup(IMG_TREE_t **rt_p, ELE_INFO_t *ele_info) {
  ELE_INFO_t *ep;

  if ((*rt_p)->l) img_tree_cleanup(&(*rt_p)->l, ele_info);
  if ((*rt_p)->r) img_tree_cleanup(&(*rt_p)->r, ele_info);

  ele_info->ele->img_no --;
  ep = partner((*rt_p)->to_image)->ele_info;
  if (ep == (*rt_p)->to_image->ele_info) {
    /* Self MSPs: mark on first encounter (stat='m'), free on second. */
    if ((*rt_p)->to_image->to_msp->stat != 'm') (*rt_p)->to_image->to_msp->stat = 'm';
    else {
      MSP_free((*rt_p)->to_image->to_msp);
    }
  } else if (!ep->ele) {
    MSP_free((*rt_p)->to_image->to_msp);
  }

  free(*rt_p);
  *rt_p = NULL;
}


void img_tree_free(IMG_TREE_t **rt_p, ELE_INFO_t *ele_info) {
  if ((*rt_p)->l) img_tree_free(&(*rt_p)->l, ele_info);
  if ((*rt_p)->r) img_tree_free(&(*rt_p)->r, ele_info);

  ele_info->ele->img_no --;
  free(*rt_p);
  *rt_p = NULL;
}


void edge_tree_cleanup(ELE_INFO_t *ele_info, EDGE_TREE_t **rt_p) {
  ELE_INFO_t *ep;

  if ((*rt_p)->l) edge_tree_cleanup(ele_info, &(*rt_p)->l);
  if ((*rt_p)->r) edge_tree_cleanup(ele_info, &(*rt_p)->r);

  if (ele_info->index == (*rt_p)->to_edge->ele1_info->index) ep = (*rt_p)->to_edge->ele2_info;
  else ep = (*rt_p)->to_edge->ele1_info;
  if (!ep->ele) EDGE_free((*rt_p)->to_edge);

  ele_info->ele->edge_no --;
  free(*rt_p);
  *rt_p = NULL;
}


void edge_tree_free(EDGE_TREE_t **rt_p) {
  if ((*rt_p)->l) edge_tree_free(&(*rt_p)->l);
  if ((*rt_p)->r) edge_tree_free(&(*rt_p)->r);

  free(*rt_p);
  *rt_p = NULL;
}


void fam_data_free(FAM_DATA_t **fd) {
  FAM_DATA_t *cur, *next;

  cur = *fd;
  while (cur != NULL) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *fd = NULL;
}


void fam_data_cleanup(FAM_DATA_t **fd) {
  FAM_DATA_t *cur, *next;

  cur = *fd;
  while (cur != NULL) {
    next = cur->next;
    fam_cleanup(&cur->to_family);
    free(cur);
    cur = next;
  }
  *fd = NULL;
}


void fam_cleanup(FAMILY_t **f) {
  if ((*f)->members != NULL) ele_data_free(&(*f)->members);
  if ((*f)->relatives != NULL) fam_data_free(&(*f)->relatives);
  free(*f);
  *f = NULL;
}


/* Unused */
void msp_data_free(MSP_DATA_t **md) {
  MSP_DATA_t *cur, *next;

  cur = *md;
  while (cur != NULL) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *md = NULL;
}


void frag_data_free(FRAG_DATA_t **fd) {
  FRAG_DATA_t *cur, *next;

  cur = *fd;
  while (cur != NULL) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *fd = NULL;
}


void frag_data_cleanup(FRAG_DATA_t **fd) {
  FRAG_DATA_t *cur, *next;

  cur = *fd;
  while (cur != NULL) {
    next = cur->next;
    free(cur->to_frag);
    free(cur);
    cur = next;
  }
  *fd = NULL;
}
