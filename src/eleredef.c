/*
 * eleredef.c  --  Stage 3: element redefinition using the syntopy algorithm
 *
 * Algorithm overview
 * ------------------
 * This is the most complex stage of the RECON pipeline.  It takes the
 * initial elements from eledef and refines their boundaries by:
 *
 * 1. Identifying "Potential Cut Points" (PCPs): positions where many
 *    full-length image endpoints cluster, suggesting the true boundary
 *    between two distinct repeat subfamilies.
 *
 * 2. Clustering PCPs into "To-Be-Determined" boundaries (TBDs) via
 *    CP_cluster() and PCP_to_TBDs().
 *
 * 3. Dissecting elements at TBD positions (dissect()), splitting MSPs
 *    that span the boundary into left and right halves.
 *
 * 4. Re-clustering the dissected images into child elements (ele_def()).
 *
 * 5. Rebuilding edges between elements -- a directed graph where an edge
 *    between elements A and B means they share full-length (or near full-
 *    length) MSP images.  Edge type 'p' = primary, 's' = secondary.
 *
 * The algorithm operates on "local networks" (clans): BFS subgraphs
 * centered on each unprocessed element, extending DEPTH hops along edges.
 * This bounds peak memory use to O(DEPTH-hop neighbourhood size).
 *
 * Multi-round convergence:
 *   The main loop repeats until no element has stat 'z', 'w', or 't'
 *   (all elements are 'v' or 'X').  New elements created by dissection
 *   may require further processing in subsequent rounds.
 *
 * Element status transitions (see ELE_INFO_t.stat in ele.h for full list):
 *   'z' -> edges_and_cps() -> 't'
 *   't' -> ele_redef()     -> 'v' (no further split) or 'w' (combo update)
 *   'w' -> combo_update()  -> 'v'
 *   'v' -> (ready for edgeredef)
 *   'X' -> dismissed (image_count dropped to 0 after dissection)
 *
 * Named constants from recon_defs.h (constants formerly #defined here):
 *   ELEREDEF_CUTOFF_SINGLE (was CUTOFF1 = 0.5)
 *   ELEREDEF_CUTOFF_DOUBLE (was CUTOFF2 = 0.9)
 *   ELEREDEF_MAX_IMAGES    (was MAX_IMG  = 1200)
 *   MIN_ELEMENT_LEN_BP     (was TOO_SHORT = 30)
 *   SPLIT_RATIO_TOLERANCE  (was FUDGE = 2)
 *   ELEMENT_BOUNDARY_MARGIN (was MARGIN = 10000)
 *   MAX_DISSECT_PASSES     (was FLURRY = 10; note: this is a distance
 *                           threshold in bp, not a count despite the name)
 *   BFS_CLAN_DEPTH         (was DEPTH = 3)
 *
 * Usage
 *   eleredef seq_list [start] [clan_ct] [-l log_level]
 *
 *   -l <level>  log verbosity: 0=silent 1=error 2=warn 3=info 4=debug
 *               (default: 3=info)
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#include <fcntl.h>
#include <sys/types.h>
#include <dirent.h>
#include <time.h>

#include "ele.h"
#include "eleredef.h"
#include "seqlist.h"
#include "treeview.h"

/* All constants now come from recon_defs.h via ele.h -> msps.h -> bolts.h.
 * Backward-compat aliases (CUTOFF1, CUTOFF2, MAX_IMG, TOO_SHORT, FUDGE,
 * MARGIN, FLURRY) are defined there and remain usable throughout this file. */

/* ---- Per-program log level storage (required by recon_log.h) ---- */
int   recon_log_level = RECON_LOG_INFO;
FILE *recon_log_fp    = NULL;

//typedef struct img_node {
//  short recorded;
//  IMAGE_t *to_image;
//  struct img_node *sib;
//  struct img_node *children;
//} IMG_NODE_t;

void report_cts();
void report_redef_stat();

ELE_DATA_t *ele_def(IMG_DATA_t **, float);
void generate_img_tree(ELEMENT_t *);
ELE_INFO_t *new_element();
void add_ele_info(ELE_INFO_t *);

void general_ele_redef(ELE_INFO_t *, IMAGE_t **);
void build_local_network(ELE_INFO_t *, ELE_DATA_t **, ELE_DATA_t **, IMAGE_t **);
void recruit(ELE_INFO_t *, EDGE_TREE_t *, ELE_DATA_t **, IMAGE_t **);
void cruise_local_net(ELE_DATA_t *, IMAGE_t **);
void local_ele_redef(ELE_INFO_t *, IMAGE_t **, int*);
void dissolve_local_network(ELE_DATA_t **);
int emptyDir(char *);
void stage_exit(int);
void dismiss_element(ELE_INFO_t *);
void remove_ele(ELE_INFO_t *);

void ele_redef(ELE_INFO_t *, IMAGE_t **);
IMG_DATA_t *img_data_sort(IMG_DATA_t *, int);
void PCP_to_TBDs(ELEMENT_t *);
BD_t *CP_cluster(CP_t *);
void CP_sort(CP_t **);
void BD_sort(BD_t **);
int span(ELEMENT_t *, int32_t);
void TBD_merge(ELEMENT_t *);

void dissect(ELE_INFO_t *);
int too_short(FRAG_t *);
MSP_t *add_msp(MSP_t *);
void register_image(IMAGE_t *, ELEMENT_t *);
void put_image(IMG_DATA_t **, IMAGE_t *);
void dump_image(IMAGE_t *);
void remove_image(IMAGE_t *);
void combo_update(ELE_INFO_t *);
void combo_edge_update(ELE_INFO_t *, EDGE_TREE_t **);
void CP_clean(CP_t **, ELE_INFO_t *);

void edges_and_cps(ELE_INFO_t *, IMAGE_t **);
int full_length(IMAGE_t *, float);
void add_CP(CP_t **, int32_t, ELE_INFO_t *);
void add_edge(ELE_INFO_t *, ELE_INFO_t *, char, int32_t, short);
void adjust_edge_tree(ELE_INFO_t *);
int charge_edge_array(EDGE_t **, EDGE_TREE_t *, int);
int consis_tree_build(IMG_NODE_t *, IMAGE_t *, int);
int print_consis_tree(IMG_NODE_t *rt);
int consis(IMAGE_t *, IMAGE_t *, float);
IMG_NODE_t **node_entry(IMG_NODE_t **);
void consis_tree_free(IMG_NODE_t *);
int find_prim(IMG_NODE_t *, float, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int32_t, int *, int32_t *, short *);

void combo_output(ELE_INFO_t *);
void obs_output(ELE_INFO_t *);
void fprint_ele_obs(FILE *, ELE_INFO_t *);

int frag_cmp(const void *, const void *);
int partner_cmp(const void *, const void *);
int CP_cmp(const void *, const void *);
int BD_cmp(const void *, const void *);
int fam_cmp(const void *, const void *);

int main (int argc, char *argv[]) {
  ELE_INFO_t *cur_ele_info;
  int i, ele_march, ei, rounds=0, start;
  char line[35], stat;
  short fu, to_march;
     clock_t start1, r, t, end;
     double cpu_time_used=0.0; 
     double ele_defTIME;
     double dissectTIME;
  IMAGE_t **img_ptr;
  FILE *ele_no, *msp_no, *edge_no, *size_list, *new_stat;
  FILE *seq_list;

  /* Strip optional "-l <level>" before positional arg parsing */
  if (recon_parse_log_flag(&argc, argv)) {
    fprintf(stderr, "error: -l requires a numeric log level argument\n");
    exit(1);
  }

  /* Validate command line */
  if (argc == 1) {
    printf("usage: eleredef seq_list [start] [clan_ct] [-l level]\n"
           "  seq_list  list of sequence names\n"
           "  start     1-based element index to start from (optional, default 1)\n"
           "  clan_ct   initial clan counter value (optional, default 0)\n"
           "  -l <level>  log verbosity: 0=silent 1=error 2=warn "
           "3=info(default) 4=debug\n");
    exit(1);
  }

  seq_list = fopen(argv[1], "r");
  if (!seq_list) {
    printf("Input file for sequence name list %s not found.  Exit.\n", argv[1]);
    exit(2);
  }
  GetSeqNames(seq_list);

  if (argc > 2) start = atoi(argv[2]) - 1;
  else start = 0;
  if (argc > 3) clan_ct = atoi(argv[3]);
  else clan_ct = 0;

  ele_no = fopen("summary/naive_ele_no", "r");
  if(!ele_no) {
    printf("Can not open naive_ele_no.  Exiting\n");
    exit(1);
  }
  msp_no = fopen("summary/redef_msp_no", "r");
  if (!msp_no) msp_no = fopen("summary/ori_msp_no", "r");
  if (!msp_no) {
    printf("Can not open msp_no.  Exiting\n");
    exit(1);
  }

  // May or may not exist.  All uses of edge_no first check for
  // its existence.
  edge_no = fopen("summary/naive_edge_no", "r");

  // May or may not exist.  All uses are checked.
  size_list = fopen("tmp/size_list", "r");

  // May or may not exist.  All uses are checked.
  new_stat = fopen("tmp2/redef_stat", "r");

  new_msps = fopen("summary/new_msps", "a");
  if ( !new_msps )
  {
    printf("Can not open summary/new_msps for writing! Exiting\n");
    exit(1);
  }
  unproc = fopen("summary/unproc", "a");
  if ( !unproc )
  {
    printf("Can not open summary/unproc for writing! Exiting\n");
    exit(1);
  }
  combo = fopen("summary/combo", "a");
  if ( !combo )
  {
    printf("Can not open summary/combo for writing! Exiting\n");
    exit(1);
  }
  obs = fopen("summary/obsolete", "a");
  if ( !obs )
  {
    printf("Can not open summary/obsolete for writing! Exiting\n");
    exit(1);
  }
  log_file = fopen("tmp2/log", "a");
  if ( !log_file )
  {
    printf("Can not open tmp2/log for writing! Exiting\n");
    exit(1);
  }
  /* Route RECON_LOG macros to the same file as the pipeline log */
  recon_log_fp = log_file;

  /* Open element database (written by eledef, updated by this stage) */
  ele_db_open();

  while (fgets(line, 15, ele_no)) {
    ele_ct = atoi(line);
  }
  fclose(ele_no);

  while (fgets(line, 15, msp_no)) {
    msp_index = atoi(line) - 1;
  }
  fclose(msp_no);

  if (edge_no) {
    while (fgets(line, 15, edge_no)) {
      edge_index = atoi(line) - 1;
    }
  } else edge_index = -1;

  // RMH: Create instances of the element structure
  //      (see ele.h).  The initial state of the
  //      an element is to have a '0' file_updated
  //      state and a 'z' status.
  ele_array_size = 2*ele_ct;
  all_ele = (ELE_INFO_t **) malloc(ele_array_size*sizeof(ELE_INFO_t *));
  for (i=0; i<ele_array_size; i++) {
    *(all_ele+i) = ele_info_init(i+1);
  }
  ele_info_data = NULL;

  /* get rid of large tandem reps */
  if (size_list && !new_stat) outthrow_big_tandems(size_list);
  fclose(unproc);

  /* set file_updated / stat for ele_info's */
  // RMH: This is not used in the normal workflow at
  //      the eleredef stage.
  if (new_stat) {
    printf("Reading stat file redef_stat...\n");
    ele_ct = 0;
    while (fgets(line, 35, new_stat)) {
      ele_ct ++;
      /*sscanf(line, "%d %c %d\n", &ei, &stat, &fu);*/
      /* for some bizzare reason, scanf doesn't work properly here. :( */
      for (i=0; i<35; i++) {
	if (line[i] == ' ') break;
      }
      ei = atoi(line);
      stat = line[i+1];
      fu = atoi(&line[i+3]);
      if (ei<=ele_array_size) {
	(*(all_ele+ei-1))->stat = stat;
	(*(all_ele+ei-1))->file_updated = fu;
      } else {
	cur_ele_info = ele_info_init(ei);
	cur_ele_info->stat = stat;
	cur_ele_info->file_updated = fu;
	add_ele_info(cur_ele_info);
      }
    }
  }

  msp_in_mem = 0;
  msp_left = 0;
  msp_ct = 0;
  edge_ct = 0;
  edge_in_mem = 0;
  edge_left = 0;
  files_read = 0;
  err_no = 0;

  /* re-define elements using the syntopy algorithm, and build edges */
  img_ptr = (IMAGE_t **) malloc(MAX_IMG*sizeof(IMAGE_t *));

  if ( ! img_ptr )
  {
    printf("eleredef: Error! Could not allocate memory for img_ptr: %zu bytes requested\n", ( MAX_IMG*sizeof(IMAGE_t *) ) );
    exit(-1);
  }

  to_march = 1;
  start1 = clock();
  while (to_march) {
    to_march = 0;
    rounds ++;
    for (i=start; i<ele_ct && i<ele_array_size; i++) {
      fprintf(log_file, "Evaluating definition of element %d\n", (*(all_ele+i))->index);
      if ((*(all_ele+i))->stat == 'z' || (*(all_ele+i))->stat == 'w' || (*(all_ele+i))->stat == 't') {
	to_march = 1;
	general_ele_redef(*(all_ele+i), img_ptr);
#if 0
	report_redef_stat();
#endif
      } /*else if ((*(all_ele+i))->stat == 'O' && !(*(all_ele+i))->file_updated) spit_out_ele(*(all_ele+i));*/
    }


    start = 0;

    // RMH: This appears to be a spill-over linked-list datastructure used after the all_eles pre-allocated
    //      array is exceeded.
    cur_ele_info = ele_info_data;
    while(cur_ele_info) {
      fprintf(log_file, "evaluating definition of element %d\n", cur_ele_info->index);
      if (cur_ele_info->stat == 'z' || cur_ele_info->stat == 'w' || cur_ele_info->stat == 't') {
	to_march = 1;
	general_ele_redef(cur_ele_info, img_ptr);
#if 0
	report_redef_stat();
#endif
      } /*else if (cur_ele_info->stat == 'O' && !cur_ele_info->file_updated) spit_out_ele(cur_ele_info);*/
      cur_ele_info = cur_ele_info->next;

    // printf("redeftime 1 %f \n", cpu_time_used);
    }
      end = clock()-start1;
       cpu_time_used = ((double) (end - start1)) / CLOCKS_PER_SEC; 
  //  printf("redeftime 2 %f \n", cpu_time_used);
  }
     cpu_time_used = ((double) (end - start1)) / CLOCKS_PER_SEC; 
//printf("redeftime 3 %f \n", cpu_time_used);
  report_cts();
  report_redef_stat();
  free(img_ptr);

#if 0
  for (i=start; i<ele_ct && i<ele_array_size; i++) {
    if ((*(all_ele+i))->stat == 'X') remove_ele(*(all_ele+i));
  }

  cur_ele_info = ele_info_data;
  while(cur_ele_info) {
    if (cur_ele_info->stat == 'X') remove_ele(cur_ele_info);
    cur_ele_info = cur_ele_info->next;
  }
#endif

  fprintf(log_file, "total numbers: %d elements, %d msps, %d edges\n", ele_ct, msp_index+1, edge_index+1);
  fprintf(log_file, "%d rounds, %d files read, %d msps seen, %d edges seen\n", rounds, files_read, msp_ct, edge_ct);
  fprintf(log_file, "%d errors, %d msps and %d edges left in memory, \n", err_no, msp_left, edge_left);
  RLOG_DBG("General_ele_redef %f , %f , %f \n", cpu_time_used, ele_defTIME, dissectTIME);
  RLOG_DBG("total numbers: %d elements, %d msps, %d edges\n", ele_ct, msp_index+1, edge_index+1);
  RLOG_DBG("%d rounds, %d files read, %d msps seen, %d edges seen\n", rounds, files_read, msp_ct, edge_ct);
  RLOG_DBG("%d errors, %d msps and %d edges left in memory, \n", err_no, msp_left, edge_left);
  fflush(log_file);
  fclose(log_file);

  ele_db_close();

  exit(0);
}



void report_cts() {
  FILE *fp;

  fp = fopen("summary/redef_ele_no", "w");
  fprintf(fp, "%d\n", ele_ct);
  fclose(fp);
}




void report_redef_stat() {
  int i, ele_march = ele_ct < ele_array_size ? ele_ct : ele_array_size;
  FILE *redef_stat, *fp;
  ELE_INFO_t *cur_ele_info = ele_info_data;

  redef_stat = fopen("tmp2/redef_stat", "w");
  for (i=0; i<ele_march; i++) {
    fprintf(redef_stat, "%d %c %d\n", (*(all_ele+i))->index, (*(all_ele+i))->stat, (*(all_ele+i))->file_updated);
  }
  while(cur_ele_info) {
    fprintf(redef_stat, "%d %c %d\n", cur_ele_info->index, cur_ele_info->stat, cur_ele_info->file_updated);
    cur_ele_info = cur_ele_info->next;
  }
  fclose(redef_stat);

  fp = fopen("summary/redef_msp_no", "w");
  fprintf(fp, "%d\n", msp_index+1);
  fclose(fp);

  fp = fopen("summary/naive_edge_no", "w");
  fprintf(fp, "%d\n", edge_index+1);
  fclose(fp);
}


/***********************************
 ***********************************
 **                               **
 **       DEFINING  ELEMENTS      **
 **                               **
 ***********************************
 ***********************************/


ELE_DATA_t *ele_def(IMG_DATA_t **img_data_p, float cutoff) {

      clock_t t;
    t = clock();


  /* this function assumes that the images are SORTED according to frag_cmp.
     defining elements is in essence partitioning images in the original list.
     cur_img_data marks the current image under inspection;
     *img_data_p and prev_img_data are used to maintain the un-participated list;
     *img_data_p marks the first (un-participated) image in the remaining list;
     prev_img_data marks the last un-participated image proceding cur_img_data,
     so that when cur_img_data is removed from the origial list and participated
     into an element, we know where the tail of the unparticipated list, used to
     be pointed by cur_img_data->next, should be appended to */

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
      /* start the definition of the element */
      ele_tmp->img_no = 1;
      ele_tmp->frag = cur_img_data->to_image->frag; /* lst and rst included */
      /* kick cur_img_data out of the original list */
      /* notice that whenever we're here, cur_img_data points to the first image in the remaining list, which is also pointed to by *img_data_p */
      *img_data_p = cur_img_data->next;
      /* updat the current image and move it to the image list of ele_tmp */
      /* notice that here it's always the first one for ele_tmp */
      cur_img_data->to_image->ele_info = ele_info_tmp;
      cur_img_data->next = NULL;
      ele_tmp->to_img_data = cur_img_data;
      /* move pointer cur_img_data to the next one */
      cur_img_data = *img_data_p;
      ritetime = 0;
      continue;
    }
    if (/*!strncmp(ele_tmp->frag.seq_name, cur_img_data->to_image->frag.seq_name, NAME_LEN)*/ ele_tmp->frag.seq_name == cur_img_data->to_image->frag.seq_name && ele_tmp->frag.rb - cur_img_data->to_image->frag.lb > 10) { /*checking possible images */
      if (sing_cov(&ele_tmp->frag, &cur_img_data->to_image->frag, cutoff)) { /* a good image */
	ele_tmp->img_no ++;
	/* move cur_img_data from the original list */
	/* notice that it's a bit more complicated than above */
	if (prev_img_data != NULL) {
	  prev_img_data->next = cur_img_data->next;
	} else {
	  *img_data_p = cur_img_data->next;
	}
	/* update the current image and move it to the image_list of ele_tmp */
	cur_img_data->to_image->ele_info = ele_info_tmp;
	cur_img_data->next = ele_tmp->to_img_data;
	ele_tmp->to_img_data = cur_img_data;

	/* updating definition of the element if necessary, and move cur_img_data to the right place */
	/* notice that it is a bit complicated than above */
	if (ele_tmp->frag.rb < cur_img_data->to_image->frag.rb) { 
	  ele_tmp->frag.rb = cur_img_data->to_image->frag.rb;
	  /* ele_tmp->frag.rst = cur_img_data->to_image->frag.rst; */
	  /* time to go back and check if previously unqualified images are now good to be participated to the updated element */
	  cur_img_data = *img_data_p;
	  prev_img_data = NULL;
	} else {
	  if (prev_img_data != NULL) cur_img_data = prev_img_data->next;
	  else cur_img_data = *img_data_p;
	}
      } else { /* current image not good for the ele_tmp, keep it in the remaining list and move on */
	prev_img_data = cur_img_data;
	cur_img_data = cur_img_data->next;
      }
    } else { /* time to finish defining of ele_tmp */
      /* add the element to the list of elements */
      ele_data_tmp = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
      ele_data_tmp->ele_info = ele_info_tmp;
      ele_data_tmp->next = ele_data;
      ele_data = ele_data_tmp;
      /*      if (!strncmp(ele_tmp->frag.seq_name, target, NAME_LEN)) printf("%d %s %d %d\n", ele_tmp->index, ele_tmp->frag.seq_name, ele_tmp->frag.lb, ele_tmp->frag.rb); */
      /* time to go back to the head of the remaining list and start a new element */
      ritetime = 1;
      cur_img_data = *img_data_p;
      prev_img_data = NULL;
    }
  }
  /* add the last element to the list of elements */
  ele_data_tmp = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
  ele_data_tmp->ele_info = ele_info_tmp;
  ele_data_tmp->next = ele_data;
  ele_data = ele_data_tmp;

 }

  ele_data_tmp = ele_data;
  while (ele_data_tmp) {
    if (ele_data_tmp->ele_info->ele->frag.lb > ele_data_tmp->ele_info->ele->frag.rb) {
      fprintf(log_file, "error:  ele %d reversed after ele_def\n", ele_data_tmp->ele_info->index);
      fflush(log_file);
      exit(3);
    }
    generate_img_tree(ele_data_tmp->ele_info->ele);
    /*ct = count_img_nodes(ele_data_tmp->ele_info->ele->to_img_tree);
    if (ct != ele_data_tmp->ele_info->ele->img_no) {
	err_no ++;
	fprintf(log_file, "error:  trouble generating image tree for redefed offsprings: ele %d %d %d\n", ele_data_tmp->ele_info->index, ele_data_tmp->ele_info->ele->img_no, ct);
	fflush(log_file);
	exit(2);
    }*/
        t = clock() - t;
    double ele_defTIME = ((double)t)/CLOCKS_PER_SEC;
 //   printf("eledef time %f", ele_defTIME);
    ele_data_tmp = ele_data_tmp->next;
  }


  return ele_data;
}




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





ELE_INFO_t *new_element() {
  ELE_INFO_t *ele_info_tmp;

  ele_ct ++;

  // RMH: It seems that the initial storage for elements is managed
  //      as a pre-allocated, initialized array of eles.  Then once
  //      that is exceeded they are allocated individually and it's
  //      treated as a linked list from then on.
  if (ele_ct <= ele_array_size) ele_info_tmp = *(all_ele+ele_ct-1);
  else {
    ele_info_tmp = ele_info_init(ele_ct);
    // This adds the element to the tail of the ele_info_data spill over..
    add_ele_info(ele_info_tmp);
  }

  ele_info_tmp->ele = ele_init(ele_ct);
  ele_info_tmp->file_updated = 1;

  return ele_info_tmp;
}




void add_ele_info(ELE_INFO_t *ele_info_tmp) {

  if (!ele_info_data) ele_info_data = ele_info_tmp;
  else ele_info_tail->next = ele_info_tmp;

  ele_info_tail = ele_info_tmp;

}





/***********************************
 ***********************************
 **                               **
 **    REDEFINING THE ELEMENTS    **
 **                               **
 ***********************************
 ***********************************/


;

/**********************************************
 * Organizing the traverse of the local graph *
 **********************************************/




// z : original stat
// t : after edges_and_cps
// v : after local_ele_redef
// X : after dismiss_element
void general_ele_redef(ELE_INFO_t *ele_info, IMAGE_t **img_ptr) {
  int i;
  ELE_DATA_t *local_net=NULL, *local_net_tail=NULL;
  /* Read in stage 1 element from "e#" ( where # is ele_info->index ). */
  RLOG_DBG("general_ele_redef: ele %d\n", ele_info->index);
  if (!ele_info->ele) ele_read_in(ele_info, 1);
  //RMH
  //  status still 'z', but file_updated is now 1
  //  update=0, l_hold=0 and only the IMG_TREE_t is built ( to_img_data is Null )
  //print_ele_info(ele_info);
  //print_ele(ele_info->ele);

  /* an element could have lost all its images during the re-evaluation of its partners */
  if (!ele_info->ele->img_no) {
    combo_update(ele_info);
    remove_ele(ele_info);
    ele_cleanup(&ele_info->ele);
  } else {
    /* redefining cur_ele and its adjacent neighbors */
    fprintf(new_msps, "new clan for ele %d\n", ele_info->index);
    fprintf(combo, "new clan for ele %d\n", ele_info->index);
    fprintf(obs, "new clan for ele %d\n", ele_info->index);
    /* set up the local network */
    clan_ct ++;
    clan_size = 0;
    clan_core_size = 0;

    fprintf(log_file, "new clan: %d for ele %d\n", clan_ct, ele_info->index);

    // RMH: Build a graph centered on ele_info and extending out 3 degrees
    //      of separation. TODO: Add details of edge labeling.
    //
    //   local_net is a linked list of ELE_INFO_t pointers (ELE_DATA_t) that
    //             is built up by this method.
    build_local_network(ele_info, &local_net, &local_net_tail, img_ptr);

    /* Debug: dump the local network graph as GML */
    if (recon_log_level >= RECON_LOG_DEBUG) print_all_eles_GML();

    fprintf(log_file, "clan size: %d, clan core size: %d\n", clan_size, clan_core_size);

    /*
     * Redefine elements in the local network:
     * queues all elements in local_net and calls local_ele_redef on each.
     */
    cruise_local_net(local_net, img_ptr);

    /* Debug: dump the updated graph */
    if (recon_log_level >= RECON_LOG_DEBUG) print_all_eles_GML();

    /* clearing up the local network */
    fflush(new_msps);
    fflush(combo);
    fflush(obs);
    dissolve_local_network(&local_net);
  }
}



void build_local_network(ELE_INFO_t *ele_info, ELE_DATA_t **net_p, ELE_DATA_t **net_tail_p, IMAGE_t **img_ptr) {
  ELEMENT_t *ele;
  ELE_DATA_t *que;
  RLOG_DBG("build_local_network: seed ele %d\n", ele_info->index);
  /* seed the network with the first element */
  que = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
  que->ele_info = ele_info;
  que->next = NULL;
  *net_p = que;
  *net_tail_p = que;

  ele_info->ele->l_hold = 1;

  /* breadth first search */
  while (que) {
    clan_size ++;
    RLOG_DBG("build_local_network: BFS queue entry ele %d\n",
             que->ele_info->ele->index);
    // RMH: DEPTH currently set in bolts.h to 3
    if (que->ele_info->ele->l_hold <= DEPTH) {
      clan_core_size ++;
      // RMH: This builds a network of element relationships centered on ele_info and extending
      //      DEPTH degress.  It also generates the list of endpoints for later clustering ( CP_t )
      //      Sets ele_info to 't'
      if (que->ele_info->stat == 'z') edges_and_cps(que->ele_info, img_ptr);
      // Read in elements referenced in initial graph and add to queue ( up to 3 levels deep )
      // The counter l_hold is incremented for all the edge partners of que->ele_info in this method
      if (que->ele_info->ele->edges) recruit(que->ele_info, que->ele_info->ele->edges, net_tail_p, img_ptr);
    }
    que = que->next;
  }
}


void recruit(ELE_INFO_t *ele_info, EDGE_TREE_t *rt, ELE_DATA_t **net_tail_p, IMAGE_t **img_ptr) {
  ELE_INFO_t *epi;
  ELE_DATA_t *member;
  RLOG_DBG("recruit: ele %d\n", ele_info->index);
  if (rt->l) recruit(ele_info, rt->l, net_tail_p, img_ptr);

  /* Extract the partner element using the edge information */
  epi = linked_ele(ele_info, rt->to_edge);
  RLOG_DBG("recruit: partner ele %d\n", epi->index);
  if (!epi->ele) ele_read_in(epi, 1);
  if (!epi->ele->l_hold) {
    epi->ele->l_hold = ele_info->ele->l_hold + 1;
    // TODO: Determine when the 'v' state is set
    if (epi->stat == 'v' && epi->ele->l_hold < DEPTH) epi->ele->l_hold = DEPTH;
    member = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
    member->ele_info = epi;
    member->next = NULL;
    (*net_tail_p)->next = member;
    *net_tail_p = member;
  }

  if (rt->r) recruit(ele_info, rt->r, net_tail_p, img_ptr);
}




void cruise_local_net(ELE_DATA_t *local_net, IMAGE_t **img_ptr) {
  ELE_DATA_t *que;
  int to_march = 1;

  while (to_march) {
    // TODO: Determine why to_march was necessary.  I would assume that
    //       local_ele_redef adds to the que, but they must have needed
    //       a way to restart from the begining of the queue to reprocess
    //       all the data again.
    to_march = 0;
    que = local_net;
    while (que) {
	if (que->ele_info->ele->l_hold <= DEPTH) {
	  local_ele_redef(que->ele_info, img_ptr, &to_march);
	}
	que = que->next;
    }
  }
}




void local_ele_redef(ELE_INFO_t *ele_info, IMAGE_t **img_ptr, int *march_p) {
    ELE_DATA_t *new_ele_data;
    ELEMENT_t *ele = ele_info->ele;

    RLOG_DBG("local_ele_redef(): ele %d stat=%c redef=%s PCP=%s\n",
             ele_info->index, ele_info->stat,
             ele->redef ? "yes" : "NULL",
             ele->PCP   ? "yes" : "NULL");

    if (ele->redef != NULL) {
      new_ele_data = ele->redef;
      while (new_ele_data != NULL) {
        local_ele_redef(new_ele_data->ele_info, img_ptr, march_p);
        new_ele_data = new_ele_data->next;
      }
    } else {
	if (ele_info->stat != 'v' && ele_info->stat != 'X') {
	  if (ele->PCP) {
  	    *march_p = 1;
	    ele_redef(ele_info, img_ptr);
	  } else ele_info->stat = 'v';
	}
    }
}





void dissolve_local_network(ELE_DATA_t **net_p) {
  ELE_DATA_t *que, *new_ele_data;
  int i, ele_left=0;
  char *command;
  /*MSP_DATA_t *md;*/

  int j, in_clan;
  FILE *fp;
  ELE_INFO_t *cur_ele_info;

  /* output results to tmp2/clan, which holds results temporily and clean up memory */
  que = *net_p;
  while (que) {
    dismiss_element(que->ele_info);
    que = que->next;
  }

  if (msp_in_mem) {
    err_no ++;
    fprintf(log_file, "error:  error in bookkeeping: %d msps total, %d seen, %d left in memory\n", msp_index+1, msp_ct, msp_in_mem);
    fflush(log_file);
    msp_left += msp_in_mem;
    /*md = all_msps;
    while (md) {
      if (md->to_msp) {
	fprint_msp(log_file, md->to_msp);
	fflush(log_file);
	MSP_free(md->to_msp);
      }
      md = md->next;
    }
    if (msp_in_mem) {
      fprintf(log_file, "error: further discrepency: %d msps left\n", msp_in_mem);
      msp_in_mem = 0;
    }*/
  }

  if (edge_in_mem) {
    err_no ++;
    fprintf(log_file, "error:  error in bookkeeping: %d edges total, %d seen, %d left in memory\n", edge_index+1, edge_ct, edge_in_mem);
    fflush(log_file);
    edge_left += edge_in_mem;
    edge_in_mem = 0;
  }

  if (err_no) exit(3);

  ele_data_free(net_p);
}


void dismiss_element(ELE_INFO_t *ele_info) {
    ELE_DATA_t *new_ele_data;

    if (ele_info->ele->redef) {
      new_ele_data = ele_info->ele->redef;
      while (new_ele_data) {
        dismiss_element(new_ele_data->ele_info);
	new_ele_data = new_ele_data->next;
      }
    } 
    /*printf("clearing ele %d\n", ele_info->index);*/
    ele_info->ele->l_hold = 0;
    if (ele_info->stat != 'X') {
      if (!ele_info->ele->img_no) {
	ele_info->stat = 'X';
	if (ele_info->ele->redef) combo_output(ele_info);
	else obs_output(ele_info);
      } else ele_write_out(ele_info, 1);
    }
    ele_cleanup(&ele_info->ele);
}




void remove_ele(ELE_INFO_t *ele_info) {
    char *command;

    command = (char *) malloc(50*sizeof(char));
    sprintf(command, "rm -f tmp2/e%d", ele_info->index);
    if (system(command)) {
      fprintf(log_file, "error removing file tmp2/e%d\n", ele_info->index);
      fflush(log_file);
      exit(6);
    }
    free(command);
}






/******************************
 * Redefining a given element *
 ******************************/

/*

Unless otherwise stated as 'v' elements evaluated under ele_redef do not have a stat -kn
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
    }
    else{
      RLOG_DBG("ele_redef: ele %d has no PCP\n", cur_ele->index);
    }
    if (cur_ele->TBD) {
      TBD_merge(cur_ele);
    }
    if (cur_ele->TBD) {
#if 0
      to_dissect = 1;\
      if (cur_ele->TBD->bd - cur_ele->frag.lb <= FLURRY || cur_ele->TBD->bd - cur_ele->frag.rb >= -FLURRY) {
	if (!cur_ele->TBD->next) {
	  to_dissect = 0;
	  if (cur_ele->TBD->bd - cur_ele->frag.lb <= FLURRY) cur_ele->frag.lb = cur_ele->TBD->bd;
	  else cur_ele->frag.rb = cur_ele->TBD->bd;
	  BD_free(&cur_ele->TBD);
	} else {
	  if (cur_ele->TBD->next->bd - cur_ele->frag.rb >= -FLURRY) {
	    if (!cur_ele->TBD->next->next) {
	      to_dissect = 0;
	      cur_ele->frag.lb = cur_ele->TBD->bd;
	      cur_ele->frag.rb = cur_ele->TBD->next->bd;
	      BD_free(&cur_ele->TBD);
	    }
	  }
	}
      }
#endif
      pbd = cur_ele->TBD;
      while (pbd) {
	if (pbd->bd-cur_ele->frag.lb>FLURRY && pbd->bd-cur_ele->frag.rb<-FLURRY) {
	  to_dissect ++;
	}
	pbd = pbd->next;
      }

      if (!to_dissect) {
	pbd = cur_ele->TBD;
	while (pbd) {
	  if (/*pbd->bd-cur_ele->frag.lb>0 &&*/ pbd->bd-cur_ele->frag.lb<=FLURRY) {
	    cur_ele->frag.lb = pbd->bd;
	  } else if (/*pbd->bd-cur_ele->frag.rb<0 &&*/ pbd->bd-cur_ele->frag.rb>=-FLURRY) {
	    cur_ele->frag.rb = pbd->bd;
	  }
	  pbd = pbd->next;
	}
	BD_free(&cur_ele->TBD);
      } else {
	/* dissect all images according to the TBDs */
	dissect(ele_info);
	/* redefine elements according to the dissected images, if any left */
	if (cur_ele->img_no) {
	  cur_ele->to_img_data = img_data_sort(cur_ele->to_img_data, cur_ele->img_no);

	  cur_ele->redef = ele_def(&cur_ele->to_img_data, CUTOFF1);

	}
	/* clear unnecessary memory, 'v'->'w' and update CPs for w's */
	combo_update(ele_info);
	/* ele_redef() for offspring elements.  we finish redefinition of all
         * offsprings before pulling in partners, because the combo may
         * have > 1 copy of the same family, which means other members' CPs
         * won't be fully updated until all offsprings are processed
         */
	new_ele_data = cur_ele->redef;
	while (new_ele_data != NULL) {
	  new_ele_data->ele_info->ele->update = 1;
	  new_ele_data->ele_info->ele->l_hold = cur_ele->l_hold;
	  if (new_ele_data->ele_info->ele->to_img_tree) {
	    edges_and_cps(new_ele_data->ele_info, img_ptr);
	    if (new_ele_data->ele_info->ele->PCP) ele_redef(new_ele_data->ele_info, img_ptr);
	    else new_ele_data->ele_info->stat = 'v';
	  } else {
	   if (new_ele_data->ele_info->ele->img_no) {
	    err_no ++;
	    fprintf(log_file, "error:  image tree missing in newly redefed offspring %d\n", new_ele_data->ele_info->index);
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
    /*img_data_free(&cur_ele->to_img_data);*/
}





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




/***************************
 * Determine possible TBDs *
 ***************************/



void PCP_to_TBDs(ELEMENT_t *ele) {
  RLOG_DBG("PCP_to_TBDs: ele %d, first PCP contributor ele %d\n",
           ele->index, ele->PCP->contributor->index);
  int s = 0, left;
  BD_t *pbd_tmp, *pbd_prev, *pbd, *pbds;
  CP_t *cp;

  /* sort the PCP list according to cp */
  CP_sort(&ele->PCP);
  /*clustering the PCPs into PBDs */
  pbds = CP_cluster(ele->PCP);
  /* identify TBP from PBDs */
  /* TBDs are removed from PBDs, what is left in PBDs are those unsuccessful ones */
  pbd_tmp = pbds;
  pbd_prev = NULL;
  while (pbd_tmp != NULL) {
    /* s is the KEY! */
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


// RMH:
//    CP_t is a linked list structure that holds
//    an integer and a pointer to an ele_info structure.
//       typedef struct cp_list {
//           int32_t cp;
//           struct ele_info *contributor;
//           struct cp_list *next;
//       } CP_t;
//
// CP: cp=1, ele_info.index = 2, ele_info.stat = t
// CP: cp=1597, ele_info.index = 2, ele_info.stat = t
// CP: cp=1597, ele_info.index = 2, ele_info.stat = t
// CP: cp=2, ele_info.index = 1, ele_info.stat = v
// CP: cp=1600, ele_info.index = 1, ele_info.stat = v
// CP: cp=1600, ele_info.index = 1, ele_info.stat = v
// CP: cp=2, ele_info.index = 4, ele_info.stat = t
// CP: cp=399, ele_info.index = 4, ele_info.stat = t
// CP: cp=399, ele_info.index = 4, ele_info.stat = t
// CP: cp=1, ele_info.index = 3, ele_info.stat = v
// CP: cp=400, ele_info.index = 3, ele_info.stat = v
// CP: cp=400, ele_info.index = 3, ele_info.stat = v
//
//  cps->cp appears to be endpoints?
//
BD_t *CP_cluster(CP_t *cps) {
  int32_t first = cps->cp, last = cps->cp, sum = 0;
  CP_t *begin = cps;
  int cpct = 0;
  BD_t *bds = NULL, *bd_tmp;

  while (cps != NULL) {
    // RMH: start
    //printf("CP: cp=%d, ele_info.index = %d, ele_info.stat = %c\n", cps->cp, cps->contributor->index, cps->contributor->stat);
    // RMH: end
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




int span(ELEMENT_t *ele, int32_t cut) {
  /* span requires PCP sorted according to cp_cmp and images sorted according to frag_cmp */
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





/**************************
 * Dissecting the element *
 **************************/


// RMH: Perhaps where splitting occurs?
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
	  if (tbd_tmp->bd > cur_img_data->to_image->frag.lb && tbd_tmp->bd < cur_img_data->to_image->frag.rb) {
	    if (tbd_tmp->bd - cur_img_data->to_image->frag.lb <= TOO_SHORT) {
	      cur_img_data->to_image->to_msp->score = (int32_t) (cur_img_data->to_image->frag.rb-tbd_tmp->bd+1.)/(cur_img_data->to_image->frag.rb-cur_img_data->to_image->frag.lb+1.)*cur_img_data->to_image->to_msp->score;
	      if (cur_img_data->to_image->to_msp->direction == 1) img_partner->frag.lb += tbd_tmp->bd - cur_img_data->to_image->frag.lb;
	      else img_partner->frag.rb -= tbd_tmp->bd - cur_img_data->to_image->frag.lb;
	      cur_img_data->to_image->frag.lb = tbd_tmp->bd;
	    } else if (tbd_tmp->bd - cur_img_data->to_image->frag.rb >= -TOO_SHORT) {
	      cur_img_data->to_image->to_msp->score = (int32_t) (tbd_tmp->bd-cur_img_data->to_image->frag.lb+1.)/(cur_img_data->to_image->frag.rb-cur_img_data->to_image->frag.lb+1.)*cur_img_data->to_image->to_msp->score;
	      if (cur_img_data->to_image->to_msp->direction == 1) img_partner->frag.rb -= cur_img_data->to_image->frag.rb - tbd_tmp->bd;
	      else img_partner->frag.lb += cur_img_data->to_image->frag.rb - tbd_tmp->bd;
	      cur_img_data->to_image->frag.rb = tbd_tmp->bd;
	    } else {
	      dissected = 1;
	      /* create a new MSP to hold the left products of the dissection */
	      msp_tmp = add_msp(cur_img_data->to_image->to_msp);
	      /* upgrade the content of the new MSP */
	      if (cur_img_data->to_image == &cur_img_data->to_image->to_msp->query) { 
		target_img = &msp_tmp->query;
		target_partner = &msp_tmp->sbjct;
	      } else {
		target_img = &msp_tmp->sbjct;
		target_partner = &msp_tmp->query;
	      }
	      msp_tmp->score = (int32_t) (tbd_tmp->bd-target_img->frag.lb+1.)/(target_img->frag.rb-target_img->frag.lb+1.)*msp_tmp->score;
	      if (msp_tmp->direction == 1) {
                // RMH: When tracking down the divide by zero problem I ended up signaling this
                //      spot in the code as the likely place where things went initially wrong.
                //      Explore further.
		target_partner->frag.rb -= target_img->frag.rb - tbd_tmp->bd;
	      }
	      else {
		target_partner->frag.lb += target_img->frag.rb - tbd_tmp->bd;
	      }
	      target_img->frag.rb = tbd_tmp->bd;
	      fprint_msp(new_msps, msp_tmp);
	      /* keep or ignore the new msp */
	      register_image(target_img, cur_ele);

	      /* create the right product */
	      /* we cheat here.  instead of generating a new MSP, we change the content
               * of the original one, then output it into new_msp when finished.  in
               * the meantime, the index of the MSP is unchanged, which points to the
               * oringinal MSP
               */
	      cur_img_data->to_image->to_msp->score -= msp_tmp->score;
	      if (cur_img_data->to_image->to_msp->direction == 1) img_partner->frag.lb += tbd_tmp->bd - cur_img_data->to_image->frag.lb + 1;
	      else img_partner->frag.rb -= tbd_tmp->bd - cur_img_data->to_image->frag.lb + 1;
	      cur_img_data->to_image->frag.lb = tbd_tmp->bd +1;
	    }
	  }
	  if (tbd_tmp->bd >= cur_img_data->to_image->frag.rb || !tbd_tmp->next) { /* end of this image */
	    if (dissected) {
		fprint_msp(new_msps, cur_img_data->to_image->to_msp);
	    }
	    if (too_short(&cur_img_data->to_image->frag) || too_short(&img_partner->frag)) {
	      if (next && next->to_image->to_msp == cur_img_data->to_image->to_msp) next = next->next;
	      dump_image(cur_img_data->to_image);
	    }
	    break;
	  }
	}
	tbd_tmp = tbd_tmp->next;
      }
             r = clock() - r;
    double dissectTIME = ((double)r)/CLOCKS_PER_SEC;
//printf("dissect time  per img %f, \n", dissectTIME);
      cur_img_data = next;

    }
    /* the reason for not generating new elements according to the TBDs and then dissect and partition the images is that there's no gaining in time, it's virtually the same thing */

}




int too_short(FRAG_t *f) {
  if (f->rb - f->lb <= TOO_SHORT) return 1;
  return 0;
}


MSP_t *add_msp(MSP_t *m) {
    MSP_t *msp_tmp;
    /*MSP_DATA_t *hanger;*/

    msp_tmp = MSP_malloc();
    /*hanger = msp_tmp->hanger;*/
    *msp_tmp = *m;
    /*msp_tmp->hanger = hanger;*/
    msp_index ++;
    msp_tmp->query.to_msp = msp_tmp;
    msp_tmp->query.index = 2*msp_index;
    msp_tmp->sbjct.to_msp = msp_tmp;
    msp_tmp->sbjct.index = 2*msp_index+1;
    msp_tmp->stat = 's';
    /*msp_in_mem ++;*/
    /*msp_ct ++;*/

    return msp_tmp;
}




void register_image(IMAGE_t *i, ELEMENT_t *ele) {
  IMAGE_t *ip = partner(i);
  int ct;

  if (too_short(&i->frag)) {
    /* rec_image(&i->ele_info->ele->ignored, i);
    rec_image(&ip->ele_info->ele->ignored, ip); */
    /* the msp here havn't been added to the corresponding elements, so we can simply free it */
    MSP_free(i->to_msp);
  } else {
      /*ct = count_img_nodes(i->ele_info->ele->to_img_tree);
      if (ct != i->ele_info->ele->img_no) {
	err_no ++;
	fprintf(log_file, "error:  image tree changed before inserting a new image: ele %d %d %d\n", i->ele_info->index, i->ele_info->ele->img_no, ct);
	fflush(log_file);
	exit(2);
      }*/
      insert_image(&i->ele_info->ele->to_img_tree, i);
      i->ele_info->ele->img_no ++;
      /*ct = count_img_nodes(i->ele_info->ele->to_img_tree);
      if (ct != i->ele_info->ele->img_no) {
	err_no ++;
        fprintf(log_file, "error:  trouble inserting a new image: ele %d %d %d\n", i->ele_info->index, i->ele_info->ele->img_no, ct);
	fflush(log_file);
        exit(2);
      }*/
      if (i->ele_info->ele->to_img_data) put_image(&i->ele_info->ele->to_img_data, i);
      /*ct = count_img_nodes(ip->ele_info->ele->to_img_tree);
      if (ct != ip->ele_info->ele->img_no) {
	err_no ++;
        fprintf(log_file, "error:  tree changed before inserting a new image: ele %d %d %d\n", ip->ele_info->index, ip->ele_info->ele->img_no, ct);
	fflush(log_file);
        exit(2);
      }*/
      insert_image(&ip->ele_info->ele->to_img_tree, ip);
      ip->ele_info->ele->img_no ++;
      /*ct = count_img_nodes(ip->ele_info->ele->to_img_tree);      
      if (ct != ip->ele_info->ele->img_no) {      
	err_no ++;
        fprintf(log_file, "error:  trouble inserting a new image: ele %d %d %d\n", ip->ele_info->index, ip->ele_info->ele->img_no, ct); 
	fflush(log_file);
        exit(2);      
      }*/      
      if (ip->ele_info->ele->to_img_data) put_image(&ip->ele_info->ele->to_img_data, ip);
  }
}





void put_image(IMG_DATA_t **idp, IMAGE_t *i) {
  IMG_DATA_t *img_data_tmp;

  img_data_tmp = (IMG_DATA_t *)malloc(sizeof(IMG_DATA_t));
  img_data_tmp->to_image = i;
  img_data_tmp->next = *idp;
  *idp = img_data_tmp;
}




void dump_image(IMAGE_t *i) {
    IMAGE_t *ip = partner(i);
    int ct;

    if (i->ele_info->ele->to_img_data) remove_image(i);
    if (ip->ele_info->ele->to_img_data) remove_image(ip);

    /*ct = count_img_nodes(i->ele_info->ele->to_img_tree);
    if (ct != i->ele_info->ele->img_no) {
	err_no ++;
	fprintf(log_file, "error:  image tree changed before deleting a node: ele %d %d %d\n", i->ele_info->index, i->ele_info->ele->img_no, ct);
	fflush(log_file);
	exit(2);
    }*/
    delete_image(&i->ele_info->ele->to_img_tree, i);
    i->ele_info->ele->img_no --;
    /*ct = count_img_nodes(i->ele_info->ele->to_img_tree);
    if (ct != i->ele_info->ele->img_no) {
	err_no ++;
        fprintf(log_file, "error:  trouble deleting an image node: ele %d %d %d\n", i->ele_info->index, i->ele_info->ele->img_no, ct);
	fflush(log_file);
        exit(2);
    }*/

    /*ct = count_img_nodes(ip->ele_info->ele->to_img_tree);
    if (ct != ip->ele_info->ele->img_no) {
	err_no ++;
        fprintf(log_file, "error:  image tree changed before deleting a node: ele %d %d %d\n", ip->ele_info->index, ip->ele_info->ele->img_no, ct);
	fflush(log_file);
        exit(2);
    }*/
    delete_image(&ip->ele_info->ele->to_img_tree, ip);
    ip->ele_info->ele->img_no --;
    /*ct = count_img_nodes(ip->ele_info->ele->to_img_tree);
    if (ct != ip->ele_info->ele->img_no) {
	err_no ++;
        fprintf(log_file, "error:  trouble deleting an image node: ele %d %d %d\n", ip->ele_info->index, ip->ele_info->ele->img_no, ct);
	fflush(log_file);
        exit(2);
    }*/
    
    MSP_free(i->to_msp);
}




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




/*****************************
 *updating dissected element *
 *****************************/





void combo_update(ELE_INFO_t *ele_info) {
  RLOG_DBG("combo_update before stat: %d\n", ele_info->stat);
  if (ele_info->ele->img_no < 0) {
    err_no ++;
    fprintf(log_file, "error:  combo ele %d has %d images\n", ele_info->index, ele_info->ele->img_no);
    fflush(log_file);
    exit(2);
  }
    ele_info->stat = 'X';
    if (ele_info->ele->edges) combo_edge_update(ele_info, &ele_info->ele->edges);
    if (ele_info->ele->edge_no) {
	err_no ++;
	fprintf(log_file, "error:  combo_ele %d, %d edge_node left\n", ele_info->index, ele_info->ele->edge_no);
	fflush(log_file);
	exit(4);
    }
    ele_info->ele->flimg_no = 0;
    if (ele_info->ele->PCP) CP_free(&ele_info->ele->PCP);
    if (ele_info->ele->redef) {
      if (ele_info->ele->to_img_data) {
	err_no ++;
	fprintf(log_file, "error re-defining ele %d, still images left\n", ele_info->index);
	fflush(log_file);
	exit(5);
      }
      /* all images and msps are relocated to offsprings, so img_tree_free() is enough */
      if (ele_info->ele->to_img_tree) img_tree_free(&ele_info->ele->to_img_tree, ele_info);
      if (ele_info->ele->img_no) {
	err_no ++;
	fprintf(log_file, "error:  combo_ele %d, %d img_node left\n", ele_info->index, ele_info->ele->img_no);
	fflush(log_file);
	ele_info->ele->img_no = 0;
	exit(2);
      }
      combo_output(ele_info);
    } else {
      if (ele_info->ele->img_no || ele_info->ele->to_img_data || ele_info->ele->to_img_tree) {
	err_no ++;
	fprintf(log_file, "error:  images not cleaned in obs ele %d\n", ele_info->index);
	fflush(log_file);
	exit(2);
      }
      else obs_output(ele_info);
    }
  RLOG_DBG("combo_update after stat: %d\n", ele_info->stat);
}





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





void CP_clean(CP_t **cps_p, ELE_INFO_t *cont) {
  /*int res=0;*/
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
      /*res ++;*/
    } else {
      cp_prev = cp_cur;
      cp_cur = cp_cur->next;
    }
  }
  /*return res;*/
}





/****************************************************
 ****************************************************
 **  functions for finding PCPs and set up edges   **
 ****************************************************
 ****************************************************/



//
// Graph building:
//
//       Nodes: elements
//       Edges: Images
//
//        E1 -----primary-----> E2
//            If there exists an alignment that is nearly
//            full length (within 10bp of both ends) for
//            either E1 or E2. Or...if more a set of partial
//            images could be seed as part of full length
//            alignment. Edge stat='p'.
//
//        E1 ----secondary---> E2
//            If there exists an non-full length image connecting
//            them.
//
void edges_and_cps(ELE_INFO_t *ele_info, IMAGE_t **img_ptr) {
    /* the core function in this part, calls everything below */
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

    RLOG_DBG("edges_and_cps: ele_info->index = %d, ele_info->ele->index = %d\n", ele_info->index, ele_info->ele->index);

    /* sort unprocessed images according to their partner element,
     * and then their left bounds to allocate proper amount of memory
     * when update is 1, it means the element is an offspring of a combo,
     * in which case all the images need to be processed; when update is 0,
     * those images whose partner is 'v' or 'w' should be omitted
     * no self-images are accepted to generate pcps or edges, updating or not
     */
    // If to_image_data is Null generate a list version from the tree version in to_img_tree
    if (!cur_ele->to_img_data) listify(cur_ele->to_img_tree, &cur_ele->to_img_data);
    cur_img_data = cur_ele->to_img_data;
    // If element update == 1 count all non-self images otherwise only count
    // 'z' status, non-self images.
    // TODO: This is redundant with the loop below. In the future we could
    //       set *(img_ptr+eff_img_ct) in this first loop and use it below
    //       should eff_img_ct be > 0.  NOTE: We should also fix the allocation
    //       of img_ptr.  It could either be malloc'd afresh each time for use
    //       here or realloc'd only when the current size is exceeded.  Either
    //       way we should verify were it should be free'd.
    eff_img_ct = 0;
    while(cur_img_data != NULL) {
      epi = partner(cur_img_data->to_image)->ele_info;
      if (epi->index != ele_info->index) {
        if (cur_ele->update || epi->stat == 'z') eff_img_ct ++;
      }
      cur_img_data = cur_img_data->next;
    }
    RLOG_DBG("edges_and_cps: img_no = %d, eff_img_ct = %d\n", cur_ele->img_no, eff_img_ct);

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

    cur_img_data = cur_ele->to_img_data;
    eff_img_ct = 0;
    // TODO: Redundant loop see above.
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
    // RMH: Sort by: difference in element identifiers, ascending
    //               msp_direction
    //               left bound
    //               right bound
    //      This considers all images from a pair of elements at a
    //      time ( sub sort starting from img_ptr and going for eff_img_ct records ).
    qsort(img_ptr, eff_img_ct, sizeof(IMAGE_t *), partner_cmp);

    /* recognizing full-length images, and put partial and secondary
     * images into a consistency tree, in which images connected are
     * consistent with each other (look at the consis() function of
     * definition of consistency)
     */
    ritetime = 1;
    for (i=0; i<eff_img_ct; i++) {
      cur_img = *(img_ptr+i);
      img_partner = partner(cur_img);

      RLOG_DBG("edges_and_cps: image cur_ele_info->index=%d, par_ele_info->index=%d\n", cur_img->ele_info->index, img_partner->ele_info->index);

      if (ritetime) { /* new ele_partner begins */
	epi = img_partner->ele_info;
	ritetime = 0;
	if (!epi->ele) ele_read_in(epi, 1);
	ele_partner = epi->ele;
	riteplace = i;
	prim = 0;
	prim_p = NULL;
	max_score = 0;
      }
      if (img_partner->ele_info->index == epi->index) { /* still the same partner element */
	/* full length */
        // RMH: Currently 0.9
        //      So unless the sequence is < 200bp it simply requires
        //      that the image is within 10bp of the endpoints of the element.
	if (full_length(cur_img, CUTOFF2)) {
          RLOG_DBG("edges_and_cps:     image is primary because full_length with current element!\n");
	  prim = 1;
          // Increment the full length image counter
	  cur_ele->flimg_no ++;
	}
	if (full_length(img_partner, CUTOFF2)) {
          RLOG_DBG("edges_and_cps:     image is primary because full_length with partner element!\n");
	  prim = 1;
          // Increment the full length image counter
	  ele_partner->flimg_no ++;
	}
        // If the image is full length for either/both of the
        // elements and is the highest score seen so far,
        // save it as prim_p
	if (prim == 1) {
	  prim = 0;
          // RMH: Track maximum scoring primary (full length)
          //      image.
	  if (cur_img->to_msp->iden > max_score) {
            RLOG_DBG("edges_and_cps:     **** image is the new high scoring primary!\n");
	    max_score = cur_img->to_msp->iden;
	    dir = cur_img->to_msp->direction;
	    prim_p =  cur_img->to_msp;
	  }
	}
        //
	if (!prim_p) {
          // RMH: No full length images found yet
          //      NOTE: Prequal flag is set
          RLOG_DBG("edges_and_cps:      Adding to consis tree (not full-length and no full-length found yet): e%d im = %s:%d-%d   e%d pt = %s:%d-%d\n", epi->index, cur_img->frag.seq_name, cur_img->frag.lb, cur_img->frag.rb, img_partner->ele_info->index, img_partner->frag.seq_name, img_partner->frag.lb, img_partner->frag.rb);
	  consis_tree_build(consis_rt, cur_img, 1);
          //print_consis_tree(consis_rt);
          //print_ascii_tree(consis_rt);
	}
      }
      if (img_partner->ele_info->index != epi->index || i == eff_img_ct-1) {
        RLOG_DBG("edges_and_cps:     Reached the end of the current element partner -- will reprocess this current image!\n");
        // reached the end of the current ele_partner
	// start a new ele_partner
	ritetime = 1;
        // RMH: Yikes...backing up the for loop counter here so he could reprocess this image again
        // after setting ritetime.
	if (img_partner->ele_info->index != epi->index) i --;
	// finish the current ele_partner
	if (prim_p) {
	  prim_p->stat = 'p';
	  prim = 1;
	} else { // if no full-length, find partial primary images
          // This identifies partial primary images by looking for fragments that share an outside edge
	  prim = find_prim(consis_rt->children, CUTOFF2, ele_info->ele->frag.lb, -1, 0, 0, 0, 0, 0, &token_mark, &max_score, &dir);
          if ( prim )
            RLOG_DBG("edges_and_cps:     Identified partial primary image!\n");
	}
	/* build edge */
	if (prim) {
          RLOG_DBG("edges_and_cps:     Adding primary edge!\n");
	  if (ele_info->index != epi->index) add_edge(ele_info, epi, 'p', max_score, dir);
	  else {
	    err_no ++;
	    fprintf(log_file, "error:  self edge seen: ele %d\n", ele_info->index);
	    fflush(log_file);
	  }
	  for (j=riteplace; j<=i; j++) {
	    if ((*(img_ptr+j))->to_msp->stat == 'p') {
	      cur_img = *(img_ptr+j);
	      img_partner = partner(cur_img);
              //RMH: start
              //  Primary edge endpoints are added to CP lists
              //printf("1: Adding CPs: %d, %d, to ele=%d\n", cur_img->frag.lb, cur_img->frag.rb, epi->index);
              //printf("2: Adding CPs: %d, %d, to ele=%d\n", img_partner->frag.lb, img_partner->frag.rb, ele_info->index);
              //RMH: end
	      add_CP(&cur_ele->PCP, cur_img->frag.lb, epi);
	      add_CP(&cur_ele->PCP, cur_img->frag.rb, epi);
	      add_CP(&ele_partner->PCP, img_partner->frag.lb, ele_info);
	      add_CP(&ele_partner->PCP, img_partner->frag.rb, ele_info);
	    }
	  }
	} else { /* no primary images, full-length or partial */
          RLOG_DBG("edges_and_cps:     Adding secondary edge!\n");
	  if (ele_info->index != epi->index) add_edge(ele_info, epi, 's', 0, 0);
          else {
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
      }else { RLOG_DBG("edges_and_cps:     image fell through all cases (NOT SURE)\n"); }
    }
//printf("printing edge tree\n");
//print_edge_tree(ele_info->ele->edges, 0);

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



int full_length(IMAGE_t *i, float cutoff) {
  if (!i->ele_info->ele) {
    err_no ++;
    fprintf(log_file, "error:  element %d not in memory\n", i->ele_info->index);
    fflush(log_file);
    // RMH: Looks like this could have recovered here, but exited instead.
    exit(3);
    ele_read_in(i->ele_info, 1);
  }
  // If the image is within 10 bp of the start/end of the
  // element to which it belongs, and the ratio of
  // the image's size to the element's size is > cutoff
  //
  // For 0.9 cutoff this is only being more stringent
  // for sequences < 200bp.  Not sure why they did that.
  if ( i->frag.lb - i->ele_info->ele->frag.lb < 10 &&
       i->frag.rb - i->ele_info->ele->frag.rb > -10 &&
       ((float)i->frag.rb - i->frag.lb) / (i->ele_info->ele->frag.rb - i->ele_info->ele->frag.lb) > cutoff) {
    return 1;
  }
  return 0;
}



void add_CP(CP_t **CP_ptr, int32_t cp, ELE_INFO_t *cont) {
  CP_t *new = (CP_t *) malloc(sizeof(CP_t));
  new->cp = cp;
  new->contributor = cont;
  new->next = *CP_ptr;
  *CP_ptr = new;
}



void add_edge(ELE_INFO_t *ele1_info, ELE_INFO_t *ele2_info, char type, int32_t score, short dir) {
  EDGE_t *new = EDGE_malloc();
  int ct;

  edge_index ++;
  new->index = edge_index;
  new->ele1_info = ele1_info;
  new->ele2_info = ele2_info;
  new->type = type;
  new->score = score;
  new->direction = dir;
  /*ct = count_edge_nodes(ele1_info->ele->edges);
  if (ct != ele1_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  trouble before inserting edge\n");
    fflush(log_file);
    exit(4);
  }*/
  ele1_info->ele->edge_no ++;
  insert_edge(&ele1_info->ele->edges, new);
  /*ct = count_edge_nodes(ele1_info->ele->edges);
  if (ct != ele1_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  trouble after inserting edge\n");
    fflush(log_file);
    exit(4);
  }*/
  if (ele1_info->index != ele2_info->index) {
    /*ct = count_edge_nodes(ele2_info->ele->edges);
    if (ct != ele2_info->ele->edge_no) {
      err_no ++;
      fprintf(log_file, "error:  trouble before inserting edge\n");
      fflush(log_file);
      exit(4);
    }*/
    ele2_info->ele->edge_no ++;
    insert_edge(&ele2_info->ele->edges, new);
    /*ct = count_edge_nodes(ele2_info->ele->edges);
    if (ct != ele2_info->ele->edge_no) {
      err_no ++;
      fprintf(log_file, "error:  trouble after inserting edge\n");
      fflush(log_file);
      exit(4);
    }*/
  }
}




void adjust_edge_tree(ELE_INFO_t *ele_info) {
  EDGE_t **edge_array;
  int ct, ct1;

  /*ct = count_edge_nodes(ele_info->ele->edges);
  if (ct != ele_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  trouble before adjusting the edge tree: ele %d %d %d\n", ele_info->index, ele_info->ele->edge_no, ct);
    fflush(log_file);
    exit(4);
  }
  ct = count_total_edges(ele_info->ele->edges);
  if (ct != ele_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  illegitimate edge in ele %d\n", ele_info->index);
    fflush(log_file);
    exit(5);
  }*/
  ct = ele_info->ele->edge_no;
  edge_array = (EDGE_t **) malloc(ct*sizeof(EDGE_t *));
  ct1 = charge_edge_array(edge_array, ele_info->ele->edges, 0);
  if (ct1 != ct) {
    err_no ++;
    fprintf(log_file, "error:  trouble charging the edge array in ele %d: %d charged, %d expected\n", ele_info->index, ct1, ct);
    fflush(log_file);
    exit(4);
  }
  edge_tree_free(&ele_info->ele->edges);
  build_edge_tree(&ele_info->ele->edges, edge_array, 0, ct-1);
  /*ct = count_edge_nodes(ele_info->ele->edges);
  if (ct != ele_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  trouble adjusting the edge tree: ele %d %d %d\n", ele_info->index, ele_info->ele->edge_no, ct);
    fflush(log_file);
    exit(4);
  }*/
  free(edge_array);
}




int charge_edge_array(EDGE_t **edge_array, EDGE_TREE_t *rt, int pos) {
  if (rt->l) pos = charge_edge_array(edge_array, rt->l, pos);

  *(edge_array+pos) = rt->to_edge;
  pos ++;

  if (rt->r) pos = charge_edge_array(edge_array, rt->r, pos);

  return pos;
}




/******************************************
 * functions for building the consis_tree *
 ******************************************/


// RMH - debug function to view consistency tree
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
    printf("c_tree: level=%d [%s:%d:%d-%d]",level, child_rt->to_image->frag.seq_name, child_rt->to_image->to_msp->direction,  child_rt->to_image->frag.lb, child_rt->to_image->frag.rb);
    sib_rt = child_rt->sib;
    while ( sib_rt != NULL )
    {
      if ( sib_rt->children == NULL )
        printf(", image %s:%d:%d-%d", sib_rt->to_image->frag.seq_name, sib_rt->to_image->to_msp->direction, sib_rt->to_image->frag.lb, sib_rt->to_image->frag.rb);
      else
        printf(", image %s:%d:%d-%d(*)", sib_rt->to_image->frag.seq_name, sib_rt->to_image->to_msp->direction, sib_rt->to_image->frag.lb, sib_rt->to_image->frag.rb);
      sib_rt = sib_rt->sib;
    }
    printf("\n");
    level++;
    child_rt = child_rt->children;
  }
  return 0;
}

// RMH: Insert images into a tree
//      Unique images ( e.g non overlapping ) get placed on
//      a distinct level in the tree whereas overlapping images
//      are placed on the same level.
int consis_tree_build(IMG_NODE_t *rt, IMAGE_t *im, int prequal) {
  int sum=0;
  IMG_NODE_t *nex_rt, *node;
  // RMH: If first call ( prequal = 1 ) or no significant overlap
  if (prequal || consis(rt->to_image, im, CUTOFF2)) {
    //if ( !prequal )
    //  printf("consis_tree_build: Found consis: im %s:%d-%d and rt %s:%d-%d\n", im->frag.seq_name, im->frag.lb, im->frag.rb, rt->to_image->frag.seq_name, rt->to_image->frag.lb, rt->to_image->frag.rb);
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

// RMH: Given two images:
//          return 1 if:
//             - Both images are mapped to the same element pair
//           and
//             - Both images are in the same direction
//           and
//             - Both image and it's complementary element
//               pair do not overlap significantly ( currently <10% )
int consis(IMAGE_t *i1, IMAGE_t *i2, float cutoff) {
  int res = 0;
  IMAGE_t *ip1 = partner(i1), *ip2 = partner(i2);
  if (i1->ele_info->index == i2->ele_info->index &&
      ip1->ele_info->index == ip2->ele_info->index &&
      i1->to_msp->direction == i2->to_msp->direction) {
    if (i1->to_msp->direction == 1) {
      // RMH: If either image or its partners start at the same
      //      position ( regardless of sequence name? ) they probably
      //      overlap by too much.  This looks like an optimization
      //      to avoid calling sing_cov.  But what if they are not
      //      the same sequence?  I guess that doesn't happen very
      //      often...still it probably should be checked.  For instance
      //      what if we use RECON to compare seqeunces representing elements
      //      of various families.  Then the sequences would be short enough
      //      to perhaps start at the same # even by chance without checking
      //      the sequence name.
      //RMH: Testing why this condition was used
      //if ((i1->frag.lb - i2->frag.lb)*(ip1->frag.lb - ip2->frag.lb) == 0)
      //  if ( i1->frag.seq_name != i2->frag.seq_name )
      //    printf("*****OOPs...this probably should have returned 1!\n");
      if ((i1->frag.lb - i2->frag.lb)*(ip1->frag.lb - ip2->frag.lb) > 0) {
        // RMH: If overlap between two images is >= cutoff * the length
        //      of either sequence. NOTE: sing_cov does check sequence
        //      names are the same.
        //      Here cutoff is 0.9.  Therefore below we return 1 if
        //      the two images ( and their pairs ) overlap less than
        //      0.1 * the length of either.
	if (!sing_cov(&i1->frag, &i2->frag, 1.0-cutoff) &&
            !sing_cov(&ip1->frag, &ip2->frag, 1.0-cutoff)) {
          // Do not overlap significantly for either element
	  res = 1;
	}
      }
    } else {
      if ((i1->frag.lb - i2->frag.lb)*(ip1->frag.rb - ip2->frag.rb) < 0) {
	if (!sing_cov(&i1->frag, &i2->frag, 1.0-cutoff) && !sing_cov(&ip1->frag, &ip2->frag, 1.0-cutoff)) {
          // Do not overlap significantly for either element
	  res = 1;
	}
      }
    }
  }
  return res;
}



IMG_NODE_t ** node_entry(IMG_NODE_t **node_pp) {
  if (*node_pp != NULL) {
    return node_entry(&(*node_pp)->sib);
  }
  return node_pp;
}




void consis_tree_free(IMG_NODE_t *rt) {
  if (rt->sib != NULL) {
    consis_tree_free(rt->sib);
  } 
  if (rt->children != NULL) {
    consis_tree_free(rt->children);
  }
  free(rt);
}




/*****************************************
 * functions for parsing the consis_tree *
 *****************************************/

/*
  I. function for finding primary images from the consis_tree
*/

#if 0
int find_prim(IMG_NODE_t *nd, float cutoff, int32_t score, int32_t hist, short which, int32_t *sc, short *d) {
  IMG_NODE_t *nex_node;
  int sum = 0;
  short level, further = 0;
  IMAGE_t *img_partner;

  if (!hist) level = 1;
  hist += nd->to_image->frag.rb - nd->to_image->frag.lb;
  
  if (level && hist) { /* if first image in a path isn't EOE, forget the path */
    if (nd->to_image->frag.lb - nd->to_image->ele_info->ele->frag.lb < MARGIN) {
      which = 1;
    }
    img_partner = partner(nd->to_image);
    if (nd->to_image->to_msp->direction == 1) {
      if (img_partner->frag.lb - img_partner->ele_info->ele->frag.lb < MARGIN) {
	which += 2;
      }
    } else {
      if (img_partner->frag.rb - img_partner->ele_info->ele->frag.rb > -MARGIN) {
	which += 2;
      }
    }
    if (!which) return 0;
  }

  if (nd->children != NULL) {
    nex_node = nd->children;
    while (nex_node != NULL) {
      sum += find_prim(nex_node, cutoff, score, hist, which, sc, d);
      nex_node = nex_node->sib;
    }
    if (sum && hist) {
      nd->to_image->to_msp->stat = 'p'; /* p stands for primary, NOT partial */
    }
  } else { /* the last image goes to EOE on the same element as the first image in the path */
    if ((which == 1 || which == 3) && nd->to_image->frag.rb - nd->to_image->ele_info->ele->frag.rb > -MARGIN  && (float) hist/(nd->to_image->ele_info->ele->frag.rb - nd->to_image->ele_info->ele->frag.lb) > cutoff) {
      further = 1;
    } else {
      img_partner = partner(nd->to_image);
      if ((which == 2 || which == 3) && (float) hist/(img_partner->ele_info->ele->frag.rb - img_partner->ele_info->ele->frag.lb) > cutoff) {
	if (nd->to_image->to_msp->direction == 1) {
	  if (img_partner->frag.rb - img_partner->ele_info->ele->frag.rb > -MARGIN) further = 1;
	} else {
	  if (img_partner->frag.lb - img_partner->ele_info->ele->frag.lb < MARGIN) further = 1;
	}
      }
    }
    if (further) {
      nd->to_image->to_msp->stat = 'p';
      score += nd->to_image->to_msp->score;
      if (score > *sc) {
	*sc = score;
	*d = nd->to_image->to_msp->direction;
      }
      sum = 1;
    } else {sum = 0;}
  }

  return sum;
}
#endif



// Initially called with
//       end1 = ele->frag.lb, and end2 = -1
//       efl1 = 0, efl2 = 0
//       al1 = 0, al2 = 0,
//       score = 0
// Making a guess here...but int32_t is probably not adequate in some instances for al1/al2 and probably
// has caused problems where overflow occurs.
int find_prim(IMG_NODE_t *nd, float cutoff, int32_t end1, int32_t end2, int32_t efl1, int32_t efl2, int32_t al1, int32_t al2, int32_t score, int *pmarkp, int32_t *sc, short *d) {
  int sum = 0, mark=0;
  int32_t skip1, skip2, len1, len2;
  IMAGE_t *ipt;

  if (nd->sib) sum += find_prim(nd->sib, cutoff, end1, end2, efl1, efl2, al1, al2, score, pmarkp, sc, d);

  ipt = partner(nd->to_image);
  // end2 is the partner element's right bound ( positive strand )
  if (end2 < 0) { /*first alignment in the group*/
    if (nd->to_image->to_msp->direction == 1) end2 = ipt->ele_info->ele->frag.lb;
    else end2 = ipt->ele_info->ele->frag.rb;
  }

  // skip1 is the image's left bound - element's left bound
  //    - So == 0 if they are the same, <0 if the image extends beyond the element
  //      and >0 if it starts later than the element's left bound.
  // skip2 is the partner element's right bound - the partner image's right bound
  //    - So == 0 if they are the same, <0 if the image extends past the bound, and
  //      >0 if it doesn't reach the end.
  skip1 = nd->to_image->frag.lb - end1;
  if (nd->to_image->to_msp->direction == 1) skip2 = ipt->frag.lb - end2;
  else skip2 = end2 - ipt->frag.rb;

  // If we *start* with a gap **AND** our partner *ends* with a gap
  // ( > 10bp ) add the gaps to the efls.
  //
  //   |-----ele1--------|        |------ele2---------|
  //     >10 |--img1-...           ....-----img2-| >10
  //
  if (skip1>10 && skip2>10) {
    efl1 += skip1;
    efl2 += skip2;
  }
  // Add the length of the actual images to the efls
  len1 = nd->to_image->frag.rb - nd->to_image->frag.lb;
  len2 = ipt->frag.rb - ipt->frag.lb;
  efl1 += len1;
  efl2 += len2;

  // Add the length of the actual images to the als
  al1 += len1;
  al2 += len2;

  // Score is simply the average length of the alignment * the % identity
  score += ((int32_t) nd->to_image->to_msp->iden)*(len1+len2)/2;

  if (nd->children) {
    end1 = nd->to_image->frag.rb;
    if (nd->to_image->to_msp->direction == 1) end2 = ipt->frag.rb;
    else end2 = ipt->frag.lb;
    sum += find_prim(nd->children, cutoff, end1, end2, efl1, efl2, al1, al2, score, &mark, sc, d);
  } else { /*last alignment in group*/
    skip1 = nd->to_image->ele_info->ele->frag.rb - nd->to_image->frag.rb;
    if (nd->to_image->to_msp->direction == 1) skip2 = ipt->ele_info->ele->frag.rb - ipt->frag.rb;
    else skip2 = ipt->frag.lb - ipt->ele_info->ele->frag.lb;
    // If we *end* with a gap **AND** our partner *starts* with a gap
    // ( > 10bp ) add the gaps to the efls.
    //
    //   |-----ele1--------|        |------ele2---------|
    //    ...--img1--| >10            >10 |-----img2-...
    //
    if (skip1>10 && skip2>10) {
      efl1 += skip1;
      efl2 += skip2;
    }


    //   |-----ele1--------|        |------ele2---------|
    //
    // Case 1:
    //     >10 |-img1-| >10           >10 |-img2-| >10
    //   efl1 = ~length(ele1)         efl1 = ~length(ele2)
    //
    // Case 2:
    //   |--img1---|  >10             >10      |-img2---|
    //   efl1 = ~length(ele1)         efl1 = ~length(ele2)
    //
    // Case 3:
    //    >10    |--img1---|        |-img2---| > 10
    //   efl1 = ~length(ele1)         efl1 = ~length(ele2)
    //
    // Case 4:  Good primary candidate
    //    >10    |--img1---|         >10      |--img2---|
    //   efl1 = ~length(img1)         efl1 = ~length(img2)
    //
    // Case 5:  Good primary candidate
    //   |--img1---|   >10          |--img2---|    >10
    //   efl1 = ~length(img1)         efl1 = ~length(img2)
    //
    /*if (1.0*al1/efl1 > cutoff || 1.0*al2/efl2 > cutoff) {*/
    RLOG_DBG("find_prim: al1=%d al2=%d efl1=%d efl2=%d:  al/ef %f, %f   ef-al %d, %d  ele:%d-%d, ptn:%d-%d\n", al1, al2, efl1, efl2, (1.0*al1/efl1), (1.0*al2/efl2), (efl1-al1), (efl2-al2), nd->to_image->ele_info->ele->frag.lb, nd->to_image->ele_info->ele->frag.rb, ipt->ele_info->ele->frag.lb, ipt->ele_info->ele->frag.rb);
    if ( (1.0*al1/efl1 > cutoff || 1.0*al2/efl2 > cutoff) && (efl1-al1 < 30 || efl2-al2 < 30) ) {
      sum = 1;
      mark = 1;
      if ( (al1+al2) == 0 )
      {
         // RMH: A divide by zero error has occured a few times at this point.
         //      I looked extensively through the code and cannot determine what
         //      is causing this -- although I have a reliable test case to
         //      explore this further.  For now this is the only workaround.
         //      If al1 & al2 are zero then simply set them to 1.
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
      nd->to_image->to_msp->stat  = 'p';
      *pmarkp = 1;
  }
  return sum;
}






/************************
 ************************
 ***                  ***
 ***      OUTPUT      ***
 ***                  ***
 ************************
 ************************/





void combo_output(ELE_INFO_t *ele_info) {
  char *ele_name;
  FILE *ele_fp;
  ELE_DATA_t *cur_ele_data;

  fprintf(combo, "%d %s %d %d \n", ele_info->index, ele_info->ele->frag.seq_name, ele_info->ele->frag.lb, ele_info->ele->frag.rb);
  fflush(combo);

#if 0
  /* no ele_read_in() needed, 'cuz where this function is called, both the combo and its offsprings are in the memory */
  ele_name = (char *) malloc(30*sizeof(char));
  sprintf(ele_name, "tmp2/clan/combos/e%d", ele_info->index);
  ele_fp = fopen(ele_name, "w");

  fprintf(ele_fp, "ele %d\n", ele_info->index);
  fprintf(ele_fp, "offsprings \n");
  cur_ele_data = ele_info->ele->redef;
  while (cur_ele_data != NULL) {
    if (cur_ele_data->ele_info->ele->direction == 1) fprintf(ele_fp, "%6d %10s %8d %8d\n", cur_ele_data->ele_info->index, cur_ele_data->ele_info->ele->frag.seq_name, cur_ele_data->ele_info->ele->frag.lb, cur_ele_data->ele_info->ele->frag.rb);
    else fprintf(ele_fp, "%6d %10s %8d %8d\n", cur_ele_data->ele_info->index, cur_ele_data->ele_info->ele->frag.seq_name, cur_ele_data->ele_info->ele->frag.rb, cur_ele_data->ele_info->ele->frag.lb);
    cur_ele_data = cur_ele_data->next;
  }

  fclose(ele_fp);
  free(ele_name);
  /*  remove_ele(ele_info);*/
#endif
}





void obs_output(ELE_INFO_t *ele_info) {
  char *ele_name;
  FILE *ele_fp;

  fprintf(obs, "%d %s %d %d \n", ele_info->index, ele_info->ele->frag.seq_name, ele_info->ele->frag.lb, ele_info->ele->frag.rb);
  fflush(obs);

#if 0
  ele_name = (char *) malloc(30*sizeof(char));
  sprintf(ele_name, "tmp2/clan/obsolete/e%d", ele_info->index);
  ele_fp = fopen(ele_name, "w");
  fprint_ele_obs(ele_fp, ele_info);
  fclose(ele_fp);
  free(ele_name);
  /*  remove_ele(ele_info);*/
#endif
}




void fprint_ele_obs(FILE *fp, ELE_INFO_t *ele_info) {
  int i;
  BD_t *cur_bd;
  ELEMENT_t *ele=ele_info->ele;

  fprintf(fp, "ele %d\n", ele_info->index);
  fprintf(fp, "%s %d %d \n", ele->frag.seq_name, ele->frag.lb,  ele->frag.rb);
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






/**************************************
 **************************************
 *                                    *
 *  comparison fucntions for qsort()  *
 *                                    *
 **************************************
 **************************************/






int frag_cmp(const void *i1, const void *i2) {
  int res = (*(IMG_DATA_t **)i1)->to_image->frag.seq_name - (*(IMG_DATA_t **)i2)->to_image->frag.seq_name; /*strncmp((*(IMG_DATA_t **)i1)->to_image->frag.seq_name, (*(IMG_DATA_t **)i2)->to_image->frag.seq_name, NAME_LEN)*/
  if (res == 0) {
    res = (*(IMG_DATA_t **)i1)->to_image->frag.lb - (*(IMG_DATA_t **)i2)->to_image->frag.lb;
    if (res == 0) {
      res = (*(IMG_DATA_t **)i1)->to_image->frag.rb - (*(IMG_DATA_t **)i2)->to_image->frag.rb;
    }
  }
  return res;
}

// RMH: Sort by: difference in element identifiers, ascending
//               msp_direction
//               left bound
//               right bound
//
int partner_cmp(const void *i1, const void *i2) {
  int res = partner(*((IMAGE_t **)i1))->ele_info->index - partner(*((IMAGE_t **)i2))->ele_info->index;
  if (res == 0) {
    res = (*((IMAGE_t **)i1))->to_msp->direction - (*((IMAGE_t **)i2))->to_msp->direction;
  }
  if (res == 0) {
    res = (*((IMAGE_t **)i1))->frag.lb - (*((IMAGE_t **)i2))->frag.lb;
  }
  if (res == 0) {
    res = (*((IMAGE_t **)i1))->frag.rb - (*((IMAGE_t **)i2))->frag.rb;
  }
  return res;
}


int CP_cmp(const void *cp1, const void *cp2) {
  return (*((CP_t **) cp1))->cp - (*((CP_t **) cp2))->cp;
}


int BD_cmp(const void *bd1, const void *bd2) {
  return (*((BD_t **) bd1))->bd - (*((BD_t **) bd2))->bd;
}


int fam_cmp(const void *fd1, const void *fd2) {
  return (*((FAM_DATA_t **) fd1))->to_family->index - (*((FAM_DATA_t **) fd2))->to_family->index;
}
