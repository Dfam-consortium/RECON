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
 * Source organisation
 *   eleredef.c      -- this file: main(), reporting, clan management
 *   redef_boundary.c -- CP/BD/PCP clustering
 *   redef_edges.c   -- edge building and consistency tree
 *   redef_dissect.c -- element construction, splitting, combo bookkeeping
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
#include "redef_boundary.h"
#include "redef_edges.h"
#include "redef_dissect.h"

/* ---- Per-program log level storage (required by recon_log.h) ---- */
int   recon_log_level = RECON_LOG_INFO;
FILE *recon_log_fp    = NULL;


/* ---- Forward declarations for functions defined in this file ---- */
void report_cts(void);
void report_redef_stat(void);

void general_ele_redef(ELE_INFO_t *, IMAGE_t **);
void build_local_network(ELE_INFO_t *, ELE_DATA_t **, ELE_DATA_t **, IMAGE_t **);
void recruit(ELE_INFO_t *, EDGE_TREE_t *, ELE_DATA_t **, IMAGE_t **);
void cruise_local_net(ELE_DATA_t *, IMAGE_t **);
void local_ele_redef(ELE_INFO_t *, IMAGE_t **, int *);
void dissolve_local_network(ELE_DATA_t **);
void dismiss_element(ELE_INFO_t *);
void remove_ele(ELE_INFO_t *);




int main (int argc, char *argv[]) {
  ELE_INFO_t *cur_ele_info;
  int i, ei, rounds=0, start;
  char line[35], stat;
  short fu, to_march;
  clock_t start1, end;
  double cpu_time_used=0.0;
  double ele_defTIME=0.0, dissectTIME=0.0;
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

  /* May or may not exist; callers check before using */
  edge_no  = fopen("summary/naive_edge_no", "r");
  size_list = fopen("tmp/size_list", "r");
  new_stat  = fopen("tmp2/redef_stat", "r");

  new_msps = fopen("summary/new_msps", "a");
  if (!new_msps) { printf("Can not open summary/new_msps for writing! Exiting\n"); exit(1); }
  unproc = fopen("summary/unproc", "a");
  if (!unproc)   { printf("Can not open summary/unproc for writing! Exiting\n");   exit(1); }
  combo  = fopen("summary/combo", "a");
  if (!combo)    { printf("Can not open summary/combo for writing! Exiting\n");    exit(1); }
  obs    = fopen("summary/obsolete", "a");
  if (!obs)      { printf("Can not open summary/obsolete for writing! Exiting\n"); exit(1); }
  log_file = fopen("tmp2/log", "a");
  if (!log_file) { printf("Can not open tmp2/log for writing! Exiting\n");         exit(1); }

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

  /* Pre-allocate element info array; new elements spill to ele_info_data list */
  ele_array_size = 2*ele_ct;
  all_ele = (ELE_INFO_t **) malloc(ele_array_size*sizeof(ELE_INFO_t *));
  for (i=0; i<ele_array_size; i++) {
    *(all_ele+i) = ele_info_init(i+1);
  }
  ele_info_data = NULL;

  /* Remove large tandem repeats from consideration */
  if (size_list && !new_stat) outthrow_big_tandems(size_list);
  fclose(unproc);

  /* Restore element state from a previous partial run (not used in normal flow) */
  if (new_stat) {
    printf("Reading stat file redef_stat...\n");
    ele_ct = 0;
    while (fgets(line, 35, new_stat)) {
      ele_ct ++;
      /* sscanf doesn't work reliably here; parse manually */
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
  msp_left   = 0;
  msp_ct     = 0;
  edge_ct    = 0;
  edge_in_mem = 0;
  edge_left  = 0;
  files_read = 0;
  err_no     = 0;

  /* Re-define elements using the syntopy algorithm, and build edges */
  img_ptr = (IMAGE_t **) malloc(MAX_IMG*sizeof(IMAGE_t *));
  if (!img_ptr) {
    printf("eleredef: Error! Could not allocate memory for img_ptr: "
           "%zu bytes requested\n", MAX_IMG*sizeof(IMAGE_t *));
    exit(-1);
  }

  to_march = 1;
  start1 = clock();
  while (to_march) {
    to_march = 0;
    rounds ++;
    /* Process the pre-allocated element array */
    for (i=start; i<ele_ct && i<ele_array_size; i++) {
      fprintf(log_file, "Evaluating definition of element %d\n",
              (*(all_ele+i))->index);
      if ((*(all_ele+i))->stat == 'z' ||
          (*(all_ele+i))->stat == 'w' ||
          (*(all_ele+i))->stat == 't') {
	to_march = 1;
	general_ele_redef(*(all_ele+i), img_ptr);
      }
    }

    start = 0;

    /* Process overflow elements from the ele_info_data linked list */
    cur_ele_info = ele_info_data;
    while(cur_ele_info) {
      fprintf(log_file, "evaluating definition of element %d\n",
              cur_ele_info->index);
      if (cur_ele_info->stat == 'z' ||
          cur_ele_info->stat == 'w' ||
          cur_ele_info->stat == 't') {
	to_march = 1;
	general_ele_redef(cur_ele_info, img_ptr);
      }
      cur_ele_info = cur_ele_info->next;
    }
    end = clock()-start1;
    cpu_time_used = ((double)(end - start1)) / CLOCKS_PER_SEC;
  }
  cpu_time_used = ((double)(end - start1)) / CLOCKS_PER_SEC;

  report_cts();
  report_redef_stat();
  free(img_ptr);

  fprintf(log_file, "total numbers: %d elements, %d msps, %d edges\n",
          ele_ct, msp_index+1, edge_index+1);
  fprintf(log_file, "%d rounds, %d files read, %d msps seen, %d edges seen\n",
          rounds, files_read, msp_ct, edge_ct);
  fprintf(log_file, "%d errors, %d msps and %d edges left in memory, \n",
          err_no, msp_left, edge_left);
  RLOG_DBG("General_ele_redef %f , %f , %f \n",
           cpu_time_used, ele_defTIME, dissectTIME);
  RLOG_DBG("total numbers: %d elements, %d msps, %d edges\n",
           ele_ct, msp_index+1, edge_index+1);
  RLOG_DBG("%d rounds, %d files read, %d msps seen, %d edges seen\n",
           rounds, files_read, msp_ct, edge_ct);
  RLOG_DBG("%d errors, %d msps and %d edges left in memory, \n",
           err_no, msp_left, edge_left);
  fflush(log_file);
  fclose(log_file);

  ele_db_close();

  exit(0);
}




/* =========================================================================
 * Reporting
 * ========================================================================= */

/* report_cts  --  write the final element count to summary/redef_ele_no */
void report_cts() {
  FILE *fp;

  fp = fopen("summary/redef_ele_no", "w");
  fprintf(fp, "%d\n", ele_ct);
  fclose(fp);
}


/* report_redef_stat  --  snapshot all element statuses to tmp2/redef_stat
 *
 * Also updates summary/redef_msp_no and summary/naive_edge_no so that a
 * subsequent eleredef restart can pick up where this run left off.
 */
void report_redef_stat() {
  int i, ele_march = ele_ct < ele_array_size ? ele_ct : ele_array_size;
  FILE *redef_stat, *fp;
  ELE_INFO_t *cur_ele_info = ele_info_data;

  redef_stat = fopen("tmp2/redef_stat", "w");
  for (i=0; i<ele_march; i++) {
    fprintf(redef_stat, "%d %c %d\n",
            (*(all_ele+i))->index,
            (*(all_ele+i))->stat,
            (*(all_ele+i))->file_updated);
  }
  while(cur_ele_info) {
    fprintf(redef_stat, "%d %c %d\n",
            cur_ele_info->index,
            cur_ele_info->stat,
            cur_ele_info->file_updated);
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




/* =========================================================================
 * Clan (local network) management
 *
 * A "clan" is the BFS subgraph centred on one unprocessed element and
 * extending DEPTH hops along edges.  Processing a clan means:
 *   1. build_local_network() -- BFS to collect all members
 *   2. cruise_local_net()    -- redefine each member in turn
 *   3. dissolve_local_network() -- write results and release memory
 * ========================================================================= */


/* general_ele_redef  --  entry point for processing one element's clan
 *
 * Status transitions driven here:
 *   'z' -> edges_and_cps() -> 't'
 *   't' -> local_ele_redef() -> 'v' or 'X'
 *   'w' -> combo_update() -> 'X'
 */
void general_ele_redef(ELE_INFO_t *ele_info, IMAGE_t **img_ptr) {
  int i;
  ELE_DATA_t *local_net=NULL, *local_net_tail=NULL;
  RLOG_DBG("general_ele_redef: ele %d\n", ele_info->index);
  if (!ele_info->ele) ele_read_in(ele_info, 1);

  if (!ele_info->ele->img_no) {
    /* element lost all images during a partner's re-evaluation */
    combo_update(ele_info);
    remove_ele(ele_info);
    ele_cleanup(&ele_info->ele);
  } else {
    fprintf(new_msps, "new clan for ele %d\n", ele_info->index);
    fprintf(combo,    "new clan for ele %d\n", ele_info->index);
    fprintf(obs,      "new clan for ele %d\n", ele_info->index);

    clan_ct ++;
    clan_size = 0;
    clan_core_size = 0;

    fprintf(log_file, "new clan: %d for ele %d\n", clan_ct, ele_info->index);

    build_local_network(ele_info, &local_net, &local_net_tail, img_ptr);

    if (recon_log_level >= RECON_LOG_DEBUG) print_all_eles_GML();

    fprintf(log_file, "clan size: %d, clan core size: %d\n",
            clan_size, clan_core_size);

    cruise_local_net(local_net, img_ptr);

    if (recon_log_level >= RECON_LOG_DEBUG) print_all_eles_GML();

    fflush(new_msps);
    fflush(combo);
    fflush(obs);
    dissolve_local_network(&local_net);
  }
}


/* build_local_network  --  BFS to collect all elements within DEPTH hops
 *
 * Seeds the queue with ele_info, then for each queued element:
 *   - calls edges_and_cps() (if status is 'z') to build its edge set
 *   - calls recruit() to enqueue all not-yet-queued neighbours
 *
 * On return, the linked list *net_p contains all clan members.
 */
void build_local_network(ELE_INFO_t *ele_info, ELE_DATA_t **net_p,
                         ELE_DATA_t **net_tail_p, IMAGE_t **img_ptr) {
  ELEMENT_t *ele;
  ELE_DATA_t *que;
  RLOG_DBG("build_local_network: seed ele %d\n", ele_info->index);

  que = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
  que->ele_info = ele_info;
  que->next = NULL;
  *net_p = que;
  *net_tail_p = que;

  ele_info->ele->l_hold = 1;

  while (que) {
    clan_size ++;
    RLOG_DBG("build_local_network: BFS queue entry ele %d\n",
             que->ele_info->ele->index);
    if (que->ele_info->ele->l_hold <= DEPTH) {
      clan_core_size ++;
      if (que->ele_info->stat == 'z') edges_and_cps(que->ele_info, img_ptr);
      if (que->ele_info->ele->edges)
        recruit(que->ele_info, que->ele_info->ele->edges, net_tail_p, img_ptr);
    }
    que = que->next;
  }
}


/* recruit  --  recursively enqueue unvisited edge partners of ele_info
 *
 * Walks the edge BST in-order and, for each partner not yet in the queue
 * (l_hold == 0), sets its depth and appends it to the tail of the queue.
 */
void recruit(ELE_INFO_t *ele_info, EDGE_TREE_t *rt,
             ELE_DATA_t **net_tail_p, IMAGE_t **img_ptr) {
  ELE_INFO_t *epi;
  ELE_DATA_t *member;
  RLOG_DBG("recruit: ele %d\n", ele_info->index);
  if (rt->l) recruit(ele_info, rt->l, net_tail_p, img_ptr);

  epi = linked_ele(ele_info, rt->to_edge);
  RLOG_DBG("recruit: partner ele %d\n", epi->index);
  if (!epi->ele) ele_read_in(epi, 1);
  if (!epi->ele->l_hold) {
    epi->ele->l_hold = ele_info->ele->l_hold + 1;
    if (epi->stat == 'v' && epi->ele->l_hold < DEPTH) epi->ele->l_hold = DEPTH;
    member = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
    member->ele_info = epi;
    member->next = NULL;
    (*net_tail_p)->next = member;
    *net_tail_p = member;
  }

  if (rt->r) recruit(ele_info, rt->r, net_tail_p, img_ptr);
}


/* cruise_local_net  --  redefine all elements in the local network
 *
 * Iterates the queue repeatedly until no element changes status (to_march).
 * Each pass calls local_ele_redef() on every element within DEPTH hops.
 */
void cruise_local_net(ELE_DATA_t *local_net, IMAGE_t **img_ptr) {
  ELE_DATA_t *que;
  int to_march = 1;

  while (to_march) {
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


/* local_ele_redef  --  redefine one element within the local network
 *
 * If the element has child (redef) elements, recurses into them first.
 * Otherwise, if the element has PCPs, calls ele_redef(); if not, marks
 * it 'v' (done).  Sets *march_p=1 if any redefinition occurred.
 */
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


/* dissolve_local_network  --  write results and release all clan memory
 *
 * Calls dismiss_element() on every member (which writes each element's
 * updated file and frees its ELEMENT_t), then verifies that no MSPs or
 * edges were left dangling in memory.
 */
void dissolve_local_network(ELE_DATA_t **net_p) {
  ELE_DATA_t *que;

  que = *net_p;
  while (que) {
    dismiss_element(que->ele_info);
    que = que->next;
  }

  if (msp_in_mem) {
    err_no ++;
    fprintf(log_file,
            "error:  error in bookkeeping: %d msps total, %d seen, "
            "%d left in memory\n",
            msp_index+1, msp_ct, msp_in_mem);
    fflush(log_file);
    msp_left += msp_in_mem;
  }

  if (edge_in_mem) {
    err_no ++;
    fprintf(log_file,
            "error:  error in bookkeeping: %d edges total, %d seen, "
            "%d left in memory\n",
            edge_index+1, edge_ct, edge_in_mem);
    fflush(log_file);
    edge_left += edge_in_mem;
    edge_in_mem = 0;
  }

  if (err_no) exit(3);

  ele_data_free(net_p);
}


/* dismiss_element  --  write out one clan element and free its ELEMENT_t
 *
 * Recursively processes child elements first.  If the element has no
 * images remaining it is marked 'X' and recorded as combo or obsolete;
 * otherwise it is written to disk by ele_write_out().
 */
void dismiss_element(ELE_INFO_t *ele_info) {
    ELE_DATA_t *new_ele_data;

    if (ele_info->ele->redef) {
      new_ele_data = ele_info->ele->redef;
      while (new_ele_data) {
        dismiss_element(new_ele_data->ele_info);
	new_ele_data = new_ele_data->next;
      }
    }
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


/* remove_ele  --  delete the on-disk file for a dismissed element */
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
