/*
 * eleredef.c - Element Redefinition
 *
 *   Given elements initially defined by single linkage clustering,
 *   redefine (breakup) into subelements at points where an
 *   aggregation of image endpoints exist.
 *
 * Author: Zhirong Bao
 * Minor modifications by: Robert Hubley, Institute for Systems Biology
 *
 * RMH Notes:
 *   Element info stat:  'z', 't', 'v', 'w', 'y', 'x' and 'X'
 *      'z' - appears to be initial state
 *      'v' - I have image-end-selection rule
 *      'w' - Secondary edges are defined
 *      'y' - Already been traversed has no neighbors
 *      't' - set at the end of edges_and_cps()
 *      'X' - appears to indicate deleted or 'dismissed' element
 *      PCP and TBD...need defining
 *          - breakup elements
 *  flimg_no = Full Length Image Number (count)
 */
#include <fcntl.h>
#include <sys/types.h>
#include <dirent.h>
#include <time.h>

#include "ele.h"
#include "eleredef.h"
#include "seqlist.h"
#include "treeview.h"

#define CUTOFF1 0.5
#define CUTOFF2 0.9
#define MAX_IMG 1200
// RMH: "The minimum length of an element, for example,
//       the maximal length that we expect the pairwise
//       alignment tool to spuriously extend by chance
//       from a true element boundary, used in the image
//       end selection rule and the element reevaluation
//       and update procedure."
#define TOO_SHORT 30
// RMH: Not sure yet.  Could be:
//      "The ratio cutoff for splitting an element at a
//       given position used in the element reevaluation
//       and update procedure."
//       It is a multiplier of the number of images that
//       span a putative split position.
// Original code had 2
#define FUDGE 2
//#define FUDGE 1.5
//#define FUDGE 1
#define MARGIN 10000
#define FLURRY 10

void report_cts();
void report_redef_stat();

ELE_DATA_t *ele_def(IMG_DATA_t **, float);
void generate_img_tree(ELEMENT_t *);
ELE_INFO_t *new_element();
void add_ele_info(ELE_INFO_t *);

void general_ele_redef(ELE_INFO_t *, IMAGE_t **);
void build_local_network(ELE_INFO_t *, ELE_DATA_t **, ELE_DATA_t **, IMAGE_t **, int);
void recruit(ELE_INFO_t *, EDGE_TREE_t *, ELE_DATA_t **, IMAGE_t **, int);
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
void test_consis_tree();
int consis(IMAGE_t *, IMAGE_t *, float);
IMG_NODE_t **node_entry(IMG_NODE_t **);
int consis_tree_free(IMG_NODE_t *);
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


   // DEBUG
   //test_consis_tree();
   //exit(10);

  /* processing command line */
  if (argc == 1) {
    printf("usage: eleredef seq_list start clan_ct\n where seq_list is the list of sequence names, start is the index of the element to start redefining, and clan_ct is the number to start counting the number of clans.  The latter two are optional.\n");
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

  img_ptr = (IMAGE_t **) malloc(MAX_IMG*sizeof(IMAGE_t *));

  if ( ! img_ptr )
  {
    printf("eleredef: Error! Could not allocate memory for img_ptr: %d bytes requested\n", ( MAX_IMG*sizeof(IMAGE_t *) ) );
    exit(-1);
  }

  //used as a counter to ensure every element gets redefined if applicable -kn
  to_march = 1;
  start1 = clock();
  /* re-define elements using the syntopy algorithm, and build edges */
  while (to_march) {
    to_march = 0;
    rounds++;
    for (i=start; i<ele_ct && i<ele_array_size; i++) {
      fprintf(log_file, "Evaluating definition of element %d\n", (*(all_ele+i))->index);
      fflush(log_file);
      // Check if the element has content -kn
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
    //      array is exceeded.  ....or like in ele_def this is a set of unused elements that need to be reconsidered
    //      once genereal_ele_redef has completed a round?
    cur_ele_info = ele_info_data;
    while(cur_ele_info) {
      fprintf(log_file, "evaluating definition of element %d\n", cur_ele_info->index);
      fflush(log_file);
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
  printf("General_ele_redef %f , %f , %f \n",cpu_time_used, ele_defTIME, dissectTIME);
  printf("total numbers: %d elements, %d msps, %d edges\n", ele_ct, msp_index+1, edge_index+1);
  printf("%d rounds, %d files read, %d msps seen, %d edges seen\n", rounds, files_read, msp_ct, edge_ct);
  printf("%d errors, %d msps and %d edges left in memory, \n", err_no, msp_left, edge_left);
  fflush(log_file);
  fclose(log_file);

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
  printf("generate_img_tree element: %d \n", ele->index);
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
  // Read in stage 1 element from "e#" ( where # is ele_info->index ).
  printf("general_ele_redef: %d\n",ele_info->index);
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
    fflush(new_msps);
    fprintf(combo, "new clan for ele %d\n", ele_info->index);
    fflush(combo);
    fprintf(obs, "new clan for ele %d\n", ele_info->index);
    fflush(obs);
    /* set up the local network */
    clan_ct ++;
    clan_size = 0;
    clan_core_size = 0;

    fprintf(log_file, "new clan: %d for ele %d\n", clan_ct, ele_info->index);
    fflush(log_file);

    // RMH: Build a graph centered on ele_info and extending out 3 degrees
    //      of separation. TODO: Add details of edge labeling.
    //
    //   local_net is a linked list of ELE_INFO_t pointers (ELE_DATA_t) that
    //             is built up by this method.
    build_local_network(ele_info, &local_net, &local_net_tail, img_ptr, DEPTH);

    // RMH: TODO: Describe the in-memory datastructure at this moment.  For
    //            instance, what does local_net look like.  What is ele->PCP store?
    //print_ele_data( local_net );
    //print_local_network(local_net);
    //print_all_eles_GML();
    print_all_eles_GV();

    fprintf(log_file, "clan size: %d, clan core size: %d\n", clan_size, clan_core_size);
    fflush(log_file);

    // redefining elements in the local network 
    //   -- queues up all elements in the local_net linked list
    //      and calls local_ele_redef on each
    cruise_local_net(local_net, img_ptr);

    //print_local_network(local_net);
    //print_all_eles_GML();
    print_all_eles_GV();

    /* clearing up the local network */
    fflush(new_msps);
    fflush(combo);
    fflush(obs);
    dissolve_local_network(&local_net);
  }
}


// Build a graph centered on ele_info and extending out max_depth degrees
//   of separation.  It also generates a list of endpoints for later clustering
//   (CP_t).
//
//      TODO: Add details of edge labeling.
//
//   NOTE: max_depth is currently hard coded to 3.
//
//   local_net is a linked list of ELE_INFO_t pointers (ELE_DATA_t) that
//             is built up by this method.
void build_local_network(ELE_INFO_t *ele_info, ELE_DATA_t **net_p, ELE_DATA_t **net_tail_p, IMAGE_t **img_ptr, int max_depth) {
  ELEMENT_t *ele;
  ELE_DATA_t *que;
  printf("build local network: %d\n", ele_info->index);

  // seed the network with the first element
  que = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
  que->ele_info = ele_info;
  que->next = NULL;
  *net_p = que;
  *net_tail_p = que;

  // In this context l_hold stores the depth of the element in the
  // network.
  ele_info->ele->l_hold = 1;

  // Process the queue of elements
  while (que) {
    clan_size ++;
    printf("queue entry: %d\n", que->ele_info->ele->index);
    if (que->ele_info->ele->l_hold <= max_depth) {
      clan_core_size ++;
      // NOTE: edges_and_cps sets ele_info to 't'
      if (que->ele_info->stat == 'z') edges_and_cps(que->ele_info, img_ptr);
      // Read in elements referenced in initial graph and add to queue ( up to 3 levels deep )
      if (que->ele_info->ele->edges) recruit(que->ele_info, que->ele_info->ele->edges, net_tail_p, img_ptr, max_depth);
    }
    que = que->next;
  }
}


// Given an element and an edge tree, recursively load the elements referenced in the edge tree in
// a depth first manner.  The elements are added to the network if they haven't been added yet.
void recruit(ELE_INFO_t *ele_info, EDGE_TREE_t *rt, ELE_DATA_t **net_tail_p, IMAGE_t **img_ptr, int max_depth) {
  ELE_INFO_t *epi;
  ELE_DATA_t *member;

  // If left child exists, recruit it
  if (rt->l) recruit(ele_info, rt->l, net_tail_p, img_ptr, max_depth);

  // Extract the partner element using the edge information
  epi = linked_ele(ele_info, rt->to_edge);
  printf("Element %d is recruiting element epi: %d\n", ele_info->index, epi->index);

  // If it's not already in-memory, read it in
  if (!epi->ele) ele_read_in(epi, 1);

  // If the element hasn't been added to the network yet, add it
  if (!epi->ele->l_hold) {

    // l_hold, in this context, is the depth of the element in the network.  If it's
    // zero then it hasn't been added to the network yet and should be set to one
    // higher than the current element.
    epi->ele->l_hold = ele_info->ele->l_hold + 1;

    // TODO: Determine when the 'v' state is set
    //   ... Set in ele_redef() or in local_ele_redef()...which appear to be later in the
    //       algorithm.
    //   ... this looks like 'v' causes the element to be terminal even if it's l_hold is
    //       less than the max_depth.
    if (epi->stat == 'v' && epi->ele->l_hold < max_depth) epi->ele->l_hold = max_depth;

    // Allocate memory for a new ELE_DATA_t and add it to the network
    member = (ELE_DATA_t *) malloc(sizeof(ELE_DATA_t));
    member->ele_info = epi;
    member->next = NULL;
    (*net_tail_p)->next = member;
    *net_tail_p = member;
  }

  // If right child exists, recruit it
  if (rt->r) recruit(ele_info, rt->r, net_tail_p, img_ptr, max_depth);
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

    // RMH: Debug
    printf("local_ele_redef(): ele_info->index=%d, ele_info->stat=%c, ", ele_info->index, ele_info->stat);
    if ( ele->redef != NULL )
      printf("ele->redef=*, ");
    else
      printf("ele->redef=Null, ");
    if ( ele->PCP )
      printf("ele->PCP=*\n");
    else
      printf("ele->PCP=Null\n");
    //

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
      printf("the name of the element without a PCP is %d \n", cur_ele->index);
    }
    // If there are multiple TBDs within 10bp(strict) of each other keep the
    // one with the higher support.
    if (cur_ele->TBD) {
      TBD_merge(cur_ele);
    }
    // PCP: Contains all the images for the element (preliminary cluster points?)
    // TBD: 
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
        // FLURRY = 10
        // If the TBD is >10bp from both the left/right edge of the element, dissect
        if (pbd->bd-cur_ele->frag.lb > FLURRY && pbd->bd-cur_ele->frag.rb < -FLURRY) {
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
      printf("PCP_to_TBD stat 1 %ld %ld \n", ele->index , ele->PCP->contributor->index );
  int s = 0, left;
  BD_t *pbd_tmp, *pbd_prev, *pbd, *pbds;
  CP_t *cp;

  /* sort the PCP list according to cp */
  CP_sort(&ele->PCP);
  // Clustering the PCPs into PBDs
  pbds = CP_cluster(ele->PCP);
  /* identify TBD from PBDs */
  /* TBDs are removed from PBDs, what is left in PBDs are those unsuccessful ones */
  pbd_tmp = pbds;
  pbd_prev = NULL;
  while (pbd_tmp != NULL) {
    /* s is the KEY! */

    // How many images span the cutpoint +/- 10 bp ( hardcoded )?
    //  
    // This calculates the number of images that span
    // +/-10bp of this putative split point (PBD),
    // multiplied by FUDGE (2).
    s = span(ele, pbd_tmp->bd);
    printf("Considering bd=%d, support=%d, span(s)=%d\n", pbd_tmp->bd, pbd_tmp->support, s);
    // Give the above info, this indicates that for a PBD to be
    // considered a TBD it must have support >= 2*count_of_images_spanning
    if (pbd_tmp->support >= s) {
        printf("Adding to TBD!\n");
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

//
// RMH:
// The PCP list ( what I assume is always being considered here )
// is a linked list of image endpoints for primary images (edges)
// linking this element with others.  For example
//
// This clusters points that are within 20bp of each other
// and do not have gaps of more than 10bp.  A BD entry which
// contains the average position is returned for each cluster
// along with the number of copies in the cluster (support).
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
      if (cps->cp - last <= 10)        begin = cps;
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


//
// Count how many images span the cut point +/- 10bp
// For some reason return this value scaled by
// "FUDGE" which is currently set to
//
//   Start before cut - 10 and end after cut + 10.
//
//
int span(ELEMENT_t *ele, int32_t cut) {
  /* span requires PCP sorted according to cp_cmp and images sorted according to frag_cmp */
  int left=0, right=0;
  IMG_DATA_t *id;

  id = ele->to_img_data;
  while (id && id->to_image->frag.lb <= cut-10) {
    left ++;
    if (id->to_image->frag.rb <= cut+10) right ++;
    id = id->next;
  }
  return (left-right)*FUDGE;
}




// If there are multiple TBDs within
// 10bp of each other keep the one with
// the higher support.
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
  printf("combo_update before stat: %f \n", ele_info->stat);
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
      printf("combo_update after stat: %f \n", ele_info->stat);
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
//            either E1 or E2. Or...if a set of partial
//            images could be seen as part of full length
//            alignment. Edge stat='p'.
//
//        E1 ----secondary---> E2
//            If there exists an non-full length image connecting
//            them.
//
void edges_and_cps(ELE_INFO_t *ele_info, IMAGE_t **img_ptr) {
    /* the core function in this part, calls everything below */
    int eff_img_ct, i;
    int j, ele_block_start;
    //short start_ele_partner;
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
    short dir = 0;

    // RMH: DEBUG TODO remove
    printf("edges_and_cps: ele_info->index = %d, ele_info->ele->index = %d\n", ele_info->index, ele_info->ele->index);

    /* sort unprocessed images according to their partner element,
     * and then their left bounds to allocate proper amount of memory.
     * When update is 1, it means the element is an offspring of a combo,
     * in which case all the images need to be processed; when update is 0,
     * those images whose partner is 'v' or 'w' should be omitted
     * no self-images are accepted to generate pcps or edges, updating or not
     */

    // If to_image_data is NULL populate it with a linked list version of the image tree (to_img_tree).
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

    // Process image list and count only images that map to other elements ( not ourselves ).
    // Unless the current element is set to update, further restrict the count to images where
    // the partner element is in the initial state ('z').
    eff_img_ct = 0;
    while(cur_img_data != NULL) {
      epi = partner(cur_img_data->to_image)->ele_info;
      //printImage(cur_img_data->to_image);
      //printf("Partner: index = %d\n", epi->index);
      if (epi->index != ele_info->index) {
        if (cur_ele->update || epi->stat == 'z') eff_img_ct++;
      }
      cur_img_data = cur_img_data->next;
    }
    // RMH: start
    printf("edges_and_cps: img_no = %d, eff_img_ct = %d\n", cur_ele->img_no, eff_img_ct);
    // RMH: end

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
              // Add pointer to image to the img_ptr list
              // TODO: Could be img_ptr[eff_img_ct] = cur_img_data->to_image; for readability
              *(img_ptr+eff_img_ct) = cur_img_data->to_image;
            eff_img_ct ++;
          }
        }
        cur_img_data = cur_img_data->next;
      }
      // RMH: Sort by: difference in **partner** element identifiers, ascending
      //               msp_direction
      //               left bound
      //               right bound
      //      This supports the processing of images in a element-pair fashion.
      qsort(img_ptr, eff_img_ct, sizeof(IMAGE_t *), partner_cmp);

      // recognizing full-length images, and put partial and secondary
      // images into a consistency tree, in which images connected are
      // consistent with each other (look at the consis() function of
      // definition of consistency)

      ele_partner = NULL;
      prim_p = NULL;
      epi = NULL;
      max_score = 0;
      // Start of the partner element image block
      ele_block_start = 0;

      for (i=0; i<eff_img_ct; i++) {
        cur_img = *(img_ptr+i);
        img_partner = partner(cur_img);

        epi = img_partner->ele_info;
        if (!epi->ele) ele_read_in(epi, 1);
        ele_partner = epi->ele;

        // RMH: DEBUG TODO remove
        printf("edges_and_cps: processing image cur_ele_info->index=%d, par_ele_info->index=%d\n", cur_img->ele_info->index, img_partner->ele_info->index );

        // full length: current element
        //   RMH: Currently CUTOFF2 is 0.9
        //        So image must be within 10bp of both element endpoints, and for elements
        //        < 200bp this becomes even more stringent.
        int is_full_length = 0;
        if (full_length(cur_img, CUTOFF2)) {
          printf("edges_and_cps:     image is primary because full_length with current element!\n");
          is_full_length = 1;
          // Increment the full length image counter
          cur_ele->flimg_no++;
        }
        // full length: partner element
        if (full_length(img_partner, CUTOFF2)) {
          printf("edges_and_cps:     image is primary because full_length with partner element!\n");
          is_full_length = 1;
          // Increment the full length image counter
          ele_partner->flimg_no++;
        }
        // If the image is full length for either/both of the
        // elements and is the highest score seen so far,
        // save it as prim_p.
        //
        // Original code (through 1.08):
        //      typedef struct msp {
        //        ...
        //        float iden;
        //        ...
        //      } MSP_t;
        //      ...
        //      int32_t max_score=0;
        //      ...
        //      if (prim == 1) {
        //        prim = 0;
        //        if (cur_img->to_msp->iden > max_score) {
        //          max_score = cur_img->to_msp->iden;
        //          dir = cur_img->to_msp->direction;
        //          prim_p = cur_img->to_msp;
        //        }
        //      }
        //
        // This appears to use identity rather than score
        // as a proxy for the quality of the alignment.  I am not
        // sure this is what he intended.  I have modified this
        // to use a score based on the identity *and* the length
        // of the alignment.  This is similar to the scoring used
        // in find_prim().
        //
        // Also note that in the original code the max_score was
        // truncating the float (identity) to an integer.
        if ( is_full_length ) {
          // RMH: Calculate a score based on the identity of the MSP in a similar
          //      fashion to find_prim()
          int32_t len_sum = (cur_img->frag.rb - cur_img->frag.lb) +
                            (img_partner->frag.rb - img_partner->frag.lb);
          int32_t ident_score = (int32_t) cur_img->to_msp->iden * (len_sum/2);
          // Original code
          //if (cur_img->to_msp->iden > max_score) {
          if (ident_score > max_score) {
            printf("edges_and_cps:     **** image is the new high scoring (identity) primary! [%f,%d,%d]\n", cur_img->to_msp->iden, ident_score, max_score);
            //max_score = cur_img->to_msp->iden;
            max_score = ident_score;
            dir = cur_img->to_msp->direction;
            prim_p =  cur_img->to_msp;
          }
        }

        // DEBUG
        if (!prim_p) {
          // RMH: No full length images found yet
          //      NOTE: Prequal flag is set
          printf("edges_and_cps:      Adding to consis tree (not full-length and no full-length found yet): e%d im = %s:%d-%d   e%d pt = %s:%d-%d\n", cur_img->ele_info->index, cur_img->frag.seq_name, cur_img->frag.lb, cur_img->frag.rb, img_partner->ele_info->index, img_partner->frag.seq_name, img_partner->frag.lb, img_partner->frag.rb);
          // NOTE: This adds all non-primary images to the tree (because 'prequal=1')
          consis_tree_build(consis_rt, cur_img, 1);
          //print_consis_tree(consis_rt);
          //print_consis_tree_GV(consis_rt,0,0);
          //print_ascii_tree(consis_rt);
        }

        // Have we reached the end of image set ( e.g a set of images that share the same partner element )?
        if ( i == eff_img_ct-1 || img_partner->ele_info->index != partner(*(img_ptr+i+1))->ele_info->index ) {
          int found_primary = 0;
          printf("edges_and_cps: This is the last image to process in a block\n");

          // Did we locate a high scoring full-length image?
          if (prim_p) {
            // We did, mark the single full-length image (MSP) to the elements as primary ('p')
            prim_p->stat = 'p';
            found_primary = 1;
          } else {
            // No full-length, find partial primary images that can be seen as part of a full-length alignment
            // NOTE: This flags all MSPs that can be seen as part of one or more full(ish)-length alignments
            //       between the elements.  This differs from the full_length initial check, which only flags
            //       one (highest score) MSP from the set of possible full-length MSPS.
            //
            //       Also, as the scoring function for the find_prim function is more lax than the full_length
            //       function, it is *still* possible to have a single MSP be marked primary here when it
            //       failed in full_length().
            //
            //    This may be a bit moot if these characteristics are not important other than to define
            //    a single edge as "primary" between the elements.
            //
            printf("Calling find_prim with 0 initial score!\n");
            // TODO: SUMMARIZE THIS ROUTINE HERE
            // Start off with first real node below root
            found_primary = find_prim(consis_rt->children, CUTOFF2, ele_info->ele->frag.lb, -1, 0, 0, 0, 0, 0, &token_mark, &max_score, &dir);
            printf("Tree after find_prim:\n");
            print_consis_tree_GV(consis_rt,0,0);
            // RMH: Debug - TODO remove
            if ( found_primary )
              printf("DEBUG: Identified partial primary image!\n");
          }

          // Build graph edge
          //   The graph is represented by edges between pairs of elements.  Pointers to an edge
          //   are stored in both element's edge binary trees by edge_index value.  E.g.
          //
          //              E1 --------------------------> E2
          //                        max_score
          //                        edge_index
          //
          //    Edge:
          //       index: 10,
          //       score: 100,
          //       direction: 1,
          //       status: 'p'
          //
          //     The use of a binary tree is at first glance a bit odd, as the inserted edges
          //     will initially be in sorted order ( edge_index is incremented as we build the
          //     graph ) -- this is equivalent to a linked list but in tree form.  Zhirong could
          //     have used an AVL or Red/Black tree here but chose instead to rebalance the tree
          //     at the end using an array and a divide and conquer approach.  Not sure why he
          //     didn't store the edges in an array to begin with.
          //
          if ( found_primary ) {
            printf("edges_and_cps:     Adding primary edge! with max_score = %d\n", max_score);

            // Add primary edge to the pair of elements
            if (ele_info->index != epi->index){
              // Store edge information in binary trees in both ele_info->edges and epi->edges.
              // Stored by edge_index which is a GLOBAL variable!
              // TODO: stomp out globals
              add_edge(ele_info, epi, 'p', max_score, dir);
            }else {
              err_no ++;
              fprintf(log_file, "error:  self edge seen: ele %d\n", ele_info->index);
              fflush(log_file);
            }

            // Record the primary image start/end points to both current and partner elements
            // E.g:
            //      images/msps:
            //          E1 100-200 ---  E2 300-400
            //          E1 50-300  ---  E2 200-440
            //          E1 1-100   ---  E3 157-250
            //
            //  Adds to the front of the linked list ( ie. unshift )
            //  E1->PCP:  [1,e3]-->[100,e3]-->[50,e2]-->[300,e2]-->[100,e2]-->[200,e2]
            //  E2->PCP:  [200,e1]-->[440,e1]-->[300,e1]-->[400,e1]
            //  E3->PCP:  [157,e1]-->[250,e1]
            //
            for (j=ele_block_start; j<=i; j++) {
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
          } else { // no primary images, full-length or partial
            // Add secondary edge to the pair of elements
            printf("edges_and_cps:     Adding secondary edge!\n");
            if (ele_info->index != epi->index) add_edge(ele_info, epi, 's', 0, 0);
            else {
              err_no ++;
              fprintf(log_file, "error:  self edge seen: ele %d\n", ele_info->index);
              fflush(log_file);
            }
          }

          // Free the consistency tree for this pair of elements
          if (consis_rt->children != NULL) {
            int released = consis_tree_free(consis_rt->children);
            printf("Freeing tree!: %d released\n", released);
            consis_rt->children = NULL;
          }

          // Reset the primary image and the max score
          ele_partner = NULL;
          prim_p = NULL;
          epi = NULL;
          max_score = 0;
          // Set the start of the next element image block
          ele_block_start = i+1;

        }// if reached the start of a new ele_partner

      } // for (i=0; i<eff_img_ct; i++)


  // DEBUG
  printf("PCP list:\n");
  print_cp_list(cur_ele->PCP);
  printf("printing edge tree\n");
  print_edge_tree(ele_info->ele->edges, 0);
  // END DEBUG

      if (eff_img_ct > MAX_IMG) {
        free(img_ptr);
      }
      free(consis_rt);
      free(token_image);
    } // if eff_img_ct

    // This appears to rebuild the edge tree by reordering the edges
    // (assuming they are already sorted by index) into a balanced binary tree.
    // Strangely it does this by taking an (ordered) tree structure like so:
    //
    //      2
    //       \
    //        10
    //          \
    //           11
    //            \
    //             42
    //
    // Converting it to an array:  2, 10, 11, 42
    //
    // and back to a tree structure but in the format of a balanced binary tree:
    //
    //         11
    //        /  \
    //      10    42
    //      /
    //     2
    printf("EDGE TREE Pre-adjustment\n");
    print_edge_tree_TEXT(ele_info->ele->edges, 0, "");
    if (ele_info->ele->edges) adjust_edge_tree(ele_info);
    printf("EDGE TREE Post-adjustment\n");
    print_edge_tree_TEXT(ele_info->ele->edges, 0, "");

    cur_ele->update = 0;
    ele_info->stat = 't';
}



// Current state:
//    Full Length = image reaches within 10bp of both ends of the
//    element *and* the ratio of the image size to the element size
//    is > than CUTOFF2 (0.9).  The CUTOFF is only applicable to
//    sequences < 200bp.
//
//    The effect of CUTOFF2 is therefore, is to prohibit anything but
//    full length images for elements < 200bp.
//
int full_length(IMAGE_t *i, float cutoff) {
  if (!i->ele_info->ele) {
    err_no ++;
    fprintf(log_file, "error:  element %d not in memory\n", i->ele_info->index);
    fflush(log_file);
    // TODO: WHY IS THIS HERE?
    // RMH: Looks like this could have recovered here, but exited instead.
    exit(3);
    ele_read_in(i->ele_info, 1);
  }
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



// USES GLOBALS!!!
void add_edge(ELE_INFO_t *ele1_info, ELE_INFO_t *ele2_info, char type, int32_t score, short dir) {
  EDGE_t *new = EDGE_malloc();
  int ct;

  // TODO: Yikes...global!
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



// This appears to rebuild the edge tree by reordering the edges
// (assuming they are already sorted by index) into a balanced binary tree.
// Strangely it does this by taking an (ordered) tree structure like so:
//
//      2
//       \
//        10
//          \
//           11
//            \
//             42
//
// Converting it to an array:  2, 10, 11, 42
//
// and back to a tree structure but in the format of a balanced binary tree:
//
//         11
//        /  \
//      10    42
//      /
//     2
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




// Depth first search of edge tree adding a pointer to each edge to the edge_array
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

// RMH: Insert images into a tree
//      Unique images ( e.g non overlapping ) get placed on
//      a distinct level in the tree whereas overlapping images
//      are placed on the same level.
//
//  NOTE:
//      The tree may contain overlapping images on different
//      strands.
//
int consis_tree_build(IMG_NODE_t *rt, IMAGE_t *im, int prequal) {
  int sum=0;
  IMG_NODE_t *nex_rt, *node;
  // Upon the first call to this function prequal is set to 1. This
  // supports having a top level tree node that is not compared to the
  // incoming image.
  //
  //  NOTE: It looks like it was setup by Zhirong go to support including
  //        overlapping elements as siblings "sib" on the same level.
  if (prequal || consis(rt->to_image, im, CUTOFF2)) {
    if ( prequal ) {
      printf("consis_tree_build(top_level, prequal=1)\n");
      printImage(im);
    }else {
      printf("Adding image because consis=1\n");
      printImage(im);
    }
    //if ( !prequal )
    //  printf("consis_tree_build: Found consis: im %s:%d-%d and rt %s:%d-%d\n", im->frag.seq_name, im->frag.lb, im->frag.rb, rt->to_image->frag.seq_name, rt->to_image->frag.lb, rt->to_image->frag.rb);
    // Take child and and if it's not null, consider it's siblings
    nex_rt = rt->children;
    while (nex_rt != NULL) {
      printf("Here we are\n");
      printImage(nex_rt->to_image);
      sum += consis_tree_build(nex_rt, im, 0);
      nex_rt = nex_rt->sib;
    }
    printf("SUM = %d\n", sum);
    if (!sum) {
      node = (IMG_NODE_t *) malloc(sizeof(IMG_NODE_t));
      node->recorded = 0;
      node->to_image = im;
      node->sib = NULL;
      node->children = NULL;
      // Take the child node and if its not null add the new node to the siblings list
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
  printf("consis: %s:%d-%d dir=%d and %s:%d-%d dir=%d \n", i1->frag.seq_name, i1->frag.lb, i1->frag.rb, i1->to_msp->direction, i2->frag.seq_name, i2->frag.lb, i2->frag.rb, i2->to_msp->direction);
  if (i1->ele_info->index == i2->ele_info->index &&
      ip1->ele_info->index == ip2->ele_info->index &&
      i1->to_msp->direction == i2->to_msp->direction) {
    if (i1->to_msp->direction == 1) {
      // This is an optimisation to avoid calling sing_cov when its not necessary.
      // The result of this test is negative if the two msps do not agree on orientation
      // and are therefore not consistent with each other.
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
      }else { printf("consis: Inconsistent orientation!\n"); }
    } else {
      // Same as above but for reverse orientation
      if ((i1->frag.lb - i2->frag.lb)*(ip1->frag.rb - ip2->frag.rb) < 0) {
        if (!sing_cov(&i1->frag, &i2->frag, 1.0-cutoff) && !sing_cov(&ip1->frag, &ip2->frag, 1.0-cutoff)) {
          // Do not overlap significantly for either element
          res = 1;
        }
      }else { printf("consis: Inconsistent orientation!\n"); }
    }
  }
  printf("consis: Returning: %d\n", res);
  return res;
}



IMG_NODE_t ** node_entry(IMG_NODE_t **node_pp) {
  if (*node_pp != NULL) {
    return node_entry(&(*node_pp)->sib);
  }
  return node_pp;
}




int consis_tree_free(IMG_NODE_t *rt) {
  int count = 0;
  if (rt->sib != NULL) {
    count += consis_tree_free(rt->sib);
  } 
  if (rt->children != NULL) {
    count += consis_tree_free(rt->children);
  }
  free(rt);
  count ++;
  return count;
}




/*****************************************
 * functions for parsing the consis_tree *
 *****************************************/

/*
  I. function for finding primary images from the consis_tree
*/

// Given an alignment consistency tree between two elements, this function identifies
// a path(s) that satisifies the following conditions:
//
// The images in the path(s) are marked as "p" (primary) and ...
//
// A consistency tree is defined as a tree where
//
// This is called when there are no full-length images found for a pair of elements.
// That is...there are no images that reach within 10bp (or lower for elements < 200bp
// see CUTOFF2) of both ends of the element.
//
// Initially called with
//       end1 = ele->frag.lb, and end2 = -1
//       efl1 = 0, efl2 = 0
//       al1 = 0, al2 = 0,
//       score = 0
// Making a guess here...but int32_t is probably not adequate in some instances for al1/al2 and probably
// has caused problems where overflow occurs.
//
//    end2 = -1 is a flag to indicate the first alignment in the group (initial call to find_prim)
//
// Hardcoded values:
//    Edge gap cutoff = 10bp
//    Len differences = efl1-al1 < 30 || efl2-al2 < 30
//
int find_prim(IMG_NODE_t *nd, float cutoff, int32_t end1, int32_t end2,
              int32_t efl1, int32_t efl2, int32_t al1, int32_t al2,
              int32_t score, int *pmarkp, int32_t *sc, short *d) {
  int sum = 0, mark=0;
  int32_t skip1, skip2, len1, len2;
  IMAGE_t *ipt;

  // Proccess siblings first ( inconsistent sequences that form alternate paths )
  if (nd->sib) sum += find_prim(nd->sib, cutoff, end1, end2, efl1, efl2, al1, al2, score, pmarkp, sc, d);

  ipt = partner(nd->to_image);
  // Reminder...these values are absolute positions provided in the msps file.  So an element/image
  // boundaries do not necessarily start at 1.
  if (end2 < 0) { //first alignment in the group
    // direction = 1 is positive strand
    // direction = -1 is negative strand
    if (nd->to_image->to_msp->direction == 1)
      end2 = ipt->ele_info->ele->frag.lb;
    else
      end2 = ipt->ele_info->ele->frag.rb;
  }

  // "end1/end2" is a bit of a misnomer.  Here this represents
  //  end1 = the left bound of the element ( ele_info->ele->frag.lb )
  //  end2 = the left bound of the partner element
  //printf("find_prim:  end1=%d, end2=%d, nd->direction = %d, efl1=%d, efl2=%d, al1=%d, al2=%d, score=%d, mark=%d, sc=%d, d=%d\n", end1, end2, nd->to_image->to_msp->direction, efl1, efl2, al1, al2, score, mark, *sc, *d);
  //printf("nd image:\n");
  //printImage(nd->to_image);

  // Determine the starting gap for each image (relative to each element):
  // skip1 is the image's left bound - element's left bound
  //    - So == 0 if they are the same, <0 if the image extends beyond the element
  //      and >0 if it starts later than the element's left bound.
  // skip2 is the partner element's left bound - the partner image's left bound
  //    - So == 0 if they are the same, <0 if the image extends beyond the partner element
  //      and >0 if it starts later than the element's left bound.
  skip1 = nd->to_image->frag.lb - end1;
  if ( nd->to_image->to_msp->direction == 1 )
    skip2 = ipt->frag.lb - end2;
  else
    skip2 = end2 - ipt->frag.rb;
  //printf("    frag.lb = %d, frag.rb = %d, ipt->frag.lb = %d, ipt->frag.rb = %d, skip1 = %d, skip2 = %d\n", nd->to_image->frag.lb, nd->to_image->frag.rb, ipt->frag.lb, ipt->frag.rb, skip1, skip2);

  // If we *start* with a gap **AND** our partner *starts* with a gap
  // ( > 10bp ) add the gaps to the efls.
  //
  //   |-----ele1--------|        |------ele2---------|
  //     >10 |--img1-...            >10 |--img2-...
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
    if (nd->to_image->to_msp->direction == 1)
      end2 = ipt->frag.rb;
    else
      end2 = ipt->frag.lb;
    //printf("find_prim: Calling children with end1=%d end2=%d\n", end1, end2);
    sum += find_prim(nd->children, cutoff, end1, end2, efl1, efl2, al1, al2, score, &mark, sc, d);
  } else {  // last alignment in path
    //printf("find_prim: Last alignment in group\n");
    skip1 = nd->to_image->ele_info->ele->frag.rb - nd->to_image->frag.rb;
    if (nd->to_image->to_msp->direction == 1)
      skip2 = ipt->ele_info->ele->frag.rb - ipt->frag.rb;
    else
      skip2 = ipt->frag.lb - ipt->ele_info->ele->frag.lb;

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

    // RMH: DEBUG TODO remove
    //printf("al1=%d al2=%d efl1=%d efl2=%d:  al/ef %f, %f   ef-al %d, %d  ele:%d-%d, ptn:%d-%d\n", al1, al2, efl1, efl2,(1.0*al1/efl1), (1.0*al2/efl2), (efl1-al1), (efl2-al2), nd->to_image->ele_info->ele->frag.lb, nd->to_image->ele_info->ele->frag.rb, ipt->ele_info->ele->frag.lb, ipt->ele_info->ele->frag.rb );
    // RMH: end

    //  Currently cutoff = 0.9
    //     So this appears to be checking that the alignment covers at least 90% of one or the other element **AND**
    //     it's less than 30bp short of the full length of one or the other element.
    //     The 30bp check was added in RECON1.03 (no reason given)
    //        Pre 1.03 code:
    //           if (1.0*al1/efl1 > cutoff || 1.0*al2/efl2 > cutoff) {
    //
    //     So what do al# and efl# represent?
    //     al# = actual length of the alignment or
    //
    //         |----------------------- element (500bp) ---------------------------------|
    //             |----alignment (42bp)----|      |----alignment (100bp)----|
    //
    //     al1 = 42 + 100 = 142
    //
    //     Instead of considering the aligned fraction of the element length it seems to be
    //     considering the aligned fraction of a (potentially) shorter element length
    //     (perhaps "effective length" or efl#).
    //
    //         |----------------------- element (500bp) ---------------------------------|
    //         ..9bp..|--alignment (42bp)--|..100bp..|--alignment (100bp)--|....249bp.....
    //
    //      The effective length is the length of the alignments + the length of any gaps
    //      exceeding 10bp.
    //
    //      efl1 = 42 + 100 + 100 + 249 = 491
    //
    //  So the ratio is :
    //
    //       al1/efl1 = 142/491 = 0.289, so not > cutoff (0.9)
    //
    if ( (1.0*al1/efl1 > cutoff || 1.0*al2/efl2 > cutoff) && (efl1-al1 < 30 || efl2-al2 < 30) ) {
      sum = 1; // Shouldn't this be += 1?
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
      printf("find_prim: Found a primary alignment path: score = %d, dir = %d\n", score, *d);
    }
  }
  // Marks all images in the path as primary ( upon recursion back up the tree )
  if (mark) {
      //printf("find_prim: Marking alignment as primary:\n");
      //printImage(nd->to_image);
      nd->to_image->to_msp->stat  = 'p';
      *pmarkp = 1;
  }
  //printf("find_prim: Returning sum = %d\n", sum);
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

// RMH: Sort by: difference in partner element identifiers, ascending
//               msp_direction
//               primary left bound
//               primary right bound
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


void test_consis_tree(){
  IMAGE_t *token_image, *new_image;
  IMG_NODE_t *consis_rt;
  MSP_t *msp_tmp;
  ELE_INFO_t *ele1, *ele2;
  ELEMENT_t *ele1_ele, *ele2_ele;

  // Top-level token image
  token_image = (IMAGE_t *) malloc(sizeof(IMAGE_t));
  token_image->frag.lb = 0;
  token_image->frag.rb = 0;
  token_image->to_msp = NULL;
  token_image->ele_info = NULL;

  // consistency tree
  consis_rt = (IMG_NODE_t *) malloc(sizeof(IMG_NODE_t));
  consis_rt->to_image = token_image;
  consis_rt->children = NULL;
  consis_rt->sib = NULL;

  // Elements
  ele1_ele = (ELEMENT_t *) malloc(sizeof(ELEMENT_t));
  ele1_ele->index = 1;
  ele1_ele->frag.lb = 1;
  ele1_ele->frag.rb = 305;
  ele1_ele->TBD = NULL;
  ele1_ele->redef = NULL;
  ele1_ele->direction = 1;

  ele1 = (ELE_INFO_t *) malloc(sizeof(ELE_INFO_t));
  ele1->index = 1;
  ele1->ele = ele1_ele;
  ele1->stat = 'z';
  ele1->next = NULL;

  ele2_ele = (ELEMENT_t *) malloc(sizeof(ELEMENT_t));
  ele1_ele->index = 2;
  ele2_ele->frag.lb = 1;
  ele2_ele->frag.rb = 1000;
  ele2_ele->TBD = NULL;
  ele2_ele->redef = NULL;
  ele2_ele->direction = 1;

  ele2 = (ELE_INFO_t *) malloc(sizeof(ELE_INFO_t));
  ele2->index = 2;
  ele2->ele = ele2_ele;
  ele2->stat = 'z';
  ele2->next = NULL;

  typedef struct {
    int direction;
    int iden;
    int sbjct_lb;
    int sbjct_rb;
    int query_lb;
    int query_rb;
  } MSPParams;

  // Sorted by direction then subject start
  MSPParams msp_params[] = {
    // This should generate:
    //   c_tree: level=0 [(null):0:0-0]
    //   c_tree: level=1 e1:(null):1:8-150 to e2:(null):1:40-200 : siblings[]
    //   c_tree: level=2 e1:(null):1:155-200 to e2:(null):1:400-500 : siblings[]
    //   c_tree: level=3 e1:(null):1:205-303 to e2:(null):1:600-700 : siblings[]
    //
    // And find a simple primary path with all three ( score = 62 )
    //
    //                      el1=1-305         ele2=1-1000
    // direction, iden, sbjct_lb, sbjct_rb, query_lb, query_rb
    // ---------  ----   -------  --------  --------  --------
    //{1,           60,      8,        150,      40,      200},
    //{1,          100,     155,       200,     400,      500},
    //{1,           40,     205,       303,     600,      700},
    //
    // This is more complex.  There is an inconsistency between the
    // ordering in e1/e2.  This generates:
    //
    //   c_tree: level=0 [(null):0:0-0]
    //   c_tree: level=1 e1:(null):1:8-150 to e2:(null):1:400-500 : siblings[(null):1:155-200(*e1:205),]
    //   c_tree: level=2 e1:(null):1:205-303 to e2:(null):1:600-700 : siblings[]
    //
    //                      el1=1-305         ele2=1-1000
    // direction, iden, sbjct_lb, sbjct_rb, query_lb, query_rb
    // ---------  ----   -------  --------  --------  --------
    {1,           30,      1,        80,     400,      500},
    {1,           60,      1,        150,     400,      500},
    {1,           30,      81,        180,     501,      599},
    {-1,          13,      85,        200,     1,        2  },
    {1,          100,    151,       200,      501,      599},
    {1,           30,      181,        305,     600,      650},
    {1,           40,     201,       305,     600,      700},

  };

  // It appears as if combining +/- strand msps in this structure
  // is not what we want to do.  They probably should be
  // handled seperately.

  for (int i = 0; i < sizeof(msp_params) / sizeof(msp_params[0]); i++) {
    MSP_t *msp_tmp = (MSP_t *) malloc(sizeof(MSP_t));
    msp_tmp->direction = msp_params[i].direction;
    msp_tmp->stat = '\0';
    msp_tmp->iden = msp_params[i].iden;
    msp_tmp->score = 0;
    msp_tmp->sbjct.frag.lb = msp_params[i].sbjct_lb;
    msp_tmp->sbjct.frag.rb = msp_params[i].sbjct_rb;
    msp_tmp->sbjct.to_msp = msp_tmp;
    msp_tmp->sbjct.ele_info = ele1;
    msp_tmp->query.frag.lb = msp_params[i].query_lb;
    msp_tmp->query.frag.rb = msp_params[i].query_rb;
    msp_tmp->query.to_msp = msp_tmp;
    msp_tmp->query.ele_info = ele2;

    // Using subject images here...
    consis_tree_build(consis_rt, &(msp_tmp->sbjct), 1);
  }

  print_consis_tree(consis_rt);

  printf("consis_rt->children = %d-%d\n", consis_rt->children->to_image->frag.lb, consis_rt->children->to_image->frag.rb);
  IMG_NODE_t *foo = consis_rt->children->sib;
  while ( foo != NULL ) {
    printf("sib = %d-%d\n", foo->to_image->frag.lb, foo->to_image->frag.rb);
    foo = foo->sib;
  }
  print_consis_tree_GV(consis_rt,0,0);
  //print_ascii_tree(consis_rt);
  int token_mark = 0;
  int max_score = 0;
  int dir = 0;
  int found_primary = find_prim(consis_rt->children, CUTOFF2, ele1_ele->frag.lb, -1, 0, 0, 0, 0, 0, &token_mark, &max_score, &dir);
  printf("AFTER\n");
  print_consis_tree_GV(consis_rt,0,0);
}

