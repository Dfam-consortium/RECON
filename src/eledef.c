/*
 * eledef.c  --  Stage 2: initial element definition by single-linkage clustering
 *
 * Algorithm overview
 * ------------------
 * This stage reads the sorted IMAGE files produced by imagespread and groups
 * overlapping images on the same sequence into "elements" using single-linkage
 * clustering (or optionally double-linkage via the "double" method).
 *
 * An element is a contiguous genomic interval that encompasses a set of
 * overlapping MSP images.  The representative interval (ele_frag) starts as
 * the first image and extends right as additional images are incorporated.
 * When an image falls outside the current element's right boundary by more
 * than MIN_OVERLAP_BP base pairs, the element is closed and a new one starts.
 *
 * The key insight is that images are processed in sorted order (by sequence
 * name, then by left boundary).  This means all images that could belong to
 * a given element appear as a contiguous run in the input -- except that some
 * may have been temporarily buffered in the "remain" list because they were
 * too short to qualify when first seen but might qualify after the element
 * boundary is extended.
 *
 * Fractional overlap criteria (sing_cov / doub_cov):
 *   sing_cov: overlap >= cutoff * max(len1, len2)  [default cutoff = 0.5]
 *   doub_cov: overlap >= cutoff * min(len1, len2)  [default cutoff = 0.9]
 *
 * Two-pass output:
 *   Pass 1 (ele_def):  clusters images -> elements, writes ele_def_res/img_prot
 *   Pass 2 (DUMBBELL): loads MSP detail from msp_file, writes ele_def_res/e<N>
 *
 * Usage
 *   eledef seq_list msp_file single|double [cutoff] [-l log_level]
 *
 *   -l <level>  log verbosity: 0=silent 1=error 2=warn 3=info 4=debug
 *               (default: 3=info)
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 *
 * Notes on macros
 * ---------------
 * Several inner-loop operations are implemented as macros for performance.
 * Each macro is documented at its definition site with its preconditions,
 * postconditions, and side effects.
 */

#include "seqlist.h"
#include "msps.h"
#include "recon_log.h"

/* freopen() inside the DUMBBELL macro is called for its side-effect;
 * the return value is intentionally ignored (failure is caught by the
 * subsequent fprintf calls which will silently no-op on a NULL stream). */
#pragma GCC diagnostic ignored "-Wunused-result"

/* ELEDEF_IMGPROT_CAP defined in recon_defs.h (pulled in via bolts.h).
 * IMG_CAP is the backward-compat alias. */

/* ---- Per-program log level storage (required by recon_log.h) ---- */
int   recon_log_level = RECON_LOG_INFO;
FILE *recon_log_fp    = NULL;


/* ============================================================
 * eledef-local data structures
 *
 * These three structures are internal to this stage and not shared
 * with other pipeline stages.
 * ============================================================ */

/*
 * MPROT_t  --  MSP prototype: tracks which element each MSP half was assigned to
 *
 * For each MSP (indexed by its position in the msp_file), this records the
 * element indices that the query (pe) and subject (se) images were assigned to
 * during the ele_def() clustering pass.
 *
 * pe  (was: pe) -- element index assigned to the query-side (primary/even) image
 * se  (was: se) -- element index assigned to the subject-side (secondary/odd) image
 *
 * Note: pe/se stand for "primary element" and "secondary element" in the original
 * code, reflecting that even image indices come from the query side and odd
 * image indices from the subject side.
 */
typedef struct msp_prototype {
  int pe, se;
} MPROT_t;

/*
 * IPROT_t  --  image prototype: lightweight image record for the DUMBBELL pass
 *
 * index      image index (even = query-side, odd = subject-side).
 * ele_index  element index this image was assigned to by ele_def().
 * to_msp     pointer to a scratch MSP_t into which the full MSP is loaded
 *            by img_charge() so that the DUMBBELL macro can write it out.
 */
typedef struct img_prototype {
  int   index, ele_index;
  MSP_t *to_msp;
} IPROT_t;

/*
 * EPROT_t  --  element prototype: summary record for each element
 *
 * flag    set to 1 after the element's header (index/frag/img_no) has been
 *         written to its ele_def_res/e<N> file, preventing duplicate headers.
 * index   element index (1-based).
 * img_no  count of images assigned to this element.
 * frag    representative genomic interval for the element.
 * next    link for the all_ep list built during ele_def().
 */
typedef struct ele_prototype {
  short  flag;
  int    index, img_no;
  FRAG_t frag;
  struct ele_prototype *next;
} EPROT_t;


/* ---- Function prototypes ---- */
void ele_def(int, FILE *, float, EPROT_t **, int *, MPROT_t **);
void img_charge(IPROT_t **, int, FILE *);
int  index_cmp(const void *, const void *);

/* toString helpers for eledef-local types */
void print_mprot(MPROT_t *m);
void print_iprot(IPROT_t *ip);
void print_eprot(EPROT_t *ep);

/* Output file handles (global so macros can access them) */
FILE *msp_no, *all_ele, *img_prot, *ele_no, *err, *size_list;


/* ============================================================
 * Macros
 *
 * These macros are used inside ele_def() for performance-critical
 * inner loops.  They are documented here with their preconditions,
 * side effects, and the variables they read/modify.
 * ============================================================ */

/*
 * INCLUDE_IMAGE
 *
 * Test whether the current "remain" list image (pointed to by 'cur') overlaps
 * the current element (ele_frag) enough to be admitted, and if so, incorporate
 * it into the element or move past it.
 *
 * Preconditions:
 *   cur      points to the IMG_DATA_t node under consideration (remain list)
 *   prev     points to the node before cur (or NULL if cur == remain)
 *   remain   head of the remain list
 *   ele_frag the current element's representative interval
 *   *ecp     the current element index
 *   all_mprot, img_prot, method, cutoff are in scope
 *
 * Postconditions (cov == 1):
 *   - cur is unlinked from the remain list and its memory freed
 *   - ele_frag.rb may be extended if cur's image extends further right
 *   - If rb is extended, cur is reset to remain (re-scan all remaining images)
 *   - img_ct is incremented
 *
 * Postconditions (cov == 0):
 *   - prev advances to cur, cur advances to cur->next (no change to list)
 */
#define INCLUDE_IMAGE \
    if (method == 1) {cov=sing_cov(&ele_frag, &cur->to_image->frag, cutoff);}\
            else {cov=doub_cov(&ele_frag, &cur->to_image->frag, cutoff);}\
            if (cov) {\
              img_ct ++;\
	      if (cur->to_image->index%2) {\
		(*(all_mprot+cur->to_image->index/2))->se = (*ecp);\
	      } else {\
		(*(all_mprot+cur->to_image->index/2))->pe = (*ecp);\
	      }\
	      fprintf(img_prot, "%d %d\n", (*ecp), cur->to_image->index);\
              fflush(img_prot);\
	      if (!cur->next) tail = prev;\
	      if (prev) prev->next = cur->next;\
	      else remain = cur->next;\
	      if (ele_frag.rb < cur->to_image->frag.rb) {\
		ele_frag.rb = cur->to_image->frag.rb;\
		free(cur->to_image);\
		free(cur);\
		cur = remain;\
		prev = NULL;\
	      } else {\
		free(cur->to_image);\
		free(cur);\
		if (prev) cur = prev->next;\
		else cur = remain;\
	      }\
	    } else {\
	      prev = cur;\
	      cur = cur->next;\
	    }


/*
 * SAVE_ELEMENT_TO_EP_LIST
 *
 * Append a new EPROT_t for the just-completed element to the ep linked list,
 * and write a summary line to the naive_eles file.
 *
 * Variables read: *ecp, img_ct, ele_frag, all_epp, ep_tail, all_ele
 *
 * Postconditions:
 *   - A new EPROT_t is malloc'd and appended to *all_epp via ep_tail.
 *   - A line is written to all_ele (summary/naive_eles).
 */
#define SAVE_ELEMENT_TO_EP_LIST \
      ep_tmp = (EPROT_t *) malloc(sizeof(EPROT_t));\
      ep_tmp->flag = 0;\
      ep_tmp->index = (*ecp);\
      ep_tmp->img_no = img_ct;\
      ep_tmp->frag = ele_frag;\
      ep_tmp->next = NULL;\
      if (*all_epp) ep_tail->next = ep_tmp;\
      else *all_epp = ep_tmp;\
      ep_tail = ep_tmp;\
      fprintf(all_ele, "%d %s %d %d\n", (*ecp), ele_frag.seq_name, ele_frag.lb, ele_frag.rb)


/*
 * INIT_ELEMENT_FROM_REMAIN_LIST
 *
 * Start a new element using the first image in the remain list.
 * Unlinks that image from remain, sets ritetime=0, and increments *ecp.
 *
 * Preconditions:  remain is non-NULL; cur points to remain.
 * Postconditions:
 *   - ele_frag is initialized from cur->to_image->frag
 *   - cur is unlinked from remain; remain advances to its next node
 *   - cur is freed; cur is reset to (new) remain
 *   - prev = NULL (no predecessor in the now-reduced remain list)
 *   - ritetime = 0
 */
#define INIT_ELEMENT_FROM_REMAIN_LIST \
	ritetime = 0;\
	(*ecp) ++;\
        img_ct = 1;\
	ele_frag = cur->to_image->frag;\
	if (cur->to_image->index%2) {\
	  (*(all_mprot+cur->to_image->index/2))->se = (*ecp);\
	} else {\
	  (*(all_mprot+cur->to_image->index/2))->pe = (*ecp);\
	}\
	fprintf(img_prot, "%d %d\n", (*ecp), cur->to_image->index);\
        fflush(img_prot);\
	remain = cur->next;\
	free(cur->to_image);\
	free(cur);\
	prev = NULL;\
	cur = remain


/*
 * EVALUATE_REMAINING_IMAGE
 *
 * Used when processing the remain list after a new element has been started
 * from remain (via INIT_ELEMENT_FROM_REMAIN_LIST).  For each node cur:
 *   - If it overlaps the current element by > MIN_OVERLAP_BP, apply INCLUDE_IMAGE.
 *   - Otherwise, save the element (it is complete) and set ritetime=1 to
 *     restart the remain-list scan with a new seed element.
 *
 * Note: the > 10 threshold is the MIN_OVERLAP_BP constant (literal 10).
 */
#define EVALUATE_REMAINING_IMAGE \
	  if (ele_frag.seq_name == cur->to_image->frag.seq_name \
              && ele_frag.rb - cur->to_image->frag.lb > MIN_OVERLAP_BP) {\
	    INCLUDE_IMAGE;\
	  } else {\
	    SAVE_ELEMENT_TO_EP_LIST;\
	    ritetime = 1;\
	    break;\
	  }


/*
 * SAVE_IMG_TO_REMAINING_LIST
 *
 * Append the current image (img) to the tail of the remain linked list.
 * remain* = img1.next->img2.next->img3.next = tail*
 *
 * Variables read: img (IMAGE_t value, copied into heap)
 * Variables modified: remain, tail
 */
#define SAVE_IMG_TO_REMAINING_LIST \
	 imgp = (IMAGE_t *) malloc(sizeof(IMAGE_t));\
	*imgp = img;\
	tmp = (IMG_DATA_t *) malloc(sizeof(IMG_DATA_t));\
	tmp->to_image = imgp;\
	tmp->next = NULL;\
	if (!remain) {\
	  remain = tmp;\
	} else {\
	  tail->next = tmp;\
	}\
	tail = tmp


/*
 * DUMBBELL
 *
 * Flush a batch of IPROT_t records to their element output files.
 *
 * This macro operates on a shadow copy of the all_iprot array (iprot_shadow,
 * of size iprot_ct) that has been sorted by image index.  For each record:
 *   1. img_charge() loads the corresponding MSP from msp_file.
 *   2. The element file (ele_def_res/e<N>) is opened (or re-opened).
 *   3. The MSP is written as an "msp ..." line.
 *
 * The flag field of EPROT_t prevents the element header from being written
 * more than once per element file.
 *
 * Variables read: iprot_shadow, iprot_ct, ep_array, all_ele, all_iprot,
 *                 msp_file, size_list
 */
#define DUMBBELL \
      img_charge(iprot_shadow, iprot_ct, msp_file);\
      for (i=0; i<iprot_ct; i++) {\
	if ((*(all_iprot+i))->index%2) {\
	  partner_index = (*(all_mprot+(*(all_iprot+i))->index/2))->pe;\
	} else {\
	  partner_index = (*(all_mprot+(*(all_iprot+i))->index/2))->se;\
	}\
	if ((*(all_iprot+i))->to_msp->score) {\
	  if (!(*(ep_array+((*(all_iprot+i))->ele_index-1)))->flag) {\
	    sprintf(ele_name, "ele_def_res/e%d", (*(all_iprot+i))->ele_index);\
	    freopen(ele_name, "w", all_ele);\
	    fprintf(all_ele, "index %d\n", (*(all_iprot+i))->ele_index);\
	    fprintf(all_ele, "frag %s %d %d\n", (*(ep_array+((*(all_iprot+i))->ele_index-1)))->frag.seq_name, (*(ep_array+((*(all_iprot+i))->ele_index-1)))->frag.lb, (*(ep_array+((*(all_iprot+i))->ele_index-1)))->frag.rb);\
            fprintf(all_ele, "img_no %d\n", (*(ep_array+((*(all_iprot+i))->ele_index-1)))->img_no);\
            fprintf(size_list, "%d %d\n", (*(all_iprot+i))->ele_index, (*(ep_array+((*(all_iprot+i))->ele_index-1)))->img_no);\
	    (*(ep_array+((*(all_iprot+i))->ele_index-1)))->flag = 1;\
	  }\
	  if (partner_index != (*(all_iprot+i))->ele_index) {\
	    if ((*(all_iprot+i))->index%2) fprintf(all_ele, "msp %d s %d %.1f %d %d %s %d %d %d %s %d %d\n", (*(all_iprot+i))->index, (*(all_iprot+i))->to_msp->score, (*(all_iprot+i))->to_msp->iden, (*(all_iprot+i))->to_msp->direction, partner_index, (*(all_iprot+i))->to_msp->query.frag.seq_name, (*(all_iprot+i))->to_msp->query.frag.lb, (*(all_iprot+i))->to_msp->query.frag.rb, (*(all_iprot+i))->ele_index, (*(all_iprot+i))->to_msp->sbjct.frag.seq_name, (*(all_iprot+i))->to_msp->sbjct.frag.lb, (*(all_iprot+i))->to_msp->sbjct.frag.rb);\
	    else fprintf(all_ele, "msp %d s %d %.1f %d %d %s %d %d %d %s %d %d\n", (*(all_iprot+i))->index, (*(all_iprot+i))->to_msp->score, (*(all_iprot+i))->to_msp->iden, (*(all_iprot+i))->to_msp->direction, (*(all_iprot+i))->ele_index, (*(all_iprot+i))->to_msp->query.frag.seq_name, (*(all_iprot+i))->to_msp->query.frag.lb, (*(all_iprot+i))->to_msp->query.frag.rb, partner_index, (*(all_iprot+i))->to_msp->sbjct.frag.seq_name, (*(all_iprot+i))->to_msp->sbjct.frag.lb, (*(all_iprot+i))->to_msp->sbjct.frag.rb);\
	  } else {\
	    if ((*(all_iprot+i))->index%2) fprintf(all_ele, "msp %d s %d %.1f %d %d %s %d %d %d %s %d %d\n", (*(all_iprot+i))->index, (*(all_iprot+i))->to_msp->score, (*(all_iprot+i))->to_msp->iden, (*(all_iprot+i))->to_msp->direction, partner_index, (*(all_iprot+i))->to_msp->query.frag.seq_name, (*(all_iprot+i))->to_msp->query.frag.lb, (*(all_iprot+i))->to_msp->query.frag.rb, (*(all_iprot+i))->ele_index, (*(all_iprot+i))->to_msp->sbjct.frag.seq_name, (*(all_iprot+i))->to_msp->sbjct.frag.lb, (*(all_iprot+i))->to_msp->sbjct.frag.rb);\
          }\
	}\
      }


int main (int argc, char *argv[]) {
  int ele_ct=0, msp_ct, i, method;
  float cutoff;
  char line[150], *m1="single", *m2="double";
  MPROT_t **all_mprot;
  EPROT_t *all_ep=NULL, **ep_array, *ep_tmp;

  IPROT_t **all_iprot, **iprot_shadow;
  int iprot_ct=0, partner_index;
  char ele_name[50];   /* name of element used as output file name */

  FILE *frags, *seq_list, *msp_file;

  /* Strip the optional "-l <level>" flag before positional arg parsing */
  if (recon_parse_log_flag(&argc, argv)) {
    fprintf(stderr, "error: -l requires a numeric log level argument\n");
    exit(1);
  }

  /* Validate command line */
  if (argc == 1) {
    printf("usage: eledef seq_list msp_file single|double [cutoff] [-l level]\n"
           "for method, choose 'single' or 'double'\n"
           "cutoff is optional\n"
           "-l <level>  log verbosity: 0=silent 1=error 2=warn "
           "3=info(default) 4=debug\n");
    exit(1);
  }

  seq_list = fopen(argv[1], "r");
  if (!seq_list) {
    printf("Input file for sequence name list %s not found.  Exit.\n", argv[1]);
    exit(2);
  }
  GetSeqNames(seq_list);

  msp_file = fopen(argv[2], "r");
  if (!msp_file) {
    printf("Input file of MSPs %s not found.  Exit.\n", argv[2]);
    exit(2);
  }

  if (!strcmp(argv[3], m1)) { method = 1; }
  else if (!strcmp(argv[3], m2)) { method = 2; }
  else {
    printf("please choose single or double to indicate the method to evaluate overlap.\n");
    printf("%s is not a valid option.  Exit.\n", argv[3]);
    exit(2);
  }

  if (argc > 4) {
    cutoff = atof(argv[4]);
  } else {
    /* Defaults from recon_defs.h */
    if (method == 1) { cutoff = ELEDEF_CUTOFF_SINGLE; }
    else             { cutoff = ELEDEF_CUTOFF_DOUBLE; }
  }

  /* Open auxiliary files */
  if (!(frags = fopen("images/images_sorted", "r"))) {
    printf("Can not open the fragment list file, exiting.\n");
    exit(1);
  }
  if (!(msp_no = fopen("summary/ori_msp_no", "r"))) {
    printf("Can not open msp_no, exiting.\n");
    exit(1);
  }
  err      = fopen("ele_def_res/errors",    "w");
  all_ele  = fopen("summary/naive_eles",    "w");
  img_prot = fopen("ele_def_res/img_prot",  "w");
  ele_no   = fopen("summary/naive_ele_no",  "w");
  size_list= fopen("ele_def_res/size_list", "w");

  /* Route RECON_LOG output to the error file */
  recon_log_fp = err;

  /* Read total MSP count written by imagespread */
  while (fgets(line, 100, msp_no)) {
    msp_ct = atoi(line);
  }

  /*
   * Allocate MPROT_t array: one entry per MSP, indexed by MSP position (0-based).
   * pe = element index for the query-side image; se = element index for subject-side.
   * Both are initialised to 0 (unassigned).
   */
  all_mprot = (MPROT_t **) malloc(msp_ct * sizeof(MPROT_t *));
  for (i = 0; i < msp_ct; i++) {
    all_mprot[i] = (MPROT_t *) malloc(sizeof(MPROT_t));
    all_mprot[i]->pe = 0;
    all_mprot[i]->se = 0;
  }

  /* ---- Pass 1: single-linkage clustering of images into elements ---- */
  ele_def(method, frags, cutoff, &all_ep, &ele_ct, all_mprot);

  fclose(ele_no);
  fclose(img_prot);
  fclose(msp_no);
  fclose(frags);

  img_prot = fopen("ele_def_res/img_prot", "r");

  /* Build an array of EPROT_t pointers indexed by element index (0-based) */
  ep_array = (EPROT_t **) malloc(ele_ct * sizeof(EPROT_t *));
  ep_tmp = all_ep;
  i = 0;
  while (ep_tmp) {
    ep_array[i] = ep_tmp;
    i++;
    ep_tmp = ep_tmp->next;
  }

  /*
   * Allocate IPROT_t arrays for the DUMBBELL pass.
   * ELEDEF_IMGPROT_CAP (was: IMG_CAP) sets the batch size.
   * When iprot_ct reaches the cap, DUMBBELL flushes the batch and resets.
   */
  all_iprot    = (IPROT_t **) malloc(ELEDEF_IMGPROT_CAP * sizeof(IPROT_t *));
  iprot_shadow = (IPROT_t **) malloc(ELEDEF_IMGPROT_CAP * sizeof(IPROT_t *));
  for (i = 0; i < ELEDEF_IMGPROT_CAP; i++) {
    all_iprot[i]          = (IPROT_t *) malloc(sizeof(IPROT_t));
    all_iprot[i]->to_msp  = (MSP_t  *) malloc(sizeof(MSP_t));
    iprot_shadow[i]       = all_iprot[i];
  }

  /* ---- Pass 2: DUMBBELL -- load MSP detail and write element files ---- */
  while (fgets(line, 100, img_prot)) {
    sscanf(line, "%d %d\n",
           &all_iprot[iprot_ct]->ele_index,
           &all_iprot[iprot_ct]->index);
    iprot_ct++;
    /* Flush when the buffer is full */
    if (iprot_ct == ELEDEF_IMGPROT_CAP) {
      DUMBBELL;
      iprot_ct = 0;
    }
  }

  if (iprot_ct) {
    for (i = 0; i < iprot_ct; i++) { iprot_shadow[i] = all_iprot[i]; }
    DUMBBELL;
  }

  /* ---- Cleanup ---- */
  for (i = 0; i < ELEDEF_IMGPROT_CAP; i++) {
    free(all_iprot[i]->to_msp);
    free(all_iprot[i]);
  }
  free(all_iprot);
  free(iprot_shadow);

  for (i = 0; i < msp_ct; i++) { free(all_mprot[i]); }
  free(all_mprot);

  for (i = 0; i < ele_ct; i++) { free(ep_array[i]); }
  free(ep_array);

  exit(0);
}


/*
 * ele_def  --  single-linkage clustering of sorted images into elements
 *
 * Reads IMAGE records from 'frags' (the sorted images/images_sorted file)
 * and clusters them into elements using the chosen overlap method.
 *
 * Algorithm
 * ---------
 * Images are processed in sorted order (by sequence name, then lb).  A
 * "remain" list holds images that were seen but did not qualify for the
 * current element (because they did not overlap enough at the time).
 *
 * For each image read from 'frags':
 *   - If ritetime (start_new_element) is set, start a new element with
 *     this image.
 *   - If the image overlaps ele_frag by > MIN_OVERLAP_BP and passes the
 *     coverage test (sing_cov or doub_cov), incorporate it.  If it extends
 *     ele_frag.rb, go back and re-check the remain list.
 *   - Otherwise, save the image to the remain list and continue.
 *   - When an image falls outside ele_frag's sequence or does not overlap
 *     at all, close the current element and process the remain list to
 *     start new elements from any unassigned images there.
 *
 * Parameters
 *   method     1 = sing_cov, 2 = doub_cov
 *   frags      FILE* of sorted IMAGE records
 *   cutoff     fractional-overlap threshold
 *   all_epp    output: linked list of EPROT_t records
 *   ecp        output: pointer to element count (was: *ecp in all callers)
 *   all_mprot  array of MPROT_t to be filled with element assignments
 *
 * Note on the "ritetime" variable (was: ritetime)
 *   ritetime == 1 means "right time to start a new element".  It is set
 *   at the start, whenever an element is closed, and whenever the remain
 *   list is exhausted.  The name is a pun on "right time" by the original
 *   author.
 */
void ele_def(int method, FILE *frags, float cutoff,
             EPROT_t **all_epp, int *ecp, MPROT_t **all_mprot) {
  int i, img_ct, cov;
  char line[100];
  FRAG_t ele_frag;
  char fragname[SEQ_NAME_MAX_LEN];
  int pos;
  IMAGE_t img, *imgp;
  IMG_DATA_t *cur, *prev=NULL, *remain=NULL, *tail, *tmp;
  short ritetime;      /* was: ritetime -- "right time" to start a new element */
  EPROT_t *ep_tail, *ep_tmp;

  ritetime = 1;   /* initially ready to start a new element */

  /*
   * Read IMAGE records from the sorted image file.  Each line contains:
   *   <img_index> <score> <seq_name> <lb> <rb>
   *
   * The score field is skipped (%*d) here; it was used only by imagespread
   * to write the file.
   *
   * Image index convention:
   *   even index (img.index % 2 == 0) => query side of the parent MSP
   *   odd  index (img.index % 2 == 1) => subject side of the parent MSP
   *   MSP position in file            => img.index / 2
   */
  while (fgets(line, 100, frags)) {
    sscanf(line, "%d %*d %s %d %d\n",
           &img.index, fragname, &img.frag.lb, &img.frag.rb);
    pos = GetSeqIndex(0, seq_count - 1, fragname);
    img.frag.seq_name = seq_name_table[pos];

    if (ritetime) {
      /* Start a new element seeded by this image */
      ritetime  = 0;
      (*ecp)++;
      img_ct    = 1;
      ele_frag  = img.frag;
      if (img.index % 2) {
        /* Odd index: subject side -- record in MPROT_t.se */
        all_mprot[img.index / 2]->se = (*ecp);
      } else {
        /* Even index: query side -- record in MPROT_t.pe */
        all_mprot[img.index / 2]->pe = (*ecp);
      }
      fprintf(img_prot, "%d %d\n", (*ecp), img.index);
      fflush(img_prot);
      continue;
    }

    /*
     * Check whether this image overlaps the current element.
     * Overlap test: same sequence AND ele_frag.rb - img.frag.lb > MIN_OVERLAP_BP.
     *
     * NOTE: seq_name comparison uses pointer equality (both point into
     * seq_name_table[]) rather than strcmp().
     */
    if (ele_frag.seq_name == img.frag.seq_name &&
        ele_frag.rb - img.frag.lb > MIN_OVERLAP_BP) {
      /* Candidate image: apply fractional coverage test */
      if (method == 1) { cov = sing_cov(&ele_frag, &img.frag, cutoff); }
      else             { cov = doub_cov(&ele_frag, &img.frag, cutoff); }

      if (cov) {
        /* Good image: admit to current element */
        img_ct++;
        if (img.index % 2) {
          all_mprot[img.index / 2]->se = (*ecp);
        } else {
          all_mprot[img.index / 2]->pe = (*ecp);
        }
        fprintf(img_prot, "%d %d\n", (*ecp), img.index);
        fflush(img_prot);

        if (ele_frag.rb < img.frag.rb) {
          /* Element boundary extended: re-scan remain list */
          ele_frag.rb = img.frag.rb;
          cur  = remain;
          prev = NULL;
          while (cur) {
            INCLUDE_IMAGE;
          }
        }
      } else {
        /* Image did not qualify: keep in remain list for later */
        SAVE_IMG_TO_REMAINING_LIST;
      }
    } else {
      /* Image falls outside this element: close element, process remain */
      SAVE_ELEMENT_TO_EP_LIST;
      SAVE_IMG_TO_REMAINING_LIST;
      ritetime = 1;

      /*
       * Process the remain list: try to start new elements from the
       * buffered images that were not incorporated into the element
       * that was just closed.
       */
      while (remain && ritetime) {
        cur = remain;
        INIT_ELEMENT_FROM_REMAIN_LIST;
        /* cur now points to the rest of the remain list */
        while (cur) {
          EVALUATE_REMAINING_IMAGE;
        }
      }
    }
  }   /* end of main image-reading loop */

  SAVE_ELEMENT_TO_EP_LIST;
  if (remain) ritetime = 1;

  /* Drain any images still in the remain list into their own elements */
  while (remain) {
    cur  = remain;
    prev = NULL;

    while (cur) {
      if (ritetime) {
        INIT_ELEMENT_FROM_REMAIN_LIST;
        continue;
      }
      EVALUATE_REMAINING_IMAGE;
    }
    if (!cur) {
      SAVE_ELEMENT_TO_EP_LIST;
      ritetime = 1;
    }
  }

  fprintf(ele_no, "%d\n", (*ecp));
}




/*
 * img_charge  --  load MSP records for a batch of IPROT_t entries
 *
 * Sorts the shadow[] array by image index (ascending), then scans the
 * msp_file sequentially, matching each line's position (pos) to the
 * indices in shadow[i]->index / 2.  When a match is found, the MSP is
 * parsed via scan_msp() and copied into shadow[i]->to_msp.
 *
 * Because the images are sorted and the MSP file is scanned linearly,
 * this is an O(N + M) pass rather than O(N * M).
 *
 * After loading all ct records, the file pointer is rewound so this
 * function can be called again for the next batch.
 *
 * Special handling:
 *   - Two images from the same MSP line (adjacent even/odd pair) are
 *     handled in the inner "same-line" block.
 *   - scan_msp() failure sets to_msp->score = 0 (used as a sentinel by
 *     the DUMBBELL macro to skip writing that MSP).
 */
void img_charge(IPROT_t **shadow, int ct, FILE *input) {
  int i=0, pos=0;
  char line[151];
  int scan_flag;
  MSP_t msp;

  qsort(shadow, ct, sizeof(IPROT_t *), index_cmp);

  while (fgets(line, 150, input)) {
    if (pos == shadow[i]->index / 2) {
      scan_flag = scan_msp(&msp, line);
      if (scan_flag) {
        shadow[i]->to_msp->score = 0;
        fprintf(err, "Wrong format in the MSP file line %d for image %d\n",
                pos, shadow[i]->index);
      } else {
        *(shadow[i]->to_msp) = msp;
      }
      i++;
      if (i == ct) break;

      /* Check if the next image is from the same MSP line */
      if (pos == shadow[i]->index / 2) {
        if (scan_flag) {
          shadow[i]->to_msp->score = 0;
          fprintf(err, "Wrong format in the MSP file line %d for image %d\n",
                  pos, shadow[i]->index);
        } else {
          *(shadow[i]->to_msp) = msp;
        }
        i++;
        if (i == ct) break;
      }
    }
    pos++;
  }
  rewind(input);
}


/*
 * index_cmp  --  qsort comparator: ascending order by IPROT_t.index
 */
int index_cmp(const void *i1, const void *i2) {
  return (*((IPROT_t **) i1))->index - (*((IPROT_t **) i2))->index;
}


/* ============================================================
 * toString-style print helpers for eledef-local types
 * ============================================================ */

void print_mprot(MPROT_t *m) {
  if (!m) { printf("MPROT_t: NULL\n"); return; }
  printf("MPROT_t: pe=%d (query-ele), se=%d (sbjct-ele)\n", m->pe, m->se);
}

void print_iprot(IPROT_t *ip) {
  if (!ip) { printf("IPROT_t: NULL\n"); return; }
  printf("IPROT_t: index=%d (%s), ele_index=%d",
         ip->index,
         (ip->index % 2 == 0) ? "query" : "sbjct",
         ip->ele_index);
  if (ip->to_msp) {
    printf(", msp->score=%d", ip->to_msp->score);
  } else {
    printf(", to_msp=NULL");
  }
  printf("\n");
}

void print_eprot(EPROT_t *ep) {
  if (!ep) { printf("EPROT_t: NULL\n"); return; }
  printf("EPROT_t: index=%d, img_no=%d, flag=%d, frag=[%s %d %d]\n",
         ep->index, ep->img_no, ep->flag,
         ep->frag.seq_name ? ep->frag.seq_name : "(null)",
         ep->frag.lb, ep->frag.rb);
}
