/*
 * eledef.c  --  Stage 1: image sorting and initial element definition
 *
 * Algorithm overview
 * ------------------
 * This stage reads the MSP file directly, creates IMAGE records in memory
 * (merged from the former imagespread stage), sorts them, and groups
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
 *   Pass 2 (dumbbell): loads MSP detail from msp_file, writes element database
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
 * NOTE: imagespread functionality (partitioning MSP images) has been merged
 * into this stage.  imagespread is retained in the source tree for reference
 * only.
 */

#include "seqlist.h"
#include "msps.h"
#include "recon_log.h"
#include "ele_db.h"

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
void ele_def(int, IMAGE_t *, int, float, EPROT_t **, int *, MPROT_t **);
void img_charge(IPROT_t **, int, FILE *);
int  index_cmp(const void *, const void *);

/* toString helpers for eledef-local types */
void print_mprot(MPROT_t *m);
void print_iprot(IPROT_t *ip);
void print_eprot(EPROT_t *ep);

/*
 * Static helper functions (previously macros).
 * These are file-local and will be inlined by the compiler under -O.
 */
static void save_element_to_ep_list(int ecp, int img_ct, FRAG_t ele_frag,
                                    EPROT_t **all_epp, EPROT_t **ep_tail_p);
static void save_img_to_remaining_list(IMAGE_t img,
                                       IMG_DATA_t **remain_p,
                                       IMG_DATA_t **tail_p);
static void include_image(IMG_DATA_t **cur_p, IMG_DATA_t **prev_p,
                          IMG_DATA_t **remain_p, IMG_DATA_t **tail_p,
                          FRAG_t *ele_frag_p, int ecp,
                          MPROT_t **all_mprot, int method, float cutoff,
                          int *img_ct_p);
static void init_element_from_remain_list(short *ritetime_p, int *ecp,
                                          int *img_ct_p, FRAG_t *ele_frag_p,
                                          MPROT_t **all_mprot,
                                          IMG_DATA_t **remain_p,
                                          IMG_DATA_t **cur_p,
                                          IMG_DATA_t **prev_p);
static int  evaluate_remaining_image(IMG_DATA_t **cur_p, IMG_DATA_t **prev_p,
                                     IMG_DATA_t **remain_p, IMG_DATA_t **tail_p,
                                     FRAG_t *ele_frag_p, int ecp,
                                     MPROT_t **all_mprot, int method, float cutoff,
                                     int *img_ct_p, EPROT_t **all_epp,
                                     EPROT_t **ep_tail_p, short *ritetime_p);
static void dumbbell_flush(IPROT_t **all_iprot, IPROT_t **iprot_shadow,
                           int iprot_ct, MPROT_t **all_mprot,
                           FILE *msp_file);

/* Output file handles (global so helper functions can access them) */
FILE *naive_eles, *img_prot, *ele_no, *err, *size_list;


/* ============================================================
 * Image sort comparator for the merged imagespread pass
 *
 * Sort order: seq_name (string asc) -> lb (asc) -> rb (desc) -> index (asc).
 * This matches the order produced by imagespread's sort(1) invocation.
 * ============================================================ */
static int image_sort_cmp(const void *a, const void *b) {
  const IMAGE_t *ia = (const IMAGE_t *)a;
  const IMAGE_t *ib = (const IMAGE_t *)b;
  int sc = strcmp(ia->frag.seq_name, ib->frag.seq_name);
  if (sc != 0) return sc;
  if (ia->frag.lb != ib->frag.lb) return ia->frag.lb - ib->frag.lb;
  if (ia->frag.rb != ib->frag.rb) return ib->frag.rb - ia->frag.rb;
  return ia->index - ib->index;
}


/* ============================================================
 * Static helper functions (formerly macros)
 *
 * These were converted from preprocessor macros to static functions.
 * Declaring them static allows the compiler to inline them under -O,
 * matching the original macro performance while gaining debuggability
 * and type safety.
 * ============================================================ */

/*
 * save_element_to_ep_list  --  append the current element to the ep list
 *
 * Allocates a new EPROT_t, fills it from the current element state, links
 * it onto the tail of *all_epp, and writes a summary line to all_ele
 * (summary/naive_eles).
 *
 * Parameters
 *   ecp        current element index (1-based)
 *   img_ct     number of images in this element
 *   ele_frag   representative genomic interval for this element
 *   all_epp    head pointer of the ep linked list (modified)
 *   ep_tail_p  tail pointer of the ep linked list (modified)
 */
static void save_element_to_ep_list(int ecp, int img_ct, FRAG_t ele_frag,
                                    EPROT_t **all_epp, EPROT_t **ep_tail_p) {
  EPROT_t *ep_tmp = (EPROT_t *) malloc(sizeof(EPROT_t));
  ep_tmp->flag   = 0;
  ep_tmp->index  = ecp;
  ep_tmp->img_no = img_ct;
  ep_tmp->frag   = ele_frag;
  ep_tmp->next   = NULL;
  if (*all_epp) (*ep_tail_p)->next = ep_tmp;
  else          *all_epp = ep_tmp;
  *ep_tail_p = ep_tmp;
  fprintf(naive_eles, "%d %s %d %d\n", ecp, ele_frag.seq_name, ele_frag.lb, ele_frag.rb);
}


/*
 * save_img_to_remaining_list  --  append an image to the remain list tail
 *
 * Copies img onto the heap and appends it to the tail of the remain list.
 * remain* = img1.next->img2.next->...->tail*
 *
 * Parameters
 *   img       IMAGE_t value to copy and append
 *   remain_p  head pointer of the remain list (modified on first append)
 *   tail_p    tail pointer of the remain list (always modified)
 */
static void save_img_to_remaining_list(IMAGE_t img,
                                       IMG_DATA_t **remain_p,
                                       IMG_DATA_t **tail_p) {
  IMAGE_t    *imgp = (IMAGE_t    *) malloc(sizeof(IMAGE_t));
  IMG_DATA_t *tmp  = (IMG_DATA_t *) malloc(sizeof(IMG_DATA_t));
  *imgp = img;
  tmp->to_image = imgp;
  tmp->next = NULL;
  if (!*remain_p) *remain_p = tmp;
  else            (*tail_p)->next = tmp;
  *tail_p = tmp;
}


/*
 * include_image  --  test and admit one remain-list image into the current element
 *
 * Applies the fractional coverage test (sing_cov or doub_cov) to the image
 * pointed to by *cur_p against the current element interval *ele_frag_p.
 *
 * If the image qualifies (cov == 1):
 *   - Unlinks *cur_p from the remain list and frees it.
 *   - Increments *img_ct_p and records the assignment in all_mprot.
 *   - Extends ele_frag_p->rb if the image reaches further right; if so,
 *     resets *cur_p to remain (re-scan from head).
 *   - Otherwise advances *cur_p past the freed node.
 *
 * If the image does not qualify (cov == 0):
 *   - Advances *prev_p and *cur_p without modifying the list.
 */
static void include_image(IMG_DATA_t **cur_p,  IMG_DATA_t **prev_p,
                          IMG_DATA_t **remain_p, IMG_DATA_t **tail_p,
                          FRAG_t *ele_frag_p, int ecp,
                          MPROT_t **all_mprot, int method, float cutoff,
                          int *img_ct_p) {
  int cov;
  IMG_DATA_t *cur = *cur_p;

  if (method == 1) cov = sing_cov(ele_frag_p, &cur->to_image->frag, cutoff);
  else             cov = doub_cov(ele_frag_p, &cur->to_image->frag, cutoff);

  if (cov) {
    (*img_ct_p)++;
    if (cur->to_image->index % 2)
      all_mprot[cur->to_image->index / 2]->se = ecp;
    else
      all_mprot[cur->to_image->index / 2]->pe = ecp;
    fprintf(img_prot, "%d %d\n", ecp, cur->to_image->index);

    /* Unlink cur from the remain list */
    if (!cur->next) *tail_p = *prev_p;
    if (*prev_p) (*prev_p)->next = cur->next;
    else         *remain_p = cur->next;

    if (ele_frag_p->rb < cur->to_image->frag.rb) {
      /* Element extended: restart remain scan from head */
      ele_frag_p->rb = cur->to_image->frag.rb;
      free(cur->to_image);
      free(cur);
      *cur_p  = *remain_p;
      *prev_p = NULL;
    } else {
      /* No extension: advance past freed node */
      free(cur->to_image);
      free(cur);
      *cur_p = (*prev_p) ? (*prev_p)->next : *remain_p;
    }
  } else {
    /* Image did not qualify: advance without modifying the list */
    *prev_p = cur;
    *cur_p  = cur->next;
  }
}


/*
 * init_element_from_remain_list  --  seed a new element from the remain list head
 *
 * Unlinks the head of the remain list (pointed to by *cur_p, which must equal
 * *remain_p), initializes a new element from that image, and advances the
 * remain pointer to the next node.
 *
 * Precondition:  *remain_p is non-NULL; *cur_p == *remain_p.
 * Postconditions:
 *   - *ritetime_p set to 0
 *   - *ecp incremented; *img_ct_p set to 1
 *   - *ele_frag_p initialized from the unlinked image
 *   - The image is recorded in all_mprot and written to img_prot
 *   - *remain_p, *cur_p advanced to the next node; *prev_p set to NULL
 */
static void init_element_from_remain_list(short *ritetime_p, int *ecp,
                                          int *img_ct_p, FRAG_t *ele_frag_p,
                                          MPROT_t **all_mprot,
                                          IMG_DATA_t **remain_p,
                                          IMG_DATA_t **cur_p,
                                          IMG_DATA_t **prev_p) {
  IMG_DATA_t *cur = *cur_p;   /* == *remain_p at entry */

  *ritetime_p = 0;
  (*ecp)++;
  *img_ct_p   = 1;
  *ele_frag_p = cur->to_image->frag;

  if (cur->to_image->index % 2)
    all_mprot[cur->to_image->index / 2]->se = *ecp;
  else
    all_mprot[cur->to_image->index / 2]->pe = *ecp;
  fprintf(img_prot, "%d %d\n", *ecp, cur->to_image->index);

  *remain_p = cur->next;
  free(cur->to_image);
  free(cur);
  *prev_p = NULL;
  *cur_p  = *remain_p;
}


/*
 * evaluate_remaining_image  --  process one remain-list image after a new element is started
 *
 * If the image overlaps the current element (same sequence, lb within rb by
 * MIN_OVERLAP_BP), delegates to include_image().
 * Otherwise, saves the current element and signals the caller to break out of
 * the inner remain-processing loop.
 *
 * Returns 1 if the caller should break (element boundary found), 0 to continue.
 */
static int evaluate_remaining_image(IMG_DATA_t **cur_p,  IMG_DATA_t **prev_p,
                                    IMG_DATA_t **remain_p, IMG_DATA_t **tail_p,
                                    FRAG_t *ele_frag_p, int ecp,
                                    MPROT_t **all_mprot, int method, float cutoff,
                                    int *img_ct_p, EPROT_t **all_epp,
                                    EPROT_t **ep_tail_p, short *ritetime_p) {
  if (ele_frag_p->seq_name == (*cur_p)->to_image->frag.seq_name &&
      ele_frag_p->rb - (*cur_p)->to_image->frag.lb > MIN_OVERLAP_BP) {
    include_image(cur_p, prev_p, remain_p, tail_p, ele_frag_p, ecp,
                  all_mprot, method, cutoff, img_ct_p);
    return 0;
  } else {
    save_element_to_ep_list(ecp, *img_ct_p, *ele_frag_p, all_epp, ep_tail_p);
    *ritetime_p = 1;
    return 1;   /* caller should break */
  }
}


/*
 * dumbbell_flush  --  flush one batch of IPROT_t records to the element database
 *
 * 1. Calls img_charge() to sort iprot_shadow by image index and load the
 *    corresponding MSP data from msp_file into each record's to_msp.
 * 2. Iterates over all_iprot[0..iprot_ct-1] (original insertion order) and
 *    writes each record to its element buffer.  When the element index changes,
 *    the previous element's buffer is flushed to ele_db.
 *
 * Note: the loop iterates over all_iprot (not iprot_shadow) to preserve the
 * original MSP-line order within each element; img_charge operates on
 * iprot_shadow but fills the shared IPROT_t objects that all_iprot also points to.
 *
 * Uses globals: size_list (FILE*), cur_ele_fp, cur_ele_buf, cur_ele_buf_size,
 *               cur_ele_index, g_ep_array.
 * The cur_ele_* state persists across batch calls -- the caller must flush
 * the last element after the final batch.
 */

/* Static state for the dumbbell pass -- persists across batch calls */
static FILE  *cur_ele_fp       = NULL;
static char  *cur_ele_buf      = NULL;
static size_t cur_ele_buf_size = 0;
static int    cur_ele_index    = -1;
static EPROT_t **g_ep_array    = NULL;  /* set by dumbbell_init */

static void dumbbell_init(EPROT_t **ep_array) {
  g_ep_array = ep_array;
}

static void dumbbell_flush(IPROT_t **all_iprot, IPROT_t **iprot_shadow,
                           int iprot_ct, MPROT_t **all_mprot,
                           FILE *msp_file) {
  int  i, partner_index;
  IPROT_t *ip;
  EPROT_t *ep;

  /* Reset shadow to insertion order before img_charge() sorts it */
  for (i = 0; i < iprot_ct; i++) iprot_shadow[i] = all_iprot[i];

  img_charge(iprot_shadow, iprot_ct, msp_file);

  for (i = 0; i < iprot_ct; i++) {
    ip = all_iprot[i];

    if (!ip->to_msp->score) continue;   /* scan_msp() failed for this image */

    /* When element changes, flush current buffer to ele_db */
    if (ip->ele_index != cur_ele_index) {
      if (cur_ele_fp) {
        fclose(cur_ele_fp);
        ele_db_write(cur_ele_index, cur_ele_buf, (int)cur_ele_buf_size);
        free(cur_ele_buf);
        cur_ele_buf      = NULL;
        cur_ele_buf_size = 0;
        cur_ele_fp       = NULL;
      }
      /* Open a new buffer for this element */
      cur_ele_index = ip->ele_index;
      ep = g_ep_array[ip->ele_index - 1];
      cur_ele_fp = open_memstream(&cur_ele_buf, &cur_ele_buf_size);
      if (!cur_ele_fp) { perror("dumbbell_flush: open_memstream"); exit(1); }
      fprintf(cur_ele_fp, "index %d\n",      ip->ele_index);
      fprintf(cur_ele_fp, "frag %s %d %d\n", ep->frag.seq_name, ep->frag.lb, ep->frag.rb);
      fprintf(cur_ele_fp, "img_no %d\n",     ep->img_no);
      fprintf(size_list,  "%d %d\n",         ip->ele_index, ep->img_no);
    }

    /* Compute partner element index */
    if (ip->index % 2)
      partner_index = all_mprot[ip->index / 2]->pe;
    else
      partner_index = all_mprot[ip->index / 2]->se;

    /*
     * Write the MSP line.
     *   odd  image: query-ele = partner_index, sbjct-ele = ip->ele_index
     *   even image: query-ele = ip->ele_index, sbjct-ele = partner_index
     * Self-MSPs only written for odd images to avoid duplicates.
     */
    if (partner_index != ip->ele_index) {
      if (ip->index % 2)
        fprintf(cur_ele_fp, "msp %d s %d %.1f %d %d %s %d %d %d %s %d %d\n",
                ip->index, ip->to_msp->score, ip->to_msp->iden, ip->to_msp->direction,
                partner_index,   ip->to_msp->query.frag.seq_name,
                ip->to_msp->query.frag.lb,  ip->to_msp->query.frag.rb,
                ip->ele_index,   ip->to_msp->sbjct.frag.seq_name,
                ip->to_msp->sbjct.frag.lb, ip->to_msp->sbjct.frag.rb);
      else
        fprintf(cur_ele_fp, "msp %d s %d %.1f %d %d %s %d %d %d %s %d %d\n",
                ip->index, ip->to_msp->score, ip->to_msp->iden, ip->to_msp->direction,
                ip->ele_index,   ip->to_msp->query.frag.seq_name,
                ip->to_msp->query.frag.lb,  ip->to_msp->query.frag.rb,
                partner_index,   ip->to_msp->sbjct.frag.seq_name,
                ip->to_msp->sbjct.frag.lb, ip->to_msp->sbjct.frag.rb);
    } else {
      if (ip->index % 2)
        fprintf(cur_ele_fp, "msp %d s %d %.1f %d %d %s %d %d %d %s %d %d\n",
                ip->index, ip->to_msp->score, ip->to_msp->iden, ip->to_msp->direction,
                partner_index,   ip->to_msp->query.frag.seq_name,
                ip->to_msp->query.frag.lb,  ip->to_msp->query.frag.rb,
                ip->ele_index,   ip->to_msp->sbjct.frag.seq_name,
                ip->to_msp->sbjct.frag.lb, ip->to_msp->sbjct.frag.rb);
    }
  }
  /* Do NOT flush cur_ele_fp here -- the same element may continue in the next batch */
}


int main (int argc, char *argv[]) {
  int ele_ct=0, msp_ct=0, i, method;
  float cutoff;
  char line[150], *m1="single", *m2="double";
  MPROT_t **all_mprot;
  EPROT_t *all_ep=NULL, **ep_array, *ep_tmp;

  IPROT_t **all_iprot, **iprot_shadow;
  int iprot_ct=0;

  FILE *seq_list, *msp_file;

  /* IMAGE array for the merged imagespread pass */
  IMAGE_t  *images = NULL;
  int       n_images = 0, img_capacity = 0;
  MSP_t     cur_msp;
  int       img_idx;
  int       pos;

  /* Check for -v (version) before any other parsing */
  recon_check_version_flag(argc, argv);

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
           "3=info(default) 4=debug\n"
           "-v          print version and exit\n");
    exit(1);
  }

  seq_list = fopen(argv[1], "r");
  if (!seq_list) {
    printf("Input file for sequence name list %s not found.  Exit.\n", argv[1]);
    exit(2);
  }
  GetSeqNames(seq_list);
  fclose(seq_list);

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
    if (method == 1) { cutoff = ELEDEF_CUTOFF_SINGLE; }
    else             { cutoff = ELEDEF_CUTOFF_DOUBLE; }
  }

  /* Open auxiliary output files */
  err       = fopen("ele_def_res/errors",    "w");
  naive_eles= fopen("summary/naive_eles",    "w");
  img_prot  = fopen("ele_def_res/img_prot",  "w");
  ele_no    = fopen("summary/naive_ele_no",  "w");
  size_list = fopen("ele_def_res/size_list", "w");

  /* Route RECON_LOG output to the error file */
  recon_log_fp = err;

  /* Open the element database */
  ele_db_open();

  /* ============================================================
   * Stage 1 (merged imagespread): read MSP file, create and sort
   * IMAGE array.
   *
   * Each MSP produces two IMAGE records:
   *   even index (2k)   = query side
   *   odd  index (2k+1) = subject side
   *
   * Sorted by: seq_name (string), lb (ascending), rb (descending),
   * index (ascending for stability).
   * ============================================================ */
  img_capacity = 8192;
  images = (IMAGE_t *) malloc(img_capacity * sizeof(IMAGE_t));
  if (!images) { perror("eledef: malloc images"); exit(1); }

  img_idx = -1;   /* starts at -1 so first ++ gives 0 */

  while (fgets(line, 150, msp_file)) {
    if (scan_msp(&cur_msp, line)) {
      RLOG_ERR("Wrong MSP format, skipping line\n");
      continue;
    }

    if (n_images + 2 > img_capacity) {
      img_capacity *= 2;
      images = (IMAGE_t *) realloc(images, img_capacity * sizeof(IMAGE_t));
      if (!images) { perror("eledef: realloc images"); exit(1); }
    }

    /* Even image: query side */
    img_idx++;
    pos = GetSeqIndex(0, seq_count - 1, cur_msp.query.frag.seq_name);
    if (pos < 0) {
      RLOG_ERR("Sequence %s not found in seq list\n", cur_msp.query.frag.seq_name);
      exit(4);
    }
    images[n_images].index        = img_idx;
    images[n_images].frag.lb      = cur_msp.query.frag.lb;
    images[n_images].frag.rb      = cur_msp.query.frag.rb;
    images[n_images].frag.seq_name = seq_name_table[pos];
    n_images++;

    /* Odd image: subject side */
    img_idx++;
    pos = GetSeqIndex(0, seq_count - 1, cur_msp.sbjct.frag.seq_name);
    if (pos < 0) {
      RLOG_ERR("Sequence %s not found in seq list\n", cur_msp.sbjct.frag.seq_name);
      exit(4);
    }
    images[n_images].index        = img_idx;
    images[n_images].frag.lb      = cur_msp.sbjct.frag.lb;
    images[n_images].frag.rb      = cur_msp.sbjct.frag.rb;
    images[n_images].frag.seq_name = seq_name_table[pos];
    n_images++;

    msp_ct++;
  }

  qsort(images, n_images, sizeof(IMAGE_t), image_sort_cmp);

  /* Write ori_msp_no -- read by eleredef to seed msp_index (the counter used
   * to assign non-colliding indices to MSPs created during element dissection).
   * Without this, new MSP indices start at 0 and collide with existing ones. */
  {
    FILE *ori_msp_no = fopen("summary/ori_msp_no", "w");
    if (ori_msp_no) { fprintf(ori_msp_no, "%d\n", msp_ct); fclose(ori_msp_no); }
  }

  /* Rewind msp_file so img_charge() can do its sequential scan in pass 2 */
  rewind(msp_file);

  /* ============================================================
   * Allocate MPROT_t array (one entry per MSP)
   * ============================================================ */
  all_mprot = (MPROT_t **) malloc(msp_ct * sizeof(MPROT_t *));
  for (i = 0; i < msp_ct; i++) {
    all_mprot[i] = (MPROT_t *) malloc(sizeof(MPROT_t));
    all_mprot[i]->pe = 0;
    all_mprot[i]->se = 0;
  }

  /* ============================================================
   * Pass 1: single-linkage clustering of images into elements
   * ============================================================ */
  ele_def(method, images, n_images, cutoff, &all_ep, &ele_ct, all_mprot);

  free(images);
  images = NULL;

  fclose(ele_no);
  fclose(img_prot);

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
   * Allocate IPROT_t arrays for the dumbbell pass.
   * ELEDEF_IMGPROT_CAP sets the batch size.
   */
  all_iprot    = (IPROT_t **) malloc(ELEDEF_IMGPROT_CAP * sizeof(IPROT_t *));
  iprot_shadow = (IPROT_t **) malloc(ELEDEF_IMGPROT_CAP * sizeof(IPROT_t *));
  for (i = 0; i < ELEDEF_IMGPROT_CAP; i++) {
    all_iprot[i]          = (IPROT_t *) malloc(sizeof(IPROT_t));
    all_iprot[i]->to_msp  = (MSP_t  *) malloc(sizeof(MSP_t));
    iprot_shadow[i]       = all_iprot[i];
  }

  /* ============================================================
   * Pass 2: dumbbell -- load MSP detail and write element records
   * ============================================================ */
  dumbbell_init(ep_array);

  while (fgets(line, 100, img_prot)) {
    sscanf(line, "%d %d\n",
           &all_iprot[iprot_ct]->ele_index,
           &all_iprot[iprot_ct]->index);
    iprot_ct++;
    if (iprot_ct == ELEDEF_IMGPROT_CAP) {
      dumbbell_flush(all_iprot, iprot_shadow, iprot_ct, all_mprot, msp_file);
      iprot_ct = 0;
    }
  }

  if (iprot_ct) {
    dumbbell_flush(all_iprot, iprot_shadow, iprot_ct, all_mprot, msp_file);
  }

  /* Flush the last element's buffer to ele_db */
  if (cur_ele_fp) {
    fclose(cur_ele_fp);
    ele_db_write(cur_ele_index, cur_ele_buf, (int)cur_ele_buf_size);
    free(cur_ele_buf);
    cur_ele_fp = NULL;
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

  ele_db_close();

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
 *   images     pre-sorted array of IMAGE_t records (from the in-memory sort)
 *   n_images   number of entries in the images array
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
void ele_def(int method, IMAGE_t *images, int n_images, float cutoff,
             EPROT_t **all_epp, int *ecp, MPROT_t **all_mprot) {
  int img_ct, cov;
  FRAG_t ele_frag;
  IMAGE_t img;
  IMG_DATA_t *cur, *prev=NULL, *remain=NULL, *tail=NULL;
  short ritetime;      /* was: ritetime -- "right time" to start a new element */
  EPROT_t *ep_tail=NULL;
  int ii;   /* loop counter over the images array */

  ritetime = 1;   /* initially ready to start a new element */

  /*
   * Iterate over the pre-sorted IMAGE array.
   *
   * Image index convention:
   *   even index (img.index % 2 == 0) => query side of the parent MSP
   *   odd  index (img.index % 2 == 1) => subject side of the parent MSP
   *   MSP position in file            => img.index / 2
   */
  for (ii = 0; ii < n_images; ii++) {
    img = images[ii];

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

        if (ele_frag.rb < img.frag.rb) {
          /* Element boundary extended: re-scan remain list */
          ele_frag.rb = img.frag.rb;
          cur  = remain;
          prev = NULL;
          while (cur) {
            include_image(&cur, &prev, &remain, &tail, &ele_frag, *ecp,
                          all_mprot, method, cutoff, &img_ct);
          }
        }
      } else {
        /* Image did not qualify: keep in remain list for later */
        save_img_to_remaining_list(img, &remain, &tail);
      }
    } else {
      /* Image falls outside this element: close element, process remain */
      save_element_to_ep_list(*ecp, img_ct, ele_frag, all_epp, &ep_tail);
      save_img_to_remaining_list(img, &remain, &tail);
      ritetime = 1;

      /*
       * Process the remain list: try to start new elements from the
       * buffered images that were not incorporated into the element
       * that was just closed.
       */
      while (remain && ritetime) {
        cur = remain;
        init_element_from_remain_list(&ritetime, ecp, &img_ct, &ele_frag,
                                      all_mprot, &remain, &cur, &prev);
        /* cur now points to the rest of the remain list */
        while (cur) {
          if (evaluate_remaining_image(&cur, &prev, &remain, &tail, &ele_frag,
                                       *ecp, all_mprot, method, cutoff, &img_ct,
                                       all_epp, &ep_tail, &ritetime)) break;
        }
      }
    }
  }   /* end of main image loop */

  save_element_to_ep_list(*ecp, img_ct, ele_frag, all_epp, &ep_tail);
  if (remain) ritetime = 1;

  /* Drain any images still in the remain list into their own elements */
  while (remain) {
    cur  = remain;
    prev = NULL;

    while (cur) {
      if (ritetime) {
        init_element_from_remain_list(&ritetime, ecp, &img_ct, &ele_frag,
                                      all_mprot, &remain, &cur, &prev);
        continue;
      }
      if (evaluate_remaining_image(&cur, &prev, &remain, &tail, &ele_frag,
                                   *ecp, all_mprot, method, cutoff, &img_ct,
                                   all_epp, &ep_tail, &ritetime)) break;
    }
    if (!cur) {
      save_element_to_ep_list(*ecp, img_ct, ele_frag, all_epp, &ep_tail);
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
