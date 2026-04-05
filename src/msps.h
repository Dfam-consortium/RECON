/*
 * msps.h  --  MSP (Maximal Scoring Pair) data structures and operations
 *
 * An MSP is a single high-scoring local alignment produced by a pairwise
 * sequence comparison tool (typically BLAST).  RECON converts a BLAST
 * output file into an "MSP file" where each line encodes one alignment:
 *
 *   score  %iden  q_start  q_end  query_name  s_start  s_end  sbjct_name
 *
 * Coordinate conventions
 * ----------------------
 * Coordinates are 1-based, fully closed (both endpoints are included).
 * A reversed alignment on the query strand is represented by having
 * q_start > q_end in the raw input; scan_msp() normalises both ranges
 * to lb < rb and tracks orientation in MSP_t.direction.
 *
 * Terminology: "image"
 * --------------------
 * Because each MSP produces *two* sequence intervals (one on the query
 * and one on the subject), each half is called an IMAGE.  An IMAGE_t
 * bundles a sequence interval (FRAG_t) with a back-pointer to its MSP
 * and a forward pointer to the element it has been assigned to.  Images
 * are the primary unit of work in the eledef and eleredef stages.
 *
 * Image index convention
 * ----------------------
 * Each MSP is assigned a sequential 0-based index when read.  Its two
 * images are given indices 2*msp_index (query side) and 2*msp_index+1
 * (subject side).  Even image index => query image; odd => subject image.
 * This encoding is used throughout eledef.c and eleredef.c.
 *
 * Data structures
 * ---------------
 *   FRAG_t        -- a sequence interval: name pointer + lb + rb
 *   IMAGE_t       -- an MSP half: index, FRAG_t, back-ptr to MSP,
 *                    forward-ptr to element
 *   MSP_t         -- a complete alignment: score, iden, direction,
 *                    plus query IMAGE_t and sbjct IMAGE_t embedded
 *                    directly (not as pointers) so that IMAGE_t.to_msp
 *                    can navigate back to the parent MSP
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#ifndef __MSPS_H__
#define __MSPS_H__

/* seqlist.h provides GetSeqIndex() and GetSeqNames(), used by scan_msp().
 * It also pulls in bolts.h (seq_count, seq_name_table) and recon_defs.h. */
#include "seqlist.h"


/* ============================================================
 * Core data structures
 * ============================================================ */

/*
 * FRAG_t  --  a half-open sequence interval
 *
 * seq_name  pointer into seq_name_table[]; used for identity (ptr==ptr)
 *           comparison rather than strcmp -- do not store allocated copies.
 * lb        left boundary (start), 1-based, always <= rb after normalisation.
 * rb        right boundary (end),  1-based, always >= lb after normalisation.
 */
typedef struct frag {
  char   *seq_name;
  int32_t lb, rb;
} FRAG_t;

/* A singly-linked list node wrapping a FRAG_t pointer. */
typedef struct frag_list {
  FRAG_t             *to_frag;
  struct frag_list   *next;
} FRAG_DATA_t;


/*
 * IMAGE_t  --  one half of an MSP (either the query or subject interval)
 *
 * index     unique image index; even = query side, odd = subject side.
 *           index / 2 gives the parent MSP's sequential position in the
 *           input file.
 * frag      the sequence interval covered by this image.
 * to_msp    back-pointer to the parent MSP_t; always non-NULL for a
 *           valid image.
 * ele_info  forward-pointer to the ELE_INFO_t of the element this image
 *           has been assigned to; NULL until clustering assigns it.
 */
typedef struct image {
  int              index;
  FRAG_t           frag;
  struct msp      *to_msp;
  struct ele_info *ele_info;
} IMAGE_t;

/* A singly-linked list node wrapping an IMAGE_t pointer. */
typedef struct image_list {
  IMAGE_t            *to_image;
  struct image_list  *next;
} IMG_DATA_t;

/* A binary-search-tree node wrapping an IMAGE_t pointer.
 * Used in ELEMENT_t.to_img_tree for O(log n) image lookup by index. */
typedef struct img_tree {
  IMAGE_t          *to_image;
  struct img_tree  *p, *l, *r;   /* parent, left child, right child */
} IMG_TREE_t;


/*
 * MSP_t  --  a complete pairwise alignment record
 *
 * stat      processing state of this MSP:
 *             'p' = highest-scoring MSP for an element pair, spanning
 *                   the full length of at least one element.
 *             's' = standard MSP (default on read).
 * score     raw alignment score from the comparison tool.
 * iden      percent identity (0.0 -- 100.0).
 * direction orientation of the alignment:
 *             +1 = both strands in the same direction (forward/forward
 *                  or reverse/reverse, depending on tool convention).
 *             -1 = strands in opposite directions.
 * query     IMAGE_t for the query sequence half.
 * sbjct     IMAGE_t for the subject sequence half.
 *
 * Note: query and sbjct are *embedded* (not pointers) so that the
 * address of msp->query or msp->sbjct can be compared with a
 * IMAGE_t pointer to determine which side of the MSP an image
 * belongs to (see partner() below).
 */
typedef struct msp {
  char    stat;
  int32_t score;
  float   iden;
  int     direction;
  IMAGE_t query, sbjct;
} MSP_t;

/* A singly-linked list node wrapping an MSP_t pointer. */
typedef struct msp_list {
  MSP_t            *to_msp;
  struct msp_list  *next;
} MSP_DATA_t;


/* ============================================================
 * Function prototypes
 * ============================================================ */

int      print_msp(MSP_t *);
int      fprint_msp(FILE *, MSP_t *);
int      scan_msp(MSP_t *, char *);
IMAGE_t *partner(IMAGE_t *);
int      doub_cov(FRAG_t *, FRAG_t *, float);
int      sing_cov(FRAG_t *, FRAG_t *, float);

/* toString-style helpers (definitions below) */
void print_frag(FRAG_t *);
void print_image(IMAGE_t *);
void print_msp_full(MSP_t *);
void print_img_data(IMG_DATA_t *);
void print_img_tree_inorder(IMG_TREE_t *, int);




#endif /* __MSPS_H__ */
