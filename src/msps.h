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


/* ============================================================
 * print_msp / fprint_msp
 *
 * Serialize an MSP_t to stdout or a FILE* in the canonical MSP
 * file format:
 *
 *   score %iden q_start q_end query_name s_start s_end sbjct_name
 *
 * For a reversed alignment (direction == -1) the query coordinates are
 * written in original (reversed) order so that the output faithfully
 * reproduces the original MSP file format: q_start > q_end signals
 * reverse orientation to downstream tools.
 *
 * Returns 0 on success, 1 on write error.
 * ============================================================ */
int print_msp(MSP_t *m) {
  if (m->direction == 1) {
    return (printf("%06d %3.1f %08d %08d %s %08d %08d %s \n",
                   m->score, m->iden,
                   m->query.frag.lb, m->query.frag.rb, m->query.frag.seq_name,
                   m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name)
            > 0) ? 0 : 1;
  } else if (m->direction == -1) {
    /* Reverse the query coordinates to signal orientation in the output */
    return (printf("%06d %3.1f %08d %08d %s %08d %08d %s \n",
                   m->score, m->iden,
                   m->query.frag.rb, m->query.frag.lb, m->query.frag.seq_name,
                   m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name)
            > 0) ? 0 : 1;
  } else {
    return 1;
  }
}


int fprint_msp(FILE *ofp, MSP_t *m) {
  if (m->direction == 1) {
    return (fprintf(ofp, "%06d %3.1f %08d %08d %s %08d %08d %s \n",
                    m->score, m->iden,
                    m->query.frag.lb, m->query.frag.rb, m->query.frag.seq_name,
                    m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name)
            > 0) ? 0 : 1;
  } else if (m->direction == -1) {
    return (fprintf(ofp, "%06d %3.1f %08d %08d %s %08d %08d %s \n",
                    m->score, m->iden,
                    m->query.frag.rb, m->query.frag.lb, m->query.frag.seq_name,
                    m->sbjct.frag.lb, m->sbjct.frag.rb, m->sbjct.frag.seq_name)
            > 0) ? 0 : 1;
  } else {
    return 1;
  }
}


/* ============================================================
 * scan_msp  --  parse one MSP file line into an MSP_t
 *
 * The raw input may encode a reversed alignment as q_start > q_end
 * and/or s_start > s_end.  scan_msp normalises both ranges to
 * lb <= rb and records the combined orientation in MSP_t.direction.
 *
 * Direction encoding:
 *   Each out-of-order range multiplies direction by -1.  So:
 *     - Both ranges in order   => direction = +1  (forward)
 *     - Query reversed only    => direction = -1  (reverse)
 *     - Subject reversed only  => direction = -1  (reverse)
 *     - Both reversed          => direction = +1  (forward, double-flip)
 *
 * Sequence-name interning:
 *   The raw name strings are looked up in seq_name_table[] via
 *   GetSeqIndex().  MSP_t stores a *pointer into that table*; never
 *   an allocated copy.  This enables O(1) pointer-equality tests.
 *
 * Returns 0 on success, non-zero if parsing fails or a name is not
 * found in seq_name_table[].
 * ============================================================ */
int scan_msp(MSP_t *m, char *line) {
  int32_t coord_swap_tmp;   /* was: bd_tmp -- temp for lb/rb swap during normalisation */
  char qname[SEQ_NAME_MAX_LEN], sname[SEQ_NAME_MAX_LEN];
  int pos;

  /* Use %d for int32_t fields (int32_t is int on this platform; the original
   * used %ld which caused warnings on 32-bit but worked on 64-bit due to
   * the sizeof(long)==sizeof(int) assumption -- now fixed to %d). */
  if (sscanf(line,
             "%d %f %d %d %s %d %d %s \n",
             &(m->score), &(m->iden),
             &(m->query.frag.lb), &(m->query.frag.rb), qname,
             &(m->sbjct.frag.lb), &(m->sbjct.frag.rb), sname) != 8) {
    return 1;
  }

  /* Look up and intern query sequence name */
  pos = GetSeqIndex(0, seq_count - 1, qname);
  if (pos < 0) return pos;
  m->query.frag.seq_name = seq_name_table[pos];

  /* Look up and intern subject sequence name */
  pos = GetSeqIndex(0, seq_count - 1, sname);
  if (pos < 0) return pos;
  m->sbjct.frag.seq_name = seq_name_table[pos];

  /* Initialise non-coordinate fields */
  m->direction = 1;
  m->stat      = 's';
  m->query.to_msp = m;
  m->sbjct.to_msp = m;

  /*
   * Normalise query range: if lb > rb the alignment is reversed on the
   * query strand.  Swap the coordinates and flip direction.
   */
  if (m->query.frag.lb > m->query.frag.rb) {
    m->direction      *= -1;
    coord_swap_tmp     = m->query.frag.lb;
    m->query.frag.lb   = m->query.frag.rb;
    m->query.frag.rb   = coord_swap_tmp;
  }

  /*
   * Normalise subject range: same logic.  A second flip would cancel
   * the first, giving direction = +1 again (double-reverse = forward).
   */
  if (m->sbjct.frag.lb > m->sbjct.frag.rb) {
    m->direction      *= -1;
    coord_swap_tmp     = m->sbjct.frag.lb;
    m->sbjct.frag.lb   = m->sbjct.frag.rb;
    m->sbjct.frag.rb   = coord_swap_tmp;
  }

  return 0;
}


/* ============================================================
 * partner  --  given one IMAGE_t of an MSP, return the other
 *
 * Because query and sbjct are embedded in MSP_t (not pointer-addressed),
 * the address of an IMAGE_t relative to its parent MSP tells us which
 * side it is.
 * ============================================================ */
IMAGE_t *partner(IMAGE_t *img) {
  if (img == &(img->to_msp->query)) return &(img->to_msp->sbjct);
  return &(img->to_msp->query);
}


/* ============================================================
 * sing_cov  --  single-coverage overlap test
 *
 * Returns 1 if f1 and f2 are on the same sequence AND the length of
 * their overlap is >= cutoff * length of *either* fragment.
 * Returns 0 otherwise.
 *
 * "Single coverage" means the overlap only needs to cover one of the
 * two fragments significantly; the other may be much longer.
 *
 * Parameters
 *   f1, f2   two FRAG_t intervals (seq_name must be an interned ptr)
 *   cutoff   fractional overlap threshold in [0, 1]
 * ============================================================ */
int sing_cov(FRAG_t *f1, FRAG_t *f2, float cutoff) {
  int32_t len1, len2, overlap_len, overlap_lb, overlap_rb;

  /* Pointer equality: both must be on the same interned sequence name */
  if (f1->seq_name == f2->seq_name) {
    len1       = f1->rb - f1->lb;
    len2       = f2->rb - f2->lb;
    overlap_lb = f1->lb > f2->lb ? f1->lb : f2->lb;
    overlap_rb = f1->rb < f2->rb ? f1->rb : f2->rb;
    overlap_len = overlap_rb - overlap_lb;
    if ((float) overlap_len / len1 >= cutoff ||
        (float) overlap_len / len2 >= cutoff) {
      return 1;
    }
  }
  return 0;
}


/* ============================================================
 * doub_cov  --  double-coverage overlap test
 *
 * Returns 1 if f1 and f2 are on the same sequence AND the overlap is
 * >= cutoff * length of *both* fragments.
 * Returns 0 otherwise.
 *
 * "Double coverage" requires both fragments to be significantly
 * covered by the overlap; used in the initial clustering step when
 * the "double" method is selected (less commonly used than sing_cov).
 * ============================================================ */
int doub_cov(FRAG_t *f1, FRAG_t *f2, float cutoff) {
  int32_t len1, len2, overlap_len, overlap_lb, overlap_rb;

  if (f1->seq_name == f2->seq_name) {
    len1       = f1->rb - f1->lb;
    len2       = f2->rb - f2->lb;
    overlap_lb = f1->lb > f2->lb ? f1->lb : f2->lb;
    overlap_rb = f1->rb < f2->rb ? f1->rb : f2->rb;
    overlap_len = overlap_rb - overlap_lb;
    if ((float) overlap_len / len1 >= cutoff &&
        (float) overlap_len / len2 >= cutoff) {
      return 1;
    }
  }
  return 0;
}


/* ============================================================
 * toString-style print helpers for MSP-layer types
 * ============================================================ */

/*
 * print_frag  --  one-line summary of a FRAG_t to stdout.
 */
void print_frag(FRAG_t *f) {
  if (!f) { printf("FRAG_t: NULL\n"); return; }
  printf("FRAG_t: seq_name=%s, lb=%d, rb=%d\n",
         f->seq_name ? f->seq_name : "(null)", f->lb, f->rb);
}


/*
 * print_image  --  one-line summary of an IMAGE_t to stdout.
 * Prints the image index, whether it is a query (even) or subject (odd)
 * image, the fragment interval, and whether ele_info has been assigned.
 * (The ele_info index is not printed here because struct ele_info is
 * defined in ele.h which includes msps.h, not the reverse; use
 * print_ele_info() from ele.h for full element detail.)
 */
void print_image(IMAGE_t *img) {
  if (!img) { printf("IMAGE_t: NULL\n"); return; }
  printf("IMAGE_t: index=%d (%s), seq=%s, lb=%d, rb=%d, ele_info=%s\n",
         img->index,
         (img->index % 2 == 0) ? "query" : "sbjct",
         img->frag.seq_name ? img->frag.seq_name : "(null)",
         img->frag.lb, img->frag.rb,
         img->ele_info ? "(assigned)" : "NULL");
}


/*
 * print_msp_full  --  multi-field summary of an MSP_t to stdout.
 * More verbose than print_msp(); intended for debugging.
 */
void print_msp_full(MSP_t *m) {
  if (!m) { printf("MSP_t: NULL\n"); return; }
  printf("MSP_t: stat=%c, score=%d, iden=%.1f, direction=%d\n",
         m->stat, m->score, m->iden, m->direction);
  printf("  query: "); print_image(&m->query);
  printf("  sbjct: "); print_image(&m->sbjct);
}


/*
 * print_img_data  --  walk an IMG_DATA_t linked list and print each image.
 */
void print_img_data(IMG_DATA_t *id) {
  int count = 0;
  while (id) {
    printf("  [%d] ", count++);
    print_image(id->to_image);
    id = id->next;
  }
}


/*
 * print_img_tree_inorder  --  in-order traversal of an IMG_TREE_t,
 * printing each node's image.  The 'level' parameter (start at 0)
 * is used for indentation.
 */
void print_img_tree_inorder(IMG_TREE_t *rt, int level) {
  int i;
  if (!rt) return;
  print_img_tree_inorder(rt->l, level + 1);
  for (i = 0; i < level; i++) printf("  ");
  print_image(rt->to_image);
  print_img_tree_inorder(rt->r, level + 1);
}


#endif /* __MSPS_H__ */
