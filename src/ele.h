/*
 * ele.h  --  Element, edge, and family data structures for the RECON pipeline
 *
 * This header is the primary shared library for stages 3-5 of the RECON
 * pipeline (eleredef, edgeredef, famdef).  It defines all data types used
 * to represent:
 *   - Elements: genomic intervals that result from clustering MSP images
 *   - Edges: similarity relationships between elements
 *   - Families: clusters of related elements
 *
 * It also provides all functions for initializing, reading, writing,
 * searching, and freeing these structures, along with a set of toString-
 * style print helpers for debugging.
 *
 * Element file format (tmp/e<N> and tmp2/e<N>)
 * --------------------------------------------
 * Each element is serialized to a text file with keyword-prefixed lines:
 *
 *   frag    <seq_name> <lb> <rb>
 *   direc   <direction>
 *   update  <update_flag>
 *   img_no  <image_count>
 *   flimg_no <full_len_image_count>
 *   edge_no <edge_count>
 *   msp     <img_id> <stat> <score> <iden> <dir> <ele1> <qseq> <qlb> <qrb> \
 *             <ele2> <sseq> <slb> <srb>
 *   edge    <idx> <type> <dir> <score> <ele1_idx> <ele2_idx>
 *   pcp     <position> <contributor_ele_idx>
 *   tbd     <position> <support>
 *   redef   <child_ele_idx>
 *
 * The file for element N is read from tmp/e<N> (original) or tmp2/e<N>
 * (updated); the ele_info->file_updated flag selects which path to use.
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#ifndef ELE_H
#define ELE_H

#include "msps.h"
#include "recon_log.h"

/*
 * TANDEM_SIZE_LIMIT and TANDEM_IMG_RATIO are defined in recon_defs.h
 * (pulled in transitively via msps.h -> bolts.h -> recon_defs.h).
 * The backward-compat aliases SIZE_LIMIT and TANDEM are also defined there.
 */


/* ============================================================
 * Data structures
 * ============================================================ */

/*
 * CP_t  --  Potential Cut Point
 *
 * A linked-list node representing a candidate element-boundary position.
 * CPs are accumulated in ELEMENT_t.PCP during edges_and_cps() from the
 * endpoints of full-length images.  They are later clustered into "to-be-
 * determined" boundaries (BD_t / TBD) by PCP_to_TBDs().
 *
 * cp           genomic coordinate of the candidate cut point.
 * contributor  the element whose image endpoint produced this CP.
 */
typedef struct cp_list {
  int32_t          cp;
  struct ele_info *contributor;
  struct cp_list  *next;
} CP_t;

/*
 * BD_t  --  To-Be-Determined boundary (clustered cut point)
 *
 * A linked-list node representing a well-supported element boundary
 * candidate produced by clustering the CP_t list in CP_cluster().
 *
 * bd       consensus genomic coordinate of the cluster (weighted mean).
 * support  number of CPs in the cluster that contributed to this boundary.
 */
typedef struct bd_list {
  int32_t       bd;
  int           support;
  struct bd_list *next;
} BD_t;

/*
 * EDGE_t  --  a similarity relationship between two elements
 *
 * Edges are created during edges_and_cps() when two elements share
 * full-length (or near-full-length) MSP images.  They are used by
 * edgeredef to filter spurious relationships and by famdef to build
 * repeat families via BFS on the edge graph.
 *
 * index       unique sequential edge identifier.
 * type        classification of the edge:
 *               'p' = primary (unfiltered, used for family building)
 *               'P' = promoted primary (winner after PPS filtering)
 *               'S' = demoted secondary (loser after PPS filtering)
 *               's' = secondary (initially weaker primary candidate)
 *               'c' = cycle/self edge (invalid, filtered out)
 * direction   relative orientation of the two elements:
 *               +1 = same strand, -1 = opposite strands.
 * score       MSP alignment score; used by edge_filt() to choose the
 *             highest-scoring PPS edge when demoting.
 * ele1_info   one endpoint of the edge (lower index by convention).
 * ele2_info   the other endpoint of the edge.
 */
typedef struct edge {
  int             index;
  char            type;
  int             direction;
  int32_t         score;
  struct ele_info *ele1_info, *ele2_info;
} EDGE_t;

/* Binary-search-tree node wrapping an EDGE_t pointer.
 * Keyed on EDGE_t.index for O(log n) lookup in find_edge(). */
typedef struct edge_tree {
  EDGE_t           *to_edge;
  struct edge_tree *p, *l, *r;   /* parent, left child, right child */
} EDGE_TREE_t;

/*
 * ELEMENT_t  --  a repeat element: a genomic interval with its images and edges
 *
 * index        sequential identifier (1-based); matches file name "e<index>".
 * frag         the representative genomic interval for this element.
 * direction    orientation of the element (+1 forward, -1 reverse).
 * update       non-zero if this element's coordinates were modified during
 *              redefinition and need to be re-evaluated by its partners.
 * l_hold       depth counter used during BFS local-network construction in
 *              build_local_network().  0 = not yet in any network;
 *              1..DEPTH = in network at this BFS depth.
 *              (was: l_hold -- "local hold")
 * img_no       current count of images assigned to this element.
 * flimg_no     count of "full-length" images (images that span >= CUTOFF2
 *              of the element's representative fragment).
 * edge_no      count of edges in the edge tree.
 * to_img_tree  balanced BST of IMAGE_t pointers, keyed on image index.
 * to_img_data  singly-linked list of IMG_DATA_t; populated lazily by
 *              listify() when a sorted traversal is needed.
 * edges        balanced BST of EDGE_t pointers, keyed on edge index.
 * PCP          linked list of Potential Cut Points.
 * TBD          linked list of To-Be-Determined boundaries (post-clustering).
 * redef        linked list of ELE_DATA_t pointing to child elements
 *              created when this element is dissected.
 */
typedef struct element {
  int             index;
  FRAG_t          frag;
  int             direction;
  short           update, l_hold;
  int             img_no;
  int             flimg_no;
  int             edge_no;
  IMG_TREE_t     *to_img_tree;
  IMG_DATA_t     *to_img_data;
  EDGE_TREE_t    *edges;
  CP_t           *PCP;
  BD_t           *TBD;
  struct ele_list *redef;
} ELEMENT_t;

/*
 * ELE_INFO_t  --  lightweight metadata wrapper for an element
 *
 * This structure acts as a persistent handle for an element across the
 * pipeline stages.  The heavy ELEMENT_t is only loaded into memory when
 * needed (ele_read_in) and freed when no longer needed (ele_cleanup).
 *
 * index        element identifier; matches the e<N> file name.
 * ele          pointer to the in-memory ELEMENT_t; NULL when not loaded.
 * stat         current processing state:
 *                'z' = initial state (not yet processed)
 *                't' = edges and CPs computed (after edges_and_cps())
 *                'v' = locally redefined; no further split needed
 *                'w' = combined (redefinition complete, no duplicates)
 *                'y' = edge-filtered; ready for family building
 *                'x' = recruited into a family (famdef)
 *                'O' = large tandem repeat; removed from processing
 *                'X' = dismissed / deleted element
 * file_updated flag: 0 => read from tmp/e<N>; 1 => read from tmp2/e<N>.
 * to_family    pointer to the FAMILY_t this element belongs to (set in famdef).
 * next         intrusive-list link for the ele_info_data overflow chain.
 */
typedef struct ele_info {
  int              index;
  ELEMENT_t       *ele;
  char             stat;
  short            file_updated;
  struct family   *to_family;
  struct ele_info *next;
} ELE_INFO_t;

/* Singly-linked list node wrapping an ELE_INFO_t pointer. */
typedef struct ele_list {
  struct ele_info *ele_info;
  struct ele_list *next;
} ELE_DATA_t;


/*
 * FAMILY_t  --  a repeat family: a connected component in the element
 * similarity graph.
 *
 * index     sequential family identifier (1-based).
 * name      reserved; currently always written as "unknown" in the output.
 * members   linked list of ELE_DATA_t of all elements in this family.
 * relatives reserved; not used in the current pipeline.
 */
typedef struct family {
  int         index;
  char        name[10];
  ELE_DATA_t *members;
  struct fam_list *relatives;
} FAMILY_t;

/* Singly-linked list node wrapping a FAMILY_t pointer. */
typedef struct fam_list {
  struct family   *to_family;
  struct fam_list *next;
} FAM_DATA_t;





ELEMENT_t *ele_init(int);
ELE_INFO_t *ele_info_init(int);
MSP_t *MSP_malloc();
void MSP_free(MSP_t *);
EDGE_t *EDGE_malloc();
void EDGE_free(EDGE_t *);

ELE_INFO_t *linked_ele(ELE_INFO_t *, EDGE_t *);
void outthrow_big_tandems(FILE *);
int int_cmp(const void *, const void *);
void spit_out_ele(ELE_INFO_t *);

ELEMENT_t *ele_read_in(ELE_INFO_t *, int);
void img_scan(ELE_INFO_t *, char *, IMAGE_t **, int *);
void edge_scan(char *, EDGE_t **, ELE_INFO_t *, int *);
ELE_INFO_t *get_ele_info(int);

void ele_write_out(ELE_INFO_t *, int);
int img_index_cmp(const void *, const void *);
int edge_index_cmp(const void *, const void *);
void write_out_msps(FILE *, IMG_TREE_t *);
void write_out_edges(FILE *, EDGE_TREE_t *, ELE_INFO_t *);

void build_img_tree(IMG_TREE_t **, IMAGE_t **, int, int);
int count_img_nodes(IMG_TREE_t *);
void insert_image(IMG_TREE_t **, IMAGE_t *);
void delete_image(IMG_TREE_t **, IMAGE_t *);
IMG_TREE_t *minimal_image(IMG_TREE_t *);
IMG_TREE_t *find_image(IMG_TREE_t *, int);
IMG_DATA_t **listify(IMG_TREE_t *, IMG_DATA_t **);

void build_edge_tree(EDGE_TREE_t **, EDGE_t **, int, int);
int count_edge_nodes(EDGE_TREE_t *);
int count_total_edges(EDGE_TREE_t *);
void insert_edge(EDGE_TREE_t **, EDGE_t *);
void delete_edge(EDGE_TREE_t **, EDGE_t *);
EDGE_TREE_t *minimal_edge(EDGE_TREE_t *);
EDGE_TREE_t *find_edge(EDGE_TREE_t *, int);

void ele_info_free(ELE_INFO_t **);
void ele_cleanup(ELEMENT_t **);
void BD_free(BD_t **);
void CP_free(CP_t **);
void ele_data_free(ELE_DATA_t **);
void img_data_free (IMG_DATA_t **);
void img_tree_free(IMG_TREE_t **, ELE_INFO_t *);
void img_tree_cleanup(IMG_TREE_t **, ELE_INFO_t *);
void msp_data_free(MSP_DATA_t **);
void edge_tree_cleanup(ELE_INFO_t *, EDGE_TREE_t **);
void edge_tree_free(EDGE_TREE_t **);
void fam_data_free(FAM_DATA_t **);
void fam_data_cleanup(FAM_DATA_t **);
void fam_cleanup(FAMILY_t **);
void frag_data_free(FRAG_DATA_t **);
void frag_data_cleanup(FRAG_DATA_t **);

/* Existing toString-style print helpers */
void print_ele(ELEMENT_t *rt);
void print_ele_info(ELE_INFO_t *rt);
void print_ele_data(ELE_DATA_t *rt);
void print_edge(EDGE_t *rt);
void print_edge_tree(EDGE_TREE_t *rt, int level);
void print_local_network(ELE_DATA_t *rt);
void print_edge_tree_GML(EDGE_TREE_t *rt, int rel_to_ele_id);
void print_all_eles_GML();

/* New toString-style print helpers (definitions at end of this file) */
void print_cp(CP_t *cp);
void print_cp_list(CP_t *cp);
void print_bd(BD_t *bd);
void print_bd_list(BD_t *bd);
void print_family(FAMILY_t *fam);
void print_fam_data(FAM_DATA_t *fd);

/* ============================================================
 * Per-program logging state
 *
 * Each .c that includes ele.h must define storage for these two
 * variables.  See recon_log.h for the RLOG_* macros that use them.
 *
 * Example (at file scope in each .c):
 *   int   recon_log_level = RECON_LOG_INFO;
 *   FILE *recon_log_fp    = NULL;
 * ============================================================ */
extern int   recon_log_level;
extern FILE *recon_log_fp;


/* ============================================================
 * Pipeline-wide global state
 *
 * These globals are defined in ele.c and declared extern here so that
 * all translation units that include this header share the same
 * instance via the linker.  Each program (eleredef, edgeredef, famdef)
 * gets its own copy at run time (separate processes).
 * ============================================================ */

/* all_ele        -- pre-allocated array of ELE_INFO_t pointers, indexed
 *                   0..ele_array_size-1.  New elements are added here until
 *                   the array is full, then via the ele_info_data overflow list. */
extern ELE_INFO_t **all_ele;

/* ele_ct         -- running count of elements seen so far in this stage.
 * ele_array_size -- initial allocation size of all_ele[]. */
extern int ele_ct, ele_array_size, fam_ct;

/* clan_size      -- total elements in the current BFS local network.
 * clan_core_size -- elements within DEPTH hops of the BFS seed. */
extern int clan_size, clan_core_size;

/* MSP memory-tracking counters (used by MSP_malloc / MSP_free) */
extern int32_t msp_in_mem;   /* currently allocated MSPs */
extern int32_t msp_left;     /* MSPs remaining after a dissolve (should be 0) */
extern int32_t msp_ct;       /* total MSPs allocated in this run */
extern int32_t msp_index;    /* highest MSP sequential index seen */

/* Edge memory-tracking counters (used by EDGE_malloc / EDGE_free) */
extern int32_t edge_index, edge_in_mem, edge_left, edge_ct;

/* Miscellaneous run-time counters */
extern int32_t files_read;   /* number of element files read from disk */
extern int32_t clan_ct;      /* number of BFS local networks processed */
extern int32_t err_no;       /* accumulated error count; non-zero triggers exit */

/* ele_info_data  -- head of the overflow linked list used when ele_ct
 *                   exceeds ele_array_size.
 * ele_info_tail  -- tail pointer for O(1) append to the overflow list. */
extern ELE_INFO_t *ele_info_data, *ele_info_tail;

/* Pipeline output files -- opened by each program's main() */
extern FILE *new_msps;    /* new MSP records created during dissection */
extern FILE *eles;        /* final elements output (famdef summary/eles) */
extern FILE *unproc;      /* elements removed as large tandems or dismissed */
extern FILE *combo;       /* combination output (combo elements) */
extern FILE *obs;         /* obsolete elements */
extern FILE *fams;        /* final families output (famdef summary/families) */
extern FILE *log_file;    /* primary pipeline progress log */


/* ============================================================
 * Flat-file element database (implemented in ele_db.c)
 * ============================================================ */
#include "ele_db.h"

/* All function implementations are in ele.c */

#endif /* ELE_H */
