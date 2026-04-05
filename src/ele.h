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
 * These globals are shared across all functions in ele.h and the
 * three .c files that include it (eleredef, edgeredef, famdef).
 * They are defined by whichever .c includes this header.
 * ============================================================ */

/* all_ele        -- pre-allocated array of ELE_INFO_t pointers, indexed
 *                   0..ele_array_size-1.  New elements are added here until
 *                   the array is full, then via the ele_info_data overflow list. */
ELE_INFO_t **all_ele;

/* ele_ct         -- running count of elements seen so far in this stage.
 * ele_array_size -- initial allocation size of all_ele[]. */
int ele_ct, ele_array_size, fam_ct;

/* clan_size      -- total elements in the current BFS local network.
 * clan_core_size -- elements within DEPTH hops of the BFS seed. */
int clan_size, clan_core_size;

/* MSP memory-tracking counters (used by MSP_malloc / MSP_free) */
int32_t msp_in_mem;   /* currently allocated MSPs */
int32_t msp_left;     /* MSPs remaining after a dissolve (should be 0) */
int32_t msp_ct;       /* total MSPs allocated in this run */
int32_t msp_index;    /* highest MSP sequential index seen */

/* Edge memory-tracking counters (used by EDGE_malloc / EDGE_free) */
int32_t edge_index, edge_in_mem, edge_left, edge_ct;

/* Miscellaneous run-time counters */
int32_t files_read;   /* number of element files read from disk */
int32_t clan_ct;      /* number of BFS local networks processed */
int32_t err_no;       /* accumulated error count; non-zero triggers exit */

/* ele_info_data  -- head of the overflow linked list used when ele_ct
 *                   exceeds ele_array_size.
 * ele_info_tail  -- tail pointer for O(1) append to the overflow list. */
ELE_INFO_t *ele_info_data, *ele_info_tail;

FAM_DATA_t *FAMs;

/* Pipeline output files -- opened by each program's main() */
FILE *err;         /* error / diagnostic messages */
FILE *new_msps;    /* new MSP records created during dissection */
FILE *eles;        /* final elements output (famdef summary/eles) */
FILE *unproc;      /* elements removed as large tandems or dismissed */
FILE *combo;       /* combination output (combo elements) */
FILE *obs;         /* obsolete elements */
FILE *fams;        /* final families output (famdef summary/families) */
FILE *log_file;    /* primary pipeline progress log */


/***************
 * DEBUG       *
 ***************/

void print_ele_data(ELE_DATA_t *rt) {
  if (!rt) return;
  print_ele_info( rt->ele_info );
  if (rt->next)
    print_ele_data(rt->next);
}

void print_ele_info(ELE_INFO_t *rt) {
  if (!rt) return;
  printf("ELE_INFO_t: index=%d, stat=%c, file_update=%d,", rt->index, rt->stat,
         rt->file_updated );
  if ( rt->to_family )
    printf(" to_family->index=%d,", rt->to_family->index );
  else
    printf(" to_family=Null,");
  if ( rt->ele )
    printf(" ele->index=%d", rt->ele->index );
  else
    printf(" ele=Null");
  printf("\n");
}

void print_ele(ELEMENT_t *rt) {
  if (!rt) return;
  printf("ELEMENT_t: index=%d, direction=%d, update=%d, l_hold=%d, img_no=%d, flimg_no=%d, edge_no=%d,",
         rt->index, rt->direction, rt->update, rt->l_hold, rt->img_no, rt->flimg_no, rt->edge_no );
  if ( rt->to_img_tree )
    printf(" to_img_tree->to_image->index=%d,", rt->to_img_tree->to_image->index );
  else
    printf(" to_img_tree=Null,");
  if ( rt->to_img_data )
    printf(" to_img_data->to_image->index=%d,", rt->to_img_data->to_image->index );
  else
    printf(" to_img_data=Null,");
  if ( rt->edges )
    printf(" edges->to_edge->index=%d,", rt->edges->to_edge->index );
  else
    printf(" edges=Null,");
  if ( rt->PCP )
    printf(" PCP->cp=%d,", rt->PCP->cp );
  else
    printf(" PCP=Null,");
  if ( rt->TBD )
    printf(" TBD->bd=%d,", rt->TBD->bd );
  else
    printf(" TBD=Null,");
  if ( rt->redef )
    printf(" redef->ele_info->index=%d", rt->redef->ele_info->index );
  else
    printf(" redef=Null");
  printf("\n");
}

void print_edge(EDGE_t *rt) {
  if (!rt) return;
  printf("EDGE_t: index=%d, type=%c, direction=%d, score=%d,", rt->index, rt->type,
         rt->direction, rt->score );
  if ( rt->ele1_info )
    printf(" ele1_info->index=%d,", rt->ele1_info->index );
  else
    printf(" ele1_info=Null,");
  if ( rt->ele2_info )
    printf(" ele2_info->index=%d,", rt->ele2_info->index );
  else
    printf(" ele2_info=Null,");
  printf("\n");
}

//
//  Given an edge tree do a dfs and print out edges in GML format
//
//  If rel_to_ele_id is < 1 then print all edges defined in the
//  the edge tree.  Otherwise only print the edges in which the
//  lowest element index is the same as rel_to_ele_id.  This is
//  used by other functions to print out a non-redundant set
//  of edges.  Although this doesn't work if there are duplicates
//  in the set.  I have seed a single node in which the following
//  two edges were in the tree: 6 -> 8, and 8 -> 6.
//
//  RMH
void print_edge_tree_GML(EDGE_TREE_t *rt, int rel_to_ele_id){
  if (!rt) return;
  if ( rel_to_ele_id < 1 ||
       ((rt->to_edge->ele1_info->index < rt->to_edge->ele2_info->index &&
        rt->to_edge->ele1_info->index == rel_to_ele_id ) ||
       (rt->to_edge->ele2_info->index < rt->to_edge->ele1_info->index &&
        rt->to_edge->ele2_info->index == rel_to_ele_id )) )
    {
    //printf("rel_to_ele_id=%d\n", rel_to_ele_id);
    printf("edge [\n");
    printf("       source %d\n", rt->to_edge->ele1_info->index );
    printf("       target %d\n", rt->to_edge->ele2_info->index );
    printf("       label \"index=%d, type=%c, direction=%d, score=%d\"\n", rt->to_edge->index, rt->to_edge->type,
         rt->to_edge->direction, rt->to_edge->score );
    printf("     ]\n");
  }
  if (rt->l)
    print_edge_tree_GML(rt->l, rel_to_ele_id);
  if (rt->r)
    print_edge_tree_GML(rt->r, rel_to_ele_id);
}


//
// Here because the all_eles array and ele_array_size
// datastructures are global and not passed into the
// functions. sigh...
//
void print_all_eles_GML() {
  int i = 0;
  printf("graph [\n");
  for (i=0; i<ele_ct && i<ele_array_size; i++) {
    printf("node [\n");
    printf("       id %d\n", (*(all_ele+i))->index);
    printf("       label \"e%d\"\n", (*(all_ele+i))->index);
    printf("     ]\n");
  }
  for (i=0; i<ele_ct && i<ele_array_size; i++)
    print_edge_tree_GML((*(all_ele+i))->ele->edges, (*(all_ele+i))->ele->index);
  printf("]\n");
}


void print_local_network(ELE_DATA_t *rt) {
  if (!rt) return;
  ELE_DATA_t *curr;
  curr=rt;
  while ( curr ) {
    //printf("loop over %d\n", curr->ele_info->ele->index);
    print_edge_tree_GML(curr->ele_info->ele->edges, curr->ele_info->ele->index);
    curr = curr->next;
  }
}

void print_edge_tree(EDGE_TREE_t *rt, int level) {
  int i;
  if (!rt) return;
  for ( i = 0; i < level; i++ )
    printf("  ");
  print_edge(rt->to_edge);
  if (rt->l)
    print_edge_tree(rt->l, level+1);
  if (rt->r)
    print_edge_tree(rt->r, level+1);
}






/***************
 * Initializer *
 ***************/



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
  /*MSP_DATA_t *new_hanger = (MSP_DATA_t *) malloc(sizeof(MSP_DATA_t *));*/

  msp_in_mem ++;
  msp_ct ++;
  /*new_hanger->to_msp = msp_tmp;
  new_hanger->next = all_msps;
  all_msps = new_hanger;

  msp_tmp->hanger = new_hanger;*/

  return msp_tmp;
}




void MSP_free(MSP_t *m) {
  msp_in_mem --;
  /*m->hanger->to_msp = NULL;*/
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




/***************
 * Small Tools *
 ***************/




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
 *
 * Local variables
 *   ei    element index
 *   ino   image count for this element
 *   p_ct  unique partner element count
 *   i_ct  total MSP (msp line) count scanned in element file
 *   p_id  most recently seen partner id (for deduplication after sort)
 */
void outthrow_big_tandems(FILE *size_list){
  int ei, ino;
  FILE *ele_file;
  char ele_name[50], line[150], head[10], *msp = "msp";
  int id, ele1, ele2;
  int *partners, i_ct, i, p_id, p_ct;
  int ratio;

  while (fgets(line, 25, size_list)) {
    sscanf(line, "%d %d\n", &ei, &ino);
    if (ino > SIZE_LIMIT) {
      sprintf(ele_name, "tmp/e%d", ei);
      ele_file = fopen(ele_name, "r");
      if (!ele_file) {
	fprintf(log_file, "Can not open ele file %s.  Exit.\n", ele_name);
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
    }
  }
}



int int_cmp(const void *i1, const void *i2) {
  return *((int *) i1) - *((int *) i2);
}



void spit_out_ele(ELE_INFO_t *ele_info) {
  /*  char *command = (char *) malloc(100*sizeof(char));*/

  ele_info->file_updated = 1;
  /* the reason for using ln instead of mv is to avoid changing files in tmp/ */
#if 0
  sprintf(command, "ln -s tmp/e%d unproc/.\n", ele_info->index);
  if (system(command)) {
    fprintf(log_file, "error linking tmp/e%d to unproc/\n", ele_info->index);
    fflush(log_file);
    exit(6);
  }
#endif
  fprintf(unproc, "%d\n", ele_info->index);
  /*  free(command);*/
}





/**********************
 * Read in an element *
 **********************/




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
 *
 * Image loading builds an in-memory balanced BST (to_img_tree) keyed on
 * image index.  Edge loading builds a parallel BST (edges) keyed on edge
 * index.  Both are sorted before building via qsort + build_*_tree().
 *
 * Sequence names in msp/frag lines are interned via GetSeqIndex() so
 * that pointer-equality tests work throughout the element processing code.
 *
 * The keyword variables (index, stat, frag, msp, edge, ...) are string
 * literals used as targets for strncmp() line-dispatch.
 */
ELEMENT_t *ele_read_in(ELE_INFO_t *ele_info, int stage) {
  char line[200], head[10], rest[150], *fn = (char *) malloc(64*sizeof(char));
  char fragname[SEQ_NAME_MAX_LEN];
  int pos;
  FILE *fp;
  IMAGE_t **img_array;
  EDGE_t **edge_array;
  int holder, img_encountered=0, edge_encountered=0, ct;
  CP_t *new_pcp;
  BD_t *new_tbd;
  ELE_DATA_t *cur_ele_data;

  char *index="index", *stat = "stat", *file_updated="file_updated", *family="family", *direc="direc", *update="update", *l_hold="l_hold", *img_no="img_no", *flimg_no="flimg_no", *edge_no="edge_no", *frag="frag", *pcp="pcp", *tbd="tbd", *redef="redef", *msp="msp", *edge="edge";

  if (ele_info->file_updated) snprintf(fn, 64, "tmp2/e%d", ele_info->index);
  else snprintf(fn, 64, "tmp/e%d", ele_info->index);
  fp = fopen(fn, "r");
  /*if (!fp) {
    snprintf(fn, 64, "tmp/e%d", ele_info->index);
    fp = fopen(fn, "r");
  }*/
  files_read ++;

  ele_info->ele = ele_init(ele_info->index);

  while (fgets(line, 200, fp)) {
    sscanf(line, "%s %*s", head);
/*    if (!strncmp(head, index, 10)) sscanf(rest, "%d\n", ele_info->index);
    if (!strncmp(head, stat, 10)) sscanf(rest, "%c\n", ele_info->stat);
    if (!strncmp(head, file_updated, 10)) sscanf(rest, "%d\n", ele_info->file_updated); */

    /* if (!strncmp(head, l_hold, 10)) {sscanf(line, "%*s %d\n", &ele_info->ele->l_hold); continue;} */
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
  free(fn);
  if (stage != 3 && img_encountered != ele_info->ele->img_no) {
    if (ele_info->file_updated) { /* we spit out elements with too many images, so some of their partners might have less images left than img_no */
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

  /*ele_info->ele->index = ele_info->index;*/
  if (edge_encountered) {
    qsort(edge_array, edge_encountered, sizeof(EDGE_t *), edge_index_cmp);
    build_edge_tree(&ele_info->ele->edges, edge_array, 0, edge_encountered-1);
    /*ct = count_edge_nodes(ele_info->ele->edges);
    if (ct != edge_encountered) {
	err_no ++;
	fprintf(log_file, "error:  trouble building edge tree: %d %d\n", edge_encountered, ct);
	fflush(log_file);
	exit(4);
    }*/
    free(edge_array);
  }
  if (stage != 3 && ele_info->ele->img_no) {
    qsort(img_array, img_encountered, sizeof(IMAGE_t *), img_index_cmp);
    build_img_tree(&ele_info->ele->to_img_tree, img_array, 0, img_encountered-1);
    /*ct = count_img_nodes(ele_info->ele->to_img_tree);
    if (ct != img_encountered) {
	err_no ++;
	fprintf(log_file, "error:  trouble building the image tree: %d %d\n", img_encountered, ct);
	fflush(log_file);
	exit(2);
    }*/
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


  /*printf("%s", line);*/
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
	/*msp_in_mem ++;*/
	/*msp_ct ++;*/
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
      /*msp_in_mem ++;*/
      /*msp_ct ++;*/
    }
  }
}





void edge_scan(char *line, EDGE_t **edge_array, ELE_INFO_t *ele_info, int *edge_encountered_ptr){
  EDGE_t *new = EDGE_malloc();
  int ele1, ele2;
  ELE_INFO_t *epi;

  sscanf(line, "edge %d %c %d %d %d %d\n", &new->index, &new->type, &new->direction, &new->score, &ele1, &ele2); 

  if (ele1 != ele2 && new->type != 'c') { /* bullet proof */
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
  }else EDGE_free(new);
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





/************************
 * Write out an element *
 ************************/





void ele_write_out(ELE_INFO_t *ele_info, int stage) {
  char *fn = (char *) malloc(50*sizeof(char));
  FILE *fp;
  CP_t *cur_pcp;
  BD_t *cur_bd;
  ELE_DATA_t *cur_ele_data;
  int ct;

  /*  sprintf(fn, "tmp2/clan/e%d", ele_info->index);*/
  sprintf(fn, "tmp2/e%d", ele_info->index);
  fp = fopen(fn, "w");

/*  fprintf(fp, "index %d\n", ele_info->index);
  fprintf(fp, "stat %s\n", ele_info->stat);
  fprintf(fp, "file_updated %d\n", ele_info->file_updated);
  if (ele_info->to_family) fprintf(fp, "family %d\n", ele_info->to_family->index); */

  fprintf(fp, "frag %s %d %d\n", ele_info->ele->frag.seq_name, ele_info->ele->frag.lb, ele_info->ele->frag.rb);
  /*  if (ele_info->ele->l_hold) fprintf(fp, "l_hold %d\n", ele_info->ele->l_hold); */
  if (ele_info->ele->direction != 1) fprintf(fp, "direc %d\n", ele_info->ele->direction);
  if (ele_info->ele->update) fprintf(fp, "update %d\n", ele_info->ele->update);
  /*ct = count_edge_nodes(ele_info->ele->edges);
  if (ct > ele_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  more edges than known: ele %d, %d, %d\n", ele_info->index, ele_info->ele->edge_no, ct);
    fflush(log_file);
    exit(4);
  } else if (ct < ele_info->ele->edge_no) {
    err_no ++;
    fprintf(log_file, "error:  less edges than known: ele %d, %d, %d\n", ele_info->index, ele_info->ele->edge_no, ct);
    fflush(log_file);
    exit(4);
  }*/
  /*ct = count_img_nodes(ele_info->ele->to_img_tree);
  if (ct > ele_info->ele->img_no) {
    err_no ++;
    fprintf(log_file, "error:  more images than known: ele %d, %d, %d\n", ele_info->index, ele_info->ele->img_no, ct);
    fflush(log_file);
    exit(2);
  }  else if (ct < ele_info->ele->img_no) {
    err_no ++;
    fprintf(log_file, "error:  less images than known: ele %d, %d, %d\n", ele_info->index, ele_info->ele->img_no, ct);
    fflush(log_file);
    exit(2);
  }*/
  fprintf(fp, "img_no %d\n", ele_info->ele->img_no);
  if (ele_info->ele->flimg_no) fprintf(fp, "flimg_no %d\n", ele_info->ele->flimg_no);
  fprintf(fp, "edge_no %d\n", ele_info->ele->edge_no);
  if (ele_info->ele->img_no) write_out_msps(fp, ele_info->ele->to_img_tree);
  if (ele_info->ele->edge_no) write_out_edges(fp, ele_info->ele->edges, ele_info);
  if (stage < 2 && ele_info->ele->PCP) {
    /*fprintf(fp, "PCPs \n");*/
    cur_pcp = ele_info->ele->PCP;
    while (cur_pcp) {
      fprintf(fp, "pcp %d %d\n", cur_pcp->cp, cur_pcp->contributor->index);
      cur_pcp = cur_pcp->next;
    }
  }
  if (ele_info->ele->TBD) {
    /*fprintf(fp, "TBDs \n");*/
    cur_bd = ele_info->ele->TBD;
    while (cur_bd) {
      fprintf(fp, "tbd %d %d\n", cur_bd->bd, cur_bd->support);
      cur_bd = cur_bd->next;
    }
  }
  if (ele_info->ele->redef) {
    /*fprintf(fp, "redefs \n");*/
    cur_ele_data = ele_info->ele->redef;
    while (cur_ele_data) {
      fprintf(fp, "redef %d\n", cur_ele_data->ele_info->index);
      cur_ele_data = cur_ele_data->next;
    }
  }
  fclose(fp);
  free(fn);
}




void write_out_msps(FILE *fp, IMG_TREE_t *rt) {
  if (rt->l) write_out_msps(fp, rt->l);

  if (rt->to_image->to_msp->query.ele_info->index != rt->to_image->to_msp->sbjct.ele_info->index || rt->to_image == &rt->to_image->to_msp->sbjct) { /* trick to report self MSPs only once */
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






/***************************
 * manipulating image tree *
 ***************************/





void build_img_tree(IMG_TREE_t **rt_ptr, IMAGE_t **img_array, int ori, int end) {
  int mid = (ori+end)/2;
  IMG_TREE_t *new = (IMG_TREE_t *) malloc(sizeof(IMG_TREE_t));
  new->to_image = *(img_array+mid);
  new->l = NULL;
  new->r = NULL;
  new->p = NULL;
  *rt_ptr = new;
  if (mid > ori) {
    build_img_tree(&new->l, img_array, ori, mid-1);
    new->l->p = new;
  }
  if (mid < end) {
    build_img_tree(&new->r, img_array, mid+1, end);
    new->r->p = new;
  }
}



int count_img_nodes(IMG_TREE_t *rt) {
  int ct=0;

  if (!rt) return 0;
  if (rt->l) ct += count_img_nodes(rt->l);
  if (rt->r) ct += count_img_nodes(rt->r);

  ct ++;
  return ct;
}




void insert_image(IMG_TREE_t **rt_ptr, IMAGE_t *i) {
  IMG_TREE_t *new, *x, *y;

  new = (IMG_TREE_t *) malloc(sizeof(IMG_TREE_t));
  new->to_image = i;
  new->p = NULL;
  new->l = NULL;
  new->r = NULL;

  y = NULL;
  x = *rt_ptr;
  while (x) {
    y = x;
    if (new->to_image->index < x->to_image->index) x = x->l;
    else x = x->r;
  }
  new->p = y;
  if (!y) *rt_ptr = new;
  else {
    if (new->to_image->index < y->to_image->index) y->l = new;
    else y->r = new;
  }
}




void delete_image(IMG_TREE_t **rt_ptr, IMAGE_t *i) {
  IMG_TREE_t *x, *y, *z;

  if (!(*rt_ptr)) {
    err_no ++;
    fprintf(log_file, "error:  Image tree does not exist: %d\n", i->ele_info->index);
    fflush(log_file);
    exit(2);
    return;
  }
  z= find_image(*rt_ptr, i->index);
  if (z) {
    if (!z->l || !z->r) y = z;
    else y = minimal_image(z->r);
    if (y->l) x = y->l;
    else x = y->r;
    if (x) x->p = y->p;
    if (!y->p) *rt_ptr = x;
    else {
      if (y == y->p->l) y->p->l = x;
      else y->p->r = x;
    }
    if (y != z) z->to_image = y->to_image;
    free(y);
  } else {
    err_no ++;
    fprintf(log_file, "error:  Can not find the image in the tree: %d %d\n", i->ele_info->index, i->index);
    fflush(log_file);
    exit(2);
    return;
  }
}



IMG_TREE_t *minimal_image(IMG_TREE_t *rt) {
  while (rt->l) {
    rt = rt->l;
  }

  return rt;
}





IMG_TREE_t *find_image(IMG_TREE_t *rt, int index) {
  while (rt && index != rt->to_image->index) {
    if (index < rt->to_image->index) rt = rt->l;
    else rt = rt->r;
  }

  return rt;
}




IMG_DATA_t **listify(IMG_TREE_t *rt, IMG_DATA_t **tail_ptr) {
  IMG_DATA_t *new, **new_tail_ptr;

  if (!rt) return NULL;
  if (rt->l) tail_ptr = listify(rt->l, tail_ptr);
  new = (IMG_DATA_t *) malloc(sizeof(IMG_DATA_t));
  new->to_image = rt->to_image;
  new->next = NULL;
  *tail_ptr = new;
  tail_ptr = &new->next;
  if (rt->r) tail_ptr = listify(rt->r, tail_ptr);
  return tail_ptr;
}





/**************************
 * Manipulating edge tree *
 **************************/




void build_edge_tree(EDGE_TREE_t **rt_ptr, EDGE_t **edge_array, int ori, int end) {
  int mid = (ori+end)/2;
  EDGE_TREE_t *new = (EDGE_TREE_t *) malloc(sizeof(EDGE_TREE_t));
  new->to_edge = *(edge_array+mid);
  new->l = NULL;
  new->r = NULL;
  new->p = NULL;
  *rt_ptr = new;
  if (mid > ori) {
    build_edge_tree(&new->l, edge_array, ori, mid-1);
    new->l->p = new;
  }
  if (mid < end) {
    build_edge_tree(&new->r, edge_array, mid+1, end);
    new->r->p = new;
  }
}




int count_edge_nodes(EDGE_TREE_t *rt) {
  int ct=0;

  if (!rt) return 0;
  if (rt->l) ct += count_edge_nodes(rt->l);
  if (rt->r) ct += count_edge_nodes(rt->r);

  if (rt->to_edge->type != 'c' && rt->to_edge->ele1_info->index != rt->to_edge->ele2_info->index) ct ++;
  return ct;
}




int count_total_edges(EDGE_TREE_t *rt) {
  int ct=0;

  if (!rt) return 0;
  if (rt->l) ct += count_total_edges(rt->l);
  if (rt->r) ct += count_total_edges(rt->r);

  ct ++;
  return ct;
}





void insert_edge(EDGE_TREE_t **rt_ptr, EDGE_t *ed) {
  EDGE_TREE_t *new, *x, *y;

  new = (EDGE_TREE_t *) malloc(sizeof(EDGE_TREE_t));
  new->to_edge = ed;
  new->p = NULL;
  new->l = NULL;
  new->r = NULL;

  y = NULL;
  x = *rt_ptr;
  while (x) {
    y = x;
    if (new->to_edge->index < x->to_edge->index) x = x->l;
    else x = x->r;
  }
  new->p = y;
  if (!y) *rt_ptr = new;
  else {
    if (new->to_edge->index < y->to_edge->index) y->l = new;
    else y->r = new;
  }
}




void delete_edge(EDGE_TREE_t **rt_ptr, EDGE_t *ed) {
  EDGE_TREE_t *x, *y, *z;

  if (!(*rt_ptr)) {
    err_no ++;
    fprintf(log_file, "error:  edge tree empty: %d %d %d\n",  ed->index, ed->ele1_info->index, ed->ele2_info->index); 
    fflush(log_file);
    exit(4);
  }
  z= find_edge(*rt_ptr, ed->index);
  if (z) {
    if (!z->l || !z->r) y = z;
    else y = minimal_edge(z->r);
    if (y->l) x = y->l;
    else x = y->r;
    if (x) x->p = y->p;
    if (!y->p) *rt_ptr = x;
    else {
      if (y == y->p->l) y->p->l = x;
      else y->p->r = x;
    }
    if (y != z) z->to_edge = y->to_edge;
    free(y);
  } else {
    err_no ++;
    fprintf(log_file, "error:  Can not find edge %d between ele %d and %d\n", ed->index, ed->ele1_info->index, ed->ele2_info->index);
    fflush(log_file);
    exit(4);
  }
}




EDGE_TREE_t *minimal_edge(EDGE_TREE_t *rt) {
  while (rt->l) {
    rt = rt->l;
  }

  return rt;
}




EDGE_TREE_t *find_edge(EDGE_TREE_t *rt, int index) {
  while (rt && index != rt->to_edge->index) {
    if (index < rt->to_edge->index) rt = rt->l;
    else rt = rt->r;
  }

  return rt;
}








/************************
 ************************
 ***                  ***
 ***   CLEANING UPS   ***
 ***                  ***
 ************************
 ************************/




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



void ele_data_free (ELE_DATA_t **e) {
  ELE_DATA_t *cur, *next;

  cur = *e;
  while (cur != NULL) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *e = NULL;
}



/* void ele_data_cleanup (ELE_DATA_t **e) {
  ELE_DATA_t *cur, *next, *dangler;

  cur = *e;
  while (cur != NULL) {
    next = cur->next;
    if (cur->ele_info->ele->redef) {
      ele_data_cleanup(&cur->ele_info->ele->redef);
    }
    ele_cleanup(&cur->ele_info->ele);
    free(cur);
    cur = next;
  }
  *e = NULL;
} */



void ele_cleanup(ELEMENT_t **e) {
  ELE_INFO_t *ele_info;
  int ct;

  if (!(*e)) return;

  if ((*e)->frag.lb > (*e)->frag.rb) {
    fprintf(log_file, "error:  ele %d reversed before clean up\n", (*e)->index);
    fflush(log_file);
    exit(3);
  }
  /* reason for finding ele_info here instead of passing it over when
  calling the function is that this way, ele_cleanup() will be more 
  independent of the ele_info structure */
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
    /*ct = count_img_nodes((*e)->to_img_tree);
    if (ct != (*e)->img_no) {
      err_no ++;
      fprintf(log_file, "error:  image keeping not synchronized: ele %d %d %d\n", (*e)->index, (*e)->img_no, ct);
      fflush(log_file);
      exit(2);
    }*/
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


void BD_free (BD_t **bds) {
  BD_t *cur, *next;

  cur = *bds;
  while (cur != NULL) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *bds = NULL;
}



void CP_free (CP_t **cps) {
  CP_t *cur, *next;

  cur = *cps;
  while (cur != NULL) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *cps = NULL;
}



void img_data_free (IMG_DATA_t **i) {
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
    /* this is a hack here.  we cann't simply free self msps 'cuz 
       that way when we run across the second image in the tree, the
       to_image pointer will be pointing to space already freed.  So
       we mark a self msp by changing its stat when we first see it, 
       and when we see it the second time, which means when a self 
       msp is alredy marked, we free it.
       
       Same thing for self edges.  only we change the type of the edge
       to 'm' */
    if ((*rt_p)->to_image->to_msp->stat != 'm') (*rt_p)->to_image->to_msp->stat = 'm';
    else {
      MSP_free((*rt_p)->to_image->to_msp);
      /*msp_in_mem --;*/
    }
  } else if (!ep->ele) {
    MSP_free((*rt_p)->to_image->to_msp);
    /*msp_in_mem --;*/
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




/* void img_data_cleanup (IMG_DATA_t **i) {
  IMG_DATA_t *cur, *next;

  cur = *i;
  while (cur != NULL) {
    next = cur->next;
    free(cur->to_image);
    free(cur);
    cur = next;
  }
  *i = NULL;
} */




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
    /*if (cur->to_msp) MSP_free(cur->to_msp);*/
    free(cur);
    cur = next;
  }
  *md = NULL;
}




/* void edge_data_free(EDGE_DATA_t **ed) {
  EDGE_DATA_t *cur, *next;

  cur = *ed;
  while (cur != NULL) {
    next = cur->next;
    free(cur);
    cur = next;
  }
  *ed = NULL;
}




void edge_data_cleanup(EDGE_DATA_t **ed) {
  EDGE_DATA_t *cur, *next;

  cur = *ed;
  while (cur != NULL) {
    next = cur->next;
    free(cur->to_edge);
    free(cur);
    cur = next;
  }
  *ed = NULL;
} */




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


/* ============================================================
 * Additional toString-style print helpers
 * ============================================================ */

/*
 * print_cp  --  one-line summary of a single CP_t node.
 */
void print_cp(CP_t *cp) {
  if (!cp) { printf("CP_t: NULL\n"); return; }
  printf("CP_t: cp=%d, contributor->index=%d (stat=%c)\n",
         cp->cp,
         cp->contributor ? cp->contributor->index : -1,
         cp->contributor ? cp->contributor->stat  : '?');
}

/*
 * print_cp_list  --  walk a CP_t linked list and print each node.
 */
void print_cp_list(CP_t *cp) {
  int count = 0;
  while (cp) {
    printf("  [%d] ", count++);
    print_cp(cp);
    cp = cp->next;
  }
}

/*
 * print_bd  --  one-line summary of a single BD_t node.
 */
void print_bd(BD_t *bd) {
  if (!bd) { printf("BD_t: NULL\n"); return; }
  printf("BD_t: bd=%d, support=%d\n", bd->bd, bd->support);
}

/*
 * print_bd_list  --  walk a BD_t linked list and print each node.
 */
void print_bd_list(BD_t *bd) {
  int count = 0;
  while (bd) {
    printf("  [%d] ", count++);
    print_bd(bd);
    bd = bd->next;
  }
}

/*
 * print_family  --  one-line summary of a FAMILY_t.
 */
void print_family(FAMILY_t *fam) {
  ELE_DATA_t *m;
  int mem_count = 0;
  if (!fam) { printf("FAMILY_t: NULL\n"); return; }
  m = fam->members;
  while (m) { mem_count++; m = m->next; }
  printf("FAMILY_t: index=%d, name='%s', member_count=%d\n",
         fam->index, fam->name, mem_count);
}

/*
 * print_fam_data  --  walk a FAM_DATA_t linked list and print each family.
 */
void print_fam_data(FAM_DATA_t *fd) {
  int count = 0;
  while (fd) {
    printf("  [%d] ", count++);
    print_family(fd->to_family);
    fd = fd->next;
  }
}
