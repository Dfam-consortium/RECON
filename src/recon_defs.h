/*
 * recon_defs.h  --  RECON pipeline-wide named constants
 *
 * All hard-coded thresholds from the original source are collected
 * here.  Each constant is annotated with its biological/algorithmic
 * meaning, the file it originated in, and -- where applicable --
 * whether it can be overridden at compile-time via -D flags in the
 * Makefile.
 *
 * Backward-compatibility aliases (matching the old macro names) are
 * provided below each new definition so that existing code continues
 * to compile unchanged.  Once every use-site has been migrated to the
 * new name the alias can be removed.
 *
 * Author: Robert Hubley, Institute for Systems Biology
 */

#ifndef RECON_DEFS_H
#define RECON_DEFS_H

/*
 * RECON_VERSION  --  release version string, injected by the Makefile via
 * -DRECON_VERSION=\"$(VERSION)\".  Falls back to "unknown" if built outside
 * the project Makefile.
 */
#ifndef RECON_VERSION
#define RECON_VERSION "unknown"
#endif


/* ============================================================
 * Sequence-list / global limits   (originally in bolts.h)
 * ============================================================ */

/*
 * SEQ_NAME_MAX_LEN  --  maximum characters (including the NUL
 * terminator) that a sequence identifier may occupy.  Sequence
 * names are read from the seq_list file into the seq_name_table[]
 * array.  Increasing this value requires a full rebuild but no
 * algorithmic change.
 *
 * Compile-time override:  make EXTRA_CFLAGS="-DSEQ_NAME_MAX_LEN=100"
 *
 * Originally: NAME_LEN 50 in bolts.h
 */
#ifndef SEQ_NAME_MAX_LEN
#define SEQ_NAME_MAX_LEN  50
#endif
/* Backward-compat alias -- remove after all uses are migrated */
#define NAME_LEN  SEQ_NAME_MAX_LEN


/*
 * BFS_CLAN_DEPTH  --  maximum depth of the breadth-first search used
 * to build the "clan" neighbourhood around a seed element during
 * element redefinition (eleredef) and edge redefinition (edgeredef).
 * A depth of 3 means the seed plus up to three rings of neighbours
 * are loaded into the local network before the redefinition pass
 * begins.
 *
 * Compile-time override:  make EXTRA_CFLAGS="-DBFS_CLAN_DEPTH=4"
 *
 * Originally: DEPTH 3 in bolts.h
 */
#ifndef BFS_CLAN_DEPTH
#define BFS_CLAN_DEPTH  3
#endif
/* Backward-compat alias */
#define DEPTH  BFS_CLAN_DEPTH


/* ============================================================
 * Element / tandem-repeat filter   (originally in ele.h)
 * ============================================================ */

/*
 * TANDEM_SIZE_LIMIT  --  elements with more images than this value
 * are subjected to the tandem-repeat filter in
 * outthrow_big_tandems().  Elements that exceed the filter are
 * labelled 'O' (orphan) and removed from further processing.
 *
 * Compile-time override:  make EXTRA_CFLAGS="-DTANDEM_SIZE_LIMIT=1000"
 *
 * Originally: SIZE_LIMIT 500 in ele.h
 */
#ifndef TANDEM_SIZE_LIMIT
#define TANDEM_SIZE_LIMIT  500
#endif
/* Backward-compat alias */
#define SIZE_LIMIT  TANDEM_SIZE_LIMIT


/*
 * TANDEM_IMG_RATIO  --  if the ratio (image_count / unique_partner_count)
 * exceeds this value an element is classified as a tandem repeat and
 * filtered out.  A high ratio indicates that the same small set of
 * partner elements accounts for a disproportionate share of the images,
 * which is characteristic of tandem repeats rather than dispersed
 * repeats.
 *
 * Compile-time override:  make EXTRA_CFLAGS="-DTANDEM_IMG_RATIO=10"
 *
 * Originally: TANDEM 5 in ele.h
 */
#ifndef TANDEM_IMG_RATIO
#define TANDEM_IMG_RATIO  5
#endif
/* Backward-compat alias */
#define TANDEM  TANDEM_IMG_RATIO


/* ============================================================
 * eledef -- initial clustering constants   (originally in eledef.c)
 * ============================================================ */

/*
 * ELEDEF_IMGPROT_CAP  --  hard capacity of the all_iprot / iprot_shadow
 * arrays in eledef.c.  When the image-to-element assignment buffer fills
 * to this limit the DUMBBELL macro flushes it and processing continues.
 * Increasing this value reduces I/O at the cost of more peak memory.
 *
 * Compile-time override:  make EXTRA_CFLAGS="-DELEDEF_IMGPROT_CAP=1000000"
 *
 * Originally: IMG_CAP 500000 in eledef.c
 */
#ifndef ELEDEF_IMGPROT_CAP
#define ELEDEF_IMGPROT_CAP  500000
#endif
/* Backward-compat alias */
#define IMG_CAP  ELEDEF_IMGPROT_CAP


/*
 * ELEDEF_CUTOFF_SINGLE  --  default fractional-overlap cutoff used by
 * sing_cov() during initial element clustering when the "single" method
 * is chosen.  An image is admitted to an element if it overlaps the
 * current element boundary by at least this fraction of *either* length.
 *
 * Originally: literal 0.5 in eledef.c main()
 * IMPORTANT: must remain a double literal (no 'f' suffix) to match the
 * precision of the original bare literal and avoid flipping boundary comparisons.
 */
#define ELEDEF_CUTOFF_SINGLE  0.5


/*
 * ELEDEF_CUTOFF_DOUBLE  --  default fractional-overlap cutoff for
 * doub_cov() when the "double" method is chosen.  Requires the overlap
 * to be at least this fraction of *both* lengths.
 *
 * Originally: literal 0.9 in eledef.c main()
 * IMPORTANT: must remain a double literal (no 'f' suffix).
 */
#define ELEDEF_CUTOFF_DOUBLE  0.9


/*
 * MIN_OVERLAP_BP  --  minimum base-pair overlap required before a
 * fractional-coverage test (sing_cov / doub_cov) is even attempted.
 * If the current element's right boundary extends fewer than this many
 * bp past an image's left boundary the image is considered non-
 * overlapping without further calculation.
 *
 * This literal value appears identically in eledef.c and eleredef.c.
 *
 * Originally: literal 10 in eledef.c ele_def() and eleredef.c ele_def()
 */
#define MIN_OVERLAP_BP  10


/* ============================================================
 * eleredef -- element redefinition constants   (originally in eleredef.c)
 * ============================================================ */

/*
 * ELEREDEF_CUTOFF_SINGLE  --  fractional-overlap cutoff for sing_cov()
 * during the element-redefinition pass (ele_def inside eleredef).
 * Same semantic as ELEDEF_CUTOFF_SINGLE but applied to a different
 * clustering context (partitioning images into child elements after a
 * dissection step).
 *
 * Originally: CUTOFF1 0.5 in eleredef.c
 * IMPORTANT: must remain a double literal (no 'f' suffix).
 */
#ifndef ELEREDEF_CUTOFF_SINGLE
#define ELEREDEF_CUTOFF_SINGLE  0.5
#endif
/* Backward-compat alias */
#define CUTOFF1  ELEREDEF_CUTOFF_SINGLE


/*
 * ELEREDEF_CUTOFF_DOUBLE  --  fractional-overlap cutoff for doub_cov()
 * used in the dissect() function to determine whether an image is
 * "full length" relative to an element.
 *
 * Originally: CUTOFF2 0.9 in eleredef.c
 * IMPORTANT: must remain a double literal (no 'f' suffix).
 */
#ifndef ELEREDEF_CUTOFF_DOUBLE
#define ELEREDEF_CUTOFF_DOUBLE  0.9
#endif
/* Backward-compat alias */
#define CUTOFF2  ELEREDEF_CUTOFF_DOUBLE


/*
 * ELEREDEF_MAX_IMAGES  --  capacity of the img_ptr working array
 * allocated in eleredef main() and passed to general_ele_redef() and
 * edges_and_cps().  Sets an upper bound on the number of images that
 * can be loaded simultaneously for a single element's local network.
 * If an element exceeds this limit the program exits with an error.
 *
 * Compile-time override:  make EXTRA_CFLAGS="-DELEREDEF_MAX_IMAGES=2400"
 *
 * Originally: MAX_IMG 1200 in eleredef.c
 */
#ifndef ELEREDEF_MAX_IMAGES
#define ELEREDEF_MAX_IMAGES  1200
#endif
/* Backward-compat alias */
#define MAX_IMG  ELEREDEF_MAX_IMAGES


/*
 * MIN_ELEMENT_LEN_BP  --  minimum acceptable element length in base
 * pairs.  Elements (or sub-images produced by dissection) shorter than
 * this are treated as spurious alignment extensions and discarded.
 * The value corresponds to the maximum length that a pairwise alignment
 * tool is expected to spuriously extend past a true element boundary by
 * chance.
 *
 * Originally: TOO_SHORT 30 in eleredef.c
 */
#define MIN_ELEMENT_LEN_BP  30
/* Backward-compat alias */
#define TOO_SHORT  MIN_ELEMENT_LEN_BP


/*
 * SPLIT_RATIO_TOLERANCE  --  multiplier applied in span() when deciding
 * whether a candidate breakpoint (potential cut point / PCP) is
 * sufficiently supported to become a "to-be-determined" boundary (TBD).
 * The support count of a clustered PCP must be >= span() * FUDGE to
 * qualify.  Increasing this value makes the algorithm more conservative
 * about splitting elements.
 *
 * Originally: FUDGE 2 in eleredef.c
 */
#define SPLIT_RATIO_TOLERANCE  2
/* Backward-compat alias */
#define FUDGE  SPLIT_RATIO_TOLERANCE


/*
 * ELEMENT_BOUNDARY_MARGIN  --  a positional guard value used to
 * constrain where breakpoints may be placed.  In the span() context it
 * prevents breakpoints within this many bp of element boundaries.
 * NOTE: In the original code this constant was defined but appears to
 * be unused in the active code paths (only inside #if 0 blocks).
 * Retained for documentation completeness.
 *
 * Originally: MARGIN 10000 in eleredef.c
 */
#define ELEMENT_BOUNDARY_MARGIN  10000
/* Backward-compat alias */
#define MARGIN  ELEMENT_BOUNDARY_MARGIN


/*
 * MAX_DISSECT_PASSES  --  proximity threshold (in bp) used in ele_redef()
 * and TBD_merge() to decide whether a TBD boundary is "close enough" to
 * an element endpoint to adjust the endpoint rather than dissect the
 * element.  Also used in TBD_merge() to merge TBDs within this distance.
 * Despite the name "FLURRY", this is a distance threshold, not a count.
 *
 * Originally: FLURRY 10 in eleredef.c
 */
#define MAX_DISSECT_PASSES  10
/* Backward-compat alias */
#define FLURRY  MAX_DISSECT_PASSES


/* ============================================================
 * edgeredef -- edge redefinition constants   (originally in edgeredef.c)
 * ============================================================ */

/*
 * EDGEREDEF_EDGE_CUTOFF  --  element size-ratio threshold used by
 * edge_filt() and edge_update() to classify an edge as a "partial
 * peripheral similarity" (PPS) edge.  If the length ratio of one
 * element to the other (smaller/larger) is below this value the edge
 * is considered a PPS edge and is demoted from primary ('p') to
 * secondary ('S') status.  Only the highest-scoring PPS edge per
 * element is promoted to 'P' (confirmed primary).
 *
 * Originally: CUTOFF3 0.7 in edgeredef.c
 * IMPORTANT: must remain a double literal (no 'f' suffix).
 */
#define EDGEREDEF_EDGE_CUTOFF  0.7
/* Backward-compat alias */
#define CUTOFF3  EDGEREDEF_EDGE_CUTOFF


/* ============================================================
 * EDGE_REPAIR_GUARD  --  suppress spurious promotions in edge_repair()
 * ============================================================
 *
 * Compiling with -DEDGE_REPAIR_GUARD activates a partner check inside
 * edge_repair() (edgeredef.c Pass 2).  Before promoting an element's best
 * 'S' edge to 'P', the guard reads each adjacent element and checks whether
 * any of them already hold a 'p' or 'P' edge pointing back at the element
 * being considered.  If such a reciprocal primary edge exists the element is
 * already reachable in the primary-edge graph and the promotion is skipped.
 *
 * Without the guard, edge_repair() promotes the best remaining 'S' (PPS)
 * edge for any element whose own primary edges were all demoted by the PPS
 * filter.  This can bridge biologically unrelated families when a demoted
 * element's highest-scoring PPS edge points into a different repeat family.
 *
 * Cost: roughly doubles the I/O of edge_repair()'s second pass (one extra
 * partner read per candidate, with early termination on the first reciprocal
 * found).  This is a small fraction of total pipeline work.
 *
 * Incompatible with ORIGINAL_BUGS: the guard is compiled out when
 * ORIGINAL_BUGS is defined, since that flag restores the original no-op
 * assignment which makes the guard moot.
 *
 * To build with the guard:
 *   make TUNABLES="-DEDGE_REPAIR_GUARD"
 */


/* ============================================================
 * ORIGINAL_BUGS  --  regression / equivalence-testing flag
 * ============================================================
 *
 * Compiling with -DORIGINAL_BUGS restores two pre-existing bugs that
 * were present in the released code, allowing output-level equivalence
 * testing against the original binaries.
 *
 * DO NOT use -DORIGINAL_BUGS in production builds.
 *
 * Bug 1 — edgeredef.c edge_repair() (type promotion no-op)
 *   Original code:  hanger->to_edge->type == 'P';
 *   Fix:            hanger->to_edge->type  = 'P';
 *   Effect: When an element's only primary edge was demoted to 'S', the
 *   repair function was supposed to promote it back to 'P'.  The
 *   comparison operator means the promotion never occurred, so affected
 *   elements remained with secondary edges only and were subsequently
 *   treated as isolated (stat 'y') leaf nodes in famdef.
 *
 * Bug 2 — ele.h ele_read_in() redef list loading (lost list head)
 *   Original code:  ele_info->ele->redef = cur_ele_data->next;
 *   Fix:            ele_info->ele->redef = cur_ele_data;
 *   Effect: The newly allocated redef ELE_DATA_t node was immediately
 *   discarded; the redef list always remained NULL.  Cross-element
 *   redef relationships (used by eleredef's recluster step) were never
 *   populated when loading element files from disk.
 *
 * To build with original bugs (for equivalence testing):
 *   make TUNABLES="-DORIGINAL_BUGS"
 */


#endif /* RECON_DEFS_H */
