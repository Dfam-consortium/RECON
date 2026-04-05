/*
 * bolts.h  --  Global constants and shared state for the RECON pipeline
 *
 * This header is the root of the include chain:
 *
 *   bolts.h  <-  seqlist.h  <-  msps.h  <-  ele.h  <-  *.c
 *
 * It defines:
 *   - Pipeline-wide named constants (via recon_defs.h)
 *   - Two global variables that hold the sorted sequence-name table
 *     built by GetSeqNames() in seqlist.h and consulted by GetSeqIndex()
 *     and scan_msp() throughout the codebase.
 *
 * Global state
 * ------------
 * These globals are declared here and defined by the single translation
 * unit that includes this header:
 *
 *   seq_count        -- total number of sequences in the input list
 *   seq_name_table   -- dynamically allocated array of seq_count char*
 *                       pointers, each pointing to a NUL-terminated
 *                       sequence name of at most SEQ_NAME_MAX_LEN chars.
 *                       The array is sorted lexically so that
 *                       GetSeqIndex() can use binary search.
 *
 * Sequence-name interning
 * -----------------------
 * FRAG_t.seq_name and IMAGE_t.frag.seq_name are *not* copies of the
 * sequence name string -- they are pointers into seq_name_table[].  This
 * means identity comparisons (ptr == ptr) rather than strcmp() can be
 * used to test whether two fragments are on the same sequence, which is
 * the pattern used throughout eledef.c and eleredef.c.
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#ifndef BOLTS_H
#define BOLTS_H

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "recon_defs.h"


/*
 * seq_count  (was: seq_no)
 * Total number of sequences in the input seq_list file, as read from
 * the first line of that file by GetSeqNames().  Used as the upper
 * bound in GetSeqIndex() binary searches and as a bound check in
 * GetSeqNames() itself.
 */
int32_t seq_count;
/* Backward-compat alias -- remove after all call sites are migrated */
#define seq_no seq_count


/*
 * seq_name_table  (was: seq_names)
 * Dynamically allocated array of char* pointers, one per sequence.
 * Each pointer addresses a malloc'd string of at most SEQ_NAME_MAX_LEN
 * characters.  Allocated by GetSeqNames(); never freed (pipeline
 * programs are short-lived).  Must be lexically sorted on input for the
 * binary-search invariant in GetSeqIndex() to hold.
 */
char **seq_name_table;
/* Backward-compat alias */
#define seq_names seq_name_table


#endif /* BOLTS_H */
