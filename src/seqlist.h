/*
 * seqlist.h  --  Sequence-name table management for the RECON pipeline
 *
 * Provides two functions used by all five RECON programs:
 *
 *   GetSeqNames()   -- reads the seq_list file and populates the global
 *                      seq_name_table[] array.
 *
 *   GetSeqIndex()   -- binary searches seq_name_table[] for a given name
 *                      and returns its 0-based index, or -1 if not found.
 *
 * The seq_list file format
 * ------------------------
 * Line 1:   an integer giving the total number of sequences that follow.
 * Lines 2+: one sequence identifier per line, either as a bare name or
 *           in FASTA header format ("> name ...").  Only the first
 *           whitespace-delimited token after any leading ">" and spaces
 *           is stored.  Lines must be in lexical (ASCII) sort order for
 *           GetSeqIndex()'s binary search to be correct.
 *
 * Dependency on bolts.h
 * ---------------------
 * seq_count and seq_name_table are declared in bolts.h; this file
 * only reads and writes those globals.
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#ifndef __SEQLIST_H__
#define __SEQLIST_H__

#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "bolts.h"


/* ---- function prototypes ---- */
void GetSeqNames(FILE *);
int  GetSeqIndex(int, int, char *);


#endif /* __SEQLIST_H__ */
