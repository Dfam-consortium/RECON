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



/*
 * GetSeqNames  --  read the seq_list file into seq_name_table[].
 *
 * Reads the leading integer from seq_list (the sequence count), then
 * reads that many name lines.  For each line:
 *   - Skips a leading '>' and any whitespace following it (FASTA headers).
 *   - Copies up to SEQ_NAME_MAX_LEN-1 chars of the first whitespace-
 *     delimited token into a freshly allocated buffer.
 *   - Stores a pointer to that buffer in seq_name_table[].
 *
 * Preconditions
 *   seq_list must be positioned at the beginning of the file.
 *   The names in seq_list must be sorted lexically; GetSeqIndex() relies
 *   on this invariant.
 *
 * Side effects
 *   Allocates seq_name_table[] and each name buffer (never freed).
 *   Sets seq_count to the declared count from the first line.
 *   Exits with code 2 on format errors or count mismatches.
 */
void GetSeqNames(FILE *seq_list) {
  char line[256];
  int seq_read_count = -1;   /* was: seq_ct  -- counts names read so far */
  char *name_start;

  /* Read the mandatory first line: the total sequence count. */
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-result"
  fgets(line, 255, seq_list);   /* return unused: EOF/error caught by atoi below */
#pragma GCC diagnostic pop
  seq_count = atoi(line);
  if (!seq_count) {
    printf("First line in the list of sequences must be an integer,\n"
           "which is the number of sequences in the list.  Exit\n");
    exit(2);
  }

  seq_name_table = (char **) malloc(seq_count * sizeof(char *));

  /* Read each sequence name into the table. */
  while (fgets(line, 255, seq_list)) {
    seq_read_count++;
    if (seq_read_count == seq_count) {
      printf("More sequence names than indicated at the beginning of "
             "the file.  Exit.\n");
      exit(2);
    }

    seq_name_table[seq_read_count] = (char *) malloc(SEQ_NAME_MAX_LEN * sizeof(char));

    /* Skip the leading '>' and any spaces that follow it, so that names
     * grepped directly from a FASTA file work without pre-processing. */
    name_start = line;
    if (line[0] == '>') {
      name_start++;
      while (isspace(*name_start)) name_start++;
    }

    /* Copy at most SEQ_NAME_MAX_LEN-1 characters, then NUL-terminate
     * at the first whitespace boundary. */
    strncpy(seq_name_table[seq_read_count], name_start, SEQ_NAME_MAX_LEN - 1);
    name_start = seq_name_table[seq_read_count];
    while (!isspace(*name_start)) name_start++;
    *name_start = '\0';
  }
}




/*
 * GetSeqIndex  --  binary search of seq_name_table[] for seq_name.
 *
 * Parameters
 *   left, right  inclusive index bounds for this search step (call with
 *                left=0, right=seq_count-1 from external callers).
 *   seq_name     the name to look up.
 *
 * Returns the 0-based index into seq_name_table[] if found, or -1 if
 * seq_name is not present.
 *
 * Precondition: seq_name_table[] must be sorted in the same lexical
 * order as produced by a standard sort of the input seq_list file.
 *
 * Note: strncmp is used (not strcmp) to guard against names longer than
 * SEQ_NAME_MAX_LEN; only the first SEQ_NAME_MAX_LEN chars are compared.
 */
int GetSeqIndex(int left, int right, char *seq_name) {
  int pos, cmp_result;

  if (left <= right) {
    pos = (left + right) / 2;
    cmp_result = strncmp(seq_name, seq_name_table[pos], SEQ_NAME_MAX_LEN);
    if (cmp_result < 0) return GetSeqIndex(left, pos - 1, seq_name);
    else if (cmp_result > 0) return GetSeqIndex(pos + 1, right, seq_name);
    else return pos;
  } else {
    return -1;
  }
}


#endif /* __SEQLIST_H__ */
