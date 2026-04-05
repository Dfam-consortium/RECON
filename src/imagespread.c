/*
 * imagespread.c  --  Stage 1: distribute MSP images across output files
 *
 * Algorithm overview
 * ------------------
 * This is the first stage of the RECON pipeline.  It reads a sorted MSP
 * file and a sequence-name list, then "spreads" the images (one per MSP
 * half) across one or more output files (images/spread1, spread2, ...).
 *
 * Each MSP produces two IMAGE records:
 *   - an even-indexed record for the query interval
 *   - an odd-indexed record for the subject interval
 *
 * The records are written to output files partitioned by the lexical rank
 * of the query or subject sequence name in the seq_name_table[].  This
 * partitioning allows the next stage (eledef) to sort each partition
 * independently, keeping sort-memory bounded for large datasets.
 *
 * Output
 *   images/spread1 ... images/spreadN  -- partitioned IMAGE records, one
 *          per line: "<img_index> <score> <seq_name> <lb> <rb>"
 *   images/errors                      -- lines that failed scan_msp()
 *   summary/ori_msp_no                 -- total MSP count for eledef
 *
 * Usage
 *   imagespread seq_name_list msp_list [num_output_files] [-l log_level]
 *
 *   -l <level>  log verbosity: 0=silent 1=error 2=warn 3=info 4=debug
 *               (default: 3=info)
 *
 * Author: Zhirong Bao
 * Modifications: Robert Hubley, Institute for Systems Biology
 */

#include "seqlist.h"
#include "msps.h"
#include "recon_log.h"
#include "minunit.h"

/* ---- Per-program log level storage (required by recon_log.h) ---- */
int   recon_log_level = RECON_LOG_INFO;
FILE *recon_log_fp    = NULL;


int GetOutputInterval(char *, int, int);


/* ---- MinUnit unit tests ---- */
MU_TEST(test_check) {
  mu_check(1==1);
}

MU_TEST_SUITE(test_suite) {
  MU_RUN_TEST(test_check);
}


int main (int argc, char *argv[]) {
  FILE *msp_file, *seq_list;       /* input files */
  FILE **output, *err, *msp_no;    /* output files */
  int num_output_files = 1;        /* was: noof -- Number Of Output Files */
  int i;
  char line[150], output_name[50];
  MSP_t cur;
  int img_ct = -1;   /* image counter; starts at -1 so first ++ gives 0 */

  /* Strip the optional "-l <level>" flag before positional arg parsing */
  if (recon_parse_log_flag(&argc, argv)) {
    fprintf(stderr, "error: -l requires a numeric log level argument\n");
    exit(1);
  }

  /* Unit test invocation: imagespread -t */
  if (argc == 2 && strcmp(argv[1], "-t") == 0) {
    MU_RUN_SUITE(test_suite);
    MU_REPORT();
    return MU_EXIT_CODE;
  }

  /* Validate command line */
  if (argc < 3) {
    printf("usage:\n"
           "  imagespread seq_name_list msp_list [num_output_files] [-l level]\n"
           "  last parameter optional with default of 1\n"
           "  -l <level>  log verbosity: 0=silent 1=error 2=warn "
           "3=info(default) 4=debug\n");
    exit(1);
  }

  /* Open input files */
  seq_list = fopen(argv[1], "r");
  if (!seq_list) {
    printf("Input file for sequence name list %s not found.  Exit.\n", argv[1]);
    exit(2);
  }
  msp_file = fopen(argv[2], "r");
  if (!msp_file) {
    printf("Input MSP file %s not found.  Exit.\n", argv[2]);
    exit(2);
  }
  if (argc == 4) { num_output_files = atoi(argv[3]); }
  if (num_output_files < 1) num_output_files = 1;

  /* Open output files */
  err = fopen("images/errors", "w");
  if (!err) {
    printf("Can not open images/errors for output.  Exit.\n");
    exit(3);
  }
  /* Route RECON_LOG output to the same error file */
  recon_log_fp = err;

  msp_no = fopen("summary/ori_msp_no", "w");
  if (!msp_no) {
    printf("Can not open images/msp_no for output.  Exit.\n");
    exit(3);
  }

  output = (FILE **) malloc(num_output_files * sizeof(FILE *));
  for (i = 0; i < num_output_files; i++) {
    sprintf(output_name, "images/spread%d", i + 1);
    output[i] = fopen(output_name, "w");
    if (!output[i]) {
      printf("Can not open %s for output.  Exit.\n", output_name);
      exit(3);
    }
  }

  /* Load the sorted sequence-name table into seq_name_table[] */
  GetSeqNames(seq_list);

  /*
   * Main loop: read each MSP, emit two IMAGE lines (one query, one subject).
   *
   * The image index (img_ct) increments by 1 for each IMAGE emitted.
   * For MSP at file position k (0-based), the query image gets index 2k
   * and the subject image gets 2k+1.  GetOutputInterval() maps each image
   * to the appropriate output-file partition based on its sequence name's
   * rank in seq_name_table[].
   */
  while (fgets(line, 150, msp_file)) {
    if (scan_msp(&cur, line)) {
      RLOG_ERR("Wrong format:\n%s", line);
      continue;
    }

    /* Query-side image (even index) */
    img_ct++;
    i = GetOutputInterval(cur.query.frag.seq_name, seq_count, num_output_files);
    if (i < 0) {
      RLOG_ERR("Sequence name %s not found in the given list of sequences.  Exit.\n",
               cur.query.frag.seq_name);
      exit(4);
    }
    fprintf(output[i], "%d %d %s %d %d \n",
            img_ct, cur.score, cur.query.frag.seq_name,
            cur.query.frag.lb, cur.query.frag.rb);

    /* Subject-side image (odd index) */
    img_ct++;
    i = GetOutputInterval(cur.sbjct.frag.seq_name, seq_count, num_output_files);
    if (i < 0) {
      RLOG_ERR("Sequence name %s not found in the given list of sequences.  Exit.\n",
               cur.query.frag.seq_name);
      exit(4);
    }
    fprintf(output[i], "%d %d %s %d %d \n",
            img_ct, cur.score, cur.sbjct.frag.seq_name,
            cur.sbjct.frag.lb, cur.sbjct.frag.rb);
  }

  /* img_ct is the index of the last image written (0-based).
   * img_ct/2 + 1 gives the number of complete MSPs processed. */
  fprintf(msp_no, "%d\n", img_ct / 2 + 1);

  exit(0);
}


/*
 * GetOutputInterval  --  map a sequence name to an output-file partition
 *
 * Computes which of the 'noof' output files should receive IMAGE records
 * for the given sequence, based on the sequence's rank in seq_name_table[].
 *
 * The partitioning divides the seq_count sequences as evenly as possible
 * among num_output_files files:
 *   partition = floor(seq_rank / (seq_count / num_output_files))
 *
 * The last few sequences that overflow the floor division are assigned
 * to the final partition (num_output_files-1).
 *
 * Parameters
 *   seq_name        sequence identifier (interned pointer into seq_name_table)
 *   seq_count_in    total number of sequences (passed explicitly to avoid
 *                   relying on the global -- was: seq_no)
 *   num_output_files number of output-file partitions
 *
 * Returns the 0-based partition index, or -1 if seq_name is not found.
 */
int GetOutputInterval(char *seq_name, int seq_count_in, int num_output_files) {
  int interval;

  interval = GetSeqIndex(0, seq_count_in - 1, seq_name) /
             (seq_count_in / num_output_files);
  if (interval < 0) {
    printf("Sequence name %s not found in the given list of sequences.  Exit.\n",
           seq_name);
    return -1;
  }
  /* Clamp to last partition for the last few sequences */
  if (interval > num_output_files - 1) interval = num_output_files - 1;

  return interval;
}
