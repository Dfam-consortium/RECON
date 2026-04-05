/*
 * recon_log.h  --  Unified logging infrastructure for the RECON pipeline
 *
 * Overview
 * --------
 * Each of the five RECON programs previously had its own ad-hoc logging:
 *   - imagespread and eledef wrote diagnostics to a local FILE *err.
 *   - eleredef, edgeredef, and famdef wrote progress traces to the
 *     global FILE *log_file defined in ele.h.
 *
 * This header provides a single set of log-level macros that all programs
 * can use.  The existing log_file and err variables are not removed; new
 * code uses the macros below and existing code is migrated incrementally.
 *
 * Integration per program
 * -----------------------
 * Each .c file must provide storage for the two extern variables declared
 * below:
 *
 *   int   recon_log_level = RECON_LOG_INFO;   // default: info-level output
 *   FILE *recon_log_fp    = NULL;             // NULL => writes to stderr
 *
 * After opening the program's primary log file, set recon_log_fp so that
 * log output goes to the same destination:
 *
 *   log_file = fopen("tmp2/log", "a");
 *   recon_log_fp = log_file;                  // unify destinations
 *
 * Log levels
 * ----------
 *   RECON_LOG_SILENT (0) -- no output at all
 *   RECON_LOG_ERROR  (1) -- error conditions that require attention
 *   RECON_LOG_WARN   (2) -- warnings (reserved for future use)
 *   RECON_LOG_INFO   (3) -- one-line progress trace per element (default)
 *   RECON_LOG_DEBUG  (4) -- verbose diagnostics; replaces #if 0 debug blocks
 *
 * CLI flag
 * --------
 * The helper recon_parse_log_flag() strips a "-l <level>" option from
 * argv before the program's normal positional-argument parsing runs.
 * Call it as the very first statement in main():
 *
 *   recon_parse_log_flag(&argc, argv);
 *
 * The extended usage string for each program should include:
 *   "  -l <level>   logging verbosity: 0=silent 1=error 2=warn "
 *                   "3=info (default) 4=debug\n"
 *
 * Author: Robert Hubley, Institute for Systems Biology
 */

#ifndef RECON_LOG_H
#define RECON_LOG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* ---- Log level constants ---- */
#define RECON_LOG_SILENT  0
#define RECON_LOG_ERROR   1
#define RECON_LOG_WARN    2
#define RECON_LOG_INFO    3
#define RECON_LOG_DEBUG   4


/*
 * Storage provided by each .c file.
 *   recon_log_level : current verbosity threshold
 *   recon_log_fp    : destination file; NULL causes fallback to stderr
 */
extern int   recon_log_level;
extern FILE *recon_log_fp;


/*
 * RECON_LOG(level, fmt, ...)
 *
 * Core logging macro.  Writes a formatted message to recon_log_fp (or
 * stderr if recon_log_fp is NULL) when 'level' is at or below the
 * current recon_log_level.  The output is flushed immediately so that
 * log entries appear even if the program exits abnormally.
 */
#define RECON_LOG(level, fmt, ...)                               \
  do {                                                           \
    if ((level) <= recon_log_level) {                            \
      FILE *_rfp = recon_log_fp ? recon_log_fp : stderr;        \
      fprintf(_rfp, fmt, ##__VA_ARGS__);                        \
      fflush(_rfp);                                              \
    }                                                            \
  } while (0)


/* Convenience wrappers for the four active levels */
#define RLOG_ERR(fmt, ...)   RECON_LOG(RECON_LOG_ERROR, fmt, ##__VA_ARGS__)
#define RLOG_WARN(fmt, ...)  RECON_LOG(RECON_LOG_WARN,  fmt, ##__VA_ARGS__)
#define RLOG_INFO(fmt, ...)  RECON_LOG(RECON_LOG_INFO,  fmt, ##__VA_ARGS__)
#define RLOG_DBG(fmt, ...)   RECON_LOG(RECON_LOG_DEBUG, fmt, ##__VA_ARGS__)


/*
 * recon_parse_log_flag  --  strip "-l <level>" from argv before the
 * program's positional-argument parsing, and set recon_log_level.
 *
 * Parameters
 *   argc_p  pointer to main's argc
 *   argv    main's argv array (modified in place if -l is found)
 *
 * Returns 0 on success.  Returns 1 if "-l" appears as the last
 * argument with no value following it (the caller should print usage
 * and exit).
 *
 * Example
 *   int main(int argc, char *argv[]) {
 *     if (recon_parse_log_flag(&argc, argv)) {
 *       fprintf(stderr, "error: -l requires a numeric argument\n");
 *       exit(1);
 *     }
 *     // remaining positional argument parsing unchanged ...
 *   }
 */
static inline int recon_parse_log_flag(int *argc_p, char **argv) {
  int i, j;
  for (i = 1; i < *argc_p; i++) {
    if (argv[i][0] == '-' && argv[i][1] == 'l' && argv[i][2] == '\0') {
      if (i + 1 >= *argc_p) return 1;          /* -l with no value */
      recon_log_level = atoi(argv[i + 1]);
      /* Shift remaining args left by 2 to remove "-l <value>" */
      for (j = i; j < *argc_p - 2; j++) argv[j] = argv[j + 2];
      *argc_p -= 2;
      return 0;
    }
  }
  return 0;   /* absent is fine; default level already set */
}


#endif /* RECON_LOG_H */
