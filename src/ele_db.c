/*
 * ele_db.c  --  Flat-file element database implementation
 *
 * See ele_db.h for the on-disk format.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "ele_db.h"

/* Header layout (all text, 23 bytes total):
 *   bytes  0-2  : "ELE"
 *   bytes  3-10 : 8-digit zero-padded element index
 *   byte   11   : status: 'A' = active, 'D' = deprecated
 *   bytes 12-21 : 10-digit zero-padded data length
 *   byte   22   : '\n'
 */
#define DB_HDR_LEN   23
#define DB_STAT_OFF  11   /* byte offset of status within the header */

static FILE *db_fp  = NULL;   /* elements.db  -- read+write */
static FILE *idx_fp = NULL;   /* elements.idx -- read+write */

static const int64_t IDX_EMPTY = -1;   /* sentinel: no record for this index */


void ele_db_open(void) {
  /* Open or create elements.db */
  db_fp = fopen("ele_store/elements.db", "r+b");
  if (!db_fp) db_fp = fopen("ele_store/elements.db", "w+b");
  if (!db_fp) { perror("ele_db_open: cannot open elements.db"); exit(1); }

  /* Open or create elements.idx */
  idx_fp = fopen("ele_store/elements.idx", "r+b");
  if (!idx_fp) idx_fp = fopen("ele_store/elements.idx", "w+b");
  if (!idx_fp) { perror("ele_db_open: cannot open elements.idx"); exit(1); }
}


void ele_db_close(void) {
  if (db_fp)  { fclose(db_fp);  db_fp  = NULL; }
  if (idx_fp) { fclose(idx_fp); idx_fp = NULL; }
}


char *ele_db_read(int ele_index, int *data_len_out) {
  int64_t offset = IDX_EMPTY;
  char    hdr[DB_HDR_LEN + 1];
  int     data_len;
  char   *buf;

  /* Look up byte offset in index */
  if (fseek(idx_fp, (int64_t)(ele_index - 1) * 8, SEEK_SET) != 0) return NULL;
  if (fread(&offset, 8, 1, idx_fp) != 1 || offset == IDX_EMPTY)   return NULL;

  /* Read and validate header */
  if (fseek(db_fp, offset, SEEK_SET) != 0)      return NULL;
  if (fread(hdr, DB_HDR_LEN, 1, db_fp) != 1)    return NULL;
  hdr[DB_HDR_LEN] = '\0';

  if (hdr[DB_STAT_OFF] == 'D') return NULL;   /* should not happen if index is correct */

  /* Parse data length */
  data_len = atoi(hdr + 12);

  /* Read data */
  buf = (char *) malloc(data_len + 1);
  if (!buf) { perror("ele_db_read: malloc"); exit(1); }
  if ((int)fread(buf, 1, data_len, db_fp) != data_len) {
    free(buf);
    return NULL;
  }
  buf[data_len] = '\0';

  if (data_len_out) *data_len_out = data_len;
  return buf;
}


void ele_db_write(int ele_index, const char *data, int data_len) {
  int64_t idx_pos   = (int64_t)(ele_index - 1) * 8;
  int64_t old_offset = IDX_EMPTY;
  int64_t new_offset;
  char    hdr[DB_HDR_LEN + 1];

  /* Read existing index entry; ignore read failure: old_offset stays IDX_EMPTY */
  if (fseek(idx_fp, idx_pos, SEEK_SET) == 0) {
    size_t nr = fread(&old_offset, 8, 1, idx_fp);
    (void)nr;
  }

  /* Deprecate the old record if one exists */
  if (old_offset != IDX_EMPTY) {
    if (fseek(db_fp, old_offset + DB_STAT_OFF, SEEK_SET) == 0) {
      fputc('D', db_fp);
      fflush(db_fp);
    }
  }

  /* Append new record at end of database */
  if (fseek(db_fp, 0, SEEK_END) != 0) { perror("ele_db_write: seek"); exit(1); }
  new_offset = ftell(db_fp);

  snprintf(hdr, sizeof(hdr), "ELE%08d%c%010d\n", ele_index, 'A', data_len);
  if (fwrite(hdr,  DB_HDR_LEN, 1, db_fp) != 1 ||
      fwrite(data, data_len,   1, db_fp) != 1) {
    perror("ele_db_write: fwrite");
    exit(1);
  }
  fflush(db_fp);

  /* Update index */
  if (fseek(idx_fp, idx_pos, SEEK_SET) != 0) { perror("ele_db_write: idx seek"); exit(1); }
  if (fwrite(&new_offset, 8, 1, idx_fp) != 1) { perror("ele_db_write: idx write"); exit(1); }
  fflush(idx_fp);
}


void ele_db_compact(void) {
  FILE    *new_db, *new_idx;
  int64_t  idx_size, n_entries, i;
  int64_t  old_offset, new_offset;
  char     hdr[DB_HDR_LEN + 1];
  int      data_len;
  char    *buf;

  /* How many index entries are there? */
  if (fseek(idx_fp, 0, SEEK_END) != 0) { perror("ele_db_compact: idx seek"); exit(1); }
  idx_size  = ftell(idx_fp);
  n_entries = idx_size / 8;

  new_db  = fopen("ele_store/elements_new.db",  "wb");
  new_idx = fopen("ele_store/elements_new.idx", "wb");
  if (!new_db || !new_idx) { perror("ele_db_compact: create"); exit(1); }

  new_offset = 0;

  for (i = 0; i < n_entries; i++) {
    int64_t this_new_offset = -1;

    /* Read index entry */
    if (fseek(idx_fp, i * 8, SEEK_SET) != 0) { perror("ele_db_compact: idx seek"); exit(1); }
    if (fread(&old_offset, 8, 1, idx_fp) != 1) old_offset = IDX_EMPTY;

    if (old_offset != IDX_EMPTY && old_offset >= 0) {
      /* Read and validate record header */
      if (fseek(db_fp, old_offset, SEEK_SET) == 0 &&
          fread(hdr, DB_HDR_LEN, 1, db_fp) == 1) {
        hdr[DB_HDR_LEN] = '\0';

        if (hdr[DB_STAT_OFF] == 'A') {
          data_len = atoi(hdr + 12);
          buf = (char *) malloc(data_len);
          if (!buf) { perror("ele_db_compact: malloc"); exit(1); }

          if ((int)fread(buf, 1, data_len, db_fp) == data_len) {
            /* Write to new database */
            this_new_offset = new_offset;
            /* Update status in header to A (already A, just rewrite cleanly) */
            fwrite(hdr, DB_HDR_LEN, 1, new_db);
            fwrite(buf, data_len,   1, new_db);
            new_offset += DB_HDR_LEN + data_len;
          }
          free(buf);
        }
      }
    }

    /* Write new index entry */
    if (fseek(new_idx, i * 8, SEEK_SET) != 0) { perror("ele_db_compact: new idx seek"); exit(1); }
    fwrite(&this_new_offset, 8, 1, new_idx);
  }

  fclose(new_db);
  fclose(new_idx);

  /* Replace originals atomically */
  fclose(db_fp);
  fclose(idx_fp);
  db_fp  = NULL;
  idx_fp = NULL;

  rename("ele_store/elements_new.db",  "ele_store/elements.db");
  rename("ele_store/elements_new.idx", "ele_store/elements.idx");
}
