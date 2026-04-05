/*
 * ele_db.h  --  Flat-file element database: single indexed store replacing
 *               the per-element tmp/e<N> / tmp2/e<N> file pairs.
 *
 * Format
 * ------
 * Two files live in ele_store/:
 *
 *   elements.db   Sequential records, each with a fixed 23-byte text header:
 *                   ELE<8-digit-index><status><10-digit-data-len>\n
 *                   <data_len bytes of text: existing element file format>
 *                 status = 'A' (active) or 'D' (deprecated).
 *                 Deprecated records are marked in-place (single-byte write
 *                 at a fixed offset) when an element is updated; the new
 *                 version is appended at the end.
 *
 *   elements.idx  Binary array of int64_t offsets, one per element (0-indexed
 *                 by ele_index-1).  Value = byte offset of the active record
 *                 in elements.db; -1 means not yet written.  The file grows
 *                 via fseek-and-write as new element indices are created.
 *
 * The database is opened once per pipeline stage (ele_db_open), updated via
 * ele_db_write, read via ele_db_read, and compacted at the end of famdef
 * (ele_db_compact) to produce a clean, deprecation-free final file.
 */
#ifndef __ELE_DB_H__
#define __ELE_DB_H__

#include <stdio.h>
#include <stdint.h>

/* Open or reopen the element database.
 * Must be called once at the start of each pipeline stage.
 * ele_store/ must already exist. */
void ele_db_open(void);

/* Close the element database file handles. */
void ele_db_close(void);

/* Read the text content of element ele_index.
 * Returns a malloc'd, NUL-terminated buffer the caller must free.
 * Returns NULL if ele_index has no record in the database. */
char *ele_db_read(int ele_index, int *data_len_out);

/* Write (or update) the element record for ele_index.
 * If a record already exists, it is deprecated in-place and the new
 * version is appended.  The index is updated atomically. */
void ele_db_write(int ele_index, const char *data, int data_len);

/* Compact the database: rewrite elements.db containing only active records
 * in index order, rebuild elements.idx, and replace the originals.
 * Called once at the end of famdef. */
void ele_db_compact(void);

#endif /* __ELE_DB_H__ */
