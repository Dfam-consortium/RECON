# RECON — De Novo Repeat Family Identification

RECON is a C pipeline that identifies repeat families from biological
sequences by clustering pairwise sequence alignments (MSPs) into elements
and grouping related elements into families.

**Original authors:** Zhirong Bao and Sean R. Eddy (Washington University
School of Medicine, 2002).  
**Continued development:** Robert Hubley, Institute for Systems Biology.

---

## AI Assistance Disclosure

This repository has been substantially refactored and optimized with the
assistance of AI coding agents (Claude, GitHub Copilot).

Correctness is validated by comparing pipeline output against the released
RECON 1.08 binary on a suite of real genomic datasets — not by manual code
review alone.  See [VALIDATION.md](VALIDATION.md) for the full description
of the equivalence-testing infrastructure.

---

## Table of Contents

1. [Building](#building)
2. [Running](#running)
3. [Input Format](#input-format)
4. [Output Format](#output-format)
5. [Algorithm Overview](#algorithm-overview)
6. [Pipeline Stages](#pipeline-stages)
7. [Source Code Organisation](#source-code-organisation)
8. [Element Database Format](#element-database-format)
9. [License](#license)

---

## Building

Run make from the project root (where this README lives):

```sh
make          # build eledef eleredef edgeredef famdef into bin/
make install  # same (binaries are written directly to bin/)
make clean    # remove compiled objects and binaries
```

The compiler can be changed by overriding `CC`:

```sh
make CC=cc
```

**Tunable constants** can be overridden on the command line:

```sh
make TUNABLES="-DELEREDEF_MAX_IMAGES=2400 -DBFS_CLAN_DEPTH=4"
```

To build with the original-bug behaviour (for equivalence testing against
RECON 1.08):

```sh
make TUNABLES="-DORIGINAL_BUGS"
```

See [VALIDATION.md](VALIDATION.md) for details on `ORIGINAL_BUGS` and the
equivalence-testing scripts.

---

## Running

The pipeline is driven by `scripts/run_recon.sh` (previously recon.pl):

```sh
scripts/run_recon.sh <bin_dir> <seq_list> <msp_file> [num_sections] [work_dir]
```

`run_recon.sh` automatically detects whether it is invoking the new
pipeline (this repo, v1.09+, where `imagespread` is merged into `eledef`)
or the old pipeline (RECON 1.08, where `imagespread` is a separate binary).

Alternatively, invoke each stage directly from the pipeline working
directory (see [Pipeline Stages](#pipeline-stages) for the directory
conventions each stage expects).

---

## Input Format

### MSP file

Pairwise local alignments in the MSP format — one alignment per line:

```
score %iden start end query_name start end subject_name
```

MSPs can be generated from BLAST output using the `MSPCollect` tool:

```sh
MSPCollect.pl blast_output.txt > alignments.msp
# for multiple BLAST outputs, concatenate:
cat run1.msp run2.msp > alignments.msp
```

### Sequence name list

The first positional argument to each stage is a file listing the input
sequence names, one per line.  The first line must be an integer giving
the total count.  Names must be sorted in lexical order.

```
3
chr1
chr2
chrX
```

The easiest way to generate this file is from the MSP file itself, since
every sequence that participated in at least one alignment is present as
either a query or subject name (columns 5 and 8):

```sh
awk '{print $5; print $8}' alignments.msp | sort -u > names.tmp
wc -l < names.tmp | cat - names.tmp > seqnames
rm names.tmp
```

Alternatively, generate it directly from the FASTA file that was searched:

```sh
grep '^>' sequences.fa | sed 's/^>//' | awk '{print $1}' | sort -u > names.tmp
wc -l < names.tmp | cat - names.tmp > seqnames
rm names.tmp
```

The FASTA-based approach is more correct if some sequences had no
alignments at all (they would be absent from the MSP file but should still
appear in the name list).

---

## Output Format

All output files are written to the `summary/` directory in the pipeline
working directory.

### `summary/eles`

One line per element that was assigned to a family:

```
family_index  element_index  strand  sequence_name  start  end
```

`strand` is `1` (same orientation as the family prototype) or `-1`
(reverse complement).

### `summary/families`

One line per family:

```
family_index  copy_count  unknown
```

### `summary/fam_no` / `summary/ele_no`

Single integers: the total number of families and the total number of
elements assigned to families, respectively.

### Intermediate files

The pipeline also produces several intermediate directories.  Unless you
need to inspect individual elements, these can be removed after a run:

| Directory | Stage that writes it | Contents |
|---|---|---|
| `ele_store/` | eledef → famdef | Flat-file element database (see [Element Database Format](#element-database-format)) |
| `ele_def_res/` | eledef | Per-stage working files (old pipeline: per-element `e<N>` text files) |
| `ele_redef_res/` | eleredef | Per-stage working files |
| `edge_redef_res/` | edgeredef | Per-stage working files; `redef_stat` used by famdef |
| `summary/` | all stages | Final output + intermediate counts |

---

## Algorithm Overview

RECON defines repeat families by asking: *which genomic intervals are
copies of the same ancestral sequence?*

The answer is built from pairwise local alignments (Maximum Scoring Pairs,
MSPs), which serve as evidence that two genomic intervals share a common
origin.  The pipeline proceeds in four stages that progressively refine a
graph of similarity relationships:

```
MSPs  →  [eledef]  →  [eleredef]  →  [edgeredef]  →  [famdef]  →  families
           images        elements       edge graph       BFS
          clusters       refined       filtered
```

### Key concepts

**Image** — one side of an MSP alignment viewed from the perspective of a
specific element.  If an MSP connects sequence regions A and B, then A is
an image from B's perspective and vice versa.

**Element** — a genomic interval produced by clustering overlapping images
using single-coverage clustering.  An element represents one copy of a
repeat at a specific position in the genome.

**Edge** — a similarity relationship between two elements.  A *primary*
edge (type `p`) exists when one or more full-length images connect the two
elements.  A *secondary* edge (`s`) exists when only partial images connect
them.

**Family** — a connected component in the element similarity graph formed
by following primary edges.  All elements in the same family are copies of
the same ancestral repeat.

---

## Pipeline Stages

### Stage 1: `eledef` — Initial element definition

`eledef` reads the MSP file, sorts alignments by sequence name and
position, and clusters overlapping images into initial elements using
single-coverage greedy clustering.

Each element is stored in the flat-file element database
(`ele_store/elements.db`).  A naive element summary is written to
`summary/naive_eles`.

**Key algorithm:** greedy single-coverage clustering — images sorted by
left bound are greedily merged into an element as long as they overlap the
current element's extent by the `ELEREDEF_CUTOFF_SINGLE` fraction.

In the original 1.08 release this stage was preceded by a separate
`imagespread` binary that partitioned the MSP file into sorted sections.
That step is now performed internally by `eledef`.

### Stage 2: `eleredef` — Element redefinition (syntopy algorithm)

This is the most complex stage.  It refines element boundaries by
identifying positions where an element should be split (Potential Cut
Points, PCPs), then recursively re-clustering the dissected images.

The algorithm operates on *local networks* (clans): BFS subgraphs centred
on each unprocessed element and extending `BFS_CLAN_DEPTH` hops along
edges.  This bounds peak memory to the size of one clan rather than the
whole dataset.

**Element status machine:**

| Status | Meaning | Next step |
|--------|---------|-----------|
| `z` | Initial; not yet processed | `edges_and_cps()` → `t` |
| `t` | Edges and PCPs computed | `ele_redef()` → `v` or `w` |
| `w` | Split into children (combo) | `combo_update()` → `X` |
| `v` | Redefinition complete | Ready for edgeredef |
| `X` | Dismissed; zero images remain | Written to `summary/obsolete` |
| `O` | Large tandem repeat; filtered out | Excluded from families |

**PCP / TBD clustering pipeline:**

1. For each element, the endpoints of all full-length images are collected
   as Potential Cut Points (PCPs) on the element's coordinate system.
2. PCPs are sorted by coordinate and grouped into Boundary (BD) candidates:
   consecutive PCPs within 20 bp of the cluster start and 10 bp of the
   running last are merged into one BD.
3. A minimum support threshold is computed for each BD using `span()`:
   the number of images that start before the cut position and end before
   it, multiplied by `SPLIT_RATIO_TOLERANCE`.
4. BDs with support ≥ threshold are promoted to TBDs; the rest are freed.
5. Adjacent TBDs within 10 bp of each other are merged.

**Dissection:**

TBDs that lie more than `MAX_DISSECT_PASSES` bp from both element
boundaries trigger dissection.  Each image spanning a TBD is split into
left and right halves; a new MSP record is created for the left half and
the original is trimmed to the right half.  Images shorter than
`MIN_ELEMENT_LEN_BP` after trimming are discarded.  The remaining images
are re-clustered into child elements.

**Edge and CP construction:**

For each element, its images are grouped by partner element.  An image
qualifies as *full-length* if it starts within 10 bp of the element's left
boundary, ends within 10 bp of its right boundary, and covers at least
`ELEREDEF_CUTOFF_DOUBLE` of the element's length.

If no single full-length image exists between two elements, a *consistency
tree* is searched for a chain of partial images that collectively tile a
full-length connection.  The consistency tree organises partial images by
their non-overlapping relationship: each level holds images consistent
(non-overlapping, same element pair, same direction) with the image at the
level above.

### Stage 3: `edgeredef` — Edge graph filtering

`edgeredef` applies a Proportional Pairwise Similarity (PPS) filter to
remove spurious edges.  For each element, primary edges are ranked by score.
Edges whose score falls below `EDGEREDEF_EDGE_CUTOFF` times the best edge
for that element are demoted to type `S`; the winner is promoted to `P`.

Elements with no remaining primary edges are assigned status `O` (outcast)
and excluded from family building.

### Stage 4: `famdef` — Family definition by BFS

`famdef` performs a breadth-first search over the element graph following
only primary (`p` / `P`) edges.  Each connected component becomes one
repeat family.  The BFS also tracks cumulative orientation along each path
from the seed, recording whether each member is on the same or opposite
strand as the family prototype.

Final output is written to `summary/eles` and `summary/families`.  At the
end of this stage the element database is compacted to remove deprecated
records.

---

## Source Code Organisation

### Original structure (RECON 1.08)

```
eledef.c        ~700 lines; initial element definition (includes imagespread logic)
eleredef.c      ~2200 lines; all of stage 2 in one file
edgeredef.c     ~400 lines; edge filtering
famdef.c        ~300 lines; BFS family assignment
ele.c           ~900 lines; shared data structure library
```

### Refactored structure (v1.09)

**`imagespread` merged into `eledef`**

The separate `imagespread` binary partitioned MSPs into sections to work
around memory limits.  Modern hardware makes this unnecessary; the
partitioning logic was folded into `eledef` and the intermediate
`images/spread*` files eliminated.

**Flat-file element database**

The original pipeline wrote each element to its own text file
(`tmp/e1`, `tmp/e2`, …), creating tens of thousands of small files.
These were replaced by a single indexed flat-file database.  See
[Element Database Format](#element-database-format).

**`eleredef.c` split into four focused files**

The 2,500-line `eleredef.c` was reorganised by logical responsibility:

| File | Lines | Responsibility |
|---|---|---|
| `eleredef.c` | ~630 | `main()`, reporting, clan BFS orchestration |
| `redef_boundary.c` | ~310 | CP/BD/PCP clustering; boundary detection |
| `redef_edges.c` | ~610 | Edge graph construction; consistency tree |
| `redef_dissect.c` | ~700 | Element construction; dissection; combo bookkeeping |

Each module has a corresponding header declaring its public interface.

**Tunable constants in `recon_defs.h`**

Magic numbers formerly `#define`d locally in each source file are now in
`recon_defs.h` with descriptive names, overridable at build time:

| `recon_defs.h` name | Original | Default |
|---|---|---|
| `ELEREDEF_CUTOFF_SINGLE` | `CUTOFF1` | 0.5 |
| `ELEREDEF_CUTOFF_DOUBLE` | `CUTOFF2` | 0.9 |
| `ELEREDEF_MAX_IMAGES` | `MAX_IMG` | 1200 |
| `MIN_ELEMENT_LEN_BP` | `TOO_SHORT` | 30 |
| `SPLIT_RATIO_TOLERANCE` | `FUDGE` | 2 |
| `MAX_DISSECT_PASSES` | `FLURRY` | 10 |
| `BFS_CLAN_DEPTH` | `DEPTH` | 3 |
| `EDGEREDEF_EDGE_CUTOFF` | `CUTOFF3` | 0.7 |

**Performance optimisations**

Two major I/O bottlenecks were fixed:

- *O(N²) `report_redef_stat()` calls:* The original `edgeredef` rewrote
  the full status of all N elements to disk after processing each element —
  O(N²) file writes on a 122,749-element dataset.  The fix calls
  `report_redef_stat()` once after all elements are processed.

- *`fflush()` in inner loops:* The original code called `fflush()` after
  every `fprintf` in per-element and per-clan loops.  Non-critical flushes
  inside tight loops were removed; end-of-clan flushes were retained.

Combined result: approximately **12× speedup** on a large real-world dataset
(20 min → 1 min 43 sec) with identical output.

---

## Element Database Format

### Motivation

The original RECON wrote each element to its own file (`tmp/e<N>`), creating
tens of thousands of small files in hot directories on large runs.  The new
database consolidates all elements into two files in `ele_store/`:

```
ele_store/
  elements.db    sequential records (text)
  elements.idx   binary offset index
```

### `elements.db` — record file

Records are appended sequentially.  Each record has a fixed 23-byte ASCII
header followed by a variable-length text body:

```
ELE<8-digit-index><status><10-digit-data-len>\n<data...>
```

| Field | Bytes | Description |
|---|---|---|
| `ELE` | 0–2 | Magic prefix |
| index | 3–10 | Zero-padded element index (1-based) |
| status | 11 | `A` = active, `D` = deprecated |
| data-len | 12–21 | Zero-padded byte count of the following data |
| `\n` | 22 | Header terminator |
| data | 23… | Element text in the existing element-file format |

**Deprecation:** when an element is updated, its old record's status byte
is overwritten with `D` in-place (a single-byte seek-and-write, no data
movement).  The new version is appended at the end.  This makes updates
O(1) regardless of record size.

**Compaction:** at the end of `famdef`, `ele_db_compact()` rewrites the
database with only active records in index order, atomically replacing the
originals via `rename()`.

### `elements.idx` — index file

A flat binary array of `int64_t` values (8 bytes each), indexed by
`ele_index - 1`.  Each entry holds the byte offset of the element's active
record in `elements.db`, or `-1` if no record has been written yet.

**Lookup:** reading element N requires one seek into the index, one 8-byte
read to get the record offset, one seek into the database, and one read of
the 23-byte header + data body.

### Element record body

The text format inside each record is unchanged from the original per-file
format:

```
frag     <seq_name> <lb> <rb>
direc    <direction>
update   <update_flag>
img_no   <image_count>
flimg_no <full_length_image_count>
edge_no  <edge_count>
msp      <img_id> <stat> <score> <iden> <dir> \
             <ele1> <qseq> <qlb> <qrb>  <ele2> <sseq> <slb> <srb>
...
edge     <idx> <type> <dir> <score> <ele1_idx> <ele2_idx>
...
pcp      <position> <contributor_ele_idx>
...
tbd      <position> <support>
...
redef    <child_ele_idx>
...
```

---

## License

Copyright (C) 2001 Washington University School of Medicine  
Copyright (C) 2011–2026 Institute for Systems Biology

Source code is freely available under the GNU General Public License.
See the `LICENSE` file for details.

See [CHANGELOG.md](CHANGELOG.md) for version history.
