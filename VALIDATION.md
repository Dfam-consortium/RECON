# Validation and Equivalence Testing

## Overview

Because much of the refactoring and optimization work in this repository was
performed with AI coding assistance rather than by manual rewrite, output
correctness is validated systematically by comparing the pipeline's results
against those of the original released binary (RECON 1.08) -- input/output
equivalance testing.  This document describes the validation infrastructure,
the test scripts, and the `-DORIGINAL_BUGS` compilation flag that preserves
the original 1.08 behaviour.  

The short term goal is to make some long waited improvements to to the codebase
while maintaining API stability and output consistency in-line with the principals
of https://rewrites.bio. The long-term goal will be to develop unit tests for 
the core algorithms so that future changes can be validated not by equivalence 
testing but rather by fine-scale characterisation of functional correctness. 

---

## The `ORIGINAL_BUGS` Compilation Flag

Two pre-existing bugs were present in the released RECON 1.08 binary.
Fixing them changes the pipeline output, which would make a bitwise
comparison against the released binary impossible even when the refactoring
itself is correct.

To separate the two concerns, both bugs are preserved behind an
`#ifdef ORIGINAL_BUGS` guard in the source:

```c
#ifdef ORIGINAL_BUGS
  /* ... original (buggy) behaviour ... */
#else
  /* ... fixed behaviour ... */
#endif
```

Building with `-DORIGINAL_BUGS` re-enables the original behaviour, producing
output that is bit-for-bit identical to RECON 1.08 on the same input.
Building without the flag (the default) applies the bug fixes.

### Bug 1 — `edgeredef.c`: off-by-one in loop bound

The original loop iterated one element past the valid range under certain
conditions, reading uninitialized memory.  The fix corrects the upper bound.

### Bug 2 — `edgeredef.c`: incorrect PPS ratio denominator

The original edge-filter ratio used the wrong element's length in the
denominator of the PPS (Proportional Pairwise Similarity) calculation,
causing some edges to be incorrectly promoted or demoted.

---

## Equivalence Testing Scripts

All scripts live in the `scripts/` directory.

### `scripts/run_recon.sh`

Runs the complete four-stage RECON pipeline in a scratch directory and
prints the working directory path on success.

```
run_recon.sh <bin_dir> <seq_list> <msp_file> [num_sections] [work_dir]
```

Supports both pipeline variants automatically:

- **New pipeline (this repo, v1.09+):** `eledef` handles image spreading
  and sorting internally; `imagespread` is not needed.
- **Old pipeline (RECON 1.08):** if `imagespread` is present alongside
  `eledef` in `bin_dir`, it is run first to produce `images/images_sorted`
  before `eledef` is invoked.

### `scripts/compare_recon_runs.sh`

Compares the output files from two pipeline runs and prints a per-file
MATCH / DIFF summary.  Files compared:

| File | Content |
|---|---|
| `summary/eles` | One line per element: family, element, strand, sequence, start, end |
| `summary/families` | One line per family: index, copy count |
| `summary/fam_no` | Total family count |
| `summary/ele_no` | Total element count in families |
| `summary/naive_eles` | Initial (pre-redefinition) element definitions |
| `summary/redef_ele_no` | Element count after redefinition |
| `summary/ori_msp_no` | Original MSP count |
| `edge_redef_res/redef_stat` | Per-element status after edge redefinition |

### `scripts/equiv_test.sh`

Builds the refactored source, runs it alongside a reference binary, and
compares the outputs.

```
equiv_test.sh [--original-bugs] <ref_bin_dir> <seq_list> <msp_file> [num_sections]
equiv_test.sh --self-compare <seq_list> <msp_file> [num_sections]
```

**Modes:**

| Mode | Reference | Test build | Purpose |
|---|---|---|---|
| `<ref_bin_dir>` | external released binary | current src/ (bug-fixed) | verify correctness of bug fixes |
| `--original-bugs <ref_bin_dir>` | external released binary | current src/ with `-DORIGINAL_BUGS` | verify refactoring does not change output |
| `--self-compare` | src/ with `-DORIGINAL_BUGS` | src/ without flag | characterise the impact of bug fixes |

**Test-set caching:**

When run in `<ref_bin_dir>` mode, `equiv_test.sh` caches both the input
files and the reference pipeline output under `test-sets/<sig>/`, where
`<sig>` is the first 12 hex characters of the MD5 of `msps.out`.  On
subsequent runs with the same MSP file the cached reference output is reused,
making the test fast even for large datasets.

Cache layout:
```
test-sets/<sig>/
  manifest                 # key=value metadata (label, ref_bin_dir, created)
  seqnames                 # copy of seq_list
  msps.out                 # copy of msp_file
  reference/
    summary/{eles,families,fam_no,ele_no,naive_eles,redef_ele_no,ori_msp_no}
    edge_redef_res/redef_stat
```

To force a fresh reference run:
```sh
rm -rf test-sets/<sig>/reference/
```

**Exit codes:** 0 = outputs match, 1 = outputs differ, 2 = build or runtime error.

### `scripts/batch_equiv_test.sh`

Runs `equiv_test.sh` over multiple `round-N/` subdirectories found under a
RepeatModeler working directory.

```
batch_equiv_test.sh [--original-bugs] <ref_bin_dir> <data_dir> [num_sections]
```

Rounds are processed in ascending numeric order starting from round-2.
A round is included if it contains both `seqnames` and `msps.out`.

---

## How to Validate a New Change

### Check that refactoring does not change output

```sh
scripts/equiv_test.sh --original-bugs /usr/local/RECON \
    path/to/seqnames path/to/msps.out
```

A PASSED result means the refactored code with `-DORIGINAL_BUGS` produces
bit-for-bit identical output to the released binary.

### Check that bug fixes change output in the expected direction

```sh
scripts/equiv_test.sh --self-compare \
    path/to/seqnames path/to/msps.out
```

A FAILED result here is expected and intentional — it shows the scope of
the output change introduced by the bug fixes.

### Run across all rounds of a RepeatModeler dataset

```sh
scripts/batch_equiv_test.sh --original-bugs /usr/local/RECON \
    /path/to/RepeatModeler/run/
```
