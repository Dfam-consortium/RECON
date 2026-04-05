#!/bin/sh
#
# batch_equiv_test.sh  --  Run equivalence tests across multiple round directories
#
# Usage:
#   batch_equiv_test.sh [--original-bugs] <ref_bin_dir> <data_dir> [num_sections]
#
#   --original-bugs  passed through to equiv_test.sh; builds with -DORIGINAL_BUGS
#                    to verify refactoring equivalence independent of bug fixes
#   ref_bin_dir      directory containing the REFERENCE (released) RECON binaries
#   data_dir         directory containing round-# subdirectories
#   num_sections     number of imagespread sections (default: 1)
#
# For each round-N subdirectory (N >= 2, ascending) that contains both
# round-N/seqnames and round-N/msps.out, calls equiv_test.sh to build
# the test binaries, run both pipelines, and compare outputs.
#
# Exit codes:
#   0  all rounds passed
#   1  one or more rounds failed or produced diffs
#   2  usage error or no testable rounds found
#

set -u

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

EXTRA_FLAGS=""
if [ $# -ge 1 ] && [ "$1" = "--original-bugs" ]; then
    EXTRA_FLAGS="--original-bugs"
    shift
fi

if [ $# -lt 2 ]; then
    echo "usage: batch_equiv_test.sh [--original-bugs] ref_bin_dir data_dir [num_sections]" >&2
    exit 2
fi

REF_BIN_DIR="$1"
DATA_DIR="$2"
NSECT="${3:-1}"

if [ ! -d "$REF_BIN_DIR" ]; then
    echo "error: ref_bin_dir not found: $REF_BIN_DIR" >&2
    exit 2
fi
if [ ! -d "$DATA_DIR" ]; then
    echo "error: data_dir not found: $DATA_DIR" >&2
    exit 2
fi

# ---- Collect testable rounds in ascending numeric order ----
rounds=""
n=2
while true; do
    rdir="$DATA_DIR/round-$n"
    if [ ! -d "$rdir" ]; then
        break
    fi
    if [ -f "$rdir/seqnames" ] && [ -f "$rdir/msps.out" ]; then
        rounds="$rounds $n"
    else
        echo "WARNING: $rdir exists but is missing seqnames or msps.out -- skipping" >&2
    fi
    n=$((n + 1))
done

if [ -z "$rounds" ]; then
    echo "error: no testable round-N directories (N>=2) found under $DATA_DIR" >&2
    exit 2
fi

echo "Found rounds:$rounds"
echo ""

# ---- Run equiv_test.sh for each round ----
passed=0
failed=0
skipped=0

for n in $rounds; do
    rdir="$DATA_DIR/round-$n"
    echo "======================================================"
    echo "  round-$n: $rdir"
    echo "======================================================"

    "$SCRIPT_DIR/equiv_test.sh" \
        $EXTRA_FLAGS \
        "$REF_BIN_DIR" \
        "$rdir/seqnames" \
        "$rdir/msps.out" \
        "$NSECT"
    status=$?

    if [ $status -eq 0 ]; then
        passed=$((passed + 1))
        echo "  round-$n: PASSED"
    elif [ $status -eq 1 ]; then
        failed=$((failed + 1))
        echo "  round-$n: FAILED (outputs differ; see equiv_test_ref/ equiv_test_fixed/ in cwd)"
    else
        skipped=$((skipped + 1))
        echo "  round-$n: ERROR (pipeline or build failure)"
    fi
    echo ""
done

# ---- Summary ----
total=$((passed + failed + skipped))
echo "======================================================"
echo "Batch result: $total rounds tested"
echo "  Passed:  $passed"
echo "  Failed:  $failed"
echo "  Errors:  $skipped"
echo "======================================================"

if [ $failed -gt 0 ] || [ $skipped -gt 0 ]; then
    exit 1
fi
exit 0
