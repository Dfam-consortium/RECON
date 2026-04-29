#!/bin/sh
#
# run_all_tests.sh  --  Run equivalence tests for every cached test set
#
# Usage:
#   run_all_tests.sh [--original-bugs] [test-sets-dir]
#
#   --original-bugs   pass through to equiv_test.sh (default: use bug-fixed build)
#   test-sets-dir     directory containing test-set subdirectories
#                     (default: test-sets/ relative to this script)
#
# Each subdirectory of test-sets-dir that contains a manifest file, seqnames,
# and msps.out is treated as a test set.  The ref_bin_dir stored in the
# manifest is used as the reference binary for that test.
#
# Prints a one-line PASS/FAIL summary per test set, then a final tally.
# Exit code: 0 if all tests pass, 1 if any fail.
#

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

ORIGINAL_BUGS_FLAG=""
if [ "${1:-}" = "--original-bugs" ]; then
    ORIGINAL_BUGS_FLAG="--original-bugs"
    shift
fi

TEST_SETS_DIR="${1:-$SCRIPT_DIR/../test-sets}"
TEST_SETS_DIR="$(cd "$TEST_SETS_DIR" && pwd)"

passed=0
failed=0
skipped=0

for ts_dir in "$TEST_SETS_DIR"/*/; do
    manifest="$ts_dir/manifest"
    seqnames="$ts_dir/seqnames"
    msps="$ts_dir/msps.out"

    # Skip directories that are not proper test sets
    if [ ! -f "$manifest" ] || [ ! -f "$seqnames" ] || [ ! -f "$msps" ]; then
        continue
    fi

    label=$(grep "^label=" "$manifest" | cut -d= -f2-)
    ref_bin_dir=$(grep "^ref_bin_dir=" "$manifest" | cut -d= -f2-)
    sig=$(basename "$ts_dir")

    if [ -z "$ref_bin_dir" ] || [ ! -d "$ref_bin_dir" ]; then
        printf "SKIP  [%s] %s  (ref_bin_dir not found: %s)\n" "$sig" "$label" "$ref_bin_dir"
        skipped=$((skipped + 1))
        continue
    fi

    out="${TMPDIR:-/tmp}/equiv_out_$$.txt"
    "$SCRIPT_DIR/equiv_test.sh" $ORIGINAL_BUGS_FLAG \
        "$ref_bin_dir" "$seqnames" "$msps" > "$out" 2>&1
    status=$?

    if [ $status -eq 0 ]; then
        printf "PASS  [%s] %s\n" "$sig" "$label"
        passed=$((passed + 1))
    else
        printf "FAIL  [%s] %s\n" "$sig" "$label"
        cat "$out"
        failed=$((failed + 1))
    fi
    rm -f "$out"
done

total=$((passed + failed + skipped))
echo ""
echo "Results: $passed passed, $failed failed, $skipped skipped ($total test sets found)"

[ $failed -eq 0 ]
