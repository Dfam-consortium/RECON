#!/bin/bash
#
# compare_recon_runs.sh  --  Diff two RECON output directories
#
# Usage:
#   compare_recon_runs.sh <run_a_dir> <run_b_dir> [label_a] [label_b]
#
#   run_a_dir   working directory from first  run_recon.sh invocation
#   run_b_dir   working directory from second run_recon.sh invocation
#   label_a     human-readable label for run A (default: "run_a")
#   label_b     human-readable label for run B (default: "run_b")
#
# Exit codes:
#   0  all compared files are identical
#   1  one or more files differ
#   2  usage error or missing directory
#
# Files compared (all key pipeline outputs):
#   summary/eles          -- final element-to-family assignments
#   summary/families      -- family member counts
#   summary/fam_no        -- total number of families
#   summary/ele_no        -- total elements assigned to families
#   summary/naive_eles    -- element count after eledef
#   summary/redef_ele_no  -- element count after eleredef
#   summary/ori_msp_no    -- MSP count from imagespread
#   tmp/redef_stat        -- per-element status after eleredef
#                           (via edge_redef_res symlink, if present)
#
# For each file, the script reports:
#   - MATCH   : files are byte-identical (ignoring trailing whitespace)
#   - DIFF    : files differ; a unified diff is shown
#   - MISSING : file absent from one or both runs
#

set -u

if [ $# -lt 2 ]; then
    echo "usage: compare_recon_runs.sh run_a_dir run_b_dir [label_a] [label_b]" >&2
    exit 2
fi

A="$1"
B="$2"
LA="${3:-run_a}"
LB="${4:-run_b}"

if [ ! -d "$A" ]; then echo "error: directory not found: $A" >&2; exit 2; fi
if [ ! -d "$B" ]; then echo "error: directory not found: $B" >&2; exit 2; fi

# Files to compare relative to the run root.
# redef_stat lives in the edge_redef_res or eleredef_res symlink target.
COMPARE_FILES="
summary/eles
summary/families
summary/fam_no
summary/ele_no
summary/naive_eles
summary/redef_ele_no
summary/ori_msp_no
edge_redef_res/redef_stat
"

any_diff=0

compare_file() {
    relpath="$1"
    fa="$A/$relpath"
    fb="$B/$relpath"
    exists_a=0; exists_b=0
    [ -f "$fa" ] && exists_a=1
    [ -f "$fb" ] && exists_b=1

    if [ $exists_a -eq 0 ] && [ $exists_b -eq 0 ]; then
        printf "%-40s  MISSING from both runs\n" "$relpath"
        return
    fi
    if [ $exists_a -eq 0 ]; then
        printf "%-40s  MISSING from %s\n" "$relpath" "$LA"
        any_diff=1
        return
    fi
    if [ $exists_b -eq 0 ]; then
        printf "%-40s  MISSING from %s\n" "$relpath" "$LB"
        any_diff=1
        return
    fi

    # Strip trailing whitespace before comparing (pipeline outputs have
    # trailing spaces in several fields).
    if diff -q \
        <(sed 's/[[:space:]]*$//' "$fa") \
        <(sed 's/[[:space:]]*$//' "$fb") \
        > /dev/null 2>&1
    then
        printf "%-40s  MATCH\n" "$relpath"
    else
        printf "%-40s  DIFF\n" "$relpath"
        any_diff=1
        echo "--- $LA/$relpath"
        echo "+++ $LB/$relpath"
        diff \
            <(sed 's/[[:space:]]*$//' "$fa") \
            <(sed 's/[[:space:]]*$//' "$fb") \
            | head -80
        echo ""
    fi
}

echo "Comparing RECON outputs:"
echo "  A ($LA): $A"
echo "  B ($LB): $B"
echo ""

for f in $COMPARE_FILES; do
    compare_file "$f"
done

echo ""
if [ $any_diff -eq 0 ]; then
    echo "Result: ALL FILES MATCH"
    exit 0
else
    echo "Result: DIFFERENCES FOUND"
    exit 1
fi
