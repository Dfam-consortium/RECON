#!/bin/sh
#
# equiv_test.sh  --  Build two RECON variants and check output equivalence
#
# Usage:
#   equiv_test.sh [--original-bugs] <ref_bin_dir> <seq_list> <msp_file> [num_sections]
#   equiv_test.sh --self-compare <seq_list> <msp_file> [num_sections]
#
# Modes:
#
#   <ref_bin_dir>          Compare against an external reference binary tree
#                          (e.g. /usr/local/RECON).  By default, builds the
#                          current src/ WITHOUT -DORIGINAL_BUGS (bug-fixed).
#                          Use --original-bugs to build WITH -DORIGINAL_BUGS
#                          to verify the refactoring alone is equivalent to
#                          the released binary.
#
#   --self-compare         Build both variants from src/ and compare them:
#                            reference = ORIGINAL_BUGS build (released behaviour)
#                            fixed     = current build     (bug-fixed behaviour)
#                          Use this to characterise the impact of bug fixes.
#                          Caching is skipped in this mode.
#
# --original-bugs          When used with <ref_bin_dir>, build src/ with
#                          -DORIGINAL_BUGS so bug fixes don't muddy the
#                          refactoring equivalence check.
#
# num_sections defaults to 1 (matching original recon.pl behaviour).
#
# Test-set caching:
#   When used with <ref_bin_dir>, each unique msps.out is stored as a
#   self-contained test set under:
#
#     test-sets/<sig>/
#       manifest              key=value metadata (label, sources, ref_bin_dir)
#       seqnames              copy of the seq_list file
#       msps.out              copy of the msp_file
#       reference/            cached reference pipeline output (key files only)
#         summary/{eles,families,fam_no,ele_no,naive_eles,redef_ele_no,ori_msp_no}
#         edge_redef_res/redef_stat
#
#   <sig> is the first 12 hex characters of the md5 of msps.out, giving a
#   stable, collision-resistant directory name independent of file paths.
#
#   The manifest records the ref_bin_dir used to populate reference/.  If a
#   subsequent call uses a different ref_bin_dir the cached reference is
#   considered stale and is regenerated.
#
#   To force a fresh reference run without changing ref_bin_dir:
#     rm -rf test-sets/<sig>/reference/
#
# Exit codes:
#   0  outputs match
#   1  outputs differ
#   2  build or runtime error
#

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJ_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TEST_SETS_DIR="$SCRIPT_DIR/../test-sets"

# ---- Parse flags ----
SELF_COMPARE=0
USE_ORIGINAL_BUGS=0

while [ $# -gt 0 ]; do
    case "$1" in
        --self-compare)   SELF_COMPARE=1;       shift ;;
        --original-bugs)  USE_ORIGINAL_BUGS=1;  shift ;;
        *)                break ;;
    esac
done

if [ $SELF_COMPARE -eq 1 ] && [ $USE_ORIGINAL_BUGS -eq 1 ]; then
    echo "error: --self-compare and --original-bugs are mutually exclusive" >&2
    exit 2
fi

if [ $SELF_COMPARE -eq 0 ]; then
    if [ $# -lt 2 ]; then
        echo "usage: equiv_test.sh [--original-bugs] <ref_bin_dir> <seq_list> <msp_file> [num_sections]" >&2
        exit 2
    fi
    REF_BIN_DIR="$(realpath "$1")"; shift
else
    if [ $# -lt 2 ]; then
        echo "usage: equiv_test.sh --self-compare <seq_list> <msp_file> [num_sections]" >&2
        exit 2
    fi
fi

SEQ_LIST="$(realpath "$1")"
MSP_FILE="$(realpath "$2")"
NSECT="${3:-1}"

# ---- Compute test-set signature and paths ----
# Use first 12 hex chars of md5(msps.out) as a stable, path-independent key.
# The human-readable label comes from the directory that originally held msps.out.
if [ $SELF_COMPARE -eq 0 ]; then
    SIG=$(md5sum "$MSP_FILE" | cut -c1-12)
    LABEL=$(basename "$(dirname "$MSP_FILE")")
    TEST_SET_DIR="$TEST_SETS_DIR/$SIG"
    MANIFEST="$TEST_SET_DIR/manifest"
    REF_CACHE="$TEST_SET_DIR/reference"
fi

# ---- Temp directory ----
TMPROOT="$(mktemp -d)"
trap 'rm -rf "$TMPROOT"' EXIT

# ---- Determine TUNABLES for the test build ----
if [ $USE_ORIGINAL_BUGS -eq 1 ]; then
    TEST_TUNABLES="-DORIGINAL_BUGS"
    TEST_LABEL="refactored (ORIGINAL_BUGS)"
else
    TEST_TUNABLES=""
    TEST_LABEL="fixed"
fi

# ---- Build test binaries ----
TEST_BIN_DIR="$TMPROOT/bin_test"
mkdir "$TEST_BIN_DIR"

echo "==> Building ${TEST_LABEL} binaries from $PROJ_DIR ..."
(cd "$PROJ_DIR" \
    && make clean > /dev/null 2>&1 \
    && make TUNABLES="$TEST_TUNABLES" INSTALL_DIR="$TEST_BIN_DIR" install 2>&1) \
    || { echo "error: test build failed" >&2; exit 2; }

# ---- In self-compare mode, also build the ORIGINAL_BUGS reference ----
if [ $SELF_COMPARE -eq 1 ]; then
    REF_BIN_DIR="$TMPROOT/bin_bugs"
    mkdir "$REF_BIN_DIR"
    echo "==> Building ORIGINAL_BUGS (reference) binaries from $PROJ_DIR ..."
    (cd "$PROJ_DIR" \
        && make clean > /dev/null 2>&1 \
        && make TUNABLES="-DORIGINAL_BUGS" INSTALL_DIR="$REF_BIN_DIR" install 2>&1) \
        || { echo "error: ORIGINAL_BUGS build failed" >&2; exit 2; }
    REF_LABEL="original-bugs"
else
    REF_LABEL="reference"
fi

# ---- Populate test-set directory (first visit) ----
if [ $SELF_COMPARE -eq 0 ]; then
    if [ ! -d "$TEST_SET_DIR" ]; then
        echo "==> Creating test set [$SIG] for \"$LABEL\" ..."
        mkdir -p "$TEST_SET_DIR"
        cp "$SEQ_LIST" "$TEST_SET_DIR/seqnames"
        cp "$MSP_FILE"  "$TEST_SET_DIR/msps.out"
        printf "label=%s\nmsps_source=%s\nseqnames_source=%s\nref_bin_dir=%s\ncreated=%s\n" \
            "$LABEL" "$MSP_FILE" "$SEQ_LIST" "$REF_BIN_DIR" "$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
            > "$MANIFEST"
    fi

    # ---- Check whether cached reference is still valid ----
    CACHED_REF_BIN=""
    [ -f "$MANIFEST" ] && CACHED_REF_BIN=$(grep "^ref_bin_dir=" "$MANIFEST" | cut -d= -f2-)

    USE_CACHE=0
    if [ -f "$REF_CACHE/summary/eles" ]; then
        if [ "$CACHED_REF_BIN" = "$REF_BIN_DIR" ]; then
            USE_CACHE=1
        else
            echo "==> Cached reference used a different ref_bin_dir; regenerating ..."
            echo "    cached: $CACHED_REF_BIN"
            echo "    current: $REF_BIN_DIR"
            rm -rf "$REF_CACHE"
        fi
    fi

    if [ $USE_CACHE -eq 1 ]; then
        echo "==> Using cached reference output [$SIG] \"$LABEL\" ..."
        RUN_REF="$REF_CACHE"
    else
        echo "==> Running reference pipeline (caching to test-sets/$SIG/) ..."
        FRESH_REF="$TMPROOT/run_ref"
        "$SCRIPT_DIR/run_recon.sh" \
            "$REF_BIN_DIR" "$TEST_SET_DIR/seqnames" "$TEST_SET_DIR/msps.out" \
            "$NSECT" "$FRESH_REF" \
            || { echo "error: reference pipeline failed" >&2; exit 2; }

        # Cache only the files that compare_recon_runs.sh checks
        mkdir -p "$REF_CACHE/summary" "$REF_CACHE/edge_redef_res"
        for f in eles families fam_no ele_no naive_eles redef_ele_no ori_msp_no; do
            src="$FRESH_REF/summary/$f"
            [ -f "$src" ] && cp "$src" "$REF_CACHE/summary/$f"
        done
        src="$FRESH_REF/edge_redef_res/redef_stat"
        [ -f "$src" ] && cp "$src" "$REF_CACHE/edge_redef_res/redef_stat"

        # Update manifest with the ref_bin_dir that generated this cache
        sed -i "s|^ref_bin_dir=.*|ref_bin_dir=$REF_BIN_DIR|" "$MANIFEST"

        RUN_REF="$REF_CACHE"
    fi
else
    # self-compare: run reference pipeline in temp (no caching)
    RUN_REF="$TMPROOT/run_ref"
    echo "==> Running ${REF_LABEL} pipeline ..."
    "$SCRIPT_DIR/run_recon.sh" \
        "$REF_BIN_DIR" "$SEQ_LIST" "$MSP_FILE" "$NSECT" "$RUN_REF" \
        || { echo "error: reference pipeline failed" >&2; exit 2; }
fi

# ---- Run test pipeline ----
# In self-compare mode use the original paths; otherwise use the cached copies.
RUN_TEST="$TMPROOT/run_test"
if [ $SELF_COMPARE -eq 1 ]; then
    TEST_SEQ="$SEQ_LIST"
    TEST_MSP="$MSP_FILE"
else
    TEST_SEQ="$TEST_SET_DIR/seqnames"
    TEST_MSP="$TEST_SET_DIR/msps.out"
fi
echo "==> Running ${TEST_LABEL} pipeline ..."
"$SCRIPT_DIR/run_recon.sh" \
    "$TEST_BIN_DIR" "$TEST_SEQ" "$TEST_MSP" \
    "$NSECT" "$RUN_TEST" \
    || { echo "error: test pipeline failed" >&2; exit 2; }

# ---- Compare outputs ----
echo ""
"$SCRIPT_DIR/compare_recon_runs.sh" "$RUN_REF" "$RUN_TEST" "$REF_LABEL" "$TEST_LABEL"
STATUS=$?

if [ $STATUS -eq 0 ]; then
    echo ""
    echo "Equivalence test PASSED."
else
    echo ""
    echo "Equivalence test FAILED.  See diffs above."
    cp -r "$RUN_TEST" ./equiv_test_result
    echo "Test output preserved in: equiv_test_result/"
    echo "Reference output is at:   $RUN_REF"
fi

exit $STATUS
