#!/bin/sh
#
# run_recon.sh  --  Run the full RECON pipeline in a working directory
#
# Usage:
#   run_recon.sh <bin_dir> <seq_list> <msp_file> [num_sections] [work_dir]
#
#   bin_dir      directory containing the RECON binaries
#   seq_list     sequence name list file
#   msp_file     MSP input file
#   num_sections passed to imagespread when using old (1.08) reference binaries;
#               ignored when running the new eledef (imagespread merged in)
#   work_dir     optional working directory; defaults to a temp dir
#
# Pipeline stages:
#   1. eledef    -- read MSPs, sort images in-memory, cluster into elements,
#                   write element database to ele_store/
#   2. eleredef  -- element redefinition / dissection
#   3. edgeredef -- edge filtering / PPS repair
#   4. famdef    -- BFS family assignment; compacts element database
#

set -e

if [ $# -lt 3 ]; then
    echo "usage: run_recon.sh bin_dir seq_list msp_file [num_sections] [work_dir]" >&2
    exit 1
fi

BIN_DIR="$1"
SEQ_LIST="$(realpath "$2")"
MSP_FILE="$(realpath "$3")"
# $4 (num_sections): used by old-pipeline imagespread; ignored by new eledef
# $5 (work_dir) or mktemp

# Extract num_sections: default 1; only set if $4 is a positive integer
if [ -n "${4:-}" ] && expr "${4:-}" : '^[0-9][0-9]*$' > /dev/null 2>&1; then
    NSECT="${4}"
    if [ -n "${5:-}" ]; then
        WORK_DIR="$(realpath "$5")"
        mkdir -p "$WORK_DIR"
    else
        WORK_DIR="$(mktemp -d)"
    fi
elif [ -n "${4:-}" ]; then
    # $4 is not a number -- treat it as work_dir
    NSECT=1
    WORK_DIR="$(realpath "$4")"
    mkdir -p "$WORK_DIR"
else
    NSECT=1
    WORK_DIR="$(mktemp -d)"
fi

# Verify binaries
for prog in eledef eleredef edgeredef famdef; do
    if [ ! -x "$BIN_DIR/$prog" ]; then
        echo "error: $BIN_DIR/$prog not found or not executable" >&2
        exit 1
    fi
done

cd "$WORK_DIR"

LOG="$WORK_DIR/run.log"

# ---- Shared output directories ----
rm -rf summary ele_def_res ele_store ele_redef_res edge_redef_res
mkdir  summary ele_def_res ele_store ele_redef_res edge_redef_res

# ---- Stage 1: eledef ----
#
# Two pipeline variants are supported:
#
#   New (this repo, 1.09+):  eledef reads MSPs directly; imagespread is
#     merged into eledef.  No images/ directory is needed.
#
#   Old (released 1.08):  imagespread must run first to partition MSPs into
#     sorted spread files, which are then merged into images/images_sorted.
#     eledef then reads images/images_sorted.  Detected by the presence of
#     an imagespread binary alongside eledef.
#
if [ -x "$BIN_DIR/imagespread" ]; then
    # Old pipeline: run imagespread + sort before eledef
    mkdir -p images
    "$BIN_DIR/imagespread" "$SEQ_LIST" "$MSP_FILE" "$NSECT" >> "$LOG" 2>&1
    i=1
    while [ "$i" -le "$NSECT" ]; do
        spread="images/spread$i"
        [ -f "$spread" ] && sort -k 3,3 -k 4n,4n -k 5nr,5nr "$spread" >> images/images_sorted
        i=$((i + 1))
    done
    rm -f images/spread*
fi

"$BIN_DIR/eledef" "$SEQ_LIST" "$MSP_FILE" single >> "$LOG" 2>&1

# ---- Stage 2: eleredef ----
# tmp  -> ele_def_res  (read: size_list, element data via ele_db)
# tmp2 -> ele_redef_res (write: redef_stat, log)
rm -f tmp tmp2
ln -s ele_def_res  tmp
ln -s ele_redef_res tmp2

"$BIN_DIR/eleredef" "$SEQ_LIST" >> "$LOG" 2>&1

rm -f tmp tmp2

# ---- Stage 3: edgeredef ----
# tmp  -> ele_redef_res  (read: redef_stat, element data via ele_db)
# tmp2 -> edge_redef_res (write: redef_stat, log)
rm -f tmp tmp2
ln -s ele_redef_res  tmp
ln -s edge_redef_res tmp2

"$BIN_DIR/edgeredef" "$SEQ_LIST" >> "$LOG" 2>&1

rm -f tmp tmp2

# ---- Stage 4: famdef (compacts element database at exit) ----
# tmp -> edge_redef_res (read: redef_stat, log2)
rm -f tmp
ln -s edge_redef_res tmp

"$BIN_DIR/famdef" "$SEQ_LIST" >> "$LOG" 2>&1

rm -f tmp

echo "$WORK_DIR"
