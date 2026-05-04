# Changelog

## v1.10
- Added a guardrail for the last-case merging functionality so that it
  doesn't attempt to force together families when there is at least one  
  backpointing primary edge.

## v1.09 (not officially released)

- Fixed bug which prohibited last-case merging of related families dating
  back to the 1.01 (original release -- possibly earlier) version.
  NOTE: This is a breaking change which may be reverted with the use of
  the `ORIGINAL_BUGS` compilation flag.
- Merged `imagespread` binary into `eledef`; eliminated `images/` directory
- Replaced per-element `tmp/e<N>` files with flat-file element database
  (`ele_store/elements.db` + `ele_store/elements.idx`)
  This reduces the huge number of small files and allows for more efficient
  random access.
- Split 2,500-line `eleredef.c` into four focused modules
  (`redef_boundary.c`, `redef_edges.c`, `redef_dissect.c`)
- Moved all tunable constants to `recon_defs.h`; added structured logging
- Fixed O(N^2) `report_redef_stat()` bottleneck in `edgeredef`; removed
  unnecessary `fflush()` calls from inner loops (~12× speedup on large runs)
- Added `scripts/` equivalence-testing infrastructure
- Added `ORIGINAL_BUGS` compilation flag to enable bit-exact comparison
  against pre RECON 1.09 versions
- Added backwards compatibility using a dummy imagespread script

## v1.08 (2014-02-05)

- Memory allocation and file operation error checks in `eleredef.c`
- Buffer overrun fix in `eledef.c` (Stephen Ficklin)

## v1.07 (2011-06-08)

- Workaround for divide-by-zero in `find_prim()` (`eleredef.c`)

## v1.06 (2008-05-20)

- 64-bit portability: replaced `long` with `int32_t` throughout
