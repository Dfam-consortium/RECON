# Makefile for RECON
#
# Build targets:  all, install, clean
# Source files live in src/; object files are placed in src/; binaries go to bin/.
#
VERSION = 1.09

## Where files are installed (default: bin/ at the project root)
## Binaries, run_recon.sh, README.md, VALIDATION.md, and LICENSE
## are all copied here by "make install".
INSTALL_DIR = bin
MANDIR      = $(HOME)/man

## Compiler
CC = gcc

## Compiler flags
## -DASSERT        Diagnostics level 1. Assertions.
## -DVERBOSE       Diagnostics level 2. Verbosity.  (unused at the moment)
## -DMEMDEBUG      Diagnostics level 3. Memory checking. (maybe slow!)
##                                      (must link w/ -ldbmalloc)
## -DLINTING       Suppress overzealous warnings from gcc -Wall and lint
#MDEFS = -DLINTING -DASSERT -DMEMDEBUG
#MDEFS = -DASSERT

## Tunable constants -- override on the command line to change pipeline behaviour
## without editing source.  All values here are the original defaults.
## Example:  make TUNABLES="-DELEREDEF_MAX_IMAGES=2400 -DBFS_CLAN_DEPTH=4"
##
## Constants and their recon_defs.h names (original name in parentheses):
##   SEQ_NAME_MAX_LEN      -- max sequence-name length      (NAME_LEN=50)
##   BFS_CLAN_DEPTH        -- BFS depth for local networks  (DEPTH=3)
##   TANDEM_SIZE_LIMIT     -- tandem-filter image threshold  (SIZE_LIMIT=500)
##   TANDEM_IMG_RATIO      -- tandem-filter ratio threshold  (TANDEM=5)
##   ELEDEF_IMGPROT_CAP    -- eledef image batch capacity    (IMG_CAP=500000)
##   ELEREDEF_MAX_IMAGES   -- eleredef image array capacity  (MAX_IMG=1200)
##   MIN_ELEMENT_LEN_BP    -- minimum element length in bp   (TOO_SHORT=30)
##   ELEREDEF_CUTOFF_SINGLE -- eleredef single-cov cutoff   (CUTOFF1=0.5)
##   ELEREDEF_CUTOFF_DOUBLE -- eleredef double-cov cutoff   (CUTOFF2=0.9)
##   EDGEREDEF_EDGE_CUTOFF  -- edge PPS filter ratio        (CUTOFF3=0.7)
##
## Equivalence / regression testing flag:
##   ORIGINAL_BUGS         -- restore two pre-existing bugs for output
##                            comparison against the released binaries.
##                            See src/recon_defs.h for details.
##                            Example:  make TUNABLES="-DORIGINAL_BUGS"
TUNABLES =

## -I src     finds all project headers (ele.h, msps.h, ...)
## -I minunit finds the minunit test framework header
CFLAGS = -O -I src -I minunit -DRECON_VERSION=\"$(VERSION)\" $(TUNABLES)

## Compression program for dist target
COMPRESS = gzip

#######
## should not need to modify below this line
#######
SHELL = /bin/sh
LIBS  = -lm

PROGS = $(INSTALL_DIR)/eledef $(INSTALL_DIR)/eleredef $(INSTALL_DIR)/edgeredef $(INSTALL_DIR)/famdef

HDRS = src/bolts.h src/seqlist.h src/msps.h src/ele.h src/ele_db.h \
       src/treeview.h src/eleredef.h src/recon_defs.h src/recon_log.h \
       src/redef_boundary.h src/redef_edges.h src/redef_dissect.h

SRC = src/eledef.c src/eleredef.c src/edgeredef.c src/famdef.c \
      src/seqlist.c src/msps.c src/ele.c src/ele_bst.c src/ele_print.c \
      src/ele_db.c src/treeview.c \
      src/redef_boundary.c src/redef_edges.c src/redef_dissect.c

DISTFILES = $(SRC) $(HDRS) Makefile

## Shared object files compiled from the library .c files.
## OBJS_DB   -- flat-file element database (needed by all four pipeline programs)
## OBJS_SEQ  -- sequence-name table (needed by all four programs)
## OBJS_MSP  -- MSP/image/frag operations (needed by all four programs)
## OBJS_ELE  -- element/edge/family operations (needed by eleredef, edgeredef, famdef)
## OBJS_TREE -- ASCII tree printer (needed only by eleredef)
## OBJS_REDEF -- split eleredef helper modules
OBJS_DB   = src/ele_db.o
OBJS_SEQ  = src/seqlist.o
OBJS_MSP  = src/msps.o
OBJS_ELE  = src/ele.o src/ele_bst.o src/ele_print.o
OBJS_TREE = src/treeview.o
OBJS_REDEF = src/redef_boundary.o src/redef_edges.o src/redef_dissect.o

OBJS_COMMON = $(OBJS_SEQ) $(OBJS_MSP)
OBJS_FULL   = $(OBJS_COMMON) $(OBJS_ELE) $(OBJS_DB)

all: $(INSTALL_DIR) $(PROGS)

$(INSTALL_DIR):
	mkdir -p $(INSTALL_DIR)

install: all
	cp scripts/run_recon.sh $(INSTALL_DIR)/run_recon.sh
	chmod +x $(INSTALL_DIR)/run_recon.sh
	cp scripts/imagespread $(INSTALL_DIR)/imagespread
	chmod +x $(INSTALL_DIR)/imagespread
	cp README.md VALIDATION.md LICENSE $(INSTALL_DIR)/
	printf '# Support for RepeatModeler pre-2.0.8\nVersion $(VERSION) (%s)\n# See README.md for program documentation.\n' \
	  "$$(date '+%b %Y')" > $(INSTALL_DIR)/00README

PROG_NAMES = eledef eleredef edgeredef famdef

clean:
	-rm -f $(PROGS) $(addprefix src/, $(PROG_NAMES)) src/*_dbg src/*.o src/*~ *~ core
	-rm -f *.Addrs *.Counts *.pixie Makefile.bak TAGS
	-rm -f $(INSTALL_DIR)/run_recon.sh $(INSTALL_DIR)/imagespread \
	        $(INSTALL_DIR)/README.md $(INSTALL_DIR)/VALIDATION.md \
	        $(INSTALL_DIR)/LICENSE $(INSTALL_DIR)/00README

tags:
	etags -t $(DISTFILES)

dist:
	@if test -d recon-$(VERSION);        then rm -rf recon-$(VERSION);    fi
	@if test -f recon-$(VERSION).tar.gz; then rm recon-$(VERSION).tar.gz; fi
	mkdir recon-$(VERSION)/
	cp $(DISTFILES) recon-$(VERSION)/
	tar cvf recon-$(VERSION).tar recon-$(VERSION)
	$(COMPRESS) recon-$(VERSION).tar


## Shared library object rules

src/seqlist.o: src/seqlist.c src/seqlist.h src/bolts.h src/recon_defs.h
	$(CC) $(CFLAGS) $(MDEFS) -c src/seqlist.c -o $@

src/msps.o: src/msps.c src/msps.h src/seqlist.h src/bolts.h src/recon_defs.h
	$(CC) $(CFLAGS) $(MDEFS) -c src/msps.c -o $@

src/ele.o: src/ele.c src/ele.h src/msps.h src/seqlist.h src/bolts.h src/recon_defs.h src/recon_log.h src/ele_db.h
	$(CC) $(CFLAGS) $(MDEFS) -c src/ele.c -o $@

src/ele_bst.o: src/ele_bst.c src/ele.h src/msps.h src/seqlist.h src/bolts.h src/recon_defs.h src/recon_log.h src/ele_db.h
	$(CC) $(CFLAGS) $(MDEFS) -c src/ele_bst.c -o $@

src/ele_print.o: src/ele_print.c src/ele.h src/msps.h src/seqlist.h src/bolts.h src/recon_defs.h src/recon_log.h src/ele_db.h
	$(CC) $(CFLAGS) $(MDEFS) -c src/ele_print.c -o $@

src/treeview.o: src/treeview.c src/treeview.h src/msps.h src/eleredef.h
	$(CC) $(CFLAGS) $(MDEFS) -c src/treeview.c -o $@

src/ele_db.o: src/ele_db.c src/ele_db.h
	$(CC) $(CFLAGS) $(MDEFS) -c src/ele_db.c -o $@

src/redef_boundary.o: src/redef_boundary.c src/redef_boundary.h src/ele.h src/msps.h src/seqlist.h src/bolts.h src/recon_defs.h src/recon_log.h src/ele_db.h
	$(CC) $(CFLAGS) $(MDEFS) -c src/redef_boundary.c -o $@

src/redef_edges.o: src/redef_edges.c src/redef_edges.h src/redef_boundary.h src/eleredef.h src/ele.h src/msps.h src/seqlist.h src/bolts.h src/recon_defs.h src/recon_log.h src/ele_db.h
	$(CC) $(CFLAGS) $(MDEFS) -c src/redef_edges.c -o $@

src/redef_dissect.o: src/redef_dissect.c src/redef_dissect.h src/redef_boundary.h src/redef_edges.h src/ele.h src/msps.h src/seqlist.h src/bolts.h src/recon_defs.h src/recon_log.h src/ele_db.h
	$(CC) $(CFLAGS) $(MDEFS) -c src/redef_dissect.c -o $@


## Executable link rules.
## Each program is compiled from its own .c file and linked against the
## shared object files that provide the library functions it uses.

$(INSTALL_DIR)/eledef: src/eledef.c $(OBJS_COMMON) $(OBJS_DB) $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o $@ src/eledef.c $(OBJS_COMMON) $(OBJS_DB) $(LIBS)

$(INSTALL_DIR)/eleredef: src/eleredef.c $(OBJS_REDEF) $(OBJS_FULL) $(OBJS_TREE) $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o $@ src/eleredef.c $(OBJS_REDEF) $(OBJS_FULL) $(OBJS_TREE) $(LIBS)

$(INSTALL_DIR)/edgeredef: src/edgeredef.c $(OBJS_FULL) $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o $@ src/edgeredef.c $(OBJS_FULL) $(LIBS)

$(INSTALL_DIR)/famdef: src/famdef.c $(OBJS_FULL) $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o $@ src/famdef.c $(OBJS_FULL) $(LIBS)
