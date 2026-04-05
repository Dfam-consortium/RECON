# Makefile for RECON
#
# Build targets:  all, install, clean, test
#
VERSION = 1.09

## where you want things installed
BINDIR  = ../bin
MANDIR  = $(HOME)/man

## your compiler
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
##                            See recon_defs.h for details.
##                            Example:  make TUNABLES="-DORIGINAL_BUGS"
TUNABLES =

CFLAGS = -O -I ../minunit $(TUNABLES)

## your compression program for dist target
COMPRESS = gzip

#######
## should not need to modify below this line
#######
SHELL = /bin/sh
LIBS  = -lm

PROGS = eledef eleredef edgeredef famdef

HDRS = bolts.h seqlist.h msps.h ele.h ele_db.h treeview.h eleredef.h \
       recon_defs.h recon_log.h

SRC = imagespread.c eledef.c eleredef.c edgeredef.c famdef.c \
      seqlist.c msps.c ele.c ele_db.c treeview.c

DISTFILES = $(SRC) $(HDRS) Makefile

## Shared object files compiled from the library .c files.
## OBJS_DB   -- flat-file element database (needed by all four pipeline programs)
## OBJS_SEQ  -- sequence-name table (needed by all four programs)
## OBJS_MSP  -- MSP/image/frag operations (needed by all four programs)
## OBJS_ELE  -- element/edge/family operations (needed by eleredef, edgeredef, famdef)
## OBJS_TREE -- ASCII tree printer (needed only by eleredef)
OBJS_DB   = ele_db.o
OBJS_SEQ  = seqlist.o
OBJS_MSP  = msps.o
OBJS_ELE  = ele.o
OBJS_TREE = treeview.o

OBJS_COMMON = $(OBJS_SEQ) $(OBJS_MSP)
OBJS_FULL   = $(OBJS_COMMON) $(OBJS_ELE) $(OBJS_DB)

all: $(PROGS)

install: $(PROGS)
	cp $(PROGS) $(BINDIR)/

clean:
	-rm -f $(PROGS) *.o *~ *.Addrs *.Counts *.pixie Makefile.bak core TAGS

tags:
	etags -t $(DISTFILES)

dist:
	@if test -d recon-$(VERSION);        then rm -rf recon-$(VERSION);    fi
	@if test -f recon-$(VERSION).tar.gz; then rm recon-$(VERSION).tar.gz; fi
	mkdir recon-$(VERSION)/
	cp $(DISTFILES) recon-$(VERSION)/
	tar cvf recon-$(VERSION).tar recon-$(VERSION)
	$(COMPRESS) recon-$(VERSION).tar

test:
	(cd Testsuite; testsuite.sh)


## Shared library object rules

seqlist.o: seqlist.c seqlist.h bolts.h recon_defs.h
	$(CC) $(CFLAGS) $(MDEFS) -c seqlist.c

msps.o: msps.c msps.h seqlist.h bolts.h recon_defs.h
	$(CC) $(CFLAGS) $(MDEFS) -c msps.c

ele.o: ele.c ele.h msps.h seqlist.h bolts.h recon_defs.h recon_log.h ele_db.h
	$(CC) $(CFLAGS) $(MDEFS) -c ele.c

treeview.o: treeview.c treeview.h msps.h eleredef.h
	$(CC) $(CFLAGS) $(MDEFS) -c treeview.c

ele_db.o: ele_db.c ele_db.h
	$(CC) $(CFLAGS) $(MDEFS) -c ele_db.c


## Executable link rules.
## Each program is compiled from its own .c file and linked against the
## shared object files that provide the library functions it uses.

eledef: eledef.c $(OBJS_COMMON) $(OBJS_DB) $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o eledef eledef.c $(OBJS_COMMON) $(OBJS_DB) $(LIBS)

eleredef: eleredef.c $(OBJS_FULL) $(OBJS_TREE) $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o eleredef eleredef.c $(OBJS_FULL) $(OBJS_TREE) $(LIBS)

edgeredef: edgeredef.c $(OBJS_FULL) $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o edgeredef edgeredef.c $(OBJS_FULL) $(LIBS)

famdef: famdef.c $(OBJS_FULL) $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o famdef famdef.c $(OBJS_FULL) $(LIBS)
