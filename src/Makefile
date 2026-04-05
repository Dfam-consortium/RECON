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

PROGS = imagespread eledef eleredef edgeredef famdef

HDRS = bolts.h seqlist.h msps.h ele.h treeview.h eleredef.h \
       recon_defs.h recon_log.h

SRC = imagespread.c eledef.c eleredef.c edgeredef.c famdef.c

DISTFILES = $(SRC) $(HDRS) Makefile

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


## Explicit build rules with full header dependencies.
## All programs are single translation units; each .c includes its headers
## directly, so the full HDRS list is listed as a dependency for each.

imagespread: imagespread.c $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o imagespread imagespread.c $(LIBS)

eledef: eledef.c $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o eledef eledef.c $(LIBS)

eleredef: eleredef.c $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o eleredef eleredef.c $(LIBS)

edgeredef: edgeredef.c $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o edgeredef edgeredef.c $(LIBS)

famdef: famdef.c $(HDRS)
	$(CC) $(CFLAGS) $(MDEFS) -o famdef famdef.c $(LIBS)
