PROGS=capmq

all: $(PROGS)

CC=gcc

HTSDIR=../htslib

INCLUDES=-I$(HTSDIR)/include -I$(HTSDIR)
LIBS=-L$(HTSDIR)/lib -L$(HTSDIR) -lhts -Wl,--rpath,$(HTSDIR)/lib -Wl,--rpath,$(HTSDIR) -lpthread -lz -lm -ldl
CFLAGS=-O3 -g -Wall -Werror

PACKAGE_VERSION = 0.3

# If building from a Git repository, replace $(PACKAGE_VERSION) with the Git
# description of the working tree: either a release tag with the same value
# as $(PACKAGE_VERSION) above, or an exact description likely based on a tag.
# $(shell), :=, etc are GNU Make-specific.  If you don't have GNU Make,
# comment out this conditional.
ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)

# Force version.h to be remade if $(PACKAGE_VERSION) has changed.
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif

# If you don't have GNU Make but are building from a Git repository, you may
# wish to replace this with a rule that always rebuilds version.h:
# version.h: force
#	echo '#define CAPMQ_VERSION "`git describe --always --dirty`"' > $@
version.h:
	echo '#define CAPMQ_VERSION "$(PACKAGE_VERSION)"' > $@

print-version:
	@echo $(PACKAGE_VERSION)

check test: capmq t_capmq
	./t_capmq

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

clean:
	-rm -f *.o $(PROGS)

capmq: $(HTSLIB) version.h capmq.o
	$(CC) -o $@ capmq.o $(CFLAGS) $(LDFLAGS) $(LIBS)

t_capmq: t_capmq.o

force:

.PHONY: force
