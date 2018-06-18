# Directories.
export BINDIR    = $(CURDIR)/bin
export INCDIR    = $(CURDIR)/include
export LIBDIR    = $(CURDIR)/lib
export LIBSRCDIR = $(CURDIR)/libsrc
export SRCDIR    = $(CURDIR)/src

# Tool chain.
export CC = gcc

export CFLAGS = -Wall -Wextra -O3
export CFLAGS += -I $(INCDIR)

# Libraries.
export LIB = libcapb
export LIBS = $(LIBDIR)/$(LIB).a -lm

all: lib
	mkdir -p bin
	cd $(SRCDIR) && $(MAKE) all

lib:
	mkdir -p lib
	cd $(LIBSRCDIR) && $(MAKE) all

clean:
	cd $(LIBSRCDIR) && $(MAKE) clean
	cd $(SRCDIR) && $(MAKE) clean
