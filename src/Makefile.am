#===============================================================================
# Makefile.am : for automake
#===============================================================================

bin_PROGRAMS = prepromali

AM_CPPFLAGS = $(INTI_CFLAGS)
AM_CFLAGS = -Wall

prepromali_SOURCES = \
	arg.c arg.h \
	getseqs.c getseqs.h \
	prepromali.c prepromali.h \
	putseqs.c putseqs.h \
	safe.c safe.h \
	seq.h

EXTRA_DIST = doxygen.cfg

CLEANFILES = $(TARGET) *.o 
DISTCLEANFILES = libtool config.cache config.log

