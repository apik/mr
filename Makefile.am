##
## Makefile.am for SM3r2m
##
## written by Andrey Pikelner
##
## Process this file with automake to produce Makefile.in
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4
SUBDIRS =         \
	tsil-1.3 \
	mr        \
	examples  \
	tests


if HAVE_MPREP
SUBDIRS+=math
endif

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = mr.pc

EXTRA_DIST = boost mr.pc.in README.md
DISTCLEANFILES = mr.pc

$(pkgconfig_DATA): config.status

