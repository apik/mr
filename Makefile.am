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

EXTRA_DIST = Eigen boost


