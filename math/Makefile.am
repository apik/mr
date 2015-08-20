AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

MATHROOT=$(dir $(MPREP))


EX_CFLAGS = $(MR_CFLAGS) -I$(MATHROOT) -I$(top_srcdir)/tsil-1.21 -I$(top_srcdir)/mr -I$(top_srcdir)/Eigen

noinst_PROGRAMS = mr

mr_SOURCES =  mrimpl.cpp

nodist_mr_SOURCES = mr.cpp
BUILT_SOURCES = mr.cpp

mr_CPPFLAGS = $(EX_CFLAGS)
mr_LDADD = $(MR_LIBS) -L$(MATHROOT)$(ALTMATHPATH) ${LDMATH} -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

mr.cpp: mr.tm
	$(MPREP) $? -o $@

EXTRA_DIST = mr.tm mr.m
