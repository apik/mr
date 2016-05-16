AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

LDADD = ../mr/libmr.la
EXCPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.3 -I$(top_srcdir)/mr $(OPENMP_CFLAGS)

EXLDDFLAGS = $(OPENMP_CFLAGS) $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.3 -ltsil -lm

noinst_PROGRAMS = ex1Xij ex2Pole2ms ex3Pole2msRun 

# CPC examples
ex1Xij_SOURCES =  ex1Xij.cpp
ex1Xij_CPPFLAGS = $(EXCPPFLAGS)
ex1Xij_LDADD = $(EXLDDFLAGS)

ex2Pole2ms_SOURCES =  ex2Pole2ms.cpp
ex2Pole2ms_CPPFLAGS = $(EXCPPFLAGS)
ex2Pole2ms_LDADD = $(EXLDDFLAGS)

ex3Pole2msRun_SOURCES =  ex3Pole2msRun.cpp
ex3Pole2msRun_CPPFLAGS = $(EXCPPFLAGS)
ex3Pole2msRun_LDADD = $(EXLDDFLAGS)

EXTRA_DIST = README
