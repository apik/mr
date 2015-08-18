AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

LDADD = ../mr/libmr.la
EX_CFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.21 -I$(top_srcdir)/mr -I$(top_srcdir)/Eigen

noinst_PROGRAMS = ex1Xij ex2Pole2ms ex3Pole2msRun 

# CPC examples
ex1Xij_SOURCES =  ex1Xij.cpp
ex1Xij_CPPFLAGS = $(EX_CFLAGS)
ex1Xij_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

ex2Pole2ms_SOURCES =  ex2Pole2ms.cpp
ex2Pole2ms_CPPFLAGS = $(EX_CFLAGS)
ex2Pole2ms_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

ex3Pole2msRun_SOURCES =  ex3Pole2msRun.cpp
ex3Pole2msRun_CPPFLAGS = $(EX_CFLAGS)
ex3Pole2msRun_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

