AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

LDADD = ../mr/libmr.la
EX_CFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.21 -I$(top_srcdir)/mr -I$(top_srcdir)/Eigen

noinst_PROGRAMS = exPole2ms exPole2run exMHcrit exMtcrit 

if HAVE_LIBGSL
noinst_PROGRAMS += gslcrit1d
endif 

# CPC examples
exPole2ms_SOURCES =  exPole2ms.cpp
exPole2ms_CPPFLAGS = $(EX_CFLAGS)
exPole2ms_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

exPole2run_SOURCES =  exPole2run.cpp
exPole2run_CPPFLAGS = $(EX_CFLAGS)
exPole2run_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

exMHcrit_SOURCES =  exMHcrit.cpp
exMHcrit_CPPFLAGS = $(EX_CFLAGS)
exMHcrit_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

exMtcrit_SOURCES =  exMtcrit.cpp
exMtcrit_CPPFLAGS = $(EX_CFLAGS)
exMtcrit_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

if HAVE_LIBGSL
gslcrit1d_SOURCES =  gslcrit1d.cpp
gslcrit1d_CPPFLAGS = $(EX_CFLAGS)
gslcrit1d_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm -lgsl -lgslcblas
endif
