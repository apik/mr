AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

LDADD = ../mr/libmr.a
EX_CFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.21 -I$(top_srcdir)/mr -I$(top_srcdir)/Eigen

noinst_PROGRAMS = testsuite

testsuite_SOURCES =  testsuite.cpp BKKS_1205.2893.cpp JKV_0105304.cpp JKV_0212319.cpp MWeqMZ.cpp mPlanck.cpp RunDec.cpp mOSmMS_QCD.cpp
testsuite_CPPFLAGS = $(EX_CFLAGS)
testsuite_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

