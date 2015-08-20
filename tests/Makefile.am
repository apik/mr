AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

LDADD = ../mr/libmr.a
EX_CFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.3 -I$(top_srcdir)/mr

noinst_PROGRAMS = testsuite

TEST_FILES = tst-BKKS_1205.2893.cpp  tst-JKV_0105304.cpp  tst-JKV_0212319.cpp  tst-mPlanck.cpp  tst-MWeqMZ.cpp  tst-RunDec.cpp

if HAVE_EIGEN
TEST_FILES += tst-mOSmMS_QCD.cpp
endif

testsuite_SOURCES =  testsuite.cpp  catch.hpp $(TEST_FILES)
testsuite_CPPFLAGS = $(EX_CFLAGS)
testsuite_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.3 -ltsil -lm

