AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

LDADD = ../mr/libmr.la

noinst_PROGRAMS = mW mZ mH mt runAlpha runSM mratio mb yukawa
if HAVE_ROOT
noinst_PROGRAMS += mHroot
endif 

mW_SOURCES =  mW.cpp
mW_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mW_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm

mZ_SOURCES =  mZ.cpp
mZ_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mZ_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm 

mH_SOURCES =  mH.cpp
mH_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mH_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm


if HAVE_ROOT
mHroot_SOURCES =  mHroot.cpp
mHroot_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr -I$(ROOTINCDIR)
mHroot_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm -L$(ROOTLIBDIR) -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread
endif

mt_SOURCES =  mt.cpp
mt_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mt_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm

yukawa_SOURCES =  yukawa.cpp
yukawa_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
yukawa_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm

mb_SOURCES =  mb.cpp
mb_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mb_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm

runAlpha_SOURCES =  runAlpha.cpp CRunDec.cpp CRunDec.h
runAlpha_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
runAlpha_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm

runSM_SOURCES =  runSM.cpp
runSM_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
runSM_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm

mratio_SOURCES =  mratio.cpp
mratio_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mratio_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm
