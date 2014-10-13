AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

LDADD = ../mr/libmr.a
EX_CFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.21 -I$(top_srcdir)/mr

noinst_PROGRAMS = mW mZ mH mt runAlpha runSM mb yukawa lambda yt gauge plotyukawa mtgl bezrukovFit
if HAVE_ROOT
noinst_PROGRAMS += mHroot
endif 

mW_SOURCES =  mW.cpp
mW_CPPFLAGS = $(EX_CFLAGS)
mW_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

mZ_SOURCES =  mZ.cpp
mZ_CPPFLAGS = $(EX_CFLAGS)
mZ_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm 

mH_SOURCES =  mH.cpp
mH_CPPFLAGS = $(EX_CFLAGS)
mH_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm


if HAVE_ROOT
mHroot_SOURCES =  mHroot.cpp
mHroot_CPPFLAGS = $(EX_CFLAGS)
mHroot_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm -L$(ROOTLIBDIR) -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread
endif

mt_SOURCES =  mt.cpp
mt_CPPFLAGS = $(EX_CFLAGS)
mt_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

mtgl_SOURCES =  mtgl.cpp
mtgl_CPPFLAGS = $(EX_CFLAGS)
mtgl_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

yukawa_SOURCES =  yukawa.cpp
yukawa_CPPFLAGS = $(EX_CFLAGS)
yukawa_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

lambda_SOURCES =  lambda.cpp
lambda_CPPFLAGS = $(EX_CFLAGS)
lambda_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

yt_SOURCES =  yt.cpp
yt_CPPFLAGS = $(EX_CFLAGS)
yt_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

bezrukovFit_SOURCES =  bezrukovFit.cpp
bezrukovFit_CPPFLAGS = $(EX_CFLAGS)
bezrukovFit_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

plotyukawa_SOURCES =  plotyukawa.cpp
plotyukawa_CPPFLAGS = $(EX_CFLAGS)
plotyukawa_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

gauge_SOURCES =  gauge.cpp
gauge_CPPFLAGS = $(EX_CFLAGS)
gauge_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

mb_SOURCES =  mb.cpp
mb_CPPFLAGS = $(EX_CFLAGS)
mb_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

runAlpha_SOURCES =  runAlpha.cpp 
runAlpha_CPPFLAGS = $(EX_CFLAGS)
runAlpha_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

runSM_SOURCES =  runSM.cpp
runSM_CPPFLAGS = $(EX_CFLAGS)
runSM_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

# mratio_SOURCES =  mratio.cpp
# mratio_CPPFLAGS = $(EX_CFLAGS) -I$(top_srcdir)/tsil-1.21 -I$(top_srcdir)/mr
# mratio_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm
