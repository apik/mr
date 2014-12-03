AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

LDADD = ../mr/libmr.la
EX_CFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.21 -I$(top_srcdir)/mr

noinst_PROGRAMS = mW mZ mH mt runAlpha runSM runGL mb yukawa lambda yt gauge plotyukawa mtgl bezrukovFit GFermi tables poleRGI
if HAVE_ROOT
noinst_PROGRAMS += mHroot interpolation
endif 

if HAVE_BOOST
noinst_PROGRAMS += errest
endif 

mW_SOURCES =  mW.cpp
mW_CPPFLAGS = $(EX_CFLAGS)
mW_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

GFermi_SOURCES =  GFermi.cpp
GFermi_CPPFLAGS = $(EX_CFLAGS)
GFermi_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

mZ_SOURCES =  mZ.cpp
mZ_CPPFLAGS = $(EX_CFLAGS)
mZ_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm 

mH_SOURCES =  mH.cpp
mH_CPPFLAGS = $(EX_CFLAGS)
mH_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm


if HAVE_ROOT
mHroot_SOURCES =  mHroot.cpp
mHroot_CPPFLAGS = $(EX_CFLAGS) -I$(ROOTINCDIR)
mHroot_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm -L$(ROOTLIBDIR) $(ROOTLIBS)

interpolation_SOURCES =  interpolation.cpp
interpolation_CPPFLAGS = $(EX_CFLAGS) -I$(ROOTINCDIR)
interpolation_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm -L$(ROOTLIBDIR) $(ROOTLIBS)
endif


mt_SOURCES =  mt.cpp
mt_CPPFLAGS = $(EX_CFLAGS)
mt_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

mtgl_SOURCES =  mtgl.cpp
mtgl_CPPFLAGS = $(EX_CFLAGS)
mtgl_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

tables_SOURCES =  tables.cpp
tables_CPPFLAGS = $(EX_CFLAGS)
tables_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

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

runGL_SOURCES =  runGL.cpp
runGL_CPPFLAGS = $(EX_CFLAGS)
runGL_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

poleRGI_SOURCES =  poleRGI.cpp
poleRGI_CPPFLAGS = $(EX_CFLAGS)
poleRGI_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm


# mratio_SOURCES =  mratio.cpp
# mratio_CPPFLAGS = $(EX_CFLAGS) -I$(top_srcdir)/tsil-1.21 -I$(top_srcdir)/mr
# mratio_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm
if HAVE_BOOST
errest_SOURCES =  errest.cpp
errest_CPPFLAGS = $(EX_CFLAGS) -I$(top_srcdir)/tsil-1.21 -I$(top_srcdir)/mr $(BOOST_CPPFLAGS)
errest_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm $(BOOST_LDFLAGS)
endif 
