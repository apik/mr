AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}


noinst_PROGRAMS = mH mW runAlpha

mH_SOURCES =  mH.cpp
mH_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mH_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm

mW_SOURCES =  mW.cpp
mW_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mW_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm


runAlpha_SOURCES =  runAlpha.cpp CRunDec.cpp CRunDec.h
runAlpha_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
runAlpha_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm
