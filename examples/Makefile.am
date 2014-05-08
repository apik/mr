AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}


noinst_PROGRAMS = mW mZ mH mt runAlpha

mW_SOURCES =  mW.cpp
mW_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mW_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm

mZ_SOURCES =  mZ.cpp
mZ_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mZ_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm ./libmodule1.la

mH_SOURCES =  mH.cpp
mH_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mH_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm

mt_SOURCES =  mt.cpp
mt_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
mt_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm

runAlpha_SOURCES =  runAlpha.cpp CRunDec.cpp CRunDec.h
runAlpha_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1 -I$(top_srcdir)/mr
runAlpha_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.1 -ltsil -lm
