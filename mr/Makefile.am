AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

lib_LIBRARIES = libmr.a

# Main sources
MAINSRC = mr.hpp tsil.hpp tsil.cpp base.hpp

# Running
RUNSRC = betaSM.hpp betaSM.cpp betaQCD.hpp CRunDec.cpp CRunDec.h


# alpha from G-Fermi and pole masses
ALPHAGF = alphaGF.hpp alphaGF.cpp alphaGF10.cpp alphaGF11.cpp alphaGF20.cpp


# 
# MS-bar masses and constants
#
BBMS = bb.hpp bb.cpp bb01.cpp bb10.cpp bb11.cpp bb20.cpp yb01.cpp yb10.cpp yb11.cpp yb20.cpp
WWMS = WW.hpp WW.cpp WW10.cpp WW11.cpp WW20.cpp yW10.cpp yW11.cpp yW20.cpp
ZZMS = ZZ.hpp ZZ.cpp ZZ10.cpp ZZ11.cpp ZZ20.cpp yZ10.cpp yZ11.cpp yZ20.cpp
HHMS = HH.hpp HH.cpp HH10.cpp HH11.cpp HH20.cpp yH10.cpp yH11.cpp yH20.cpp
TTMS = tt.hpp tt.cpp tt01.cpp tt02.cpp tt03.cpp tt10.cpp tt11.cpp tt20.cpp yt01.cpp yt10.cpp yt11.cpp yt20.cpp

# 
# Gauge-less limit
# 
HHMSGL = HHgl10.cpp HHgl11.cpp HHgl20.cpp lamgl10.cpp lamgl11.cpp lamgl20.cpp
BBMSGL = bbgl10.cpp bbgl11.cpp bbgl20.cpp ybgl10.cpp ybgl11.cpp ybgl20.cpp
TTMSGL = ttgl10.cpp ttgl11.cpp ttgl20.cpp ytgl10.cpp ytgl11.cpp ytgl20.cpp

#
# delta-r in MS-bar  scheme with OS input
# 
DRMSOS = dr.cpp dr10.cpp dr11.cpp dr20.cpp drgl10.cpp drgl11.cpp drgl20.cpp

# W pole mass from MS input
WWMS2OS = WW2MS10.cpp WW2MS11.cpp WW2MS20.cpp
# Z pole mass from MS input
ZZMS2OS = ZZ2MS10.cpp ZZ2MS11.cpp ZZ2MS20.cpp
# H pole mass from MS input
HHMS2OS = HH2MS10.cpp HH2MS11.cpp HH2MS20.cpp


libmr_a_SOURCES =    $(MAINSRC) $(ALPHAGF) $(RUNSRC) $(BBMS) $(WWMS) $(ZZMS) $(HHMS) $(TTMS)\
 $(HHMSGL) $(BBMSGL) $(TTMSGL) $(DRMSOS)\
 $(WWMS2OS) $(ZZMS2OS) $(HHMS2OS)

libmr_a_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.21

# lib_LTLIBRARIES = libmr.la

# libmr_la_SOURCES =    $(MAINSRC) $(RUNSRC) $(BBMS) $(WWMS) $(ZZMS) $(HHMS) $(TTMS)\
#  $(HHMSGL) $(BBMSGL) $(TTMSGL) $(DRMSOS)\
#  $(WWMS2OS) $(ZZMS2OS) $(HHMS2OS)


# libmr_la_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.21




# mr_LDADD = $(MR_LIBS) -L$(top_srcdir)/tsil-1.1 -ltsil -lm

