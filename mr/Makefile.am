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
BBMS = bb.hpp bb.cpp mbbbar01.cpp mbbbar02.cpp mbbbar03.cpp mbbbar10.cpp mbbbar11.cpp mbbbar20.cpp yubb01.cpp yubb10.cpp yubb11.cpp yubb20.cpp
WWMS = WW.hpp WW.cpp mWWbar10.cpp mWWbar11.cpp mWWbar20.cpp yuWW10.cpp yuWW11.cpp yuWW20.cpp
ZZMS = ZZ.hpp ZZ.cpp mZZbar10.cpp mZZbar11.cpp mZZbar20.cpp yuZZ10.cpp yuZZ11.cpp yuZZ20.cpp
HHMS = HH.hpp HH.cpp mHHbar10.cpp mHHbar11.cpp mHHbar20.cpp yuHH10.cpp yuHH11.cpp yuHH20.cpp
TTMS = tt.hpp tt.cpp mttbar01.cpp mttbar02.cpp mttbar03.cpp mttbar10.cpp mttbar11.cpp mttbar20.cpp yutt01.cpp yutt10.cpp yutt11.cpp yutt20.cpp

# 
# Gauge-less limit
# 
HHMSGL = mHHbarGL10.cpp mHHbarGL11.cpp mHHbarGL20.cpp yuHHGL10.cpp yuHHGL11.cpp yuHHGL20.cpp
BBMSGL = mbbbarGL10.cpp mbbbarGL11.cpp mbbbarGL20.cpp yubbGL10.cpp yubbGL11.cpp yubbGL20.cpp
TTMSGL = mttbarGL10.cpp mttbarGL11.cpp mttbarGL20.cpp yuttGL10.cpp yuttGL11.cpp yuttGL20.cpp

#
# delta-r in MS-bar  scheme with OS input
# 
DRMSOS = dr.cpp dr10.cpp dr11.cpp dr20.cpp drgl10.cpp drgl11.cpp drgl20.cpp

# W pole mass from MS input
WWMS2OS = mWWos10.cpp mWWos11.cpp mWWos20.cpp
# Z pole mass from MS input
ZZMS2OS = mZZos10.cpp mZZos11.cpp mZZos20.cpp
# H pole mass from MS input
HHMS2OS = mHHos10.cpp mHHos11.cpp mHHos20.cpp


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

