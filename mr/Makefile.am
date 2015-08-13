AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

lib_LIBRARIES = libmr.a

# Main sources
MAINSRC = mr.hpp tsil.hpp tsil.cpp base.hpp mOS2mMSQCD.cpp pole2run.hpp pole2run.cpp p2ms.hpp p2ms.cpp tools.hpp tools.cpp

MAINHPP =  constants.hpp sminput.hpp timer.hpp  tsil.hpp  dr.hpp operators.hpp logger.hpp

# Running
RUNSRC = betaSM.hpp betaSM.cpp betaQCD.hpp betaQEDQCD.hpp betaQEDQCD.cpp betaQCD.cpp


# alpha from G-Fermi and pole masses
ALPHAGF = alphaGF.hpp alphaGF.cpp aGF/alphaGF10.cpp aGF/alphaGF11.cpp aGF/alphaGF20.cpp


# 
# MS-bar masses and constants
#
BBMS = bb.hpp bb.cpp mM/mbbbar01.cpp mM/mbbbar10.cpp mM/mbbbar11.cpp mM/mbbbar20.cpp yu/yubb01.cpp yu/yubb10.cpp yu/yubb11.cpp yu/yubb20.cpp
WWMS = WW.hpp WW.cpp mM/mWWbar10.cpp mM/mWWbar11.cpp mM/mWWbar20.cpp yu/yuWW10.cpp yu/yuWW11.cpp yu/yuWW20.cpp
ZZMS = ZZ.hpp ZZ.cpp mM/mZZbar10.cpp mM/mZZbar11.cpp mM/mZZbar20.cpp yu/yuZZ10.cpp yu/yuZZ11.cpp yu/yuZZ20.cpp
HHMS = HH.hpp HH.cpp mM/mHHbar10.cpp mM/mHHbar11.cpp mM/mHHbar20.cpp yu/yuHH10.cpp yu/yuHH11.cpp yu/yuHH20.cpp
TTMS = tt.hpp tt.cpp mM/mttbar01.cpp mM/mttbar10.cpp mM/mttbar11.cpp mM/mttbar20.cpp yu/yutt01.cpp yu/yutt10.cpp yu/yutt11.cpp yu/yutt20.cpp

# 
# Gauge-less limit
# 
HHMSGL = gl/mHHbarGL10.cpp gl/mHHbarGL11.cpp gl/mHHbarGL20.cpp yugl/yuHHGL10.cpp yugl/yuHHGL11.cpp yugl/yuHHGL20.cpp
BBMSGL = gl/mbbbarGL10.cpp gl/mbbbarGL11.cpp gl/mbbbarGL20.cpp yugl/yubbGL10.cpp yugl/yubbGL11.cpp yugl/yubbGL20.cpp
TTMSGL = gl/mttbarGL10.cpp gl/mttbarGL11.cpp gl/mttbarGL20.cpp yugl/yuttGL10.cpp yugl/yuttGL11.cpp yugl/yuttGL20.cpp 

#
# delta-r in MS-bar  scheme with OS input
# 
DRMSOS = dr.cpp dr/dr10.cpp dr/dr11.cpp dr/dr20.cpp dr/drgl10.cpp dr/drgl11.cpp dr/drgl20.cpp

#
# delta-r in MS-bar  scheme with MS-bar input
# 
DRMSMS = dr/drbar10.cpp dr/drbar11.cpp dr/drbar20.cpp

# W pole mass from MS input
WWMS2OS = Mm/mWWos10.cpp Mm/mWWos11.cpp Mm/mWWos20.cpp
# Z pole mass from MS input
ZZMS2OS = Mm/mZZos10.cpp Mm/mZZos11.cpp Mm/mZZos20.cpp
# H pole mass from MS input
HHMS2OS = Mm/mHHos10.cpp Mm/mHHos11.cpp Mm/mHHos20.cpp
# t pole mass from MS input
TTMS2OS = Mm/mttos01.cpp Mm/mttos10.cpp Mm/mttos11.cpp Mm/mttos20.cpp
# b pole mass from MS input
BBMS2OS = Mm/mbbos01.cpp Mm/mbbos10.cpp Mm/mbbos11.cpp Mm/mbbos20.cpp


libmr_a_SOURCES =    $(MAINSRC) $(MAINHPP) $(ALPHAGF) $(RUNSRC) $(BBMS) $(WWMS) $(ZZMS) $(HHMS) $(TTMS)\
 $(HHMSGL) $(BBMSGL) $(TTMSGL) $(DRMSOS) $(DRMSMS)\
 $(WWMS2OS) $(ZZMS2OS) $(HHMS2OS) $(TTMS2OS) $(BBMS2OS)

libmr_a_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir) -I$(top_srcdir)/tsil-1.21 -I$(top_srcdir)/Eigen

# lib_LTLIBRARIES = libmr.la

# libmr_la_SOURCES =    $(MAINSRC) $(RUNSRC) $(BBMS) $(WWMS) $(ZZMS) $(HHMS) $(TTMS)\
#  $(HHMSGL) $(BBMSGL) $(TTMSGL) $(DRMSOS)\
#  $(WWMS2OS) $(ZZMS2OS) $(HHMS2OS)


# libmr_la_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.21




# mr_LDADD = $(MR_LIBS) -L$(top_srcdir)/tsil-1.1 -ltsil -lm

