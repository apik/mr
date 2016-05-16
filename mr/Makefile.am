AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

lib_LIBRARIES = libmr.a

# Main sources
MAINSRC = mr.hpp tsil.hpp tsil.cpp base.hpp tdecl.hpp mOS2mMSQCD.cpp p2ms.hpp p2ms.cpp

MAINHPP =  constants.hpp sminput.hpp timer.hpp  tsil.hpp  dr.hpp operators.hpp logger.hpp

# Running
RUNSRC = smRGE.hpp smRGE.cpp alphas.hpp alphas.cpp


# alpha from G-Fermi and pole masses
ALPHAGF = alphaGF.hpp alphaGF.cpp aGF/alphaGF10.cpp aGF/alphaGF11.cpp aGF/alphaGF20.cpp


# 
# MS-bar masses and constants
#
BBMS = bb.hpp bb.cpp mMSmOS/mbbbar01.cpp mMSmOS/mbbbar10.cpp mMSmOS/mbbbar11.cpp mMSmOS/mbbbar20.cpp yu/yubb01.cpp yu/yubb10.cpp yu/yubb11.cpp yu/yubb20.cpp
WWMS = WW.hpp WW.cpp mMSmOS/mWWbar10.cpp mMSmOS/mWWbar11.cpp mMSmOS/mWWbar20.cpp yu/yuWW10.cpp yu/yuWW11.cpp yu/yuWW20.cpp
ZZMS = ZZ.hpp ZZ.cpp mMSmOS/mZZbar10.cpp mMSmOS/mZZbar11.cpp mMSmOS/mZZbar20.cpp yu/yuZZ10.cpp yu/yuZZ11.cpp yu/yuZZ20.cpp
HHMS = HH.hpp HH.cpp mMSmOS/mHHbar10.cpp mMSmOS/mHHbar11.cpp mMSmOS/mHHbar20.cpp yu/yuHH10.cpp yu/yuHH11.cpp yu/yuHH20.cpp
TTMS = tt.hpp tt.cpp mMSmOS/mttbar01.cpp mMSmOS/mttbar10.cpp mMSmOS/mttbar11.cpp mMSmOS/mttbar20.cpp yu/yutt01.cpp yu/yutt10.cpp yu/yutt11.cpp yu/yutt20.cpp

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
WWMS2OS = mOSmMS/mWWos10.cpp mOSmMS/mWWos11.cpp mOSmMS/mWWos20.cpp
# Z pole mass from MS input
ZZMS2OS = mOSmMS/mZZos10.cpp mOSmMS/mZZos11.cpp mOSmMS/mZZos20.cpp
# H pole mass from MS input
HHMS2OS = mOSmMS/mHHos10.cpp mOSmMS/mHHos11.cpp mOSmMS/mHHos20.cpp
# t pole mass from MS input
TTMS2OS = mOSmMS/mttos01.cpp mOSmMS/mttos10.cpp mOSmMS/mttos11.cpp mOSmMS/mttos20.cpp
# b pole mass from MS input
BBMS2OS = mOSmMS/mbbos01.cpp mOSmMS/mbbos10.cpp mOSmMS/mbbos11.cpp mOSmMS/mbbos20.cpp


libmr_a_SOURCES =    $(MAINSRC) $(MAINHPP) $(ALPHAGF) $(RUNSRC) $(BBMS) $(WWMS) $(ZZMS) $(HHMS) $(TTMS)\
 $(HHMSGL) $(BBMSGL) $(TTMSGL) $(DRMSOS) $(DRMSMS)\
 $(WWMS2OS) $(ZZMS2OS) $(HHMS2OS) $(TTMS2OS) $(BBMS2OS)

libmr_a_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.3 $(OPENMP_CFLAGS)

EXTRA_DIST = aGF/README dr/README gl/README mMSmOS/README mOSmMS/README yu/README yugl/README

pkginclude_HEADERS = \
HH.hpp \
WW.hpp \
ZZ.hpp \
alphaGF.hpp \
alphas.hpp \
base.hpp \
bb.hpp \
constants.hpp \
dr.hpp \
logger.hpp \
mr.hpp \
operators.hpp \
p2ms.hpp \
smRGE.hpp \
sminput.hpp \
tdecl.hpp \
timer.hpp \
tsil.hpp \
tt.hpp



