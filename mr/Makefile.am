AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

lib_LIBRARIES = libmr.a

libmr_a_SOURCES =  betaSM.hpp  mr.hpp betaSM.cpp betaQCD.hpp \
WW.hpp WW.cpp WW10.cpp WW11.cpp WW20.cpp yW10.cpp yW11.cpp yW20.cpp\
ZZ.hpp ZZ.cpp ZZ10.cpp ZZ11.cpp ZZ20.cpp yZ10.cpp yZ11.cpp yZ20.cpp\
HH.hpp HH.cpp HH10.cpp HH11.cpp HH20.cpp lam10.cpp lam11.cpp lam20.cpp\
bb.hpp bb.cpp bb01.cpp bb10.cpp bb11.cpp bb20.cpp yb01.cpp yb10.cpp yb11.cpp yb20.cpp \
tt.hpp tt.cpp tt01.cpp tt10.cpp tt11.cpp tt20.cpp yt01.cpp yt10.cpp yt11.cpp yt20.cpp \
tsil.hpp tsil.cpp 

libmr_a_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1

# lib_LTLIBRARIES = libmr.la

# libmr_la_SOURCES =  betaSM.hpp  mr.hpp betaSM.cpp betaQCD.hpp \
# WW.hpp WW.cpp WW10.cpp WW11.cpp WW20.cpp \
# ZZ.hpp ZZ.cpp ZZ10.cpp ZZ11.cpp ZZ20.cpp \
# HH.hpp HH.cpp HH10.cpp HH11.cpp HH20.cpp \
# bb.hpp bb.cpp bb01.cpp bb10.cpp bb11.cpp bb20.cpp yb01.cpp yb10.cpp yb11.cpp yb20.cpp \
# tt.hpp tt.cpp tt01.cpp tt10.cpp tt11.cpp tt20.cpp yt01.cpp yt10.cpp yt11.cpp yt20.cpp \
# tsil.hpp tsil.cpp 

# libmr_la_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1




# mr_LDADD = $(MR_LIBS) -L$(top_srcdir)/tsil-1.1 -ltsil -lm

# sampler_SOURCES = sampler.hpp sampler.cpp
# sampler_CPPFLAGS = $(BOXSAMPLER_CFLAGS)
# sampler_LDADD = $(BOXSAMPLER_LIBS) -lsector_decomposition
