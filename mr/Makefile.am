AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

lib_LIBRARIES = libmr.a

libmr_a_SOURCES =  betaSM.hpp  mr.hpp betaSM.cpp betaQCD.hpp \
WW.hpp WW.cpp WW10.cpp WW11.cpp WW20.cpp \
HH.hpp HH.cpp HH10.cpp HH11.cpp HH20.cpp \
tt.hpp tt.cpp tt01.cpp tt10.cpp tt11.cpp tt20.cpp \
tsil.hpp tsil.cpp 

libmr_a_CPPFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.1
# mr_LDADD = $(MR_LIBS) -L$(top_srcdir)/tsil-1.1 -ltsil -lm

# sampler_SOURCES = sampler.hpp sampler.cpp
# sampler_CPPFLAGS = $(BOXSAMPLER_CFLAGS)
# sampler_LDADD = $(BOXSAMPLER_LIBS) -lsector_decomposition
