AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

noinst_PROGRAMS = SM3r2m

SM3r2m_SOURCES =  betaSM.hpp  SM3r2m.cpp  SM3r2m.hpp betaSM.cpp \
tt.hpp tt.cpp \
WW.hpp WW.cpp WW10.cpp WW11.cpp WW20.cpp \
HH.hpp HH.cpp HH10.cpp HH11.cpp HH20.cpp \
topal.cpp topas.cpp topalas.cpp topalal.cpp LambertW.cc LambertW.h \
mtopQCDEW2l.cpp  mtopQCDEW2l.hpp \
tsil.hpp tsil.cpp 
SM3r2m_CPPFLAGS = $(SM3R2M_CFLAGS) -I$(top_srcdir)/tsil-1.1
SM3r2m_LDADD = $(SM3R2M_LIBS) -L$(top_srcdir)/tsil-1.1 -ltsil -lm

# sampler_SOURCES = sampler.hpp sampler.cpp
# sampler_CPPFLAGS = $(BOXSAMPLER_CFLAGS)
# sampler_LDADD = $(BOXSAMPLER_LIBS) -lsector_decomposition
