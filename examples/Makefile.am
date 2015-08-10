AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

LDADD = ../mr/libmr.la
EX_CFLAGS = $(MR_CFLAGS) -I$(top_srcdir)/tsil-1.21 -I$(top_srcdir)/mr -I$(top_srcdir)/Eigen

noinst_PROGRAMS = exPole2ms exPole2run exMHcrit exMtcrit exTables exTadpoles exPhaseDiagram exPole2fixed runAlpha exLandauPoles

# if HAVE_ROOT
# noinst_PROGRAMS += mHroot interpolation
# endif 

# # if HAVE_BOOST
# noinst_PROGRAMS += errest
# # endif 

# CPC examples
exPole2ms_SOURCES =  exPole2ms.cpp
exPole2ms_CPPFLAGS = $(EX_CFLAGS)
exPole2ms_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

exPole2run_SOURCES =  exPole2run.cpp
exPole2run_CPPFLAGS = $(EX_CFLAGS)
exPole2run_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

exPole2fixed_SOURCES =  exPole2fixed.cpp
exPole2fixed_CPPFLAGS = $(EX_CFLAGS)
exPole2fixed_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

exPhaseDiagram_SOURCES =  exPhaseDiagram.cpp
exPhaseDiagram_CPPFLAGS = $(EX_CFLAGS)
exPhaseDiagram_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

exMHcrit_SOURCES =  exMHcrit.cpp
exMHcrit_CPPFLAGS = $(EX_CFLAGS)
exMHcrit_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

exMtcrit_SOURCES =  exMtcrit.cpp
exMtcrit_CPPFLAGS = $(EX_CFLAGS)
exMtcrit_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

exTables_SOURCES =  exTables.cpp
exTables_CPPFLAGS = $(EX_CFLAGS)
exTables_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

exTadpoles_SOURCES =  exTadpoles.cpp
exTadpoles_CPPFLAGS = $(EX_CFLAGS)
exTadpoles_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

exLandauPoles_SOURCES =  exLandauPoles.cpp
exLandauPoles_CPPFLAGS = $(EX_CFLAGS)
exLandauPoles_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm

runAlpha_SOURCES =  runAlpha.cpp CRunDec.cpp
runAlpha_CPPFLAGS = $(EX_CFLAGS)
runAlpha_LDADD = $(MR_LIBS) -L$(top_srcdir)/mr -lmr -L$(top_srcdir)/tsil-1.21 -ltsil -lm
