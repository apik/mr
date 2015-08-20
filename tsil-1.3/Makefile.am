AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS} -I m4

lib_LIBRARIES = libtsil.a
####################################################################
# 
#    Makefile for TSIL (Two-loop Self-energy Integral Library)
#
######################## DATA SIZE #################################
# Uncomment one of these to choose the floating point data size (if
# unset, LONG will be used):
#
 TSIL_SIZE	= -DTSIL_SIZE_LONG    # long double/long double complex
#
# TSIL_SIZE	= -DTSIL_SIZE_DOUBLE  # double/double complex
#
########################## COMPILER ################################
# Uncomment the choices appropriate for your computing environment.
# (Not guaranteed to be complete. Please send info to the authors if
# you succeed in making things work on a system requiring different
# flags.)
#
# Intel C compiler:
# =================
# CC	 	= icc
# TSIL_OPT 	= -O3 -unroll -w
#
# GNU C Compiler:
# ===============
 CC		= gcc
 TSIL_OPT	= -O3 -funroll-loops 
#
################### INSTALLATION DIRECTORIES #######################
#
# TSIL_LIBDIR	= /usr/lib
# TSIL_INCDIR	= /usr/include
#
################### TEST SUITE MODE ################################
#
# TSIL_TESTOPT	= -DTSIL_TEST_STU
# TSIL_TESTOPT	= -DTSIL_TEST_ST
#
####################################################################
####################################################################
#
#  Users should not need to modify anything below this point!!!
#
####################################################################
####################################################################

libtsil_a_SOURCES =   initialize.c \
		    evaluate.c \
		    generic.c \
		    rk6.c \
		    rk5.c \
		    functions.c \
		    dilog.c \
		    trilog.c \
		    dBds.c \
		    dSds.c \
		    dTds.c \
		    dUds.c \
		    dsMds.c \
		    initB.c \
		    initS.c \
		    initT.c \
		    initU.c \
		    initV.c \
		    initM.c \
		    special.c \
		    analyticAB.c \
                    analyticI.c \
		    analyticS.c \
		    analyticT.c \
		    analyticTbar.c \
		    analyticU.c \
	   	    analyticV.c \
	   	    analyticM.c \
	   	    setTbar.c \
	    	    setV.c \
	    	    setbold.c \
                    internal.h tsil.h tsil_funcs.h tsil_global.h tsil_names.h tsil_params.h


libtsil_a_CFLAGS = $(TSIL_RESCALE) $(TSIL_SIZE) $(TSIL_OPT)

EXTRA_DIST = CHANGELOG.txt  ISSUES.txt  LICENSE.txt  README.txt  foo.txt

pkginclude_HEADERS = tsil.h 
