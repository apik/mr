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

CFLAGS	= $(TSIL_RESCALE) $(TSIL_SIZE) $(TSIL_OPT)
LDFLAGS	= 
LIBS	= -L. -ltsil -lm
ARFLAGS = -rus

TSIL_SRCS = initialize.c \
	    evaluate.c \
	    fevaluate.c \
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
	    setbold.c

TSIL_OBJS = $(TSIL_SRCS:.c=.o)

MAIN_OBJ = basecalc.o

TEST_SRC = testprog.c

TEST_OBJ = $(TEST_SRC:.c=.o)

all:	libtsil.a tsil tsil.tst

libtsil.a: $(TSIL_OBJS)
	$(AR) $(ARFLAGS) $@ $?

tsil:	$(MAIN_OBJ) libtsil.a
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(MAIN_OBJ) $(LIBS)

$(TEST_OBJ):	$(TEST_SRC)
	$(CC) -c $(CFLAGS) $(TSIL_TESTOPT) $(TEST_SRC)

tsil.tst:	$(TEST_OBJ) libtsil.a
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(TEST_OBJ) $(LIBS)

install: libtsil.a
	cp libtsil.a $(TSIL_LIBDIR)
	cp tsil.h $(TSIL_INCDIR)

tidy:
	rm -rf core* a.out *~

clean:
	rm -rf core* a.out tsil tsil.tst *~ *.o *.T libtsil.a


# ===== Header dependencies: =====

generic.o initialize.o: tsil_params.h

$(TSIL_OBJS):	internal.h tsil.h tsil_funcs.h tsil_global.h

tsil_main.o:	tsil.h

$(TEST_OBJ):	tsil_testparams.h internal.h tsil.h
