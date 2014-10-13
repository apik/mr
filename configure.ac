AC_INIT([Standard Model 3-loop running and 2-loop matching], [0.1], [pikelner@theor.jinr.ru],
             [mr], [http://theor.jinr.ru/~pikelner/])
AC_PREREQ([2.59])
AC_CONFIG_MACRO_DIR([m4])
dnl AM_INIT_AUTOMAKE([1.8 -Wall no-define])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_HEADERS([config.h])
: ${CXXFLAGS=""}
AC_PROG_CXX
AC_PROG_CC

AC_DEFINE([HAVE_ROOT], [], [CERN ROOT installation found])
dnl ROOT_PATH(, AC_DEFINE([HAVE_ROOT]))
AM_CONDITIONAL(HAVE_ROOT, false)
ROOT_PATH(, AM_CONDITIONAL(HAVE_ROOT, true))

dnl AM_CONDITIONAL(HAVE_ROOT, test x$ROOTLIBDIR != x)
AC_PROG_RANLIB
dnl LT_INIT

dnl PKG_CHECK_MODULES(SM3R2M, ginac >= 1.3.7)
AC_CONFIG_FILES(Makefile tsil-1.21/Makefile mr/Makefile examples/Makefile)
AC_OUTPUT
