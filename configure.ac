AC_INIT([Standard Model 3-loop running and 2-loop matching], [0.1], [pikelner@theor.jinr.ru],
             [mr], [http://theor.jinr.ru/~pikelner/])
AC_PREREQ([2.59])
dnl AM_INIT_AUTOMAKE([1.8 -Wall no-define])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_HEADERS([config.h])
: ${CXXFLAGS="-fopenmp"}
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
dnl PKG_CHECK_MODULES(SM3R2M, ginac >= 1.3.7)
AC_CONFIG_FILES(Makefile tsil-1.1/Makefile mr/Makefile examples/Makefile)
AC_OUTPUT
