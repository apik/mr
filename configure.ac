dnl AC_INIT([diana], 2.3.6, https://github.com/apik/diana/)
dnl AC_CONFIG_SRCDIR(src/init.c)




AC_INIT([Standard Model 3-loop running and 2-loop matching], [0.1], [pikelner at theor.jinr.ru],
             [SM3r2m], [http://theor.jinr.ru/~pikelner/])
AC_PREREQ([2.59])
dnl AM_INIT_AUTOMAKE([1.8 -Wall no-define])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_HEADERS([config.h])
: ${CXXFLAGS="-fopenmp"}
AC_PROG_CXX
AC_PROG_CC
dnl PKG_CHECK_MODULES(SM3R2M, ginac >= 1.3.7)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT(src/Makefile)
