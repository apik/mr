AC_INIT([Standard Model 3-loop running and 2-loop matching], [0.1], [pikelner@theor.jinr.ru],
             [mr], [http://theor.jinr.ru/~pikelner/])
AC_PREREQ([2.59])
AC_CONFIG_MACRO_DIR([m4])
dnl AM_INIT_AUTOMAKE([1.8 -Wall no-define])
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_HEADERS([config.h])
: ${CXXFLAGS="-O3"}
AC_PROG_CXX
AC_PROG_CC

AM_CONDITIONAL(HAVE_BOOST, false)
AX_BOOST_BASE([1.30], AM_CONDITIONAL(HAVE_BOOST, true), [])

dnl AC_DEFINE([HAVE_ROOT], [], [CERN ROOT installation found])
dnl ROOT_PATH(, AC_DEFINE([HAVE_ROOT]))
AM_CONDITIONAL(HAVE_ROOT, false)
ROOT_PATH(, AM_CONDITIONAL(HAVE_ROOT, true))


AC_MSG_NOTICE([-------------- Testing math.])

AC_ARG_WITH([mcc-path],
            [AS_HELP_STRING([--with-mcc-path],
            [build Mathematica MathLink interface @<:@default=auto@:>@. Set this to full path to mcc tool])],
            [test "x$with_mcc_path" == "xyes" && with_mcc_path=auto],
            [with_mcc_path=auto])

AC_MSG_NOTICE([test "x$with_mcc_path" == "xyes" && with_mcc_path=auto])

# ------------------------- # -------------------------
if test "x$with_mcc_path" == "xauto" ; then
AC_MSG_CHECKING([for mcc in system path])
AC_PATH_PROGS_FEATURE_CHECK([MCC], [mcc mcc9.0 mcc8.0 mcc7.0 mcc6.0 ],
[[
  cat - <<_ACEOF > conftest.cpp
    #include <mathlink.h>
    int main () { return 0; }
_ACEOF
  rm -f conftest$ac_exeext
  $ac_path_MCC conftest.cpp -o conftest$ac_exeext > conftest.err
  if test -s conftest$ac_exeext; then
    grep -v '^ *+' conftest.err >conftest.er1
    cat conftest.er1 >&5
    mv -f conftest.er1 conftest.err
  fi
  test ! -s conftest.err && test -s conftest$ac_exeext \
    && ac_cv_path_MCC=$ac_path_MCC ac_path_MCC_found=:
  rm -f conftest.cpp conftest$ac_exeext conftest.err
]], [with_mcc_path=])
if test "x$ac_cv_path_MCC" != "x" ; then
  AC_SUBST([MCC],[$ac_cv_path_MCC])
  AC_MSG_RESULT($ac_cv_path_MCC)
else
  AC_MSG_RESULT(no)     
fi                      
if test `uname -s` == "Linux" ; then
  case `uname -m` in
    ia64)   SystemIDList="Linux-IA64";;
    i?86)   SystemIDList="Linux";;
    x86_64) SystemIDList="Linux-x86-64";;
    *)      SystemIDList="Unknown";;
  esac
  math_p1a=/opt/products/mathematica
  math_p1b=/usr/local/Wolfram/Mathematica
  math_p2=SystemFiles/Links/MathLink/DeveloperKit/$SystemIDList/CompilerAdditions
  math_ps=$math_p1a/9.0/$math_p2$PATH_SEPARATOR$math_p1b/9.0/$math_p2$PATH_SEPARATOR$math_p1a/8.0/$math_p2$PATH_SEPARATOR$math_p1b/8.0/$math_p2$PATH_SEPARATOR$math_p1a/7.0/$math_p2$PATH_SEPARATOR$math_p1b/7.0/$math_p2$PATH_SEPARATOR$math_p1a/6.0/$math_p2$PATH_SEPARATOR$math_p1b/6.0/$math_p2
  AC_PATH_PROGS(MPREP,[mprep],[],[$math_ps])
  if test "x$ac_cv_path_MCC" == "x" ; then
    AC_PATH_PROGS(MCC,[mcc mcc9.0 mcc8.0 mcc7.0 mcc6.0],[],[$math_ps])
  fi
fi
# -------------------------
else
  if test "x$with_mcc_path" != "xno" ; then
    if test -d $with_mcc_path ; then
      ac_cv_path_MCC=$with_mcc_path/mcc
      ac_cv_path_MPREP=$with_mcc_path/mprep
    else
      ac_cv_path_MCC=$with_mcc_path
      ac_cv_path_MPREP=$with_mcc_path/mprep
    fi
    if test "x$ac_cv_path_MCC" != "x" && test -x $ac_cv_path_MCC ; then
      AC_SUBST([MCC],[$ac_cv_path_MCC])
    fi
    if test "x$ac_cv_path_MPREP" != "x" && test -x $ac_cv_path_MPREP ; then
      AC_SUBST([MPREP],[$ac_cv_path_MPREP])
    fi
  fi
fi
# ------------------------- # ------------------------- # -------------------------
dnl AM_CONDITIONAL(HAVE_MCC,   test "x$with_mcc_path" != "xno" && test "x$ac_cv_path_MCC"   != "x")
AM_CONDITIONAL(HAVE_MPREP, test "x$with_mcc_path" != "xno" && test "x$ac_cv_path_MPREP" != "x")

AC_MSG_NOTICE([$MPREP])

AC_CANONICAL_HOST       
case $host_os in        
  darwin* )   
          LDMATH="-lMLi4 -lc++ -framework Foundation"
          AC_SUBST(LDMATH)  
          ALTMATHPATH="AlternativeLibraries"
          AC_SUBST(ALTMATHPATH)  
          ;;
  linux*)
          LDMATH="-lML64i3 -lm -lpthread -lrt -lstdc++"
          AC_SUBST(LDMATH)  
          ALTMATHPATH=""
          AC_SUBST(ALTMATHPATH)  
          ;;
   *BSD*)
          LDMATH="-lML64i3 -lm -lpthread -lrt -lstdc++"
          AC_SUBST(LDMATH)  
          ALTMATHPATH=""
          AC_SUBST(ALTMATHPATH)  
          ;;
    *)

    AC_MSG_ERROR([Mathlink on your platform is not currently supported])
    ;;
esac

AC_MSG_NOTICE([$host_os])

AC_PROG_RANLIB
dnl LT_INIT

dnl PKG_CHECK_MODULES(SM3R2M, ginac >= 1.3.7)
AC_CONFIG_FILES(Makefile tsil-1.21/Makefile mr/Makefile examples/Makefile tests/Makefile math/Makefile Minuit/Makefile)
AC_OUTPUT
