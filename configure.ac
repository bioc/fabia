AC_INIT([fabia], 1.5.3) dnl package name, version


dnl Find the compiler and compiler flags to use
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
echo "could not determine R_HOME"
exit 1
fi

CC=`"${R_HOME}/bin/R" CMD config CC`
CPP=`"${R_HOME}/bin/R" CMD config CPP`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
LDFLAGS=`"${R_HOME}/bin/R" CMD config LDFLAGS`
AC_PROG_CC
AC_PROG_CPP

dnl and do substitution in the src/Makevars.in 
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT

