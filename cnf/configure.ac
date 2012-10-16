#############################################################################
##
#W  configure.ac                                            Laurent Bartholdi
##
#H   @(#)$Id$
##
#Y Copyright (C) 2009, Laurent Bartholdi
##
#############################################################################

AC_PREREQ(2.68)
LT_PREREQ([2.4.2])
AC_INIT(fr,,laurent.bartholdi@gmail.com)
AC_CONFIG_SRCDIR([src/fr_dll.c])
AC_CONFIG_AUX_DIR([cnf])
AC_CONFIG_MACRO_DIR([cnf/m4])
LT_INIT([disable-static dlopen win32-dll])

AC_PREFIX_DEFAULT([${PWD}])

# Checks for programs.
AC_PROG_CC
AX_CC_MAXOPT

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

# Locates GAP
AC_FIND_GAP

# Finds if DOS/MacOSX install
AC_CYGWIN
AM_CONDITIONAL([SYS_IS_CYGWIN], [test "$CYGWIN" = "yes"])
if test "$CYGWIN" = "yes"; then
  AC_DEFINE(SYS_IS_CYGWIN32, 1, are we on CYGWIN?)
else
  AC_DEFINE(SYS_IS_CYGWIN32, 0, are we on CYGWIN?)
fi

case "$host" in
  *-darwin* )
  AC_DEFINE(SYS_IS_DARWIN, 1, are we on DARWIN?)
  ;;
  * )
  AC_DEFINE(SYS_IS_DARWIN, 0, are we on DARWIN?)
  ;;
esac

EXTERN="\$(CURDIR)/bin/$TARGET/extern"

################################################################
# Check for gsl library

LIB_TARGET=""

GSLDIR="yes"
GSLINCLUDE=""
GSLLIB=""
AC_ARG_WITH(gsl,
 [  --with-gsl=<path>|yes|extern
    Location at which the GSL library was installed.
    If the argument is omitted, the library is assumed to be reachable
    under the standard search path (/usr, /usr/local,...).  Otherwise
    you must give the <path> to the directory which contains the
    library. The special value "extern", which is the default, asks Float
    to compile a version of gsl in the subdirectory extern/.
 ],
 [GSLDIR="$withval"]
)

AC_ARG_WITH(gsl-include,
 [  --with-gsl-include=<location>
    Location at which the GSL include files were installed.],
 [GSLINCLUDE="$withval"]
)

AC_ARG_WITH(gsl-lib,
 [  --with-gsl-lib=<location>
    Location at which the GSL library files were installed.],
 [GSLLIB="$withval"]
)

if test "$GSLDIR" = extern; then
    LIB_TARGET="$LIB_TARGET gsllib"
else
    if test "$GSLDIR" != yes; then
        if test "$GSLINCLUDE" == ""; then GSLINCLUDE="$GSLDIR/include"; fi
	if test "$GSLLIB" == ""; then GSLLIB="$GSLDIR/lib"; fi
	CPPFLAGS="$CPPFLAGS -I$GSLINCLUDE"
    	GAP_CPPFLAGS="$GAP_CPPFLAGS -I$GSLINCLUDE"
        LIBS="$LIBS -L$GSLLIB"
    fi
    AC_CHECK_HEADER(gsl/gsl_vector.h,[],[AC_MSG_ERROR([library gsl not found. Specify its location using --with-gsl])],[])
    AC_CHECK_LIB([gsl],[gsl_multiroot_fsolver_set],,
	[AC_ERROR([The GSL library could not be found. It is needed for IMG calculations.])],[-lgslcblas])
    AC_CHECK_LIB([gslcblas],[cblas_ctrmv],,
        [AC_ERROR([The GSL CBlas library could not be found. It is needed for IMG calculations.])])
fi

if test "$GSLLIB" != ""; then
    LIBS="$GAP_LDFLAGS -L$GSLLIB -Wl,-rpath,$GSLLIB"
fi

LIBS="$LIBS -lgsl -lgslcblas"

################################################################
# external programs configuration

AC_PATH_PROGS(DOT,[dot],[],[$PATH$PATH_SEPARATOR/usr/local/graphviz/bin])

if test -z "$DOT"; then
    AC_WARN([Could not find 'dot' (debian package graphviz)... you won't be able to draw automata])
fi

AC_PATH_PROG(DISP,[display])

if test -z "$DISP"; then
    AC_WARN([Could not find 'display' (debian package imagemagick)... you won't be able to draw automata])
fi

AC_PATH_PROG(APPLETVIEWER,[appletviewer])

if test -z "$APPLETVIEWER"; then
    AC_WARN([Could not find 'appletviewer' (debian package java-6-sdk)... you won't be able to draw spiders])
fi

# Check for java compiler
AC_CHECK_PROGS(JAVAC,[javac])

if test -z "$JAVAC"; then
    AC_WARN([Could not find java compiler... you won't be able to draw spiders])
else
    JAVABUILD="java/javaplot.class java/javaview.jar"
fi
AC_SUBST(JAVABUILD)
AC_SUBST(LIBS)
AC_SUBST(LIB_TARGET)

################################################################
# generate files

mkdir -p bin/$GAPARCH
CONFIG_STATUS=bin/$GAPARCH/config.status

AC_CONFIG_FILES([$GAP_MAKEFILE:cnf/Makefile.in])

if test "$GAP_MAKEFILE" != Makefile; then
    ln -sf "$GAP_MAKEFILE" Makefile
fi

AC_CONFIG_SUBDIRS(hurwitz)

AC_OUTPUT
