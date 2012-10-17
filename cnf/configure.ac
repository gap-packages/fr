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

AC_SYS_IS_CYGWIN
AC_SYS_IS_DARWIN

################################################################
# external libraries configuration
EXTERN="\$(CURDIR)/bin/$TARGET/extern"
LIB_TARGET=""

AC_CHECK_LEVMAR

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
