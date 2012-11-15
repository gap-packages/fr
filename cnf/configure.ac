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
AC_INIT(fr,,laurent.bartholdi@gmail.com)
AC_CONFIG_AUX_DIR([cnf])
AC_CONFIG_MACRO_DIR([cnf/m4])
AC_PREFIX_DEFAULT([${PWD}])

# Locates GAP
AC_FIND_GAP

AC_SYS_IS_CYGWIN
AC_SYS_IS_DARWIN

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

AC_OUTPUT
