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
AC_PROG_CC
AC_FIND_GAP

AC_SYS_IS_CYGWIN
AC_SYS_IS_DARWIN

################################################################
# external programs configuration

AC_PATH_PROGS(DOT,[dot],[],[$PATH$PATH_SEPARATOR/usr/local/graphviz/bin])

if test -z "$DOT"; then
    AC_WARN([Could not find 'dot' (debian package graphviz)... you won't be able to draw automata])
fi

AC_PATH_PROG(PSVIEWER,[display])
if test -n "$PSVIEWER"; then
    PSVIEWER="$PSVIEWER -flatten -";
else
    AC_PATH_PROG(PSVIEWER,[gv])
    if test -n "$PSVIEWER"; then
        PSVIEWER="$PSVIEWER -";
    else
	PSVIEWER=false;
    	AC_PATH_PROG(SVGVIEWER,[rsvg-view])

    	if test -n "$SVGVIEWER"; then
            SVGVIEWER="$SVGVIEWER --stdin";
	else
	    AC_WARN([Could not find 'display' (debian package imagemagick) or 'gv' (debian package gv) or 'rsvg-view' (debian package librsvg2-bin)... you won't be able to draw automata])
    	fi
    fi
fi

################################################################
# generate files

mkdir -p bin/$GAPARCH
CONFIG_STATUS=bin/$GAPARCH/config.status

AC_CONFIG_FILES([$GAP_MAKEFILE:cnf/Makefile.in bin/$GAPARCH/files.g:cnf/files.g.in])

if test "$GAP_MAKEFILE" != Makefile; then
    ln -sf "$GAP_MAKEFILE" Makefile
fi

AC_OUTPUT
