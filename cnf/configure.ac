#############################################################################
##
#W  configure.ac                                            Laurent Bartholdi
##
#H   @(#)$Id$
##
#Y Copyright (C) 2009, Laurent Bartholdi
##
#############################################################################

AC_PREREQ(2.00)
AC_INIT(fr,,laurent.bartholdi@gmail.com)
AC_CONFIG_SRCDIR([src/fr_dll.c])
AC_CONFIG_AUX_DIR(cnf)

# Checks for programs.
AC_PROG_CC

# Check for -fno-stack-protector, because we link within GAP
AC_CACHE_CHECK([whether $CC accepts -fno-stack-protector],
    [ns_cv_cc__nostackprotector],
    [save_CFLAGS=$CFLAGS
     CFLAGS="$CFLAGS -fno-stack-protector"
     AC_LINK_IFELSE([AC_LANG_PROGRAM([], [])],
                    [ns_cv_cc__nostackprotector=yes],
                    [ns_cv_cc__nostackprotector=no])
     CFLAGS=$save_CFLAGS])
if test $ns_cv_cc__nostackprotector = yes; then
     CFLAGS="$CFLAGS -fno-stack-protector"
fi

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([float.h stdlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE

################################################################
# check for GAP
AC_ARG_VAR([GAPDIR], [Location of the GAP root directory, e.g. ../..])

if test -z "$GAPDIR"; then
    for dir in ../.. /Applications/gap4r5 /usr/local/src/gap /usr/local/gap; do
	if test -f $dir/sysinfo.gap; then
	    GAPDIR=$dir
	    break
	fi
    done
fi

AC_ARG_VAR([CONFIGNAME],[Name of GAP build configuration])
if test -n "$CONFIGNAME"; then
    SYSINFO="$GAPDIR/sysinfo.gap-$CONFIGNAME"
    MAKEFILE="Makefile-$CONFIGNAME"
    GAPPROG="$GAPDIR/bin/gap-$CONFIGNAME.sh"
else
    SYSINFO="$GAPDIR/sysinfo.gap"
    MAKEFILE="Makefile"
    GAPPROG="$GAPDIR/bin/gap.sh"
fi

if ! test -f "$SYSINFO"; then
    AC_ERROR([Could not locate the GAP root directory;
        specify its location with './configure GAPDIR=DIR [CONFIGNAME=NAME)]'])
fi

GAPDIR=`cd $GAPDIR && pwd` # make path absolute
AC_SUBST(GAPDIR)

. "$SYSINFO"
TARGET="$GAParch"

echo checking target... "$TARGET"

XTARGET="`cnf/config.guess`-$CC-`echo $TARGET | sed 's/.*-//'`"
if test "$XTARGET" != "$GAParch_system"; then
   AC_WARN([The guessed target $XTARGET is not the gap target $GAParch_system. Cross your fingers])
fi
AC_SUBST(TARGET)

echo checking gap executable... "$GAPPROG"

if ! test -e "$GAPPROG"; then
    AC_WARN([Could not find GAP executable $GAPPROG; I won't compile the documentation])
fi
AC_SUBST(GAPPROG)

GAC="$GAPDIR/bin/$TARGET/gac"

echo checking gac compiler... $GAC

if ! test -e "$GAC"; then
    AC_ERROR([Could not find GAP compiler $GAC])
fi
AC_SUBST(GAC)

eval `grep '^c_[[a-z_]]*=' $GAC`

AC_SUBST(c_compiler)
AC_SUBST(c_options)
AC_SUBST(c_linker)
AC_SUBST(c_link_options)
AC_SUBST(c_libs)
AC_SUBST(c_dyn_options)
AC_SUBST(c_dyn_linker)
AC_SUBST(c_dyn_linking)
AC_SUBST(c_dynlibs)
AC_SUBST(c_addlibs)

EXTERN="\$(CURDIR)/bin/$TARGET/extern"

################################################################
# gmp configuration

eval `grep '^GMP_[[A-Z]]*=' "$GAPDIR/bin/$TARGET/config.log"`

AC_ARG_WITH(gmp,
 [  --with-gmp=<location>
    Location at which the GMP library, needed for MPFR, was installed.
    If the argument is omitted, then MPFR will use the version of gmp
    included in the GAP distribution.
 ],
 [GMP_CFLAGS="-I$withval/include"; GMP_LDFLAGS="-L$withval/lib"; GMP_LIBS="-lgmp"]
)

AC_ARG_WITH(gmp-include,
 [  --with-gmp-include=<location>
    Location at which the GMP include files were installed.],
 [GMP_CFLAGS="-I$withval"]
)

AC_ARG_WITH(gmp-lib,
 [  --with-gmp-lib=<location>
    Location at which the GMP library files were installed.],
 [GMP_LDFLAGS="-L$withval"; GMP_LIBS="-lgmp"]
)

# check gmp presence
CPPFLAGS="$CPPFLAGS $GMP_CFLAGS"
GACFLAGS="$GACFLAGS -p $GMP_CFLAGS"
AC_CHECK_HEADER(gmp.h,[],[AC_MSG_ERROR([library gmp not found. Specify its location using --with-gmp])],[])

AC_SUBST(GMP_CFLAGS)
AC_SUBST(GMP_LDFLAGS)
AC_SUBST(GMP_LIBS)

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
    	GACFLAGS="$GACFLAGS -p -I$GSLINCLUDE"
        LIBS="$LIBS -L$GSLLIB"
    fi
    AC_CHECK_HEADER(gsl/gsl_vector.h,[],[AC_MSG_ERROR([library gsl not found. Specify its location using --with-gsl])],[])
    AC_CHECK_LIB([gsl],[gsl_multiroot_fsolver_set],,
	[AC_ERROR([The GSL library could not be found. It is needed for IMG calculations.])],[-lgslcblas])
    AC_CHECK_LIB([gslcblas],[cblas_ctrmv],,
        [AC_ERROR([The GSL CBlas library could not be found. It is needed for IMG calculations.])])
fi

if test "$GSLLIB" != ""; then
    GACFLAGS="$GACFLAGS -L -L$GSLLIB -L -Wl,-rpath,$GSLLIB"
fi

GACFLAGS="$GACFLAGS -L -lgsl -L -lgslcblas"

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
AC_SUBST(GACFLAGS)
AC_SUBST(LIBS)
AC_SUBST(LIB_TARGET)

################################################################
# generate files

WITHGMP=""
INCLGMP=""
LINKGMP=""
if test "$GMPDIR" != ""; then
   WITHGMP="$WITHGMP --with-gmp=$GMPDIR"
fi
if test "$GMPINCLUDE" != ""; then
   WITHGMP="$WITHGMP --with-gmp-include=$GMPINCLUDE"
   INCLGMP="-I$GMPINCLUDE"
fi
if test "$GMPLIB" != ""; then
   WITHGMP="$WITHGMP --with-gmp-lib=$GMPLIB"
   LINKGMP="-L$GMPLIB"
fi

AC_SUBST(WITHGMP)
AC_SUBST(INCLGMP)
AC_SUBST(LINKGMP)

mkdir -p bin/$TARGET
CONFIG_STATUS=bin/$TARGET/config.status

AC_CONFIG_FILES([$MAKEFILE:cnf/Makefile.in])

if test "$MAKEFILE" != Makefile; then
    ln -sf "$MAKEFILE" Makefile
fi

AC_CONFIG_SUBDIRS(hurwitz)

AC_OUTPUT
