# autoconf for the Hurwitz problem programs

AC_PREREQ(2.61)
AC_INIT(fr/hurwitz, 0.0, laurent.bartholdi@gmail.com)
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_AUX_DIR(../cnf)
AC_CONFIG_MACRO_DIR(../cnf/m4)

AC_PROG_CC
AC_PROG_CXX
AC_CANONICAL_TARGET

################################################################
# check for GAP
AC_FIND_GAP

mkdir -p ../bin/$GAPARCH

################################################################
# check for cholmod library

AC_ARG_WITH(cholmod,
 [  --with-cholmod=<location>
    Location at which the cholmod library, needed for layout, was installed.],
 [CHOLMODINCLUDE="-I$withval/include"; CHOLMODLIB="$withval/lib"]
)

AC_ARG_WITH(cholmod-include,
 [  --with-cholmod-include=<location>
    Location at which the cholmod include files were installed.],
 [CHOLMODINCLUDE="-I$withval"]
)

AC_ARG_WITH(cholmod-lib,
 [  --with-cholmod-lib=<location>
    Location at which the cholmod library files were installed.
 ],
 [CHOLMODLIB="$withval"]
)

AC_CHECK_HEADER(suitesparse/cholmod.h,[CHOLMODINCLUDE="-I/usr/include/suitesparse"])

if test "$CHOLMODINCLUDE" != ""; then
    CPPFLAGS="$CPPFLAGS $CHOLMODINCLUDE"
fi

AC_CHECK_HEADER(cholmod.h,,AC_MSG_ERROR([cholmod.h not found. Specify its location using --with-cholmod.]))

if test "$CHOLMODLIB" != ""; then
    LDFLAGS="$LDFLAGS -L$CHOLMODLIB"
fi
AC_CHECK_LIB(cholmod,cholmod_allocate_triplet,,AC_MSG_ERROR([libcholmod not found. Specify its location using --with-cholmod.]))

# check for dogleg library

AC_ARG_WITH(libdogleg,
 [  --with-libdogleg=<location>
    Location at which the libdogleg library, needed for layout, was installed.],
 [LIBDOGLEGINCLUDE="-I$withval/include"; LIBDOGLEGDIR="$withval/lib"]
)

AC_ARG_WITH(libdogleg-include,
 [  --with-libdogleg-include=<location>
    Location at which the libdogleg include files were installed.],
 [LIBDOGLEGINCLUDE="-I$withval"]
)

AC_ARG_WITH(libdogleg-lib,
 [  --with-libdogleg-lib=<location>
    Location at which the libdogleg library files were installed.
 ],
 [LIBDOGLEGLIB="$withval"]
)

if test "$LIBDOGLEGINCLUDE" != ""; then
    CPPFLAGS="$CPPFLAGS $LIBDOGLEGINCLUDE"
fi
AC_CHECK_HEADER(dogleg.h,,AC_MSG_ERROR([dogleg.h not found. Specify its location using --with-libdogleg.
The package may be downloaded from https://github.com/Oblong/libdogleg.git]))

if test "$LIBDOGLEGLIB" != ""; then
    LDFLAGS="$LDFLAGS -L$LIBDOGLEGLIB"
fi
AC_CHECK_LIB(dogleg,dogleg_optimize,,AC_MSG_ERROR([libdogleg not found. Specify its location using --with-libdogleg.]))

# check for levmar library

AC_CHECK_LEVMAR

# check for standard C++ stuff

AC_CHECK_HEADERS(malloc/malloc.h)
AC_CHECK_HEADERS(malloc.h)

AC_CHECK_LIB(m,sincos,AC_DEFINE(HAVE_SINCOS,1,do we have sincos?))

AC_SYS_IS_CYGWIN
AC_SYS_IS_DARWIN

# hack to specify correct location of DLL, for Linux

LDLAYOUT="-lcholmod -ldogleg"

if test "$CHOLMODLIB" != ""; then
    LDLAYOUT="-LCHOLMODLIB $LDLAYOUT"
fi
if test "$LIBDOGLEGLIB" != ""; then
    LDLAYOUT="-L$LIBDOGLEGLIB -Wl,-rpath,$LIBDOGLEGLIB $LDLAYOUT"
fi

AC_SUBST(CHOLMODINCLUDE)
AC_SUBST(LIBDOGLEGINCLUDE)

AC_SUBST(LDLAYOUT)

AC_OUTPUT(Makefile)
