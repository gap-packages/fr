# autoconf for the Hurwitz problem programs

AC_PREREQ(2.61)
AC_INIT(fr/hurwitz, 0.0, laurent.bartholdi@gmail.com)
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_AUX_DIR(../cnf)

AC_PROG_CXX
AC_CANONICAL_TARGET

if ! test -d ../bin/$TARGET; then
    AC_ERROR([$TARGET should be the directory under bin/ to contain the executables; specify it with TARGET=...])
fi

# check for cholmod library

AC_ARG_WITH(cholmod,
 [  --with-cholmod=<location>
    Location at which the cholmod library, needed for layout, was installed.
 ],
 [CHOLMODINCLUDE="-I$withval/include"; CHOLMODLIB="$withval/lib"]
)

AC_ARG_WITH(cholmod-include,
 [  --with-cholmod-include=<location>
    Location at which the cholmod include files were installed.],
 [CHOLMODINCLUDE="-I$withval"]
)

AC_ARG_WITH(cholmod-lib,
 [  --with-cholmod-lib=<location>
    Location at which the cholmod library files were installed.],
 [CHOLMODLIB="$withval"]
)

if test "$CHOLMODINCLUDE" == ""; then
    AC_CHECK_HEADER(cholmod.h,,[CHOLMODINCLUDE="-I/usr/include/sparsesuite"])
fi

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
    Location at which the libdogleg library, needed for layout, was installed.
 ],
 [LIBDOGLEGINCLUDE="-I$withval/include"; LIBDOGLEGDIR="$withval/lib"]
)

AC_ARG_WITH(libdogleg-include,
 [  --with-libdogleg-include=<location>
    Location at which the libdogleg include files were installed.],
 [LIBDOGLEGINCLUDE="-I$withval"]
)

AC_ARG_WITH(libdogleg-lib,
 [  --with-libdogleg-lib=<location>
    Location at which the libdogleg library files were installed.],
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

AC_ARG_WITH(levmar,
 [  --with-levmar=<location>
    Location at which the levmar library, needed for layout, was installed.
 ],
 [LEVMARINCLUDE="-I$withval/include"; LEVMARLIB="$withval/lib"]
)

AC_ARG_WITH(levmar-include,
 [  --with-levmar-include=<location>
    Location at which the levmar include files were installed.],
 [LEVMARINCLUDE="-I$withval"]
)

AC_ARG_WITH(levmar-lib,
 [  --with-levmar-lib=<location>
    Location at which the levmar library files were installed.],
 [LEVMARLIB="$withval"]
)

if test "$LEVMARINCLUDE" != ""; then
    CPPFLAGS="$CPPFLAGS $LEVMARINCLUDE"
fi
AC_CHECK_HEADER(levmar.h,,AC_MSG_ERROR([levmar.h not found. Specify its location using --with-levmar.
The package may be downloaded from http://www.ics.forth.gr/~lourakis/levmar/]))

if test "$LEVMARLIB" != ""; then
    LDFLAGS="$LDFLAGS -L$LEVMARLIB"
fi
AC_CHECK_LIB(levmar,dlevmar_dif,,AC_MSG_ERROR([liblevmar not found. Specify its location using --with-levmar.]),[-llapack -lblas])

# check for standard C++ stuff

AC_CHECK_HEADERS(malloc/malloc.h)
AC_CHECK_HEADERS(malloc.h)

AC_CHECK_LIB(m,sincos,AC_DEFINE(HAVE_SINCOS,1,do we have sincos?))

case "$host" in
    *-darwin*)
        AC_DEFINE(SYS_IS_DARWIN, 1, are we on DARWIN?)
        SYS_IS_DARWIN=1
        ;;
    *)
        AC_DEFINE(SYS_IS_DARWIN, 0, are we on DARWIN?)
        SYS_IS_DARWIN=0
        ;;
esac

# hack to specify correct location of DLL, for Linux

LDLAYOUT="-lcholmod -ldogleg"
LDHSOLVE="-llevmar -llapack -lblas"

if test "$CHOLMODLIB" != ""; then
    LDLAYOUT="-LCHOLMODLIB $LDLAYOUT"
fi
if test "$LIBDOGLEGLIB" != ""; then
    LDLAYOUT="-L$LIBDOGLEGLIB -Wl,-rpath,$LIBDOGLEGLIB $LDLAYOUT"
fi

if test "$LEVMARLIB" != ""; then
    LDHSOLVE="-L$LEVMARLIB $LDHSOLVE"
fi

AC_SUBST(CHOLMODINCLUDE)
AC_SUBST(LIBDOGLEGINCLUDE)
AC_SUBST(LEVMARINCLUDE)

AC_SUBST(LDLAYOUT)
AC_SUBST(LDHSOLVE)

AC_SUBST(TARGET)

AC_OUTPUT(Makefile)
