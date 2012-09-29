# autoconf for the Hurwitz problem programs

AC_PREREQ(2.61)
AC_INIT(fr/hurwitz, 0.0, laurent.bartholdi@gmail.com)
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_AUX_DIR(../cnf)

AC_PROG_CXX
AC_CANONICAL_TARGET

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

if test "$SYS_IS_DARWIN" == 1; then
    LDDOGLEG="-ldogleg"
    LDSUITESPARSE="-lcholmod -lblas -lcolamd -lamd -llapack -lsuitesparseconfig"
else
    LDDOGLEG="-Wl,-rpath,/usr/local/src/libdogleg -ldogleg"
    LDSUITESPARSE="-lcholmod"
fi

LDLEVMAR="-llevmar -llapack -lblas"

AC_SUBST(LDDOGLEG)
AC_SUBST(LDSUITESPARSE)
AC_SUBST(LDLEVMAR)

AC_OUTPUT(Makefile)
