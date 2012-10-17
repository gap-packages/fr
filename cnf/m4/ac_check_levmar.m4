# check for levmar library
# sets the LEVMAR_CFLAGS, LEVMAR_LDFLAGS and LEVMAR_LIBS appropriately

AC_DEFUN([AC_CHECK_LEVMAR],[

AC_ARG_WITH(levmar,
 [  --with-levmar=<location>
    Location at which the levmar library, needed for layout, was installed.],
 [LEVMAR_CFLAGS="-I$withval/include"; LEVMAR_LDFLAGS="-L$withval/lib"]
)

AC_ARG_WITH(levmar-include,
 [  --with-levmar-include=<location>
    Location at which the levmar include files were installed.],
 [LEVMAR_CFLAGS="-I$withval"]
)

AC_ARG_WITH(levmar-lib,
 [  --with-levmar-lib=<location>
    Location at which the levmar library files were installed.
 ],
 [LEVMAR_LDFLAGS="-L$withval"]
)

LEVMAR_LIBS="-llevmar -llapack -lblas"

AC_LANG_PUSH([C])

lm_CFLAGS=$CFLAGS
CFLAGS="$CFLAGS $LEVMAR_CFLAGS"
AC_CHECK_HEADER(levmar.h,,AC_MSG_ERROR([levmar.h not found. Specify its location using --with-levmar.
The package may be downloaded from http://www.ics.forth.gr/~lourakis/levmar/]))
CFLAGS=$lm_CFLAGS

lm_LDFLAGS=$LDFLAGS
LDFLAGS="$LDFLAGS $LEVMAR_LDFLAGS"
AC_CHECK_LIB(levmar,dlevmar_dif,,AC_MSG_ERROR([liblevmar not found. Specify its location using --with-levmar.]),[-llapack -lblas])
LDFLAGS=$lm_LDFLAGS

AC_LANG_POP([C])

AC_SUBST(LEVMAR_CFLAGS)
AC_SUBST(LEVMAR_LDFLAGS)
AC_SUBST(LEVMAR_LIBS)
])
