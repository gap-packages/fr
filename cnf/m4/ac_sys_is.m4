# Finds if DOS/MacOSX install
AC_DEFUN([AC_SYS_IS_CYGWIN],[
AC_CYGWIN
AM_CONDITIONAL([SYS_IS_CYGWIN], [test "$CYGWIN" = "yes"])
if test "$CYGWIN" = "yes"; then
  AC_DEFINE(SYS_IS_CYGWIN32, 1, are we on CYGWIN?)
else
  AC_DEFINE(SYS_IS_CYGWIN32, 0, are we on CYGWIN?)
fi
])

AC_DEFUN([AC_SYS_IS_DARWIN],[
case "$host" in
  *-darwin* )
  AC_DEFINE(SYS_IS_DARWIN, 1, are we on DARWIN?)
  ;;
  * )
  AC_DEFINE(SYS_IS_DARWIN, 0, are we on DARWIN?)
  ;;
esac
])
