#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(biozim, 0.1, mail@sebastianbauer.info)
AC_CONFIG_SRCDIR([src/biozim.cpp])
AM_INIT_AUTOMAKE
AC_CONFIG_HEADER([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lcmph':
AC_CHECK_LIB([cmph], [main])
# FIXME: Replace `main' with a function in `-lm':
AC_CHECK_LIB([m], [main])
# FIXME: Replace `main' with a function in `-lsundials_cvode':
AC_CHECK_LIB([sundials_cvode], [main])
# FIXME: Replace `main' with a function in `-lsundials_nvecserial':
AC_CHECK_LIB([sundials_nvecserial], [main])

AC_CHECK_LIB(gmp, __gmpz_init)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stddef.h stdint.h stdlib.h string.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_TIME
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_FUNC_STRTOD
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([gettimeofday memset pow strchr strdup strtol])

AC_PROG_SWIG

dnl ******************************
dnl openmp checking
dnl ******************************

AX_OPENMP
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CFLAGS $OPENMP_CXXFLAGS"
LDFLAGS="$LDFLAGS $OPENMP_CFLAGS"

# TODO: Add better sundials stuff
SUNDIALS_CFLAGS=`sundials-config -m cvode -t s -l c -s cflags`
SUNDIALS_LIBS=`sundials-config -m cvode -t s -l c -s libs`

AC_SUBST(SUNDIALS_CFLAGS)
AC_SUBST(SUNDIALS_LIBS)

dnl ******************************
dnl libsbml checking
dnl ******************************

PKG_CHECK_MODULES(LIBSBML, [libsbml >= 3.3.2])

dnl ******************************
dnl liblua5.1 checking
dnl ******************************

AC_ARG_WITH(lua5.1,
		AC_HELP_STRING([--with-lua5.1],[include support for lua5.1 (default is auto)]),
		[with_lua=$withval],[with_lua=auto])

if test "x$with_lua" != "xno" ; then
	PKG_CHECK_MODULES(LUA5_1, [lua5.1], enable_lua=yes, enable_lua=no)
	if test "x$enable_lua" != "xyes" ; then
		if test "x$with_lua" == "xyes"; then
			dnl For the error message
			PKG_CHECK_MODULES(LUA5_1, [lua5.1])
		fi
	fi
else
	dnl asked to ignore lua
	enable_lua=no
fi
if test "x$enable_lua" == "xyes" ; then
	AC_DEFINE(HAVE_LIBLUA, 1, [Enable support for liblua 5.1])
fi

dnl ******************************

AC_CONFIG_FILES([Makefile])
AC_CONFIG_FILES([src/Makefile])
AC_OUTPUT

echo "
Configuration: 
        Source code location:     ${srcdir}
        Destination path prefix:  ${prefix}
        C Compiler:               ${CC}
        C++ Compiler:             ${CXX}
        CPPFLAGS:                 ${CPPFLAGS}
        CXXFLAGS:                 ${CXXFLAGS}
        CFLAGS:                   ${CFLAGS}
        LDFLAGS:                  ${LDFLAGS}

        Enable Lua 5.1 support:   ${enable_lua}    
"