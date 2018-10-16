#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT([dcm_to_4dfp],
	m4_esyscmd([(git describe --always --tag || echo unknown )| tr -d 'v\n']),
	[karchie@wustl.edu])
AC_CONFIG_SRCDIR([dcm_to_4dfp.c])
AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AM_PROG_CC_C_O

# Checks for libraries.
AC_ARG_WITH([TRX], [AS_HELP_STRING([--with-TRX],
		[give location of TRX code])],
		[AC_SUBST([TRX],[$with_TRX])],
		[AC_MSG_FAILURE([Must use --with-TRX to specify location of TRX directory])])
AC_CHECK_FILE(${with_TRX}/rec.o, [true],
	[AC_MSG_ERROR([cannot find required TRX object rec.o])])
AC_CHECK_FILE(${with_TRX}/Getifh.o, [true],
	[AC_MSG_ERROR([cannot find require TRX object Getifh.o])])
AC_CHECK_FILE(${with_TRX}/endianio.o, [true],
	[AC_MSG_ERROR([cannot find required TRX object endianio.o])])

# Checks for header files.
AC_HEADER_DIRENT
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([fcntl.h limits.h malloc.h netdb.h netinet/in.h stddef.h stdlib.h string.h sys/file.h sys/param.h sys/socket.h sys/time.h sys/timeb.h unistd.h regex.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_HEADER_TIME
AC_STRUCT_TM

AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(int)

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STAT
AC_FUNC_VPRINTF
AC_CHECK_FUNCS([gettimeofday memset mkdir re_comp regcomp sqrt strchr strerror strrchr strstr uname])

AC_SEARCH_LIBS(sqrt,m)

# needed for unit tests created by 'make check'
AC_SEARCH_LIBS([CU_basic_run_tests],[cunit])

if test -d .git ; then
AC_SUBST(GITPREREQ, ".git/HEAD .git/index")
fi

AC_CONFIG_FILES([Makefile include/Makefile])
AC_OUTPUT
