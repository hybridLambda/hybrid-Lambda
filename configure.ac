#AC_INIT([hybrid-Lambda],[beta-0.3],[sha.joe.zhu@gmail.com])
AC_INIT([hybrid-Lambda],esyscmd([sh -c "cat version | tr -d '\n'"]),[sha.joe.zhu@gmail.com])
#AC_INIT([hybrid-Lambda],esyscmd([sh -c "date \"+%Y%m%d\" | tr -d '\n'"]),[joezhu@well.ox.ac.uk])

#dnl AC_CONFIG_SRCDIR([simu.c])
#AC_CONFIG_HEADER([config.h])
# tar-ustar: make dist can include files with long names
AM_INIT_AUTOMAKE([1.9 tar-ustar subdir-objects]) 
#AM_INIT_AUTOMAKE([subdir-objects])
m4_ifndef([AM_SILENT_RULES], [m4_define([AM_SILENT_RULES],[])])
AM_SILENT_RULES([yes])

AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for libraries.
#AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
AC_LANG(C++) 
#AC_CHECK_HEADERS([iostream string sstream fstream vector iomanip valarray math.h stdio.h ctime])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

#AC_CHECK_FUNCS([floor])
#AC_CHECK_FUNCS([pow])
AC_HEADER_STDBOOL

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

