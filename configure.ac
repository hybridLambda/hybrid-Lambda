AC_INIT([hybrid-Lambda],[beta-0.3],[joezhu@well.ox.ac.uk])
#AC_INIT([hybrid-Lambda],esyscmd([sh -c "date \"+%Y%m%d\" | tr -d '\n'"]),[joezhu@well.ox.ac.uk])

#dnl AC_CONFIG_SRCDIR([simu.c])
#AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE([subdir-objects])
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


#
# Platform specific setup
#
#############################
AC_CANONICAL_HOST
# Check for which host we are on and setup a few things
# specifically based on the host
case $host_os in
  darwin* )
   AC_PROG_CXX([clang])
    export CXX=clang++
    export CXXFLAGS=-I/opt/local/include -stdlib=libc++ 
    AC_CHECK_HEADERS([boost/math/special_functions/gamma.hpp],[],[AC_MSG_ERROR(['boost' header files not found. Install from the port: sudo port install boost])])
         ;;
  linux*)
  AC_PROG_CXX
    AC_CHECK_HEADERS([boost/math/special_functions/gamma.hpp],[],[AC_MSG_ERROR(['boost' header files not found. Try downloading from :http://www.boost.org. Or: sudo apt-get install libboost-all-dev])])        
		;;
   *BSD*)
        # Do something specific for BSD
        ;;
    *)
        #Default Case
        AC_MSG_ERROR([Your platform is not currently supported])
        ;;
esac


# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC

AC_CHECK_FUNCS([floor])
AC_CHECK_FUNCS([pow])
AC_HEADER_STDBOOL

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

