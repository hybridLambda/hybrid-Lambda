#AC_INIT([hybrid-Lambda],esyscmd([sh -c "date \"+%Y%m%d\" | tr -d '\n'"]),[joe.zhu@pg.canterbury.ac.nz])
AC_INIT([hybrid-Lambda],[beta-0.2],[joezhu@well.ox.ac.uk])
dnl AC_CONFIG_SRCDIR([simu.c])
#AC_CONFIG_SRCDIR([main.cpp])
#AC_CONFIG_HEADER([config.h])

AM_INIT_AUTOMAKE
AC_PROG_INSTALL
AC_PREREQ
AC_CANONICAL_HOST

# Checks for programs.
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for libraries.
#AC_CHECK_LIB([m], [log])

# Checks for header files.
AC_HEADER_STDC
#AC_CHECK_HEADERS([float.h stdlib.h string.h unistd.h])
AC_LANG(C++) 
AC_CHECK_HEADERS([iostream string sstream fstream vector iomanip valarray math.h stdio.h ctime])

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
   export CXXFLAGS=-I/opt/local/include
    AC_CHECK_HEADERS([boost/math/special_functions/gamma.hpp],[],[AC_MSG_ERROR(['boost' header files not found. Install from the port: sudo port install boost])])
         ;;
  linux*)
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

DX_INIT_DOXYGEN(hybrid-Lambda, Doxyfile, html/)
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)


# Check for presence of pdfLaTeX
AC_CHECK_PROG(PDFLATEX, pdflatex, pdflatex)
AC_CHECK_PROG(BIBTEX, bibtex,bibtex)
if test -z "$PDFLATEX"; then
  AC_MSG_WARN([Unable to create PDF version of the user manual.])
fi

AM_CONDITIONAL([HAVE_PDFLATEX], test -n "$PDFLATEX")

AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile])
AC_OUTPUT

