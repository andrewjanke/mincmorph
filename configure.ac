# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

# The arguments are package name, and package version.
AC_INIT([mincmorph],[1.0])

# The argument is any source file distributed with the package.
AC_CONFIG_SRCDIR([mincmorph.c])

AM_INIT_AUTOMAKE

# The argument is the name of the generated header file.
# It is recommended to leave it as "config.h".
AC_CONFIG_HEADERS([config.h])

# This macro eases building the package; see m4/README.
smr_WITH_BUILD_PATH


# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL


# Checks for libraries.  See m4/README.
mni_REQUIRE_VOLUMEIO


# Finish up by writing output files.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
