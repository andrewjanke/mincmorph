# Require autoconf 2.57 or newer.
AC_PREREQ([2.57])

# The arguments are package name, and package version.
AC_INIT([mincmorph],[1.5.1],[Andrew L Janke <a.janke@gmail.com>])
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

# for clean MINC2.0 volume_io
AC_DEFINE_UNQUOTED(VIO_PREFIX_NAMES, 1, [Play nice with the other kids volume_io])

# Finish up by writing output files.
AC_CONFIG_FILES([Makefile])
AC_OUTPUT(epm-header)
