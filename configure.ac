dnl This file is part of smithlab_cpp
dnl
dnl Copyright (C) 2010-2019: Andrew D. Smith and Meng Zhou
dnl
dnl Author: Andrew D. Smith and Meng Zhou
dnl
dnl This is free software: you can redistribute it and/or modify it
dnl under the terms of the GNU General Public License as published by
dnl the Free Software Foundation, either version 3 of the License, or
dnl (at your option) any later version.
dnl
dnl This software is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
dnl General Public License for more details.

AC_INIT([libsmithlab_cpp], [1.0], [andrewds@usc.edu],
        [], [https://github.com/smithlabcode])
dnl the config.h is not currently #included in the source, and only
dnl used to keep command lines short.
AC_CONFIG_HEADERS([config.h])
AC_PREREQ([2.63])  dnl 4-argument AC_CHECK_HEADER
AM_INIT_AUTOMAKE([foreign])

dnl Less output to the terminal. Use "--disable-silent-rules" when
dnl running ./configure to reverse this.
AM_SILENT_RULES([yes])

AC_CONFIG_SRCDIR([smithlab_utils.cpp])
AC_CONFIG_MACRO_DIR([m4])

dnl check for c++11
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

AM_PROG_AR
dnl Turn off shared libraries by default
LT_INIT([disable-shared])

dnl the line breaks in thie message are intended
hts_fail_msg="

Failed to locate HTSLib on your system. Please use the LDFLAGS and
CPPFLAGS variables to specify the directories where the HTSLib library
and headers can be found.
"

AC_ARG_ENABLE([hts],
  [AS_HELP_STRING([--enable-hts], [Enable HTSLib @<:@yes@:>@])],
  [enable_hts=yes], [enable_hts=no])
AS_IF([test "x$enable_hts" = "xyes"],
  [AC_CHECK_LIB([hts], [hts_version], [], [AC_MSG_FAILURE([$hts_fail_msg])])])
AM_CONDITIONAL([ENABLE_HTS], [test "x$enable_hts" = "xyes"])

dnl check for required libraries
AC_CHECK_LIB([z],[zlibVersion])
AC_CHECK_LIB([gslcblas],[cblas_dgemm])
AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
