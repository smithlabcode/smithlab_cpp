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
AC_PREREQ([2.63])  dnl 4-argument AC_CHECK_HEADER
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_CONFIG_SRCDIR([smithlab_os.cpp])
AC_CONFIG_MACRO_DIR([m4])

dnl check for c++11
AC_LANG(C++)
AC_PROG_CXX
AX_CXX_COMPILE_STDCXX_11([noext], [mandatory])

AM_PROG_AR
LT_INIT

dnl check for HTSLib
AX_WITH_HTSLIB
AS_IF([test "x$ax_cv_htslib" != "xyes"],
      [AC_MSG_WARN([HTSlib not found.

HTSLib is optional for using libsmithlab_cpp. Please use
--with-htslib=DIR to specify htslib installation if you have it
installed locally.
])])
AM_CONDITIONAL([HAS_HTSLIB], [test "x$ax_cv_htslib" = "xyes"])

dnl check for required libraries
AC_CHECK_LIB([z],[zlibVersion])
AC_CHECK_LIB([gslcblas],[main])
AC_CHECK_LIB([gsl],[main])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
