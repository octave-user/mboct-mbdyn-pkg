dnl  Copyright (C) 2018(-2020) Reinhard <octave-user@a1.net>

dnl This program is free software; you can redistribute it and/or modify
dnl it under the terms of the GNU General Public License as published by
dnl the Free Software Foundation; either version 3 of the License, or
dnl (at your option) any later version.

dnl This program is distributed in the hope that it will be useful,
dnl but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl GNU General Public License for more details.

dnl You should have received a copy of the GNU General Public License
dnl along with this program; If not, see <http://www.gnu.org/licenses/>.

AC_PREREQ([2.69])
AC_INIT([mboct-mbdyn-pkg], [0.1.0])
AC_CONFIG_HEADERS([config.h])

AC_CONFIG_MACRO_DIRS([m4])

AM_INIT_AUTOMAKE

AC_CHECK_PROG(MKOCTFILE, mkoctfile, mkoctfile)

if test -z "$MKOCTFILE"; then
  AC_MSG_ERROR([mkoctfile not found], 1);
fi

AC_CHECK_PROG(OCTAVE_CONFIG, octave-config, octave-config)

if test -z "$OCTAVE_CONFIG"; then
  AC_MSG_ERROR([octave-config not found], 1);
fi

OCTINCLUDEDIR="${OCTINCLUDEDIR:-`${MKOCTFILE} -p INCFLAGS`}"
OCTLIBDIR="${OCTLIBDIR:-`${MKOCTFILE} -p OCTLIBDIR`}"
OCTLIBS="${OCTLIBS:-`${MKOCTFILE} -p LIBOCTAVE` `${MKOCTFILE} -p LIBOCTINTERP`}"
OCTLDFLAGS="${OCTLDFLAGS:--L${OCTLIBDIR} -Wl,-rpath=${OCTLIBDIR} -Wl,--no-undefined}"
CXX="${CXX:-`${MKOCTFILE} -p CXX`}"
CXXFLAGS="${CXXFLAGS} `${MKOCTFILE} -p CXXFLAGS`"

AC_PROG_CXX
AC_LANG(C++)

CPPFLAGS="${OCTINCLUDEDIR} ${CPPFLAGS}"
AC_CHECK_HEADERS(octave/oct.h, [], AC_MSG_ERROR([octave header files not found], 1))

TLDFLAGS="${LDFLAGS}"
TLIBS="${LIBS}"

dnl Force mkoctfile to link to Octave's libraries in order to avoid unresolved symbols at runtime
LDFLAGS="${LDFLAGS} ${OCTLDFLAGS}"
LIBS="${LIBS} ${OCTLIBS}"

AC_MSG_CHECKING([for octave libraries])
AC_TRY_LINK([
#include <octave/oct.h>
],[
octave_value_list ovl;
ovl.append(octave_value(1.5));
],[use_octlibs=yes],[use_octlibs=no])

AC_MSG_RESULT([${use_octlibs}])

if ! test "${use_octlibs}" = "yes"; then
  LDFLAGS="${TLDFLAGS}"
  LIBS="${TLIBS}"
fi

AC_CHECK_HEADERS(algorithm cassert cmath cstring fstream iostream map set sstream string, [], AC_MSG_ERROR([standard C++ headers not found], 1))

OF_OCTAVE_LIST_ALT_SYMS([

[dnl
  [is_bool_type],
  [islogical],
  [[octave_value ().islogical ();]],
  [OV_ISLOGICAL],
  [],
  []
]

],[oct-alt-includes.h])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
