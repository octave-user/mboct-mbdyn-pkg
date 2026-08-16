## Copyright (C) 2014(-2024) Reinhard <octave-user@a1.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} @var{test_type} = mbdyn_testsuite_test_type()
## Returns the desired test type.
## If @var{test_type} == "full" then the full test suite will be executed.
## However, if @var{test_type} == "quick" then only quick tests will be performed.
## @end deftypefn

function test_type = mbdyn_testsuite_test_type()
  if (nargout > 1 || nargin > 0)
    print_usage();
  endif
  test_type = getenv("MBD_TEST_TYPE");
  if (isempty(test_type))
    test_type = "quick";
  endif
endfunction
