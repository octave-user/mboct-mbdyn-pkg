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
## @deftypefn {Function File} @var{num_threads} = mbdyn_solver_num_threads_default()
## Return the default number of threads used by all tests
## @end deftypefn

function num_threads = mbdyn_solver_num_threads_default()
  if (nargout > 1 || nargin > 0)
    print_usage();
  endif
  num_threads = int32(str2num(getenv("MBD_NUM_THREADS")));
  if (isempty(num_threads) || num_threads < 1)
    num_threads = 1;
  endif
endfunction
