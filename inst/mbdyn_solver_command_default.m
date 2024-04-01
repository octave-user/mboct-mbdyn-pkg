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
## @deftypefn {Function File} @var{num_threads} = mbdyn_solver_command_default()
## Return the default mbdyn command (e.g. "mbdyn -C")
## @end deftypefn

function cmd = mbdyn_solver_command_default()
  if (nargout > 1 || nargin > 0)
    print_usage();
  endif

  persistent mbdyn_command = getenv("MBOCT_MBDYN_PKG_MBDYN_SOLVER_COMMAND");

  if (isempty(mbdyn_command))
    mbdyn_command = "mbdyn -C";
  endif

  persistent cnt = int32(0);

  cmd = sprintf(mbdyn_command, ++cnt);
endfunction
