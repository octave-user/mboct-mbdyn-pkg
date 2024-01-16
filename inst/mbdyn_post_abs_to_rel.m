## Copyright (C) 2011(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} mbdyn_post_abs_to_rel(@var{ref_node}, @var{abs_mov_file}, @var{rel_mov_file}, @var{ref_only}, @var{input_mode}, @var{output_mode})
##
## Converts absolute node positions an velocities stored in file "<@var{abs_mov_file}>" to relative positions and velocities with respect to node <@var{ref_node}>.
## The result will be stored in file <@var{rel_mov_file}>.
##
## @end deftypefn

function mbdyn_post_abs_to_rel(ref_node, abs_mov_file, rel_mov_file, ref_only, input_mode, output_mode)
  if (nargin < 3 || nargin > 6)
    print_usage();
  endif

  if (nargin < 4)
    ref_only = 0;
  endif

  if (nargin < 5)
    input_mode = "euler123";
  endif

  if (nargin < 6)
    output_mode = "euler123";
  endif

  if (~isscalar(ref_node))
    print_usage();
  endif

  if (~ischar(abs_mov_file))
    print_usage();
  endif

  if (~ischar(rel_mov_file))
    print_usage();
  endif

  mbdyn_path_init();

  command = sprintf("exec awk -f abs2rel.awk -v RefNode=%d -v RefOnly=%d -v InputMode=%s -v OutputMode=%s '%s.mov' > '%s.mov'", ...
                    ref_node, ...
                    ref_only, ...
                    input_mode, ...
                    output_mode, ...
                    mbdyn_convert_path(abs_mov_file), ...
                    mbdyn_convert_path(rel_mov_file));

  rc = shell(command);

  if (0 ~= rc)
    error("command \"%s\" failed with status %d", command, rc);
  endif
endfunction
