## Copyright (C) 2014(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{node_id}, @var{y}, @var{ydot}] = mbdyn_post_load_output_abs(@var{output_file})
##
## Load output data from abstract nodes from MBDyn output file @var{output_file}.
##
## @var{node_id} @dots{} Array of node numbers.
##
## @var{y} @dots{} Cell array of abstract node values.
##
## @var{ydot} @dots{} Cell array of abstract node time derivatives.
##
## @end deftypefn

function [node_id, y, ydot] = mbdyn_post_load_output_abs(abs_filename, filter_node_id, append_rows)
  if (nargin < 1 || nargin > 3 || nargout > 3)
    print_usage();
  endif

  if (nargin < 2)
    filter_node_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  abs_filename = mbdyn_post_output_filename(abs_filename, ".abs");

  [node_id, data] = mbdyn_post_load_output(abs_filename, 2, filter_node_id, append_rows);

  for i=1:numel(data)
    y{i} = data{i}(:,1);

    if (nargout >= 3)
      if (columns(data{i}) >= 2)
        ydot{i} = data{i}(:, 2);
      else
        ydot{i} = [];
      endif
    endif
  endfor
endfunction
