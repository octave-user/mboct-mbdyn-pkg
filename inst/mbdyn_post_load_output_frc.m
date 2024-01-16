## Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{force_id}, @var{force_node_id}, @var{force}, @var{arm}] = mbdyn_post_load_output_frc(@var{mbdyn_filename}, @var{filter_force_id}, @var{append_rows})
##
## Loads the MBDyn .frc file named "<@var{mbdyn_filename}.frc>".
##
## @var{filter_node_id} @dots{} Load only the data of the corresponding nodes if present.
## If empty, the data of all nodes is loaded.
##
## @var{append_rows} @dots{} Hint for memory reallocation.
##
## @var{force_id} @dots{} Structural force element number
##
## @var{force_node_id} @dots{} Structural node number
##
## @var{force} @dots{} Structural force values
##
## @var{arm} @dots{} The arm of the force, in the global frame (i.e. referred to point @{0, 0, 0@} and oriented as the global frame).
##
## @end deftypefn

function [force_id, force_node_id, force, arm] = mbdyn_post_load_output_frc(mbdyn_filename, filter_force_id, append_rows, auto_resize_rows)
  if (nargin < 1 || nargout > 4)
    print_usage();
  endif

  if (nargin < 2)
    filter_force_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  if (nargin < 4)
    auto_resize_rows = true;
  endif

  frc_filename = mbdyn_post_output_filename(mbdyn_filename, ".frc");

  column_count = 0;

  if (nargout >= 2)
    ++column_count;
  endif

  if (nargout >= 3)
    column_count += 3;
  endif

  if (nargout >= 4)
    column_count += 3;
  endif

  [force_id, data] = mbdyn_post_load_output(frc_filename, column_count, filter_force_id, append_rows, 4, 1, auto_resize_rows);

  if (numel(data) > 0)
    for i=1:numel(data)
      if (nargout >= 2)
        if (columns(data{i}) >= 8)
          force_node_id{i} = data{i}(:, [1, 8]);
        else
          force_node_id{i} = data{i}(:, 1);
        endif
      endif

      if (nargout >= 3)
        if (columns(data{i}) >= 11)
          force{i} = data{i}(:, [2:4, 9:11]);
        elseif (columns(data{i}) >= 4)
          force{i} = data{i}(:, 2:4);
        elseif (columns(data{i}) >= 2)
          force{i} = data{i}(:, 2:end);
        else
          force{i} = [];
        endif
      endif

      if (nargout >= 4)
        if (columns(data{i}) >= 14)
          arm{i} = data{i}(:, [5:7, 12:14]);
        elseif (columns(data{i}) >= 7)
          arm{i} = data{i}(:, 5:7);
        else
          arm{i} = [];
        endif
      endif
    endfor
  else
    force = {};
    force_node_id = {};
    arm = {};
  endif
endfunction
