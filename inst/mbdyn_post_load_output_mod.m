## Copyright (C) 2018(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{elem_id}, @var{q}, @var{qdot}, @var{qddot}] = mbdyn_post_load_output_mod(@var{mbdyn_filename}, @var{num_steps})
## Load modal element data from MBDyn output file "<@var{mbdyn_filename}>.mod".
##
## @var{elem_id} @dots{} Modal joint number.
##
## @var{q} @dots{} Modal displacements.
##
## @var{qdot} @dots{} Modal velocities.
##
## @var{qddot} @dots{} Modal accelerations.
##
## @var{num_steps} @dots{} Number of output time steps returned from mbdyn_post_load_output_out.
##
## @end deftypefn

function [elem_id, q, qdot, qddot] = mbdyn_post_load_output_mod(mbdyn_filename, num_steps)
  if (nargin ~= 2 || nargout > 4)
    print_usage();
  endif

  mod_filename = mbdyn_post_output_filename(mbdyn_filename, ".mod");

  data = load("-ascii", mod_filename);

  data_elem_id = floor(data(:, 1));
  elem_id = unique(data_elem_id);

  for i=1:numel(elem_id)
    elem_idx = find(data_elem_id == elem_id(i));

    num_modes = numel(elem_idx) / num_steps;

    if (nargout >= 2)
      q{i} = reshape(data(elem_idx, 2), num_modes, num_steps);
    endif

    if (nargout >= 3)
      qdot{i} = reshape(data(elem_idx, 3), num_modes, num_steps);
    endif

    if (nargout >= 4)
      qddot{i} = reshape(data(elem_idx, 4), num_modes, num_steps);
    endif
  endfor
endfunction
