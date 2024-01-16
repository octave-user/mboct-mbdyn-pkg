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
## @deftypefn {Function File} @var{R} = mbdyn_post_angles_to_rotation_mat(@var{node_id}, @var{res}, @var{log_dat})
## Convert rotation angles to rotation matrices for all structural node numbers in array <@var{node_id}>.
##
## @var{node_id} @dots{} Array of node numbers.
##
## @var{res} @dots{} Data structure containing the output from mbdyn_post_load_output_struct.
##
## @var{log_dat} @dots{} Return value from mbdyn_post_load_log.
##
## @end deftypefn

function R = mbdyn_post_angles_to_rotation_mat(node_id, res, log_dat, idx_t)
  if (nargin < 3 || nargin > 4 || nargout > 1)
    print_usage();
  endif

  if (nargin < 4)
    idx_t = 1:numel(res.t);
  endif

  R = cell(1, numel(node_id));

  for i=1:numel(node_id)
    node_idx = find(node_id(i) == res.node_id);
    node_idx_log = find([log_dat.nodes.label] == node_id(i));

    if (numel(node_idx_log) ~= 1)
      error("node_id=%d not found in log_dat", node_id(i));
    endif

    orientation_description = log_dat.nodes(node_idx_log).orientation_description;

    if (isempty(node_idx))
      warning("node_id=%d not found in result", node_id(i));
      Rij = repmat(log_dat.nodes(node_idx_log).R0, [1, 1, numel(idx_t)]);
    else
      Phi = res.trajectory{node_idx}(idx_t, 4:6).';

      switch(orientation_description)
        case "euler123"
          Rij = euler123_to_rotation_matrix(Phi);
        case "euler313"
          Rij  = euler313_to_rotation_matrix(Phi);
        case "euler321"
          Rij  = euler321_to_rotation_matrix(Phi);
        case "phi"
          Rij  = rotation_vector_to_rotation_matrix(Phi);
        otherwise
          error("orientation_description = \"%s\" not supported", log_dat.nodes(node_idx).orientation_description);
      endswitch
    endif

    R{i} = Rij;
  endfor
endfunction
