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
    
    if (numel(node_idx) ~= 1)
      error("node_id=%d not found in result", node_id(i));
    endif

    node_idx_log = find([log_dat.nodes.label] == node_id(i));

    if (numel(node_idx_log) ~= 1)
      error("node_id=%d not found in log_dat", node_id(i));
    endif
    
    orientation_description = log_dat.nodes(node_idx_log).orientation_description;
    
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

    R{i} = Rij;
  endfor
endfunction

%!test
%! state = rand("state");
%! unwind_protect
%! rand("seed", 0);
%! N = 10;
%! orient = {"euler123", "euler321", "phi"};
%! func = {@rotation_matrix_to_euler123, ...
%!         @rotation_matrix_to_euler321, ...
%!         @rotation_matrix_to_rotation_vector};
%! res.node_id = 1:numel(orient);
%! for i=1:numel(orient)
%!   res.t = 1:N;
%!   res.trajectory{i} = [zeros(N, 3), (2 * rand(N, 3) - 1) * 0.5 * pi];
%!   log_dat.nodes(i).label = i;
%!   log_dat.nodes(i).orientation_description = orient{i};
%! endfor
%! R = mbdyn_post_angles_to_rotation_mat(res.node_id, res, log_dat);
%! for i=1:numel(R)
%!   Phi = feval(func{i}, R{i}).';
%!   assert(Phi, res.trajectory{i}(:, 4:6), sqrt(eps) * pi);
%! endfor
%! unwind_protect_cleanup
%! rand("state", state);
%! end_unwind_protect

%!test
%! state = rand("state");
%! unwind_protect
%! rand("seed", 0);
%! N = 10;
%! orient = {"euler123", "euler321", "phi"};
%! func = {@rotation_matrix_to_euler123, ...
%!         @rotation_matrix_to_euler321, ...
%!         @rotation_matrix_to_rotation_vector};
%! res.node_id = 1:numel(orient);
%! for i=1:numel(orient)
%!   res.t = 1:N;
%!   res.trajectory{i} = [zeros(N, 3), (2 * rand(N, 3) - 1) * 0.5 * pi];
%!   log_dat.nodes(i).label = i;
%!   log_dat.nodes(i).orientation_description = orient{i};
%! endfor
%! R = mbdyn_post_angles_to_rotation_mat(res.node_id, res, log_dat, 1:2:N);
%! for i=1:numel(R)
%!   Phi = feval(func{i}, R{i}).';
%!   assert(Phi, res.trajectory{i}(1:2:N, 4:6), sqrt(eps) * pi);
%! endfor
%! unwind_protect_cleanup
%! rand("state", state);
%! end_unwind_protect
