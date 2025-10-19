## Copyright (C) 2023(-2025) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{sol}] = mbdyn_post_load_output_sol(@var{output_file})
##
## Loads data for solid and beam elements from an MBDyn output file "<@var{output_file}>".
##
## @var{output_file} @dots{} Name of the MBDyn output file without extension.
##
## @var{mesh} @dots{} Finite element mesh data structure.
##
## @var{sol} @dots{} Finite element solution data structure.
##
## @seealso{fem_post_sol_external}
## @end deftypefn

function [mesh, sol] = mbdyn_post_load_output_sol(output_file)
  if (nargin ~= 1 || nargout ~= 2 || ~ischar(output_file))
    print_usage();
  endif

  [mesh, elem_types] = mbdyn_post_load_mesh_sol(output_file);

  mov_file = [output_file, ".mov"];
  sol_file = [output_file, ".sol"];

  [info, err, msg] = stat(mov_file);

  if (err ~= 0)
    error("file \"%s\" not found: %s", mov_file, msg);
  endif

  have_sol_file = false;

  [info, err, msg] = stat(sol_file);

  if (err == 0)
    have_sol_file = true;
  endif

  [sol.t, dt, niter, reserr, solerr, solconv, outputflag] = mbdyn_post_load_output_out(output_file, 1024, false);

  sol.t = sol.t(outputflag);

  [node_id, node_data] = mbdyn_post_load_output(mov_file, 6, [], numel(sol.t), 3, 1, false);

  sol.def = zeros(rows(mesh.nodes), 6, numel(sol.t));

  for i=1:numel(node_data)
    switch (mesh.node_type(mesh.node_id_transl(node_id(i))))
      case 123
        idx_def = 1:6;
        scale_def = [ones(1, 3), repmat(pi / 180, 1, 3)];
      case 100
        idx_def = 1:6;
        scale_def = ones(1, 6);
      case 999
        idx_def = 1:3;
        scale_def = ones(1, 3);
      otherwise
        error("invalid node type %d", mesh.node_type(mesh.node_id_transl(node_id(i))));
    endswitch

    for j=1:numel(idx_def)
      sol.def(mesh.node_id_transl(node_id(i)), j, :) = node_data{i}(:, idx_def(j)) * scale_def(j) - mesh.nodes(mesh.node_id_transl(node_id(i)), j);
    endfor
  endfor

  inum_elem_solid = 0;
  sol.stress.tau = struct();
  sol.strain.epsilon = struct();
  elem_stress_sum = elem_strain_sum = zeros(rows(mesh.nodes), 6, numel(sol.t));
  elem_stress_cnt = zeros(rows(mesh.nodes), 1);

  for l=1:numel(elem_types)
    switch (elem_types(l).elem_type)
      case {"beam2", "beam3"}
        continue;
    endswitch

    if (~isfield(mesh.elements, elem_types(l).elem_type))
      continue;
    endif

    elem_nodes = getfield(mesh.elements, elem_types(l).elem_type);

    if (have_sol_file)
      [elem_id, stress_strain] = mbdyn_post_load_output([output_file, ".sol"], columns(elem_nodes) * 6 * 2, inum_elem_solid + (1:rows(elem_nodes)), numel(sol.t));

      elem_stress = elem_strain = zeros(rows(elem_nodes), columns(elem_nodes), 6, numel(sol.t));

      for i=1:numel(elem_id)
        for j=1:columns(elem_nodes)
          for k=1:6
            elem_strain(elem_id(i) - inum_elem_solid, j, k, 1:rows(stress_strain{i})) = stress_strain{i}(:, (j - 1) * 12 + k);
            elem_stress(elem_id(i) - inum_elem_solid, j, k, 1:rows(stress_strain{i})) = stress_strain{i}(:, (j - 1) * 12 + 6 + k);
          endfor
        endfor
      endfor

      sol.strain.epsilon = setfield(sol.strain.epsilon, elem_types(l).elem_type, elem_strain);
      sol.stress.tau = setfield(sol.stress.tau, elem_types(l).elem_type, elem_stress);

      for k=1:columns(elem_nodes)
        for i=1:6
          for j=1:numel(sol.t)
            elem_strain_sum(elem_nodes(:, k), i, j) += elem_strain(:, k, i, j);
            elem_stress_sum(elem_nodes(:, k), i, j) += elem_stress(:, k, i, j);
          endfor
        endfor
        ++elem_stress_cnt(getfield(mesh.elements,elem_types(l).elem_type)(:, k));
      endfor

      inum_elem_solid += rows(elem_nodes);
    endif
  endfor

  if (have_sol_file)
    sol.stress.taum = sol.strain.epsilonm = struct();

    for l=1:numel(elem_types)
      if (~isfield(mesh.elements, elem_types(l).elem_type))
        continue;
      endif

      elem_nodes = getfield(mesh.elements, elem_types(l).elem_type);
      elem_stressm = elem_strainm = zeros(rows(elem_nodes), columns(elem_nodes), 6, numel(sol.t));

      for k=1:columns(elem_nodes)
        for i=1:6
          for j=1:numel(sol.t)
            elem_strainm(:, k, i, j) = elem_strain_sum(elem_nodes(:, k), i, j) ./ elem_stress_cnt(elem_nodes(:, k));
            elem_stressm(:, k, i, j) = elem_stress_sum(elem_nodes(:, k), i, j) ./ elem_stress_cnt(elem_nodes(:, k));
          endfor
        endfor
      endfor

      sol.strain.epsilonm = setfield(sol.strain.epsilonm, elem_types(l).elem_type, elem_strainm);
      sol.stress.taum = setfield(sol.stress.taum, elem_types(l).elem_type, elem_stressm);
    endfor
  endif
endfunction
