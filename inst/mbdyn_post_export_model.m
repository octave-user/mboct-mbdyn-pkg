## Copyright (C) 2020(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{sol}] = mbdyn_post_export_model(@var{res}, @var{idx_t})
##
## Convert the output data from MBDyn to Finite Element model data suiteable for processing with mboct-fem-pkg
##
## @end deftypefn

function [mesh, sol] = mbdyn_post_export_model(res, idx_t)
  if (nargin < 1 || nargin > 2 || nargout > 2)
    print_usage();
  endif

  if (nargin < 2)
    idx_t = 1:numel(res.t);
  endif

  R = mbdyn_post_angles_to_rotation_mat([res.log_dat.nodes.label], res, res.log_dat, idx_t);

  elem_types_fe = {"beam2", "beam2", "beam3"};
  elem_types_mb = {"joints", "beams2", "beams3"};
  elem_node_num = int32([2, 2, 3]);
  elem_node_order = {[1:2], [1:2], [1,3,2]};

  num_aux_nodes = int32(0);

  for i=1:numel(elem_types_mb)
    num_aux_nodes += elem_node_num(i) * numel(getfield(res.log_dat, elem_types_mb{i}));
  endfor

  mesh.nodes = zeros(numel(res.log_dat.nodes) + num_aux_nodes, 6);

  mesh.nodes(1:numel(res.log_dat.nodes), :) = [[res.log_dat.nodes.X0].', [res.log_dat.nodes.Phi0].'];

  iaux_node = zeros(1, numel(elem_types_mb), "int32");

  mesh.elements = struct();

  for k=1:numel(elem_types_mb)
    elem_k = getfield(res.log_dat, elem_types_mb{k});
    elnodes = struct("nodes", cell(1, numel(elem_k)));

    for i=1:numel(elem_k)
      Xi = zeros(3, numel(elem_k(i).nodes));

      for j=1:numel(elem_k(i).nodes)
        node_idx = find([res.log_dat.nodes.label] == elem_k(i).nodes(j).label);

        if (numel(node_idx) ~= 1)
          error("invalid node label %d", elem_k(i).nodes(j).label);
        endif

        Xi(:, j) = res.log_dat.nodes(node_idx).X0 + res.log_dat.nodes(node_idx).R0 * elem_k(i).nodes(j).offset(:);
      endfor

      node_id = numel(res.log_dat.nodes) + sum(iaux_node(1:k)) + (1:columns(Xi));
      mesh.nodes(node_id, 1:3) = Xi.';
      elnodes(i).nodes = node_id(elem_node_order{k});
      iaux_node(k) += columns(Xi);
    endfor

    if (isfield(mesh.elements, elem_types_fe{k}))
      elnodesprev = getfield(mesh.elements, elem_types_fe{k});
      elnodesprev(end + (1:numel(elnodes))) = elnodes;
    else
      elnodesprev = elnodes;
    endif

    mesh.elements = setfield(mesh.elements, elem_types_fe{k}, elnodesprev);
  endfor

  iaux_node = zeros(1, numel(elem_types_mb), "int32");

  for k=1:numel(elem_types_mb)
    elem_k = getfield(res.log_dat, elem_types_mb{k});
    elnodes = struct("nodes", cell(1, 2 * numel(elem_k)));
    ielem = int32(0);

    for i=1:numel(elem_k)
      node_id = numel(res.log_dat.nodes) + sum(iaux_node(1:k)) + (1:numel(elem_k(i).nodes));

      for j=[1, numel(elem_k(i).nodes)]
        node_idx = find([res.log_dat.nodes.label] == elem_k(i).nodes(j).label);

        if (numel(node_idx) ~= 1)
          error("invalid node label %d", elem_k(i).nodes(j).label);
        endif

        elnodes(++ielem).nodes = int32([node_id(j), node_idx]);
      endfor

      iaux_node(k) += numel(elem_k(i).nodes);
    endfor

    mesh.elements.beam2(end + (1:numel(elnodes))) = elnodes;
  endfor

  if (isfield(mesh.elements, "beam3"))
    beam3degen = struct("nodes", cell(1, 2 * numel(mesh.elements.beam3)));
    ielem = int32(0);
    for i=1:numel(mesh.elements.beam3)
      beam3degen(++ielem).nodes = mesh.elements.beam3(i).nodes([1, 3]);
      beam3degen(++ielem).nodes = mesh.elements.beam3(i).nodes([3, 2]);
    endfor
    mesh.elements.beam2(end + (1:numel(beam3degen))) = beam3degen;
    mesh.elements = rmfield(mesh.elements, "beam3");
  endif

  sol.t = res.t(idx_t);
  sol.def = zeros(rows(mesh.nodes), columns(mesh.nodes), numel(idx_t));

  for i=1:numel(res.log_dat.nodes)
    node_idx = find(res.node_id == res.log_dat.nodes(i).label);
    if (~isempty(node_idx))
      if (numel(node_idx) ~= 1)
        error("invalid node label %d", res.log_dat.nodes(i).label);
      endif

      for j=1:3
        sol.def(node_idx, j, :) = res.trajectory{node_idx}(idx_t, j) - res.log_dat.nodes(i).X0(j);
        sol.def(node_idx, j + 3, :) = res.trajectory{node_idx}(idx_t, j + 3) - res.log_dat.nodes(i).Phi0(j);
      endfor
    endif
  endfor

  iaux_node = zeros(1, numel(elem_types_mb), "int32");

  for k=1:numel(elem_types_mb)
    elem_k = getfield(res.log_dat, elem_types_mb{k});

    for i=1:numel(elem_k)
      node_id = numel(res.log_dat.nodes) + sum(iaux_node(1:k)) + (1:numel(elem_k(i).nodes));

      for j=1:numel(elem_k(i).nodes)
        node_idx = find(res.node_id == elem_k(i).nodes(j).label);

        if (~isempty(node_idx))
          if (numel(node_idx) ~= 1)
            error("invalid node label %d", res.log_dat.nodes(i).label);
          endif

          Xj = res.trajectory{node_idx}(idx_t, 1:3).';
          Phij = res.trajectory{node_idx}(idx_t, 4:6).';
          Rj = R{node_idx};

          for l=1:3
            sol.def(node_id(j), l, :) = Xj(l, :) - mesh.nodes(node_id(j), l);
            sol.def(node_id(j), l + 3, :) = Phij(l, :) - mesh.nodes(node_id(j), l + 3);
            for m=1:3
              sol.def(node_id(j), l, :) += Rj(l, m, :) * elem_k(i).nodes(j).offset(m);
            endfor
          endfor
        endif
      endfor

      iaux_node(k) += numel(elem_k(i).nodes);
    endfor
  endfor
endfunction
