## Copyright (C) 2025(-2025) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}, @var{sol}] = mbdyn_post_load_output_eig_sol(@var{output_file})
## @deftypefnx {} [@dots{}] = mbdyn_post_load_output_eig_sol(@var{output_file}, @var{options})
## @deftypefnx {} [@dots{}] = mbdyn_post_load_output_eig_sol(@var{output_file}, @var{options}, @var{index})
##
## Loads modal data for solid and beam elements from an MBDyn output file "<@var{output_file}>".
##
## @var{output_file} @dots{} Name of the MBDyn output file without extension.
##
## @var{options} @dots{} Options passed to mbdyn_post_load_output_eig
##
## @var{index} @dots{} In case of multiple analyses, load the results for the corresponding @var{index}.
##
## @var{mesh} @dots{} Finite element mesh data structure
##
## @var{sol} @dots{} Finite element solution data structure
##
## @seealso{mbdyn_post_load_output_eig}
## @end deftypefn

function [mesh, sol] = mbdyn_post_load_output_eig_sol(output_file, options, index)
  if (nargin ~= 1 || nargout ~= 2 || ~ischar(output_file))
    print_usage();
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (nargin < 3)
    index = 0;
  endif

  [mesh] = mbdyn_post_load_mesh_sol(output_file);

  modal = mbdyn_post_load_output_eig(output_file, options, index);

  sol.f = modal.f;

  sol.def = zeros(rows(mesh.nodes), 6, numel(sol.f));

  node_idx_trans = zeros(numel(modal.labels), 1, "int32");

  for i=1:numel(modal.labels)
    node_idx_trans(i) = find(modal.labels(i) == mesh.node_id);
  endfor

  for i=1:numel(sol.f)
    for j=1:3
      sol.def(node_idx_trans, j, i) = modal.VR(modal.idx + j, i);
    endfor
  endfor
endfunction
