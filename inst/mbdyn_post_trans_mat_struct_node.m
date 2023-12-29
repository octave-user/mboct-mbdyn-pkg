## Copyright (C) 2014(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{T}] = mbdyn_post_trans_mat_struct_node(@var{dof_info}, @var{number_dofs}, @var{node_label}, @var{type}, @var{offset}, @var{direction})
##
## Assemble a transformation matrix for a structural node needed in mbdyn_post_frequency_response.
##
## @var{dof_info} @dots{} Structural node degrees of freedom returned as the second argument from mbdyn_post_load_log_node.
##
## @var{number_dofs} @dots{} Number of degrees of freedom (size of the jacobian matrix).
##
## @var{node_label} @dots{} Label of the node where the force and/or torque is applied.
##
## @var{type} @dots{} Must be either "force" for the force transformation matrix or "displacement" for the displacement transformation matrix.
##
## @var{offset} @dots{} Optional the offset between force and node. In this case only a eccentric force and no torque are the input quantity.
##
## @var{direction} @dots{} Optional the direction of the force. In this case only the magnitude of the force is the input quantity.
##
## @var{T} @dots{} Transformation matrix for a force and/or torque: @var{p} = @var{T}.' * @var{p}_@var{i}.
##
## @seealso{mbdyn_post_frequency_response}
## @end deftypefn

function T = mbdyn_post_trans_mat_struct_node(dof_info, number_dofs, node_label, type, offset, direction, component)
  if (nargin < 4 || nargin > 7 || nargout > 1)
    print_usage();
  endif

  if (nargin < 5)
    offset = [];
  endif

  if (nargin < 6)
    direction = [];
  endif

  if (nargin < 7)
    component = [];
  endif

  node_idx = find(dof_info.struct_node_labels == node_label);

  if (length(node_idx) == 0)
    error("node_id %d not found!", node_label);
  endif

  if (numel(component))
    if (numel(component) > 3)
      error("numel(component) must be maximum three");
    endif

    dof = int32(component);

    if (any(dof < 1 | dof > 6))
      error("dof must be in range [1:6]");
    endif
  else
    dof = int32(1:6);
  endif

  switch (type)
    case "displacement"
      dof += dof_info.struct_node_dofs( node_idx );
    case "force"
      dof += dof_info.struct_node_force_index( node_idx );
    otherwise
      error("type=\"%s\" is invalid!",type);
  endswitch

  T = sparse(1:numel(dof), dof, ones(1, numel(dof)), numel(dof), number_dofs);

  if (numel(offset))
    if (rows(T) ~= 6)
      error("component must not be used in combination with offset");
    endif

    T = sparse([eye(3), -skew(offset)]) * T;
  endif

  if (numel(direction))
    if (rows(T) ~= 3)
      error("component must not be used in combination with direction");
    endif

    direction /= norm(direction);
    T = sparse(direction.') * T;
  endif
endfunction
