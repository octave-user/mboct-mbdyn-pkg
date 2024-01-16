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
## @deftypefn {Function File} [@var{nodes}, @var{dof_info}] = mbdyn_post_load_log_node (@var{mbdyn_filename})
##
## Parses the MBDyn log-file "<@var{mbdyn_filename}>.log" and returns information about the structural node degrees of freedom.
##
## @var{mbdyn_filename} @dots{} MBDyn output file without extension
##
## @var{nodes} @dots{} Struct array with information about all nodes including dummy nodes in ascending order
##
## @var{nodes}(@var{i}).X0 @dots{} Initial position of the structural node @var{i}
##
## @var{nodes}(@var{i}).Phi0 @dots{} Angles of rotation of the structural node @var{i}
##
## @var{nodes}(@var{i}).R0 @dots{} Initial rotation matrix of the structural node @var{i}
##
## @var{nodes}(@var{i}).orientation_description @dots{} Orientation_description of the node: one of  ("euler123", "euler321", "euler313","phi","mat")
##
## @var{dof_info}.struct_node_dofs @dots{} Index of the first degree of freedom minus one of the node in the global vector of degrees of freedom (e.g X_i = X(struct_node_dofs(i) + 1:3); g_i = X(struct_node_dofs(i) + 4:6))
##
## @var{dof_info}.struct_node_force_index @dots{} Index of the first equation of the force and moment equilibrium in the global residual vector (e.g F_i = f(struct_node_force_index(i) + 1:3); M_i = f(struct_node_force_index(i) + 4:6))
##
## @var{dof_info}.struct_node_labels @dots{} Labels of the corresponding structural nodes
##
## @var{X}(struct_node_dofs(i) + 1:6) @dots{} Is the vector of degrees of freedom of a static or dynamic structural node with label @var{struct_node_labels}(i)
##
##
## @end deftypefn

function [nodes, dof_info] = mbdyn_post_load_log_node(mbdyn_filename)
  if (nargin < 1 || nargout > 2)
    print_usage();
  endif

  options.nodes = true;
  options.dof_info = (nargout >= 2);
  options.vars = false;
  options.beams2 = false;
  options.beams3 = false;

  data = mbdyn_post_load_log(mbdyn_filename, options);

  nodes = data.nodes;

  if (nargout >= 2)
    dof_info = data.dof_info;
  endif
endfunction
