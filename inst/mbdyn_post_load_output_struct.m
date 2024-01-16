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
## @deftypefn {Function File} [@var{t}, @var{trajectory}, @var{deformation}, @var{velocity}, @var{acceleration}, @var{node_id}, @var{force}, @var{force_id}, @var{force_node_id}, @var{orientation_description}] = mbdyn_post_load_output_struct(@var{mbdyn_filename}, @var{filter_node_id}, @var{filter_force_id})
##
## Loads data from MBDyn .mov, .out and .frc files ("<@var{mbdyn_filename}.mov>", "<@var{mbdyn_filename}.out>" and "<@var{mbdyn_filename}.frc>").
##
## @var{filter_node_id} @dots{} Load only the data of the corresponding nodes if present.
## If empty the data of all nodes is loaded.
##
## @var{t} @dots{} Simulation time.
##
## @var{trajectory} @dots{} Matrix of position and orientation of structural nodes.
##
## @var{deformation} @dots{} The difference between actual position/orientation and initial position/orientation.
##
## @var{velocity} @dots{} Matrix of velocity and angular velocity of structural nodes.
##
## @var{acceleration} @dots{} Matrix of acceleration and angular acceleration.
##
## @var{node_id} @dots{} Structural node numbers.
##
## @var{force} @dots{} Structural force values.
##
## @var{force_id} @dots{} Structural force element number.
##
## @var{force_node_id} @dots{} Structural node number.
##
## @var{orientation_description} @dots{} Cell array of strings. One of ("euler123", "euler321", "euler313", "phi")
##
## @end deftypefn

function [t, trajectory, deformation, velocity, acceleration, node_id, force,  force_id, force_node_id, orientation_description, t_trc, OutputFlag] = mbdyn_post_load_output_struct(mbdyn_filename, filter_node_id, filter_force_id, auto_resize_rows)

  if (nargin < 1 || nargout > 12)
    print_usage();
  endif

  if (nargin < 2)
    filter_node_id = [];
  endif

  if (nargin < 3)
    filter_force_id = [];
  endif

  if (nargin < 4)
    auto_resize_rows = false;
  endif

  [t_trc, TStep, NIter, ResErr, SolErr, SolConv, OutputFlag] = mbdyn_post_load_output_out(mbdyn_post_output_filename(mbdyn_filename), 1024, auto_resize_rows);

  t = t_trc(find(OutputFlag));

  if (nargout >= 2)
    [node_id, trajectory, velocity, acceleration, orientation_description, deformation] = mbdyn_post_load_output_mov(mbdyn_filename, filter_node_id, numel(t), 6, auto_resize_rows);
  endif

  if (nargout >= 7)
    if (2 == exist(mbdyn_post_output_filename(mbdyn_filename, ".frc"), "file"))
      [force_id, force_node_id, force] = mbdyn_post_load_output_frc(mbdyn_filename, filter_force_id, numel(t), auto_resize_rows);
    else
      force_id = [];
      force_node_id = [];
      force = {};
    endif
  endif
endfunction
