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
## @deftypefn {Function File} @var{eccentric_acceleration} = mbdyn_post_rigid_body_kinematics(@var{res}, @var{log_dat}, @var{response})
##
## Compute position, velocity and acceleration of a point rigidly attached to a structural node.
##
## @var{response}(@var{i}).node_label @dots{} Label of the node where the response <@var{i}> is measured.
##
## @var{response}(@var{i}).offset @dots{} Offset between the point, where the response <@var{i}> is measured, and the node
## with respect to the reference frame of the node.
##
## @var{response}(@var{i}).direction @dots{} Measurement direction of the response <@var{i}>
## with respect to the reference frame of the node.
##
## @end deftypefn

function [eccentric_acceleration, eccentric_velocity, eccentric_position] = mbdyn_post_rigid_body_kinematics(res, log_dat, response)
  if (nargin ~= 3 || nargout > 3)
    print_usage();
  endif

  for i=1:length(response)
    node_idx_i = find(response(i).node_label == res.node_id);

    if (length(node_idx_i) == 0)
      error("node_label = %d not found in node_id!", response(i).node_label);
    endif

    if (isfield(response(i), "offset") && length(response(i).offset) ~= 0)
      oi = response(i).offset;
    else
      error("offset required for response(%d)", i);
    endif

    if (isfield(response(i), "direction") && length(response(i).direction) ~= 0)
      ni = response(i).direction;
      ni /= norm(ni);
    endif

    ## nonlinear case
    R1 = mbdyn_post_angles_to_rotation_mat(response(i).node_label, res, log_dat){1};

    if (isfield(response(i), "direction") && length(response(i).direction) > 0)
      eccentric_acceleration{i} = zeros(length(res.t), 1);
    else
      eccentric_acceleration{i} = zeros(length(res.t), 3);
    endif

    if (nargout >= 2)
      eccentric_velocity{i} = eccentric_acceleration{i};
    endif

    if (nargout >= 3)
      eccentric_position{i} = eccentric_acceleration{i};
    endif

    node_idx_i_log = find([log_dat.nodes.label] == response(i).node_label);

    if (numel(node_idx_i_log) ~= 1)
      error("node %d not found in log_dat", response(i).node_id);
    endif

    X1_0 = log_dat.nodes(node_idx_i_log).X0;
    R1_0 = log_dat.nodes(node_idx_i_log).R0;
    Xi_0 = X1_0 + R1_0 * oi;

    if (isfield(response(i), "direction") && numel(response(i).direction) > 0)
      R1_0_ni = R1_0 * ni;
    endif

    for j=1:numel(res.t)
      X1 = res.trajectory{node_idx_i}(j, 1:3).';
      XP1 = res.velocity{node_idx_i}(j, 1:3).';
      XPP1 = res.acceleration{node_idx_i}(j, 1:3).';
      omega1 = res.velocity{node_idx_i}(j, 4:6).';
      omegaP1 = res.acceleration{node_idx_i}(j, 4:6).';
      R1_oi = R1(:, :, j) * oi;
      XPPi = XPP1 - cross(R1_oi, omegaP1) - cross(omega1, cross(R1_oi, omega1));

      if (nargout >= 2)
        XPi = XP1 + cross(omega1, R1_oi);
      endif

      if (nargout >= 3)
        Xi = X1 + R1_oi;
      endif

      if (isfield(response(i), "direction") && numel(response(i).direction) > 0)
        R1_ni = R1(:, :, j) * ni;

        eccentric_acceleration{i}(j) = R1_ni.' * XPPi;

        if (nargout >= 2)
          eccentric_velocity{i}(j) = R1_ni.' * XPi;
        endif

        if (nargout >= 3)
          eccentric_position{i}(j) = R1_ni.' * Xi - R1_0_ni.' * Xi_0;
        endif
      else
        eccentric_acceleration{i}(j, :) = XPPi.';

        if (nargout >= 2)
          eccentric_velocity{i}(j, :) = XPi.';
        endif

        if (nargout >= 3)
          eccentric_position{i}(j, :) = (Xi - Xi_0).';
        endif
      endif
    endfor
  endfor
endfunction
