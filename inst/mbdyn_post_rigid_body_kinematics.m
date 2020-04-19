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

%!test
%! r1 = 15e-3;
%! omega0 = 100;
%! Phi0 = 30 * pi / 180;
%! t = linspace(0, abs(2 * pi / omega0), 100);
%! omegaP = 10;
%! Phi1 = Phi0 + omega0 * t + 0.5 * omegaP * t.^2;
%! omega1 = omega0 + omegaP * t;
%! a_rad = -r1 * omega1.^2;
%! a_tan = repmat(r1 * omegaP, 1, length(t));
%! v_tan = r1 * omega1;
%! n1 = [1; 0; 0];
%! o1 = [r1; 0; 0];
%! x1 = [r1 * cos(Phi1);
%!       r1 * sin(Phi1);
%!       zeros(1, length(t))];
%! res.t = t.';
%! res.trajectory{1} = [zeros(5, length(t));
%!                      Phi1].';
%! res.velocity{1} = [zeros(5, length(t));
%!                    omega1].';
%! res.acceleration{1} = [zeros(5, length(t));
%!                        repmat(omegaP, 1, length(t))].';
%! res.node_id = 1234;
%! response(1).node_label = res.node_id;
%! response(1).offset = o1;
%! response(1).direction = n1;
%! response(2).node_label = res.node_id;
%! response(2).offset = o1;
%! response(2).direction = [0; 1; 0];
%! response(3).node_label = res.node_id;
%! response(3).offset = o1;
%! response(3).direction = -n1;
%! response(4).node_label = res.node_id;
%! response(4).offset = o1;
%! response(4).direction = -[0; 1; 0];

%! log_dat.nodes(1).label = res.node_id;
%! log_dat.nodes(1).orientation_description = "euler123";
%! log_dat.nodes(1).X0 = zeros(3, 1);
%! log_dat.nodes(1).R0 = eye(3);
%! [a, v, s] = mbdyn_post_rigid_body_kinematics(res, log_dat, response);
%! tol = eps^0.8;
%! assert(a{1}, a_rad.', tol * norm(a_rad));
%! assert(v{1}, zeros(length(res.t), 1), tol * norm(v_tan));
%! assert(a{2}, a_tan.', tol * norm(a_tan));
%! assert(v{2}, v_tan.', tol * norm(v_tan));
%! assert(a{3}, -a_rad.', tol * norm(a_rad));
%! assert(v{3}, zeros(length(res.t), 1), tol * norm(v_tan));
%! assert(a{4}, -a_tan.', tol * norm(a_tan));
%! assert(v{4}, -v_tan.', tol * norm(v_tan));

%!test
%! t = linspace(0, 10, 100);
%! x0 = 1;
%! v0 = 10;
%! a0 = 20;
%! a = repmat(a0, 1, length(t));
%! v = v0 + a0 * t;
%! x = x0 + v0 * t + 0.5 * a0 * t.^2;

%! n1 = [5; 7; 9];
%! n1 /= norm(n1);
%! o1 = [1; 2; 3];
%! res.t = t.';
%! res.trajectory{1} = [n1 * x;
%!                      zeros(3, length(t))].';
%! res.velocity{1} = [n1 * v;
%!                    zeros(3, length(t))].';
%! res.acceleration{1} = [n1 * a;
%!                        zeros(3, length(t))].';
%! res.node_id = 1234;
%! response(1).node_label = res.node_id;
%! response(1).offset = o1;
%! response(1).direction = n1;

%! log_dat.nodes(1).label = res.node_id;
%! log_dat.nodes(1).orientation_description = "euler123";
%! log_dat.nodes(1).X0 = zeros(3, 1);
%! log_dat.nodes(1).R0 = eye(3);
%! [a1, v1, s1] = mbdyn_post_rigid_body_kinematics(res, log_dat, response);
%! tol = eps;
%! assert(a1{1}, a.', tol * norm(a));
%! assert(v1{1}, v.', tol * norm(v));

%!test
%! r1 = 15e-3;
%! omega0 = 100;
%! Phi0 = 30 * pi / 180;
%! t = linspace(0, abs(2 * pi / omega0), 100);
%! omegaP = 500;
%! Phi1 = Phi0 + omega0 * t + 0.5 * omegaP * t.^2;
%! omega1 = omega0 + omegaP * t;
%! o1 = [r1; 0; 0];
%! x1 = r1 * [cos(Phi1);
%!            sin(Phi1);
%!            zeros(1, length(t))];
%! xP1 = r1 * [-sin(Phi1) .* omega1;
%!              cos(Phi1) .* omega1;
%!              zeros(1, length(t))];
%! xPP1 = r1 * [-cos(Phi1) .* omega1.^2 - sin(Phi1) * omegaP;
%!              -sin(Phi1) .* omega1.^2 + cos(Phi1) * omegaP;
%!              zeros(1, length(t))];
%! res.t = t.';
%! res.trajectory{1} = [zeros(5, length(t));
%!                      Phi1].';
%! res.velocity{1} = [zeros(5, length(t));
%!                    omega1].';
%! res.acceleration{1} = [zeros(5, length(t));
%!                        repmat(omegaP, 1, length(t))].';
%! res.node_id = 1234;
%! response(1).node_label = res.node_id;
%! response(1).offset = o1;

%! log_dat.nodes(1).label = res.node_id;
%! log_dat.nodes(1).orientation_description = "euler123";
%! log_dat.nodes(1).X0 = x1(:, 1)*2;
%! log_dat.nodes(1).R0 = euler123_to_rotation_matrix([0; 0; 2*Phi0]);
%! [a, v, s] = mbdyn_post_rigid_body_kinematics(res, log_dat, response);
%! tol = eps^0.8;
%! assert(s{1}, (x1 - log_dat.nodes(1).X0 - log_dat.nodes(1).R0 * o1).', tol * norm(x1));
%! assert(v{1}, xP1.', tol * norm(xP1));
%! assert(a{1}, xPP1.', tol * norm(xPP1));

%!demo
%! r1 = 15e-3;
%! omega0 = 100;
%! Phi0 = 30 * pi / 180;
%! t = linspace(0, abs(2 * pi / omega0), 100);
%! omegaP = 10;
%! Phi1 = Phi0 + omega0 * t + 0.5 * omegaP * t.^2;
%! omega1 = omega0 + omegaP * t;
%! a_rad = -r1 * omega1.^2;
%! a_tan = repmat(r1 * omegaP, 1, length(t));
%! v_tan = r1 * omega1;
%! n1 = [1; 0; 0];
%! o1 = [r1; 0; 0];
%! x1 = [r1 * cos(Phi1);
%!       r1 * sin(Phi1);
%!       zeros(1, length(t))];
%! res.t = t.';
%! res.trajectory{1} = [zeros(5, length(t));
%!                      Phi1].';
%! res.velocity{1} = [zeros(5, length(t));
%!                    omega1].';
%! res.acceleration{1} = [zeros(5, length(t));
%!                        repmat(omegaP, 1, length(t))].';
%! res.node_id = 1234;
%! response(1).node_label = res.node_id;
%! response(1).offset = o1;
%! response(1).direction = n1;
%! response(2).node_label = res.node_id;
%! response(2).offset = o1;
%! response(2).direction = [0; 1; 0];
%! response(3).node_label = res.node_id;
%! response(3).offset = o1;
%! response(3).direction = -n1;
%! response(4).node_label = res.node_id;
%! response(4).offset = o1;
%! response(4).direction = -[0; 1; 0];

%! log_dat.nodes(1).label = res.node_id;
%! log_dat.nodes(1).orientation_description = "euler123";
%! log_dat.nodes(1).X0 = zeros(3, 1);
%! log_dat.nodes(1).R0 = eye(3);
%! [a, v, s] = mbdyn_post_rigid_body_kinematics(res, log_dat, response);
%! tol = eps^0.8;
%! assert(a{1}, a_rad.', tol * norm(a_rad));
%! assert(v{1}, zeros(length(res.t), 1), tol * norm(v_tan));
%! assert(a{2}, a_tan.', tol * norm(a_tan));
%! assert(v{2}, v_tan.', tol * norm(v_tan));
%! assert(a{3}, -a_rad.', tol * norm(a_rad));
%! assert(v{3}, zeros(length(res.t), 1), tol * norm(v_tan));
%! assert(a{4}, -a_tan.', tol * norm(a_tan));
%! assert(v{4}, -v_tan.', tol * norm(v_tan));
