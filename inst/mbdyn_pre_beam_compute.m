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
## @deftypefn {Function File} @var{beam} = mbdyn_pre_beam_compute(@var{X}, @var{N}, @var{interpolation_points})
##
## Compute nodal and gauss point positions and orientations for a curved beam model defined by grid points @var{X}.
##
## @var{X} @dots{} Vertices of the curve to be interpolated.
##
## @var{N} @dots{} The number of beam elements to be generated.
##
## @var{interpolation_points} @dots{} The number of points used to interpolate the length of the beam.
##
## @end deftypefn

function beam = mbdyn_pre_beam_compute(X, N, interpolation_points, options)
  if (nargin < 2 || nargin > 4 || nargout > 1)
    print_usage();
  endif

  if (nargin < 3)
    interpolation_points = 1;
  endif

  if (rows(X) ~= 3)
    error("X must be a 3xN vector");
  endif

  if (columns(X) < 2)
    error("the number of columns of X is not sufficient");
  endif

  if (~isscalar(N) || N < 1 || floor(N) ~= N)
    error("N must be an integer >= 1");
  endif

  N = double(N);

  if (~isscalar(interpolation_points)
      || interpolation_points < 1
      || floor(interpolation_points) ~= interpolation_points)
    error("interpolation_points must be an integer >= 1");
  endif

  if (nargin < 4)
    options = struct();
  endif

  if (~isfield(options, "smooth_curvature"))
    options.smooth_curvature = false;
  endif

  interpolation_points = double(interpolation_points);

  beam.X = X;

  pkg load nurbs;

  knots = [0,0,linspace(0,1,columns(X)-1),1,1];

  beam.crv = nrbmak(X,knots);

  [beam.dcrv, beam.dcrv2] = nrbderiv(beam.crv);

  beam.ti = linspace(0, 1, columns(X) * interpolation_points);

  beam.si = mbdyn_curved_beam_length_vector(beam, beam.ti);

  beam.sn = linspace(beam.si(1), beam.si(end), 2 * N + 1);

  beam.sg = zeros(1, 2 * N);

  beam.sg(1:2:end) = beam.si(end) / (2 * N) * (2 * (1:N) - 1 - 1 / sqrt(3));
  beam.sg(2:2:end) = beam.si(end) / (2 * N) * (2 * (1:N) - 1 + 1 / sqrt(3));

  [beam.Xn, beam.Rn] = mbdyn_curved_beam_interpolation(beam, beam.sn);

  [beam.Xg, beam.Rg] = mbdyn_curved_beam_interpolation(beam, beam.sg);

  if (options.smooth_curvature)
    beam.cosPhi = ones(1, size(beam.Rn, 3));

    for i=1:size(beam.Rn, 3) - 1
      R1 = beam.Rn(:, :, i);
      R2 = beam.Rg(:, :, i);
      R3 = beam.Rn(:, :, i + 1);

      [R2, cosPhi2] = mbdyn_pre_beam_orient_smooth(R1, R2);
      [R3, cosPhi3] = mbdyn_pre_beam_orient_smooth(R2, R3);

      beam.Rg(:, :, i) = R2;
      beam.Rn(:, :, i + 1) = R3;
      beam.cosPhi(i + 1) = cosPhi3;
    endfor
  endif

  for i=1:N
    beam.beams(i).nidx = [2 * i - 1, 2 * i, 2 * i + 1];
    beam.beams(i).gidx = [2 * i - 1, 2 * i];
  endfor

  for i=1:columns(beam.sn)
    ds = 0;

    if (i > 1)
      ds += (beam.sn(i) - beam.sn(i - 1)) / 2;
    endif

    if (i < columns(beam.sn))
      ds += (beam.sn(i + 1) - beam.sn(i)) / 2;
    endif

    beam.bodies(i).ds = ds;
  endfor
endfunction

function [X, jac, hess] = mbdyn_curved_beam_nurbs_interpolation(beam, t)
  [X, jac, hess] = nrbdeval(beam.crv, beam.dcrv, beam.dcrv2, t);
endfunction

function [X, R] = mbdyn_curved_beam_position_orientation(beam, t)
  [X, jac, hess] = mbdyn_curved_beam_nurbs_interpolation(beam, t);

  if (iscell(jac))
    jac = jac{:};
  endif

  if (iscell(hess))
    hess = hess{:};
  endif

  R = zeros(3, 3, columns(X));

  for i=1:columns(X)
    e1 = jac(:,i);

    e1 /= norm(e1);

    vi = hess(:,i);

    e3 = cross(e1, vi);

    if (norm(e3) < sqrt(eps))
      vi = mbdyn_curved_beam_compute_e2(e1);
      e3 = cross(e1, vi);
    endif

    e2 = cross(e3, e1);
    e2 /= norm(e2);

    e3 /= norm(e3);

    R(:, :, i) = [e1, e2, e3];

    mbdyn_pre_beam_check_rotation_matrix(R(:, :, i));
  endfor
endfunction

function [e2, theta] = mbdyn_curved_beam_compute_e2(e1)
  e1 /= norm(e1);

  theta = [0;
           asin(-e1(3));
           atan2(e1(2), e1(1)) ];

  e2 = [-sin(theta(3));
        cos(theta(3));
        0];
endfunction

function [X, R] = mbdyn_curved_beam_interpolation(beam, s)
  t = mbdyn_curved_beam_parameters(beam, s);
  [X, R] = mbdyn_curved_beam_position_orientation(beam, t);
endfunction

function t = mbdyn_curved_beam_parameters(beam, s)
  t = interp1(beam.si, beam.ti, s, 'spline');
endfunction

function s = mbdyn_curved_beam_length_vector(beam, t)
  s = cumtrapz(t, mbdyn_curved_beam_segment_length_ds(beam, t));
endfunction

function s = mbdyn_curved_beam_segment_length(beam, t0=0, t1=1)
                                % best performance with quadv
                                % quadl ... 9.89776s
                                % quadv ... 1.75s
                                % quadcc ... 2.11s
  s = quadv(@(t) mbdyn_curved_beam_segment_length_ds(beam, t), t0, t1);
endfunction

function ds = mbdyn_curved_beam_segment_length_ds(beam, t)
  [X, jac] = nrbdeval(beam.crv, beam.dcrv, beam.dcrv2, t);

  if (iscell(jac))
    jac = jac{:};
  endif

  ds = sqrt(sum(jac.^2,1));
endfunction

function [R2s, cosPhi] = mbdyn_pre_beam_orient_smooth(R1, R2)
  e1 = R2(:, 1);
  e3 = cross(R2(:, 1), R1(:, 2));
  e2 = cross(R1(:, 3), R2(:, 1));

  e3a = cross(e1, e2);
  e2a = cross(e3a, e1);

  e2b = cross(e3, e1);
  e3b = cross(e1, e2b);

  cosPhi2a = e2a.' * e2 / (norm(e2a) * norm(e2));
  cosPhi3b = e3b.' * e3 / (norm(e3b) * norm(e3));

  if (cosPhi2a > cosPhi3b)
    cosPhi = cosPhi2a;
    R2s = [e1, e2a, e3a];
  else
    cosPhi = cosPhi3b;
    R2s = [e1, e2b, e3b];
  endif

  R2s *= diag(1 ./ norm(R2s, "cols"));

  mbdyn_pre_beam_check_rotation_matrix(R2s);
endfunction
