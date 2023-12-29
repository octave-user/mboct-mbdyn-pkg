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

function beam = mbdyn_pre_beam_compute(X, N, interpolation_points)
  if (nargin < 2 || nargin > 3 || nargout > 1)
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

  ds = sqrt(sum(jac.^2,1));
endfunction

%!test
%! N = 4;

%! Theta2 = 88*pi/180;
%! Theta3 = -10*pi/180;

%! X = [ linspace(1,7,N);
%!       zeros(1,N);
%!       zeros(1,N) ];

%! R = [cos(Theta2)*cos(Theta3),-sin(Theta3),sin(Theta2)*cos(Theta3);
%!      cos(Theta2)*sin(Theta3),cos(Theta3),sin(Theta2)*sin(Theta3);
%!      -sin(Theta2),0,cos(Theta2)];

%! n = [ 1, 1;
%!       0, 0;
%!       0, 0 ];

%! beam = mbdyn_pre_beam_compute(R * X, N);
%! assert_simple(R, beam.Rn(:,:,1), sqrt(eps));

%!test
%! N = 100;
%! X = [ 1,2,3,4;
%!       0,3,0,0;
%!       0,4,0,0 ];
%! beam = mbdyn_pre_beam_compute(X, N, 40);
%! norm_dXn = norm(beam.Xn(:,2:end) - beam.Xn(:,1:end-1), 2, 'cols');
%! f = max(abs(1 - norm_dXn/mean(norm_dXn)));
%! assert_simple(f < 1e-2);

%!test
%! close all;
%! N = 4;

%! Theta2 = 88*pi/180;
%! Theta3 = -10*pi/180;

%! X = [ linspace(1,7,N);
%!       zeros(1,N);
%!       zeros(1,N) ];

%! R = [cos(Theta2)*cos(Theta3),-sin(Theta3),sin(Theta2)*cos(Theta3);
%!      cos(Theta2)*sin(Theta3),cos(Theta3),sin(Theta2)*sin(Theta3);
%!      -sin(Theta2),0,cos(Theta2)];

%! n = [ 1, 1;
%!       0, 0;
%!       0, 0 ];
%! beam = mbdyn_pre_beam_compute(R*X,N);
%! figure("visible", "off");
%! mbdyn_pre_beam_plot(beam,struct("s",0.001,"Rn",false,"Rg",false));
%! title("straight beam");
%! assert_simple(R,beam.Rn(:,:,1),sqrt(eps));

%!test
%! f_print_input_file = false;
%! F = 10;
%! d = 1e-3;
%! A = d^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = d^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! D = 10e-3;
%! L = 1e-3;
%! n = 1;
%! N = n * 7;
%! Phi = linspace(0, 2 * pi * n, n * 36);
%! X = [0.5 * D * cos(Phi);
%!      0.5 * D * sin(Phi);
%!      linspace(0, L, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 10);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         nonlinear solver: nox, modified, 10,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 12,\n");
%!     fputs(fd, "             krylov subspace size, 12;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "joint: 2, total pin joint,\n");
%!     fprintf(fd, " %d,\n", columns(beam.Xn));
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fputs(fd, "     position constraint, active, active, inactive, null,\n");
%!     fputs(fd, "     orientation constraint, active, active, active, null;\n");
%!     fprintf(fd, "force: 1, absolute, %d, position, reference, %d, null, 0., 0., -1., mult, time, F;\n", columns(beam.Xn), columns(beam.Xn));
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   if (f_print_input_file)
%!     spawn_wait(spawn("cat", {fname}));
%!   endif
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert_simple(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert_simple(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert_simple(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert_simple(bodies(i).node, int32(i));
%!   endfor
%!   R = G * d^4 / (8 * n * D^3);
%!   wref = -F / R;
%!   tol = 1e-2;
%!   assert_simple(deformation{end}(end, 3), wref, tol * abs(wref));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! f_print_input_file = false;
%! F = 10;
%! d = 1e-3;
%! A = d^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = d^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! D = 10e-3;
%! L = 1e-3;
%! n = 1;
%! N = n * 7;
%! Phi = linspace(0, 2 * pi * n, n * 36);
%! X = [0.5 * D * cos(Phi);
%!      0.5 * D * sin(Phi);
%!      linspace(0, L, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 10);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         nonlinear solver: nox, modified, 10,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 12,\n");
%!     fputs(fd, "             krylov subspace size, 12;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "joint: 2, total pin joint,\n");
%!     fprintf(fd, " %d,\n", columns(beam.Xn));
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fputs(fd, "     position constraint, active, active, inactive, null,\n");
%!     fputs(fd, "     orientation constraint, active, active, active, null;\n");
%!     fprintf(fd, "force: 1, absolute, %d, position, reference, %d, null, 0., 0., -1., mult, time, F;\n", columns(beam.Xn), columns(beam.Xn));
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   if (f_print_input_file)
%!     spawn_wait(spawn("cat", {fname}));
%!   endif
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert_simple(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert_simple(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert_simple(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert_simple(bodies(i).node, int32(i));
%!   endfor
%!   R = G * d^4 / (8 * n * D^3);
%!   wref = -F / R;
%!   tol = 1e-2;
%!   assert_simple(deformation{end}(end, 3), wref, tol * abs(wref));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! f_print_input_file = false;
%! F = 10;
%! d = 1e-3;
%! A = d^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = d^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! D = 10e-3;
%! L = 1e-3;
%! n = 1;
%! N = n * 7;
%! Phi = linspace(0, 2 * pi * n, n * 36);
%! X = [0.5 * D * cos(Phi);
%!      0.5 * D * sin(Phi);
%!      linspace(0, L, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 10);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "joint: 2, total pin joint,\n");
%!     fprintf(fd, " %d,\n", columns(beam.Xn));
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fputs(fd, "     position constraint, active, active, inactive, null,\n");
%!     fputs(fd, "     orientation constraint, active, active, active, null;\n");
%!     fprintf(fd, "force: 1, absolute, %d, position, reference, %d, null, 0., 0., -1., mult, time, F;\n", columns(beam.Xn), columns(beam.Xn));
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   if (f_print_input_file)
%!     spawn_wait(spawn("cat", {fname}));
%!   endif
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert_simple(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert_simple(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert_simple(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert_simple(bodies(i).node, int32(i));
%!   endfor
%!   R = G * d^4 / (8 * n * D^3);
%!   wref = -F / R;
%!   tol = 1e-2;
%!   assert_simple(deformation{end}(end, 3), wref, tol * abs(wref));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! close all;
%! N = 10;
%! X = [ 1,2,3,4;
%!       0,3,0,0;
%!       0,4,0,0 ];
%! beam = mbdyn_pre_beam_compute(X,N,50);

%! figure("visible", "off");
%! plot(beam.ti,beam.si,'-;s(t);1');
%! xlabel('t []');
%! ylabel('s [m]');
%! title('curve length versus parameter');
%! grid on;
%! grid minor on;

%! norm_dXn = norm(beam.Xn(:,2:end) - beam.Xn(:,1:end-1),2,'cols');
%! f = max(abs(1-norm_dXn/mean(norm_dXn)));
%! figure("visible","off");
%! stem(norm_dXn,'o-;norm(dXn);1');
%! xlabel('node #');
%! ylabel('norm(dXn) [m]');
%! assert_simple(f < 0.18);
%! mbdyn_pre_beam_plot(beam,struct("s",0.1,"X",true,"Rn",true,"Rg",true));

%!test
%! close all;
%! r = 0.5;
%! v = 0.5 * r;
%! h = 2;
%! n = 1;
%! N = 4;
%! k = 200;
%! Phi = linspace(0,2*pi*n,k*N+1);
%! X = [ r*cos(Phi);
%!       r*sin(Phi);
%!       h * Phi./(2*pi) ];
%!
%! Xt = [ -r*sin(Phi);
%!         r*cos(Phi);
%!         repmat(h / (2 * pi),1,length(Phi)) ];
%!
%! t = 2*pi*n/N*Xt(:,[1,end]);
%! beam = mbdyn_pre_beam_compute(X,N);
%! figure("visible", "off");
%! hold on;
%! plot3(X(1,1:k:k*N+1),X(2,1:k:k*N+1),X(3,1:k:k*N+1),'x;curve;3');
%! set(plot3(X(1,:),X(2,:),X(3,:),'-;curve;3'),'linewidth',1);
%! set(gca(),'dataaspectratio',[1,1,1]);
%! xlabel('x');
%! ylabel('y');
%! zlabel('z');
%! grid on;
%! grid minor on;
%! title('helix');
%! mbdyn_pre_beam_plot(beam,struct("s",0.08,"Rn",true,"Rg",true));
%! hold on;
%! plot3(X(1,:),X(2,:),X(3,:),'--;curve;0');

%!demo
%! f_print_input_file = false;
%! F = 8;
%! d = 1.3e-3;
%! A = d^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = d^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! D = 10e-3;
%! L = 5e-3;
%! n = 5;
%! N = n * 21;
%! Phi = linspace(0, 2 * pi * n, n * 3600);
%! X = [0.5 * D * cos(Phi);
%!      0.5 * D * sin(Phi);
%!      linspace(0, L, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 10);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");

%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "joint: 2, total pin joint,\n");
%!     fprintf(fd, " %d,\n", columns(beam.Xn));
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fputs(fd, "     position constraint, active, active, inactive, null,\n");
%!     fputs(fd, "     orientation constraint, active, active, active, null;\n");
%!     fprintf(fd, "force: 1, absolute, %d, position, reference, %d, null, 0., 0., -1., mult, time, F;\n", columns(beam.Xn), columns(beam.Xn));
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   if (f_print_input_file)
%!     spawn_wait(spawn("cat", {fname}));
%!   endif
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options.output_file);
%!   res.log_dat = mbdyn_post_load_log(fname);
%!   res.bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(res.log_dat.nodes)
%!     assert_simple(res.log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert_simple(res.log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(res.log_dat.beams3)
%!     for j=1:3
%!       assert_simple(res.log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(res.bodies)
%!     assert_simple(res.bodies(i).node, int32(i));
%!   endfor
%!   R = G * d^4 / (8 * n * D^3);
%!   wref = -F / R;
%!   tol = 1e-2;
%!   assert_simple(res.deformation{end}(end, 3), wref, tol * abs(wref));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
