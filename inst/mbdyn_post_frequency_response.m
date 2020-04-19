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
## @deftypefn {Function File} [@var{F}] = mbdyn_post_frequency_response(@var{modal}, @var{dof_info}, @var{excitation}, @var{response}, @var{omega}, @var{p}, @var{options})
## Compute the frequency response function of the linearized equations of motion of a MBDyn model.
##
## @var{F} @dots{} Complex frequency response function
##
## @var{modal} @dots{} Data structure returned from mbdyn_post_load_output_eig()
##
## @var{dof_info} @dots{} Data structure returned form mbdyn_post_load_log_node()
##
## @var{excitation}.node_label @dots{} Label of the node that receives the force @var{p}
##
## @var{excitation}.offset @dots{} Offset between force and node
##
## @var{excitation}.direction @dots{} Direction of the force @var{p}
##
## @var{response}(@var{i}).node_label @dots{} Label of the node where the response @var{i} is measured
##
## @var{response}(@var{i}).offset @dots{} Offset between the point where the response @var{i} is measured and the node
##
## @var{response}(@var{i}).direction @dots{} Measurement direction of the response @var{i}
##
## @var{omega} @dots{} Angular velocity of the excitation force @var{p}
##
## @var{p} @dots{} Complex excitation force
##
## @var{options}.singular @dots{} If true, it is assumed that df/dy + j * omega * df/dy_dot is singular
##
## @var{options}.matrix_type @dots{} "wre" means F(omega_idx, response_idx, excitation_idx)
##
##                                   "rew" means F(response_idx, excitation_idx, omega_idx)
##
## @end deftypefn

function F = mbdyn_post_frequency_response(modal, dof_info, excitation, response, omega, p, options)
  if (nargin < 6 || nargin > 7 || nargout > 1)
    print_usage();
  endif

  if (nargin < 7)
    options = struct("singular", false, "matrix_format", "wre");
  endif

  if (~isfield(options,"singular"))
    options.singular = false;
  endif

  if (~isfield(options, "number_of_processors"))
    options.number_of_processors = 1;
  endif

  if (~isfield(options, "matrix_format"))
    options.matrix_format = "wre";
  endif

  for i=1:length(excitation)
    if (~isfield(excitation, "offset"))
      excitation(i).offset = [];
    endif

    if (~isfield(excitation, "direction"))
      excitation(i).direction = [];
    endif

    if (~isfield(excitation, "component"))
      excitation(i).component = [];
    endif
  endfor

  for i=1:length(response)
    if (~isfield(response, "offset"))
      response(i).offset = [];
    endif

    if (~isfield(response, "direction"))
      response(i).direction = [];
    endif

    if (~isfield(response, "component"))
      response(i).component = [];
    endif
  endfor

  Jac1 = modal.Aplus;
  Jac2 = modal.Aminus;
  dCoef1 = modal.dCoef;
  dCoef2 = -modal.dCoef;

  df_dy = (Jac2 - Jac1) / (dCoef1 - dCoef2);
  df_dy_dot = -dCoef1 * df_dy - Jac1;

  number_dofs = length(df_dy);

  Tp = sparse([], [], [], 0, number_dofs);

  for i=1:length(excitation)
    Tp_i = mbdyn_post_trans_mat_struct_node(dof_info, number_dofs, excitation(i).node_label, "force", excitation(i).offset, excitation(i).direction, excitation(i).component);
    Tp((end + 1):(end + rows(Tp_i)), :) = Tp_i;
  endfor

  Tu = sparse([], [], [], 0, number_dofs);

  for i=1:length(response)
    Tu_i = mbdyn_post_trans_mat_struct_node(dof_info, number_dofs, response(i).node_label, "displacement", response(i).offset, response(i).direction, response(i).component);
    Tu((end + 1):(end + rows(Tu_i)), :) = Tu_i;
  endfor

  N_response = rows(Tu);

  if (columns(p) == 1)
    N_excitation = rows(Tp);
  else
    N_excitation = 1;
  endif

  switch(options.matrix_format)
    case "wre"
      F = zeros(length(omega), N_response, N_excitation);
    case "rew"
      F = zeros(N_response, N_excitation, length(omega));
  endswitch

  if (options.number_of_processors == 1)
    for i=1:length(omega)
      y = mbdyn_post_frequency_response_helper(i, df_dy, df_dy_dot, Tp, p, Tu, omega, options);

      switch(options.matrix_format)
        case "wre"
          F(i, :, :) = y;
        case "rew"
          F(:, :, i) = y;
      endswitch
    endfor
  else
    options_par.number_of_processors = options.number_of_processors;
    options_par.number_of_parameters = length(omega);
    y = run_parallel(options_par, @mbdyn_post_frequency_response_helper, df_dy, df_dy_dot, Tp, p, Tu, omega, options);

    for i=1:length(omega)
      switch(options.matrix_format)
        case "wre"
          F(i, :, :) = y{i};
        case "rew"
          F(:, :, i) = y{i};
      endswitch
    endfor
  endif
endfunction

%!test
%! ##  Frequency response of a cantilever beam
%! ##  Robert Gash, Klaus Knothe
%! ##  Strukturdynamik Band 2
%! ##  Page 25 figure 2.10, equation 2.38b
%!
%! f_print_input_file = false;
%! f_plot_response = false;
%! if (f_plot_response)
%!  close("all");
%! endif
%! N = 40;
%! M = 200;
%! l = 4;
%! X = [linspace(0, l, N);
%!      zeros(2, N)];
%! g = 9.81;
%! D = 10e-3;
%! A = D^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = D^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! B = E * Iy;
%! mu = rho * A;
%! P0 = (0.3 * l) / l^3 * 3 * B;
%! omega1 = sqrt(B / (mu * l^4));
%! wstat = P0 * l^3 / (3 * B);
%! omega_krit = [3.516, 22.035, 61.697] * omega1;
%! omega = linspace(1e-5, 1, M) * 0.5 * omega_krit(end);
%! DELTA = sqrt(omega ./ omega1) / l;
%! V = 3 ./ (DELTA * l).^3 .* (sin(DELTA * l) .* cosh(DELTA * l) - cos(DELTA * l) .* sinh(DELTA * l)) ./ (1 .+ cos(DELTA * l) .* cosh(DELTA * l));
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_frequency_response_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
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
%!     fputs(fd, "         final time: 1e-2;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         eigenanalysis: 0,\n");
%!     fputs(fd, "          suffix format, \"%02d\",\n");
%!     fputs(fd, "          output sparse matrices,\n");
%!     fputs(fd, "          output eigenvectors,\n");
%!     fputs(fd, "          results output precision, 16,\n");
%!     fputs(fd, "          parameter, 0.01, use arpack, 100, 200, 0;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 1;\n");
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
%!   options.logfile = [fname, ".stdout"];
%!   options.verbose = false;
%!   info = mbdyn_solver_run(fname, options);
%!   log_dat = mbdyn_post_load_log(fname);
%!   modal = mbdyn_post_load_output_eig(options.output_file);
%!   excitation.node_label = columns(beam.Xn);
%!   excitation.offset = [0; 0; 0];
%!   excitation.direction = [0; 0; 1];
%!   response.node_label = columns(beam.Xn);
%!   response.offset = [0; 0; 0];
%!   response.direction = [0; 0; 1];
%!   options.matrix_type = "wre";
%!   F = mbdyn_post_frequency_response(modal, log_dat.dof_info, excitation, response, omega, P0, options);
%!   if (f_plot_response)
%!     figure("visible", "off");
%!     subplot(2, 1, 1);
%!     hold on;
%!     semilogy(omega / (2 * pi), abs(F), "-;abs(F);1");
%!     semilogy(omega / (2 * pi), abs(V * wstat), "-;abs(V);0");
%!     grid minor on;
%!     xlabel("f [Hz]");
%!     ylabel("abs(F) [m]");
%!     title("frequency response magnitude");
%!     subplot(2, 1, 2);
%!     hold on;
%!     plot(omega / (2 * pi), 180 / pi * arg(F), "-;arg(F);1");
%!     plot(omega / (2 * pi), 180 / pi * arg(V * wstat), "-;arg(V);0");
%!     grid minor on;
%!     xlabel("f [Hz]");
%!     ylabel("arg(F) [deg]");
%!     title("frequency response phase");
%!   endif
%!   tol = 1e-2;
%!   assert(sum(abs(F - V.' * wstat).^2) / sum(abs(V * wstat).^2) < tol);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! ##  Frequency response of a cantilever beam
%! ##  Robert Gash, Klaus Knothe
%! ##  Strukturdynamik Band 2
%! ##  Page 25 figure 2.10, equation 2.38b
%!
%! f_print_input_file = false;
%! f_plot_response = true;
%! if (f_plot_response)
%!  close("all");
%! endif
%! N = 40;
%! M = 200;
%! l = 4;
%! X = [linspace(0, l, N);
%!      zeros(2, N)];
%! g = 9.81;
%! D = 10e-3;
%! A = D^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = D^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! B = E * Iy;
%! mu = rho * A;
%! P0 = (0.3 * l) / l^3 * 3 * B;
%! omega1 = sqrt(B / (mu * l^4));
%! wstat = P0 * l^3 / (3 * B);
%! omega_krit = [3.516, 22.035, 61.697] * omega1;
%! omega = linspace(1e-5, 1, M) * 0.5 * omega_krit(end);
%! DELTA = sqrt(omega ./ omega1) / l;
%! V = 3 ./ (DELTA * l).^3 .* (sin(DELTA * l) .* cosh(DELTA * l) - cos(DELTA * l) .* sinh(DELTA * l)) ./ (1 .+ cos(DELTA * l) .* cosh(DELTA * l));
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_frequency_response_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
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
%!     fputs(fd, "         final time: 1e-2;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         eigenanalysis: 0,\n");
%!     fputs(fd, "          suffix format, \"%02d\",\n");
%!     fputs(fd, "          output sparse matrices,\n");
%!     fputs(fd, "          output eigenvectors,\n");
%!     fputs(fd, "          results output precision, 16,\n");
%!     fputs(fd, "          parameter, 0.01, use arpack, 100, 200, 0;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 1;\n");
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
%!   options.logfile = [fname, ".stdout"];
%!   options.verbose = false;
%!   info = mbdyn_solver_run(fname, options);
%!   log_dat = mbdyn_post_load_log(fname);
%!   modal = mbdyn_post_load_output_eig(options.output_file);
%!   excitation.node_label = columns(beam.Xn);
%!   excitation.offset = [0; 0; 0];
%!   excitation.direction = [0; 0; 1];
%!   response.node_label = columns(beam.Xn);
%!   response.offset = [0; 0; 0];
%!   response.direction = [0; 0; 1];
%!   options.matrix_type = "wre";
%!   F = mbdyn_post_frequency_response(modal, log_dat.dof_info, excitation, response, omega, P0, options);
%!   if (f_plot_response)
%!     figure("visible", "off");
%!     subplot(2, 1, 1);
%!     hold on;
%!     semilogy(omega / (2 * pi), abs(F), "-;abs(F);1");
%!     semilogy(omega / (2 * pi), abs(V * wstat), "-;abs(V);0");
%!     grid minor on;
%!     xlabel("f [Hz]");
%!     ylabel("abs(F) [m]");
%!     title("frequency response magnitude");
%!     subplot(2, 1, 2);
%!     hold on;
%!     plot(omega / (2 * pi), 180 / pi * arg(F), "-;arg(F);1");
%!     plot(omega / (2 * pi), 180 / pi * arg(V * wstat), "-;arg(V);0");
%!     grid minor on;
%!     xlabel("f [Hz]");
%!     ylabel("arg(F) [deg]");
%!     title("frequency response phase");
%!   endif
%!   tol = 1e-2;
%!   assert(sum(abs(F - V.' * wstat).^2) / sum(abs(V * wstat).^2) < tol);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
