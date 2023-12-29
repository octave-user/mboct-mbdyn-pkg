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
## @deftypefn {Function File} mbdyn_post_plot_model(@var{output_file}, @var{res}, @var{idx_t})
##
## Plot the deformed shape of an MBDyn model (currently only beam elements are plotted).
##
## @end deftypefn

function mbdyn_post_plot_model(output_file, res, idx_t)
  if (nargin < 2 || nargin > 3 || nargout > 0)
    print_usage();
  endif

  if (nargin < 3)
    idx_t = 1:numel(res.t);
  endif

  R = mbdyn_post_angles_to_rotation_mat([res.log_dat.nodes.label], res, res.log_dat);

  X = [res.log_dat.nodes.X0];

  lim = zeros(3, 2);

  for i=1:3
    lim(i, 1) = min(X(i, :));
    lim(i, 2) = max(X(i, :));
  endfor

  node_ids = [res.log_dat.nodes.label];

  for k=1:numel(idx_t)
    for i=1:numel(res.log_dat.beams3)
      Xi = zeros(3, numel(res.log_dat.beams3(i).nodes));

      for j=1:numel(res.log_dat.beams3(i).nodes)
        node_idx = find(node_ids == res.log_dat.beams3(i).nodes(j).label);
        Xi(:, j) = res.trajectory{node_idx}(idx_t(k), 1:3).' + R{node_idx}(:, :, idx_t(k)) * res.log_dat.beams3(i).nodes(j).offset(:);
      endfor

      for j=1:3
        lim(j, 1) = min(lim(j, 1), min(Xi(j, :)));
        lim(j, 2) = max(lim(j, 2), max(Xi(j, :)));
      endfor
    endfor
  endfor

  hfig = 0;
  
  unwind_protect
    hfig = figure("visible", "off");

    for k=1:numel(idx_t)
      clf;
      hold on;

      for i=1:numel(res.log_dat.beams3)
        Xi = zeros(3, numel(res.log_dat.beams3(i).nodes));

        for j=1:numel(res.log_dat.beams3(i).nodes)
          node_idx = find(node_ids == res.log_dat.beams3(i).nodes(j).label);
          Xi(:, j) = res.trajectory{node_idx}(idx_t(k), 1:3).' + R{node_idx}(:, :, idx_t(k)) * res.log_dat.beams3(i).nodes(j).offset(:);
        endfor

        line("xdata", Xi(1, :), "ydata", Xi(2, :), "zdata", Xi(3, :));
      endfor

      daspect(ones(1,3));
      view(30, 30);
      
      if (lim(1, 2) > lim(1, 1))
        xlim(lim(1, :));
      endif

      if (lim(2, 2) > lim(2, 1))
        ylim(lim(2, :));
      endif

      if (lim(3, 2) > lim(3, 1))
        zlim(lim(3, :));
      endif
      
      xlabel("x");
      ylabel("y");
      zlabel("z");
      print(hfig, "-dpng", sprintf("%s_%03d.png", output_file, k));
    endfor

    pid = spawn("ffmpeg", {"-i", [output_file, "_%03d.png"], "-y", [output_file, ".mpg"]});
    
    status = spawn_wait(pid);

    if (status ~= 0)
      warning("ffmpeg failed with status %d", status);
    endif
  unwind_protect_cleanup
    if (hfig ~= 0)
      close(hfig);
    endif
  end_unwind_protect
endfunction

%!test
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 10;
%! l = 1000e-3;
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
%! X = [linspace(0, l, N);
%!      zeros(2, N)];
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_plot_model_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real g = %g;\n", g);
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
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     gravity;\n");
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
%!     fputs(fd, " gravity: uniform, 0., 0., -1., mult, time, g;\n");
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
%!   options.verbose = true;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [res.t, res.trajectory, res.deformation, res.velocity, res.acceleration, res.node_id] = mbdyn_post_load_output_struct(options.output_file);
%!   res.log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
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
%!   for i=1:numel(bodies)
%!     assert_simple(bodies(i).node, int32(i));
%!   endfor
%!   w = zeros(1, numel(res.deformation));
%!   for i=1:numel(w)
%!     w(i) = res.deformation{i}(end, 3);
%!   endfor
%!   z = zeros(1, numel(res.log_dat.nodes));
%!   for i=1:numel(res.log_dat.nodes)
%!     z(i) = l - res.log_dat.nodes(i).X0(1);
%!   endfor
%!   wref = -rho * A * g * l^4 / (24 * E * Iy) * (3 - 4 * z / l + (z / l).^4);
%!   if (f_plot_deformation)
%!     figure("visible", "off");
%!     hold on;
%!     plot(z, wref, "-;wref;0");
%!     plot(z, w, "-;w;1");
%!     grid on;
%!     grid minor on;
%!     title("deformation of a cantilever beam");
%!     xlabel("z [m]");
%!     ylabel("w [m]");
%!   endif
%!   gtk = graphics_toolkit();
%!   unwind_protect
%!     graphics_toolkit("gnuplot");
%!     mbdyn_post_plot_model([fname, "_video"], res);
%!   unwind_protect_cleanup
%!     graphics_toolkit(gtk);
%!   end_unwind_protect
%!   tol = 1e-2;
%!   assert_simple(w, wref, tol * max(abs(wref)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect


%!test
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 10;
%! l = 1000e-3;
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
%! X = [linspace(0, l, N);
%!      zeros(2, N)];
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_plot_model_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real g = %g;\n", g);
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
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     gravity;\n");
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
%!     fputs(fd, " gravity: uniform, 0., 0., -1., mult, time, g;\n");
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
%!   options.verbose = true;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [res.t, res.trajectory, res.deformation, res.velocity, res.acceleration, res.node_id] = mbdyn_post_load_output_struct(options.output_file);
%!   res.log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
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
%!   for i=1:numel(bodies)
%!     assert_simple(bodies(i).node, int32(i));
%!   endfor
%!   w = zeros(1, numel(res.deformation));
%!   for i=1:numel(w)
%!     w(i) = res.deformation{i}(end, 3);
%!   endfor
%!   z = zeros(1, numel(res.log_dat.nodes));
%!   for i=1:numel(res.log_dat.nodes)
%!     z(i) = l - res.log_dat.nodes(i).X0(1);
%!   endfor
%!   wref = -rho * A * g * l^4 / (24 * E * Iy) * (3 - 4 * z / l + (z / l).^4);
%!   if (f_plot_deformation)
%!     figure("visible", "off");
%!     hold on;
%!     plot(z, wref, "-;wref;0");
%!     plot(z, w, "-;w;1");
%!     grid on;
%!     grid minor on;
%!     title("deformation of a cantilever beam");
%!     xlabel("z [m]");
%!     ylabel("w [m]");
%!   endif
%!   gtk = graphics_toolkit();
%!   unwind_protect
%!     graphics_toolkit("gnuplot");
%!     mbdyn_post_plot_model([fname, "_video"], res);
%!   unwind_protect_cleanup
%!     graphics_toolkit(gtk);
%!   end_unwind_protect
%!   tol = 1e-2;
%!   assert_simple(w, wref, tol * max(abs(wref)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 10;
%! l = 1000e-3;
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
%! X = [linspace(0, l, N);
%!      zeros(2, N)];
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_plot_model_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real g = %g;\n", g);
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
%!     fputs(fd, "         max iterations: 20;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, "         output: iterations;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     gravity;\n");
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
%!     fputs(fd, " gravity: uniform, 0., 0., -1., mult, time, g;\n");
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
%!   options.verbose = true;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [res.t, res.trajectory, res.deformation, res.velocity, res.acceleration, res.node_id] = mbdyn_post_load_output_struct(options.output_file);
%!   res.log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
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
%!   for i=1:numel(bodies)
%!     assert_simple(bodies(i).node, int32(i));
%!   endfor
%!   w = zeros(1, numel(res.deformation));
%!   for i=1:numel(w)
%!     w(i) = res.deformation{i}(end, 3);
%!   endfor
%!   z = zeros(1, numel(res.log_dat.nodes));
%!   for i=1:numel(res.log_dat.nodes)
%!     z(i) = l - res.log_dat.nodes(i).X0(1);
%!   endfor
%!   wref = -rho * A * g * l^4 / (24 * E * Iy) * (3 - 4 * z / l + (z / l).^4);
%!   if (f_plot_deformation)
%!     figure("visible", "off");
%!     hold on;
%!     plot(z, wref, "-;wref;0");
%!     plot(z, w, "-;w;1");
%!     grid on;
%!     grid minor on;
%!     title("deformation of a cantilever beam");
%!     xlabel("z [m]");
%!     ylabel("w [m]");
%!   endif
%!   gtk = graphics_toolkit();
%!   unwind_protect
%!     graphics_toolkit("gnuplot");
%!     mbdyn_post_plot_model([fname, "_video"], res);
%!   unwind_protect_cleanup
%!     graphics_toolkit(gtk);
%!   end_unwind_protect
%!   tol = 1e-2;
%!   assert_simple(w, wref, tol * max(abs(wref)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
