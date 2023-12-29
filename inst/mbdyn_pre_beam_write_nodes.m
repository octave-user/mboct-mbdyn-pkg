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
## @deftypefn {Function File} mbdyn_pre_beam_write_nodes(@var{beam}, @var{output_file}, @var{options})
##
## Generate an MBDyn input file <@var{output_file}> containing all nodes of beam model <@var{beam}>.
##
## @var{beam} @dots{} Return value from mbdyn_pre_beam_compute.
##
## @var{output_file} @dots{} If <@var{output_file}> is a string, open a new file using <@var{output_file}> as name.
## If <@var{output_file}> it is a file descriptor, write the output to that file descriptor.
##
## @var{open_mode} @dots{} If <@var{output_file}> is a string, pass <@var{open_mode}> to fopen (e.g. "wt", "at").
##
## @var{options}.first_reference_frame_number @dots{} String or integer number of the first reference frame of the beam.
##
## @var{options}.first_node_number @dots{} String or integer number of the first node of the beam.
##
## @var{options}.id_offset @dots{} Integer offset for all body id's.
##
## @var{options}.start_node @dots{} If the beginning of the beam should be connected to an existing node, then <@var{start_node}> is
## the string name of an integer variable holding the existing node number.
##
## @var{options}.end_node @dots{} If the end of the beam should be connected to an existing node, then <@var{end_node}> is the
## string name of an integer variable holding the existing node number.
##
## @var{options}.node_type @dots{} Possible values are "static" and "dynamic".
##
## @end deftypefn

function mbdyn_pre_beam_write_nodes(beam, output_file, options)
  if (nargin < 2)
    print_usage();
  endif

  if (nargin < 3)
    options = struct();
  endif

  if (~isfield(options, "first_reference_frame_number"))
    options.first_reference_frame_number = 1;
  endif

  if (~isfield(options, "first_node_number"))
    options.first_node_number = 1;
  endif

  if (~isfield(options,"node_type"))
    options.node_type = "dynamic";
  endif

  if (~isfield(options, "id_offset"))
    options.id_offset = 0;
  endif

  if (~isfield(options, "open_mode"))
    options.open_mode = "wt";
  endif

  if (~ischar(options.first_node_number))
    options.first_node_number = sprintf("%d", options.first_node_number);
  endif

  if (~ischar(options.first_reference_frame_number))
    options.first_reference_frame_number = sprintf("%d", options.first_reference_frame_number);
  endif

  if (~isfield(options, "struct_node_relative_position"))
    options.struct_node_relative_position = "null";
  endif

  if (~isfield(options, "struct_node_relative_orientation"))
    options.struct_node_relative_orientation = "eye";
  endif

  if (~isfield(options, "struct_node_relative_velocity"))
    options.struct_node_relative_velocity = "null";
  endif

  if (~isfield(options, "struct_node_relative_angular_velocity"))
    options.struct_node_relative_angular_velocity = "null";
  endif

  fout = -1;
  owns_fd = false;

  unwind_protect
    if (ischar(output_file))
      owns_fd = true;

      [fout, msg] = fopen(output_file, options.open_mode);

      if (fout == -1)
        error("could not open file \"%s\": %s", output_file,msg);
      endif
    else
      fout = output_file;
    endif

    start_node = 1;
    end_node = columns(beam.Xn);

    if (isfield(options, "start_node"))
      ++start_node;
    endif

    if (isfield(options, "end_node"))
      --end_node;
    endif

    fprintf(fout, "\n# curved beam node [%d:%d]\n", start_node, end_node);

    for i=start_node:end_node
      fprintf(fout, "\n# curved beam: node #%d\n", i);
      fprintf(fout, "structural: %s + %d + %d - 1, %s,\n", options.first_node_number, i, options.id_offset, options.node_type);
      fprintf(fout, "\tposition, reference, %s + %d + %d - 1, %s,\n", options.first_reference_frame_number, i, options.id_offset, options.struct_node_relative_position);
      fprintf(fout, "\torientation, reference, %s + %d + %d - 1, %s,\n", options.first_reference_frame_number, i, options.id_offset, options.struct_node_relative_orientation);
      fprintf(fout, "\tvelocity, reference, %s + %d + %d - 1, %s,\n", options.first_reference_frame_number, i, options.id_offset, options.struct_node_relative_velocity);
      fprintf(fout, "\tangular velocity, reference, %s + %d + %d - 1, %s", options.first_reference_frame_number, i, options.id_offset, options.struct_node_relative_angular_velocity);

      if (isfield(options,"output_flag"))
        fprintf(fout,",\n");
        fprintf(fout, "\toutput, %s", options.output_flag);
      endif

      fprintf(fout, ";\n");
    endfor

  unwind_protect_cleanup
    if (owns_fd && fout ~= -1)
      fclose(fout);
    endif
  end_unwind_protect
endfunction

%!test
%! f_print_input_file = false;
%! options.verbose = false;
%! N = 3;
%! SI_unit_kilogram = 1e3;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! t1 = 1 / SI_unit_second;
%! R = 1000e-3 / SI_unit_meter;
%! F = 50000 / SI_unit_newton;
%! D = 40e-3 / SI_unit_meter;
%! A = D^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = D^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6 / SI_unit_pascal;
%! G = 81500e6 / SI_unit_pascal;
%! beta = 1e-6 / SI_unit_second^-1;
%! rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%! Phi = linspace(0, pi / 2, N);
%! fmin = 1 / SI_unit_second^-1;
%! fmax = 1000 / SI_unit_second^-1;
%! X = [-R * cos(Phi);
%!      R * sin(Phi);
%!      zeros(1, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! autodiff = {false, true};
%! damping = {false, true};
%! offset = {false, true};
%! orientation = {false, true};
%! balance = {false, true};
%! solvers = {"lapack", "arpack"};
%! cpu.total = zeros(1, numel(autodiff));
%! res = cell(numel(autodiff), numel(damping), numel(offset), numel(orientation), numel(balance), numel(solvers));
%! unwind_protect
%!   for i=1:numel(autodiff)
%!     for j=1:numel(damping)
%!       for k=1:numel(offset)
%!         for l=1:numel(orientation)
%!           for m=1:numel(balance)
%!             for n=1:numel(solvers)
%!           switch (offset{k})
%!             case true
%!               options.struct_node_relative_position = "1., 1., 1.";
%!             case false
%!               options.struct_node_relative_position = "null";
%!           endswitch
%!           switch (orientation{k})
%!             case true
%!               options.struct_node_relative_orientation = "euler123, pi / 2., pi / 4., pi / 8.";
%!             case false
%!               options.struct_node_relative_orientation = "eye";
%!           endswitch
%!           fd = -1;
%!           unwind_protect
%!             [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_write_beams_XXXXXX"));
%!             if (fd == -1)
%!               error("failed to open temporary file");
%!             endif
%!             fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!             fprintf(fd, " set: real t1 = %g;\n", t1);
%!             fprintf(fd, " set: real F = %g;\n", F);
%!             fprintf(fd, " set: real R = %g;\n", R);
%!             fprintf(fd, " set: real A = %g;\n", A);
%!             fprintf(fd, " set: real Ay = %g;\n", Ay);
%!             fprintf(fd, " set: real Az = %g;\n", Az);
%!             fprintf(fd, " set: real Iy = %g;\n", Iy);
%!             fprintf(fd, " set: real Iz = %g;\n", Iz);
%!             fprintf(fd, " set: real It = %g;\n", It);
%!             fprintf(fd, " set: real E = %g;\n", E);
%!             fprintf(fd, " set: real G = %g;\n", G);
%!             fprintf(fd, " set: real rho = %g;\n", rho);
%!             switch (damping{j})
%!               case true
%!                 fprintf(fd, " set: real beta = %g;\n", beta);
%!             endswitch
%!             fprintf(fd, " set: real fmin = %g;\n", fmin);
%!             fprintf(fd, " set: real fmax = %g;\n", fmax);
%!             fputs(fd, " begin: data;\n");
%!             fputs(fd, "         problem: initial value;\n");
%!             fputs(fd, " end: data;\n");
%!             fputs(fd, " begin: initial value;\n");
%!             fputs(fd, "         initial time: 0;\n");
%!             fputs(fd, "         final time: t1;\n");
%!             fputs(fd, "         time step: 0.005 * t1;\n");
%!             switch (autodiff{i})
%!             case true
%!               fputs(fd, "         linear solver: umfpack, grad, pivot factor, 0.1, scale, iterative, always,  tolerance, 1e-12, max iterations, 100;\n");
%!             otherwise
%!               fputs(fd, "         linear solver: umfpack, map, pivot factor, 0.1, scale, iterative, always,  tolerance, 1e-12, max iterations, 100;\n");
%!             endswitch
%!             fputs(fd, "         method: implicit euler;\n");
%!             fputs(fd, "         max iterations: 100;\n");
%!             fputs(fd, "         tolerance: 1e-9, test, sepnorm, 1e-14, test, norm;\n");
%!             if (options.verbose)
%!               fputs(fd, "         output: iterations, solver condition number, stat, yes;\n");
%!             endif
%!             fputs(fd, "         threads: assembly, 1;\n");
%!             fputs(fd, "         derivatives max iterations: 10;\n");
%!             fputs(fd, "         derivatives coefficient: auto;\n");
%!             fputs(fd, "         nonlinear solver: nox, modified, 30,\n");
%!             fputs(fd, "             keep jacobian matrix,\n");
%!             fputs(fd, "             use preconditioner as solver, yes,\n");
%!             fputs(fd, "             inner iterations before assembly, 6,\n");
%!             fputs(fd, "             jacobian operator, newton,\n");
%!             fputs(fd, "             solver, line search based,\n");
%!             fputs(fd, "             forcing term, type 2,\n");
%!             fputs(fd, "             direction, newton,\n");
%!             fputs(fd, "             minimum step, 1e-11,\n");
%!             fputs(fd, "             recovery step, 1e-6,\n");
%!             fputs(fd, "             recovery step type, constant,\n");
%!             fputs(fd, "             linear solver, gmres,\n");
%!             fputs(fd, "             linear solver tolerance, 1e-12,\n");
%!             fputs(fd, "             print convergence info, no,\n");
%!             fputs(fd, "             linear solver max iterations, 12,\n");
%!             fputs(fd, "             krylov subspace size, 12;\n");
%!             fputs(fd, "         eigenanalysis: t1,\n");
%!             fputs(fd, "          parameter, 1. / (2. * pi * fmin),\n");
%!             fputs(fd, "          suffix format, \"%02d\",\n");
%!             fputs(fd, "          output eigenvectors,\n");
%!             fputs(fd, "          output geometry,\n");
%!             fputs(fd, "          lower frequency limit, fmin,\n");
%!             fputs(fd, "          upper frequency limit, fmax,\n");
%!             fputs(fd, "          results output precision, 16,\n");
%!             switch (solvers{n})
%!             case "lapack"
%!               fputs(fd, "          use lapack");
%!               if (balance{m})
%!                 fputs(fd, ", balance, permute");
%!               endif
%!               fputs(fd, ";\n");
%!             case "arpack"
%!               fputs(fd, "mode, largest imaginary part,\n\tuse arpack, 20, 90, 0;\n");
%!             endswitch
%!               fputs(fd, " end: initial value;\n");
%!               fputs(fd, " begin: control data;\n");
%!               switch (autodiff{i})
%!               case true
%!                 fputs(fd, "     use automatic differentiation;\n");
%!               endswitch
%!             fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn) + 1);
%!             fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!             fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!             fputs(fd, "       joints: 1;\n");
%!             fputs(fd, "       forces: 1;\n");
%!             fputs(fd, "       output precision: 16;\n");
%!             fputs(fd, "       default output: none;\n");
%!             fputs(fd, " end: control data;\n");
%!             switch (damping{j})
%!               case false
%!                 fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!               case true
%!                 fputs(fd, " constitutive law: 1, 6, linear viscoelastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz, proportional, beta;\n");
%!             endswitch
%!             mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!             fputs(fd, " begin: nodes;\n");
%!             mbdyn_pre_beam_write_nodes(beam, fd, options);
%!             fprintf(fd, "structural: %d, dummy, %d, offset, reference, %d, null, reference, %d, eye, output, yes;\n", columns(beam.Xn) + 1, columns(beam.Xn), columns(beam.Xn), columns(beam.Xn));
%!             fputs(fd, " end: nodes;\n");
%!             fputs(fd, " begin: elements;\n");
%!             mbdyn_pre_beam_write_bodies(beam, fd, options);
%!             mbdyn_pre_beam_write_beams(beam, fd, options);
%!             fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!             fprintf(fd, "force: 1, absolute, %d, position, reference, %d, null, component, null, string, \"Time / t1 * F\", null;\n", columns(beam.Xn), columns(beam.Xn));
%!             fputs(fd, " end: elements;\n");
%!           unwind_protect_cleanup
%!             if (fd ~= -1)
%!               fclose(fd);
%!             endif
%!           end_unwind_protect
%!           options.output_file = fname;
%!           if (f_print_input_file)
%!             spawn_wait(spawn("nl", {fname}));
%!           endif
%!           if (~options.verbose)
%!             options.logfile = [fname, ".stdout"];
%!           endif
%!           start = tic();
%!           mbdyn_solver_run(fname, options);
%!           cpu.total(i) += toc(start);
%!           res{i, j, k, l, m, n}.log_dat = mbdyn_post_load_log(fname);
%!           res{i, j, k, l, m, n}.modal = mbdyn_post_load_output_eig(fname);
%!           [res{i, j, k, l, m, n}.t, res{i, j, k, l, m, n}.trajectory, res{i, j, k, l, m, n}.deformation, res{i, j, k, l, m, n}.velocity, res{i, j, k, l, m, n}.acceleration, res{i, j, k, l, m, n}.node_id, res{i, j, k, l, m, n}.force,  res{i, j, k, l, m, n}.force_id, res{i, j, k, l, m, n}.force_node_id, res{i, j, k, l, m, n}.orientation_description] = mbdyn_post_load_output_struct(fname);
%!           res{i, j, k, l, m, n}.R = mbdyn_post_angles_to_rotation_mat(columns(beam.Xn) + 1, res{i, j, k, l, m, n}, res{i, j, k, l, m, n}.log_dat);
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%!   endfor
%!   endfor
%!   tolR = 2e-6;
%!   tolU = 2e-6;
%!   tolf = 2e-5;
%!   tolD = 3e-2;
%!   for j=1:numel(damping)
%!     for i=1:numel(autodiff)
%!       for k=1:numel(offset)
%!         for l=1:numel(orientation)
%!           for m=1:numel(balance)
%!             for n=1:numel(solvers)
%!               U1 = res{i, j, k, l, m, n}.trajectory{end};
%!               R1 = res{i, j, k, l, m, n}.R{end};
%!               f1 = res{i, j, k, l, m, n}.modal.f;
%!               D1 = res{i, j, k, l, m, n}.modal.D;
%!               U2 = res{1, j, 1, 1, 1, 1}.trajectory{end};
%!               R2 = res{1, j, 1, 1, 1, 1}.R{end};
%!               f2 = res{1, j, 1, 1, 1, 1}.modal.f;
%!               D2 = res{1, j, 1, 1, 1, 1}.modal.D;
%!               M = min([10, numel(f1), numel(f2)]);
%!               df = f1(1:M) ./ f2(1:M) - 1;
%!               dD = (D1(1:M) - D2(1:M)) * damping{j};
%!               D1ref = 2 * pi * f1 * beta / 2;
%!               if (options.verbose)
%!                 fprintf(stderr, "[%d, %d, %d, %d, %d, %d]:\n", i, j, k, l, m, n);
%!                 fprintf(stderr, "%5.2f %3.2f%% %5.2f %3.2f%% %4.2f%% %4.2f%%\n", [f1(1:M).'; 100 * D1(1:M).'; f2(1:M).'; 100 * D2(1:M).'; 100 * df.'; 100 * dD.' / max(abs(D1ref(1:M))).']);
%!               endif
%!               try
%!                 assert_simple(max(abs(df)) < tolf);
%!                 assert_simple(max(abs(dD)) < tolD * max(abs(D1ref(1:M))));
%!                 assert_simple(max(max(abs(U1(:,1:3) - U2(:, 1:3)))) < tolU * max(max(abs(U2))));
%!                 for o=1:size(R1, 3)
%!                   assert_simple(max(max(max(abs(R1(:, :, o).' * R2(:, :, o) - eye(3))))) < tolR);
%!                   assert_simple(max(max(max(abs(R2(:, :, o).' * R1(:, :, o) - eye(3))))) < tolR);
%!                 endfor
%!               catch
%!                 fprintf(stderr, "mbdyn_pre_beam_write_nodes:test1[%d, %d, %d, %d, %d, %d]:\n", i, j, k, l, m, n);
%!                 fprintf(stderr, "%5.2f %3.2f%% %5.2f %3.2f%% %4.2f%% %4.2f%%\n", [f1(1:M).'; 100 * D1(1:M).'; f2(1:M).'; 100 * D2(1:M).'; 100 * df.'; 100 * dD.' / max(abs(D1ref(1:M))).']);
%!                 rethrow(lasterror());
%!               end_try_catch
%!              endfor
%!           endfor
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
