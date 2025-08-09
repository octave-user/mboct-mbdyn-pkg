## mbdyn_pre_beam_write_beams.tst:02
%!test
%! try
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 30;
%! R = 1000e-3;
%! F = 150;
%! D = 50e-3;
%! A = D^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = D^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! Phi = linspace(0, pi / 2, N);
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
%! fd = -1;
%! %unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_beam_write_beams_XXXXXX"));
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
%!     fputs(fd, "         linear solver: naive, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-4;\n");
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
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "       forces: 1;\n");
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
%!     fprintf(fd, " force: 1, absolute, %d, position, reference, node, null, 0., -1., 0., mult, time, F;\n", columns(beam.Xn));
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
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
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
%!   Xn = [log_dat.nodes.X0];
%!   Phin = atan2(Xn(2, :), -Xn(1, :));
%!   wref = F * R^3 / (2 * E * Iz) * Phin .* sin(Phin);
%!   uref = F * R^3 / (2 * E * Iz) * (sin(Phin) - Phin .* cos(Phin));
%!   Uref = [wref .* cos(Phin) + uref .* sin(Phin);
%!           -wref .* sin(Phin) + uref .* cos(Phin);
%!           zeros(1, numel(Phin))];
%!   U = zeros(3, numel(Phin));
%!   for i=1:columns(U)
%!     U(:, i) = deformation{i}(end, 1:3).';
%!   endfor
%!   tol = 1e-3;
%!   assert_simple(U, Uref, tol * norm(Uref));
%!   wref_1 = F * R^3 / (2 * E * Iz);
%!   w_1 = deformation{end}(end, 1);
%!   tol = 1e-3;
%!   assert_simple(w_1, wref_1, tol * max(abs(wref_1)));
%! %unwind_protect_cleanup
%!   if (fd ~= -1)
%!     [~] = unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       [~] = unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! %end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
