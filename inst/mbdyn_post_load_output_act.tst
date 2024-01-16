%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_act_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real L = 0.5;\n");
%!     fputs(fd, " set: real F1 = 10;\n");
%!     fputs(fd, " set: real D = 50e-3;\n");
%!     fputs(fd, " set: real A = D^2 * pi / 4.;\n");
%!     fputs(fd, " set: real Ay = 9. / 10. * A;\n");
%!     fputs(fd, " set: real Az = Ay;\n");
%!     fputs(fd, " set: real Iy = D^4 * pi / 64.;\n");
%!     fputs(fd, " set: real Iz = Iy;\n");
%!     fputs(fd, " set: real It = Iy + Iz;\n");
%!     fputs(fd, " set: real E = 210000e6;\n");
%!     fputs(fd, " set: real G = 81500e6;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-1;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-5;\n");
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
%!     fputs(fd, "     structural nodes: 3;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     beams: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 2, static,\n");
%!     fputs(fd, "                 0.5 * L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 3, static,\n");
%!     fputs(fd, "                 L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "         beam3: 1, 1, null, 2, null, 3, null, eye,\n");
%!     fputs(fd, "                linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz,\n");
%!     fputs(fd, "                same, same;\n");
%!     fputs(fd, "         force: 1, absolute, 3, position, null, 0., 0., -1, mult, time, F1;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t] = mbdyn_post_load_output_out(fname, 1024, false);
%!   log_dat = mbdyn_post_load_log(fname);
%!   [beam_id, F_I, M_I, F_II, M_II] = mbdyn_post_load_output_act(fname, 1, numel(t));
%!   tol = eps^0.4;
%!   assert_simple(F_I{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert_simple(F_II{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert_simple(M_I{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 + 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
%!   assert_simple(M_II{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 - 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
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
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_act_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real L = 0.5;\n");
%!     fputs(fd, " set: real F1 = 10;\n");
%!     fputs(fd, " set: real D = 50e-3;\n");
%!     fputs(fd, " set: real A = D^2 * pi / 4.;\n");
%!     fputs(fd, " set: real Ay = 9. / 10. * A;\n");
%!     fputs(fd, " set: real Az = Ay;\n");
%!     fputs(fd, " set: real Iy = D^4 * pi / 64.;\n");
%!     fputs(fd, " set: real Iz = Iy;\n");
%!     fputs(fd, " set: real It = Iy + Iz;\n");
%!     fputs(fd, " set: real E = 210000e6;\n");
%!     fputs(fd, " set: real G = 81500e6;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-1;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-5;\n");
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
%!     fputs(fd, "     structural nodes: 3;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     beams: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 2, static,\n");
%!     fputs(fd, "                 0.5 * L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 3, static,\n");
%!     fputs(fd, "                 L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "         beam3: 1, 1, null, 2, null, 3, null, eye,\n");
%!     fputs(fd, "                linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz,\n");
%!     fputs(fd, "                same, same;\n");
%!     fputs(fd, "         force: 1, absolute, 3, position, null, 0., 0., -1, mult, time, F1;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t] = mbdyn_post_load_output_out(fname, 1024, false);
%!   log_dat = mbdyn_post_load_log(fname);
%!   [beam_id, F_I, M_I, F_II, M_II] = mbdyn_post_load_output_act(fname, 1, numel(t));
%!   tol = eps^0.4;
%!   assert_simple(F_I{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert_simple(F_II{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert_simple(M_I{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 + 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
%!   assert_simple(M_II{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 - 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
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
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_act_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real L = 0.5;\n");
%!     fputs(fd, " set: real F1 = 10;\n");
%!     fputs(fd, " set: real D = 50e-3;\n");
%!     fputs(fd, " set: real A = D^2 * pi / 4.;\n");
%!     fputs(fd, " set: real Ay = 9. / 10. * A;\n");
%!     fputs(fd, " set: real Az = Ay;\n");
%!     fputs(fd, " set: real Iy = D^4 * pi / 64.;\n");
%!     fputs(fd, " set: real Iz = Iy;\n");
%!     fputs(fd, " set: real It = Iy + Iz;\n");
%!     fputs(fd, " set: real E = 210000e6;\n");
%!     fputs(fd, " set: real G = 81500e6;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-1;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-5;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 3;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     beams: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 2, static,\n");
%!     fputs(fd, "                 0.5 * L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 3, static,\n");
%!     fputs(fd, "                 L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "         beam3: 1, 1, null, 2, null, 3, null, eye,\n");
%!     fputs(fd, "                linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz,\n");
%!     fputs(fd, "                same, same;\n");
%!     fputs(fd, "         force: 1, absolute, 3, position, null, 0., 0., -1, mult, time, F1;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t] = mbdyn_post_load_output_out(fname, 1024, false);
%!   log_dat = mbdyn_post_load_log(fname);
%!   [beam_id, F_I, M_I, F_II, M_II] = mbdyn_post_load_output_act(fname, 1, numel(t));
%!   tol = eps^0.4;
%!   assert_simple(F_I{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert_simple(F_II{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert_simple(M_I{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 + 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
%!   assert_simple(M_II{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 - 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%!demo
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_act_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real L = 0.5;\n");
%!     fputs(fd, " set: real F1 = 10;\n");
%!     fputs(fd, " set: real D = 50e-3;\n");
%!     fputs(fd, " set: real A = D^2 * pi / 4.;\n");
%!     fputs(fd, " set: real Ay = 9. / 10. * A;\n");
%!     fputs(fd, " set: real Az = Ay;\n");
%!     fputs(fd, " set: real Iy = D^4 * pi / 64.;\n");
%!     fputs(fd, " set: real Iz = Iy;\n");
%!     fputs(fd, " set: real It = Iy + Iz;\n");
%!     fputs(fd, " set: real E = 210000e6;\n");
%!     fputs(fd, " set: real G = 81500e6;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-1;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-5;\n");
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
%!     fputs(fd, "     structural nodes: 3;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     beams: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 2, static,\n");
%!     fputs(fd, "                 0.5 * L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 3, static,\n");
%!     fputs(fd, "                 L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "         beam3: 1, 1, null, 2, null, 3, null, eye,\n");
%!     fputs(fd, "                linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz,\n");
%!     fputs(fd, "                same, same;\n");
%!     fputs(fd, "         force: 1, absolute, 3, position, null, 0., 0., -1, mult, time, F1;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t] = mbdyn_post_load_output_out(fname, 1024, false);
%!   log_dat = mbdyn_post_load_log(fname);
%!   [beam_id, F_I, M_I, F_II, M_II] = mbdyn_post_load_output_act(fname, 1, numel(t));
%!   tol = eps^0.4;
%!   assert_simple(F_I{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert_simple(F_II{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert_simple(M_I{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 + 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
%!   assert_simple(M_II{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 - 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
