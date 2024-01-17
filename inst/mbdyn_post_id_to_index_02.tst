## mbdyn_post_id_to_index.tst:02
%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_mov_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer node_id_1 = 1001;\n");
%!     fputs(fd, " set: integer body_id_1 = 2001;\n");
%!     fputs(fd, " set: integer force_id_1 = 3001;\n");
%!     fputs(fd, " set: real m1 = 2.5;\n");
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real J11 = 1.;\n");
%!     fputs(fd, " set: real J22 = 1.;\n");
%!     fputs(fd, " set: real J33 = 1.;\n");
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
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: node_id_1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: body_id_1, node_id_1, m1, null, diag, J11, J22, J33;\n");
%!     fputs(fd, "         force: force_id_1, absolute, node_id_1, position, null, 1., 0., 0, F1;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = true;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation, velocity, acceleration, res.node_id, force, res.force_id, force_node_id] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   log_dat.vars = mbdyn_post_id_to_index(res, log_dat.vars);
%!   g = -log_dat.vars.g;
%!   F = force{log_dat.vars.force_idx_1}(:, 1);
%!   m = log_dat.vars.m1;
%!   assert_simple(res.node_id, int32(log_dat.vars.node_id_1));
%!   assert_simple(t, (0:0.1:1.1).', 1e-5);
%!   assert_simple(trajectory{log_dat.vars.node_idx_1}, ...
%!          [0.5 * F / m .* t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g .* t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert_simple(velocity{log_dat.vars.node_idx_1}, ...
%!         [t .* F / m, ...
%!          zeros(numel(t), 1), g .* t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert_simple(acceleration{log_dat.vars.node_idx_1}, ...
%!         [F / m, ...
%!          zeros(numel(t), 1), ...
%!          repmat(g, numel(t), 1), ...
%!          zeros(numel(t), 3)], ...
%!          1e-5);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
