## mbdyn_post_node_id_to_node_index.tst:03
%!test
%! try
%! fd = -1;
%! %unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_node_id_to_node_index_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
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
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         set: integer node_id_1 = 1000;\n");
%!     fputs(fd, "         structural: node_id_1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, node_id_1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         force: 1, absolute, node_id_1, position, null, 1., 0., 0, 100.;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., 9.81;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, dt, niter, reserr, solerr, solconv, output_flag] = mbdyn_post_load_output_out(fname, 1024, false);
%!   [node_id, trajectory, velocity, acceleration, orientation_description, deformation] = mbdyn_post_load_output_mov(options.output_file, [log_dat.vars.node_id_1], numel(t), 6, false);
%!   g = -9.81;
%!   F = 100;
%!   m = 1;
%!   assert_simple(all(output_flag));
%!   assert_simple(node_id, int32(log_dat.vars.node_id_1));
%!   assert_simple(t, (0:0.1:1.1).', 1e-5);
%!   log_dat.vars = mbdyn_post_node_id_to_node_index(log_dat.vars.node_id_1, log_dat.vars);
%!   assert_simple(trajectory{log_dat.vars.node_idx_1}, ...
%!          [0.5 * F / m * t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g * t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert_simple(velocity{log_dat.vars.node_idx_1}, ...
%!         [t * F / m, ...
%!          zeros(numel(t), 1), g * t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert_simple(acceleration{log_dat.vars.node_idx_1}, ...
%!         [F / m * ones(numel(t), 1), ...
%!          zeros(numel(t), 1), ...
%!          repmat(g, numel(t), 1), ...
%!          zeros(numel(t), 3)], ...
%!          1e-5);
%! %unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! %end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
