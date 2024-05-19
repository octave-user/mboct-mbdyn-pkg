## mbdyn_post_abs_to_rel.tst:03
%!test
%! try
%! ## TEST3
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 200;\n");
%!     fputs(fd, " set: real F3 = 300;\n");
%!     fputs(fd, " set: real m1 = 1;\n");
%!     fputs(fd, " set: real m2 = 2;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 0.1;\n");
%!     fputs(fd, "         time step: 1e-5;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, "         output: iterations;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     output meter: closest next, 0, forever, const, 1e-2;\n");
%!     fputs(fd, "     output precision: 16;\n");
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " reference: 1,\n");
%!     fputs(fd, "     position, reference, global, 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "     orientation, reference, global, euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "     velocity, reference, global, 1., 2., 3.,\n");
%!     fputs(fd, "     angular velocity, reference, global, 0.1, 0.2, 0.3;\n");
%!     fputs(fd, " reference: 2,\n");
%!     fputs(fd, "     position, reference, 1, 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "     orientation, reference, 1, euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "     velocity, reference, 1, null,\n");
%!     fputs(fd, "     angular velocity, reference, 1, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, eye,");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 reference, 2, eye,");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., mult, time, const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., mult, time, const, F2;\n");
%!     fputs(fd, "         joint: 1, total joint,\n");
%!     fputs(fd, "                   1,\n");
%!     fputs(fd, "                    position, reference, 1, null,\n");
%!     fputs(fd, "                    position orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    rotation orientation, reference, 1, eye,\n");
%!     fputs(fd, "                   2,\n");
%!     fputs(fd, "                    position, reference, 1, null,\n");
%!     fputs(fd, "                    position orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    rotation orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    position constraint, active, active, active, null,\n");
%!     fputs(fd, "                    orientation constraint, active, active, active, null;\n");
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
%!   mbdyn_post_abs_to_rel(2, options.output_file, [options.output_file, "_rel"], false);
%!   copyfile([options.output_file, ".log"], [options.output_file, "_rel.log"]);
%!   copyfile([options.output_file, ".out"], [options.output_file, "_rel.out"]);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory_abs, deformation_abs, velocity_abs, acceleration_abs, node_id] = mbdyn_post_load_output_struct(options.output_file, 1:2);
%!   [t, trajectory_rel, deformation_rel, velocity_rel, acceleration_rel, node_id] = mbdyn_post_load_output_struct([options.output_file, "_rel"], 1:2);
%!   tol = eps^0.4;
%!   for i=1:numel(node_id)
%!     assert_simple(trajectory_rel{i}, repmat(trajectory_rel{i}(1,:), numel(t), 1), tol * max(max(abs(trajectory_abs{i}))));
%!     assert_simple(velocity_rel{i}, zeros(numel(t), 6), tol * max(max(abs(velocity_abs{i}))));
%!     ## FIXME: bug in abs2rel.awk?
%!     ## assert_simple(acceleration_rel{i}, zeros(numel(t), 6), tol * max(max(abs(acceleration_abs{i}))));
%!   endfor
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
