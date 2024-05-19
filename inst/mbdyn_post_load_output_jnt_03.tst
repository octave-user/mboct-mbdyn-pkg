## mbdyn_post_load_output_jnt.tst:03
%!test
%! try
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_jnt_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 200;\n");
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
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static,\n");
%!     fputs(fd, "                 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "                 euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 2, static,\n");
%!     fputs(fd, "                 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "                 euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., const, F2;\n");
%!     fputs(fd, "         joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "         joint: 2, clamp, 2, node, node;\n");
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
%!   t = mbdyn_post_load_output_out(options.output_file);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   forces = mbdyn_post_load_log_force(options.output_file);
%!   [joint_id, Fl, Fg] = mbdyn_post_load_output_jnt(options.output_file, 1:2, numel(t));
%!   X1 = log_dat.nodes(1).X0;
%!   X2 = log_dat.nodes(2).X0;
%!   R1 = log_dat.nodes(1).R0;
%!   R2 = log_dat.nodes(2).R0;
%!   F1 = repmat([log_dat.vars.F1, 0, 0], numel(t), 1);
%!   F2 = repmat([0, log_dat.vars.F2, 0], numel(t), 1);
%!   l1 = forces(1).arm1;
%!   l2 = forces(2).arm1;
%!   M1 = cross(repmat((R1 * l1).', numel(t), 1), F1);
%!   M2 = cross(repmat((R2 * l2).', numel(t), 1), F2);
%!   tol = eps^0.3;
%!   assert_simple(Fg{1}(:, 1:3), F1, tol * norm(F1));
%!   assert_simple(Fg{1}(:, 4:6), M1, tol * norm(M1));
%!   assert_simple(Fl{1}(:, 1:3), F1 * R1, tol * norm(F1));
%!   assert_simple(Fl{1}(:, 4:6), M1 * R1, tol * norm(M1));
%!   assert_simple(Fg{2}(:, 1:3), F2, tol * norm(F2));
%!   assert_simple(Fg{2}(:, 4:6), M2, tol * norm(M2));
%!   assert_simple(Fl{2}(:, 1:3), F2 * R2, tol * norm(F2));
%!   assert_simple(Fl{2}(:, 4:6), M2 * R2, tol * norm(M2));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
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
