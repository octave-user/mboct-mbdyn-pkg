## mbdyn_post_load_log.tst:03
%!test
%! try
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_log_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 20;\n");
%!     fputs(fd, " set: real m1 = 1;\n");
%!     fputs(fd, " set: real m2 = 2;\n");
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
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1., 2., 3.,\n");
%!     fputs(fd, "                 euler123, 0.1, 0.2, 0.3,\n");
%!     fputs(fd, "                 10., 20., 30.,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         body: 2, 2, m2, null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, null, 1., 0., 0, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, null, 0., 1., 0, F2;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
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
%!   log_dat = mbdyn_post_load_log(fname);
%!   assert_simple(log_dat.nodes(1).label, int32(1));
%!   assert_simple(log_dat.nodes(1).X0, zeros(3, 1));
%!   assert_simple(log_dat.nodes(1).R0, eye(3));
%!   assert_simple(log_dat.nodes(2).label, int32(2));
%!   assert_simple(log_dat.nodes(2).X0, [1; 2; 3]);
%!   assert_simple(log_dat.nodes(2).R0, euler123_to_rotation_matrix([0.1; 0.2; 0.3]), 1e-6);
%!   assert_simple(log_dat.vars.m1, 1);
%!   assert_simple(log_dat.vars.F1, 100);
%!   assert_simple(log_dat.vars.g, 9.81);
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
