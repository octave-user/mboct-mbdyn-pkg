## mbdyn_post_inertia_compute.tst:03
%!test
%! try
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_inertia_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real m1 = 123;\n");
%!     fputs(fd, " set: real m2 = 0.456;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1e-3;\n");
%!     fputs(fd, "         time step: 1e-4;\n");
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
%!     fputs(fd, "     output precision: 16;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 0.5, 0.6, 0.7,\n");
%!     fputs(fd, "                 euler123, 0.3 * pi, 0.7 * pi, -1.4 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "                 euler123, pi / 4., -pi / 3., pi / 2.,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, reference, global, null, diag, 0.456, 0.789, 0.123, inertial, reference, global, eye;\n");
%!     fputs(fd, "         body: 2, 2, m2, reference, global, null, diag, 45, 78, 123, inertial, reference, global, eye;\n");
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
%!   bodies = mbdyn_post_load_log_body(fname);
%!   log_dat = mbdyn_post_load_log(fname);
%!   body_groups(1).labels = [1, 2];
%!   body_groups(1).name = "body 1 + body 2";
%!   inertia = mbdyn_post_inertia_compute(body_groups, bodies, log_dat.nodes);
%!   mbdyn_post_inertia_print(inertia, [fname, "inertia.dat"]);
%!   tol = 1e-5;
%!   assert_simple(inertia(1).dm, 123.456, tol * 123.456);
%!   assert_simple(inertia(1).Xgc, zeros(3, 1), tol);
%!   assert_simple(inertia(1).J, diag([45.456, 78.789, 123.123]), tol * 123.123);
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
