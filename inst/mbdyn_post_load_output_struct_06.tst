## mbdyn_post_load_output_struct.tst:06
%!test
%! try
%! f_plot = false;
%! fd = -1;
%! %unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_struct_XXXXXX"));
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
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, null, 1., 0., 0, 100.;\n");
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
%!   [t, trajectory, deformation, velocity, acceleration, node_id, force,  force_id, force_node_id, orientation_description, t_trc, output_flag] = mbdyn_post_load_output_struct(fname, [1], [1], false);
%!   g = -9.81;
%!   F = 100;
%!   m = 1;
%!   assert_simple(all(output_flag));
%!   assert_simple(node_id, int32(1));
%!   assert_simple(t, (0:0.1:1.1).', 1e-5);
%!   assert_simple(trajectory{1}, ...
%!          [0.5 * F / m * t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g * t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert_simple(velocity{1}, ...
%!         [t * F / m, ...
%!          zeros(numel(t), 1), g * t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert_simple(acceleration{1}, ...
%!         [F / m * ones(numel(t), 1), ...
%!          zeros(numel(t), 1), ...
%!          repmat(g, numel(t), 1), ...
%!          zeros(numel(t), 3)], ...
%!          1e-5);
%!  assert_simple(force{1}(:, 1:3), [repmat(F, numel(t), 1), zeros(numel(t), 2)]);
%!  if (f_plot)
%!  figure("visible", "off");
%!  subplot(3, 1, 1);
%!  hold on;
%!  for i=1:3
%!    plot(t, trajectory{1}(:, i), sprintf("-;x%d;%d", i, i));
%!  endfor
%!  xlabel("t [s]");
%!  ylabel("x [m]");
%!  grid on;
%!  grid minor on;
%!  title("trajectory");
%!  subplot(3, 1, 2);
%!  hold on;
%!  for i=1:3
%!    plot(t, velocity{1}(:, i), sprintf("-;v%d;%d", i, i));
%!  endfor
%!  xlabel("t [s]");
%!  ylabel("v [m/s]");
%!  grid on;
%!  grid minor on;
%!  title("velocity");
%!  subplot(3, 1, 3);
%!  hold on;
%!  for i=1:3
%!    plot(t, acceleration{1}(:, i), sprintf("-;a%d;%d", i, i));
%!  endfor
%!  xlabel("t [s]");
%!  ylabel("a [m/s^2]");
%!  grid on;
%!  grid minor on;
%!  title("acceleration");
%!  endif
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
