## mbdyn_post_offset_ref_frame_node.tst:03
%!test
%! try
%! fd = -1;
%! %unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_offset_ref_frame_nodes_XXXXXX"));
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
%!     fputs(fd, "     default output: reference frames;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         set: integer ref_id_0 = 1001;\n");
%!     fputs(fd, "         set: integer ref_id_1 = 1002;\n");
%!     fputs(fd, "         set: integer ref_id_2 = 1003;\n");
%!     fputs(fd, "         set: integer node_id_1 = 2001;\n");
%!     fputs(fd, "         reference: ref_id_0,\n");
%!     fputs(fd, "                 position, reference, global, 0.4, 0.2, 0.3,\n");
%!     fputs(fd, "                 orientation, reference, global, euler123, pi / 2, pi / 4, pi,\n");
%!     fputs(fd, "                 velocity, reference, global, null,\n");
%!     fputs(fd, "                 angular velocity, reference, global, null;\n");
%!     fputs(fd, "         reference: ref_id_1,\n");
%!     fputs(fd, "                 position, reference, ref_id_0, 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_0, euler123, 0., 0., pi / 2,\n");
%!     fputs(fd, "                 velocity, reference, ref_id_0, null,\n");
%!     fputs(fd, "                 angular velocity, reference, ref_id_0, null;\n");
%!     fputs(fd, "         reference: ref_id_2,\n");
%!     fputs(fd, "                 position, reference, ref_id_1, -0.1, 0.5, 2.2,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_1, euler123, pi / 3, -pi/2, -pi / 3,\n");
%!     fputs(fd, "                 velocity, reference, ref_id_1, null,\n");
%!     fputs(fd, "                 angular velocity, reference, ref_id_1, null;\n");
%!     fputs(fd, "         structural: node_id_1, dynamic,\n");
%!     fputs(fd, "                 position, reference, ref_id_0, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_0, eye,\n");
%!     fputs(fd, "                 velocity, reference, ref_id_0, null,\n");
%!     fputs(fd, "                 angular velocity, reference, ref_id_0, null,\n");
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
%! # options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [ref.ref_id, ref.position, ref.orientation, ref.velocity, ref.angular_velocity] = mbdyn_post_load_output_rfm(options.output_file);
%!   X = mbdyn_post_offset_ref_frame_node(ref, log_dat.nodes, log_dat.vars, "ref_id_2", "node_id_1", "node");
%!   assert_simple(X, [0.6; 2.1; 5.5], 1e-6);
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
