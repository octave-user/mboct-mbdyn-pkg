## mbdyn_post_abs_to_rel.tst:11
%!test
%! try
%! ## TEST11
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!  fputs(fd, "set: real r = 10.;\n");
%!  fputs(fd, "set: real m = 1.;\n");
%!  fputs(fd, "set: real gx = 0;\n");
%!  fputs(fd, "set: real gy = -9.81;\n");
%!  fputs(fd, "set: real gz = 0.;\n");
%!  fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!  fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!  fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!  fputs(fd, "set: real initial_time = 0;\n");
%!  fputs(fd, "set: real final_time = 0.1;\n");
%!  fputs(fd, "set: integer N = 720;\n");
%!  fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!  fputs(fd, "begin: data;\n");
%!  fputs(fd, "	problem: initial value; # the default\n");
%!  fputs(fd, "end: data;\n");
%!  fputs(fd, "begin: initial value;\n");
%!  fputs(fd, "    threads: assembly, 1;\n");
%!  fputs(fd, "	initial time: initial_time;\n");
%!  fputs(fd, "	final time: final_time;\n");
%!  fputs(fd, "	time step: time_step;\n");
%!  fputs(fd, "    max iterations: 100;\n");
%!  fputs(fd, "    method: ms, 0.;\n");
%!  fputs(fd, "    nonlinear solver: nox, modified, 1000,\n");
%!  fputs(fd, "        keep jacobian matrix,\n");
%!  fputs(fd, "        inner iterations before assembly, 6,\n");
%!  fputs(fd, "        jacobian operator, newton,\n");
%!  fputs(fd, "        solver, trust region based,\n");
%!  fputs(fd, "        forcing term, type 2,\n");
%!  fputs(fd, "        direction, newton,\n");
%!  fputs(fd, "        verbose, yes,\n");
%!  fputs(fd, "        print convergence info, no,\n");
%!  fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!  fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!  fputs(fd, "        linear solver, gmres,\n");
%!  fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!  fputs(fd, "        linear solver max iterations, 10,\n");
%!  fputs(fd, "        krylov subspace size, 10,\n");
%!  fputs(fd, "        minimum step, 1e-12,\n");
%!  fputs(fd, "        use preconditioner as solver, no;\n");
%!  fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!  fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!  fputs(fd, "        derivatives max iterations: 3000;\n");
%!  fputs(fd, "        derivatives coefficient: 1e-4;\n");
%!  fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!  fputs(fd, "end: initial value;\n");
%!  fputs(fd, "begin: control data;\n");
%!  fputs(fd, "       tolerance: 1e-10;\n");
%!  fputs(fd, "       max iterations: 100;\n");
%!  fputs(fd, "	structural nodes: 1;\n");
%!  fputs(fd, "	rigid bodies: 1;\n");
%!  fputs(fd, "	joints: 1;\n");
%!  fputs(fd, "    forces: 1;\n");
%!  fputs(fd, "    gravity;\n");
%!  fputs(fd, "end: control data;\n");
%!  fputs(fd, "reference: 1,\n");
%!  fputs(fd, "    position, reference, global, 0.05, 0.02, 0.03,\n");
%!  fputs(fd, "    orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!  fputs(fd, "    velocity, reference, global, null,\n");
%!  fputs(fd, "    angular velocity, reference, global, null;\n");
%!  fputs(fd, "reference: 2,\n");
%!  fputs(fd, "   position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!  fputs(fd, "   orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!  fputs(fd, "   velocity, reference, 1, null,\n");
%!  fputs(fd, "   angular velocity, reference, 1, 0., 0., omega0;\n");
%!  fputs(fd, "begin: nodes;\n");
%!  fputs(fd, "	structural: 1, dynamic,\n");
%!  fputs(fd, "		position, reference, 2, 0.07, 0.08, 0.02,\n");
%!  fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!  fputs(fd, "		velocity, reference, 2, null,\n");
%!  fputs(fd, "		angular velocity, reference, 2, null;\n");
%!  fputs(fd, "end: nodes;\n");
%!  fputs(fd, "begin: elements;\n");
%!  fputs(fd, "	body: 1, 1, \n");
%!  fputs(fd, "		# mass\n");
%!  fputs(fd, "		m,\n");
%!  fputs(fd, "		# center of mass\n");
%!  fputs(fd, "		null,\n");
%!  fputs(fd, "		# inertia matrix\n");
%!  fputs(fd, "		diag, .1, 1., 100.;\n");
%!  fputs(fd, "	joint: 1, total pin joint,\n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		# relative position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "		# relative orientation matrix\n");
%!  fputs(fd, "		position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!  fputs(fd, "		# absolute position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "                position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, eye,                \n");
%!  fputs(fd, "                position constraint,\n");
%!  fputs(fd, "                active, active, active, null,\n");
%!  fputs(fd, "                orientation constraint,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                inactive,\n");
%!  fputs(fd, "                null;\n");
%!  fputs(fd, "	couple: 1, absolute, \n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		position, reference, node, null, \n");
%!  fputs(fd, "		component, 100., \n");
%!  fputs(fd, "			   100., \n");
%!  fputs(fd, "			   1000.;\n");
%!  fputs(fd, "    gravity: uniform, component,\n");
%!  fputs(fd, "        gx, gy, gz;\n");
%!  fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   if (~options.verbose)
%! # options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
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
