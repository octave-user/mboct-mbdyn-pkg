## mbdyn_post_load_output_ine.tst:03
%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_ine_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real m = 1.234;\n");
%!     fputs(fd, " set: real J11 = 0.011;\n");
%!     fputs(fd, " set: real J22 = 0.022;\n");
%!     fputs(fd, " set: real J33 = 0.033;\n");
%!     fputs(fd, " set: real omega1 = 1.123;\n");
%!     fputs(fd, " set: real omega2 = 2.234;\n");
%!     fputs(fd, " set: real omega3 = 3.456;\n");
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
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 omega1, omega2, omega3,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m, null, diag, J11, J22, J33;\n");
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
%!   [t] = mbdyn_post_load_output_out(fname, 1024, false);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   g = log_dat.vars.g;
%!   m = log_dat.vars.m;
%!   J11 = log_dat.vars.J11;
%!   J22 = log_dat.vars.J22;
%!   J33 = log_dat.vars.J33;
%!   J = diag([J11, J22, J33]);
%!   omega = [log_dat.vars.omega1; log_dat.vars.omega2; log_dat.vars.omega3];
%!   vref = [0; 0; -g] * t.';
%!   betaref = m * vref;
%!   gammaref = J * omega;
%!   betadotref = [0; 0; -m * g];
%!   gammadotref = zeros(3, 1);
%!   [node_id, beta, gamma, beta_dot, gamma_dot] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta, gamma, beta_dot] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta, gamma] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   assert_simple(node_id, int32(1));
%!   tol = eps^0.4;
%!   assert_simple(beta{1}, betaref.', tol * max(max(abs(betaref))));
%!   assert_simple(gamma{1}, repmat(gammaref.', numel(t), 1), tol * norm(gammaref));
%!   assert_simple(beta_dot{1}, repmat(betadotref.', numel(t), 1), tol * norm(betadotref));
%!   assert_simple(gamma_dot{1}, repmat(gammadotref.', numel(t), 1), tol * max([1,norm(gammadotref)]));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
