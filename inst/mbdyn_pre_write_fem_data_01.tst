## mbdyn_pre_write_fem_data.tst:01
%!test
%! try
%! ## TEST1
%! fd = -1;
%! omega1 = [  0,  0,  0,  0,   0];
%! omega2 = [0.5, 50,  5, 10, 100];
%! omega3 = [0.2, 10, 20, 50, 500];
%! m = 0.1;
%! s1 = 0;
%! s2 = 0;
%! s3 = 0.05;
%! J1 = 5e-3;
%! J2 = 5e-3;
%! J3 = 10e-3;
%! h3 = 0.075;
%! nsteps = 14400;
%! nrev = 1;
%! n0 = 0.1;
%! tolres = 1e-9;
%! silent = true;
%! for i=1:numel(omega1)
%! %unwind_protect
%!     unwind_protect
%!       [fd, fname1] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_write_fem_data1_XXXXXX"));
%!       if (fd == -1)
%!         error("failed to open temporary file");
%!       endif
%!       fprintf(fd, " set: real omega1 = %e;\n", omega1(i));
%!       fprintf(fd, " set: real omega2 = %e;\n", omega2(i));
%!       fprintf(fd, " set: real omega3 = %e;\n", omega3(i));
%!       fputs(fd, " set: real omega = sqrt(omega1^2 + omega2^2 + omega3^2);\n");
%!       fprintf(fd, " set: real m = %e;\n", m);
%!       fprintf(fd, " set: real s1 = %e;\n", s1);
%!       fprintf(fd, " set: real s2 = %e;\n", s2);
%!       fprintf(fd, " set: real s3 = %e;\n", s3);
%!       fprintf(fd, " set: real J1 = %e;\n", J1);
%!       fprintf(fd, " set: real J2 = %e;\n", J2);
%!       fprintf(fd, " set: real J3 = %e;\n", J3);
%!       fprintf(fd, " set: real h3 = %e;\n", h3);
%!       fprintf(fd, " set: real nsteps = %e;\n", nsteps);
%!       fprintf(fd, " set: real nrev = %e;\n", nrev);
%!       fprintf(fd, " set: real n0 = %e;\n", n0);
%!       fprintf(fd, " set: real tolres = %e;\n", tolres);
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "         problem: initial value;\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "         initial time: 0;\n");
%!       fputs(fd, "         final time: 2. * pi * nrev / abs(omega);\n");
%!       fputs(fd, "         time step: 2. * pi / (nsteps * abs(omega));\n");
%!       fputs(fd, "         linear solver: naive, colamd, scale, row max column max, always;\n");
%!       fputs(fd, "         method: ms, 0.6;\n");
%!       fputs(fd, "         max iterations: 10;\n");
%!       fputs(fd, "         tolerance: tolres;\n");
%!       fputs(fd, "         threads: assembly, 1;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: auto;\n");
%!       fputs(fd, "         nonlinear solver: nox, modified, 10,\n");
%!       fputs(fd, "             keep jacobian matrix,\n");
%!       fputs(fd, "             inner iterations before assembly, 6,\n");
%!       fputs(fd, "             jacobian operator, newton krylov,\n");
%!       fputs(fd, "             solver, line search based,\n");
%!       fputs(fd, "             forcing term, type 2,\n");
%!       fputs(fd, "             direction, newton,\n");
%!       fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!       fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!       fputs(fd, "             linear solver, gmres,\n");
%!       fputs(fd, "             linear solver max iterations, 12,\n");
%!       fputs(fd, "             krylov subspace size, 12;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "     use automatic differentiation;\n");
%!       fputs(fd, "     max iterations: 0;\n");
%!       fputs(fd, "     structural nodes: 1;\n");
%!       fputs(fd, "     joints: 1;\n");
%!       fputs(fd, "     rigid bodies: 1;\n");
%!       fputs(fd, "     gravity;\n");
%!       fputs(fd, "     output meter: closest next, 2 * pi * n0 / abs(omega), forever, 2 * pi / (100 * abs(omega));\n");
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fputs(fd, "         structural: 1, dynamic, 0., 0., h3, eye, reference, global, h3 * omega2, -h3 * omega1, 0., reference, global, omega1, omega2, omega3, accelerations, yes;\n");
%!       fputs(fd, " end: nodes;\n");
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "         body: 1, 1, m, reference, node, s1, s2, s3, diag, J1, J2, J3;\n");
%!       fputs(fd, "         joint: 1, total pin joint, 1,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position constraint, active, active, active, null,\n");
%!       fputs(fd, "            orientation constraint, inactive, inactive, inactive, null;\n");
%!       fputs(fd, "            gravity: uniform, 0., 0., -1., 9.81;\n");
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!     end_unwind_protect
%!     options.output_file = fname1;
%!     options.verbose = false;
%!     if (silent)
%!       options.logfile = [fname1, ".stdout"];
%!     endif
%!     mbdyn_solver_run(fname1, options);
%!     [t1, trajectory1, deformation1, velocity1, acceleration1, node_id1] = mbdyn_post_load_output_struct(fname1);
%! %unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname1);
%!       files = dir([fname1, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%! %end_unwind_protect
%!   fd = -1;
%! %unwind_protect
%!     unwind_protect
%!       [fd, fname2] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_write_fem_data2_XXXXXX"));
%!       if (fd == -1)
%!         error("failed to open temporary file");
%!       endif
%!       fprintf(fd, " set: real omega1 = %e;\n", omega1(i));
%!       fprintf(fd, " set: real omega2 = %e;\n", omega2(i));
%!       fprintf(fd, " set: real omega3 = %e;\n", omega3(i));
%!       fputs(fd, " set: real omega = sqrt(omega1^2 + omega2^2 + omega3^2);\n");
%!       fprintf(fd, " set: real m = %e;\n", m);
%!       fprintf(fd, " set: real s1 = %e;\n", s1);
%!       fprintf(fd, " set: real s2 = %e;\n", s2);
%!       fprintf(fd, " set: real s3 = %e;\n", s3);
%!       fprintf(fd, " set: real J1 = %e;\n", J1);
%!       fprintf(fd, " set: real J2 = %e;\n", J2);
%!       fprintf(fd, " set: real J3 = %e;\n", J3);
%!       fprintf(fd, " set: real h3 = %e;\n", h3);
%!       fprintf(fd, " set: real nsteps = %e;\n", nsteps);
%!       fprintf(fd, " set: real nrev = %e;\n", nrev);
%!       fprintf(fd, " set: real n0 = %e;\n", n0);
%!       fprintf(fd, " set: real tolres = %e;\n", tolres);
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "         problem: initial value;\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "         initial time: 0;\n");
%!       fputs(fd, "         final time: 2. * pi * nrev / abs(omega);\n");
%!       fputs(fd, "         time step: 2. * pi / (nsteps * abs(omega));\n");
%!       fputs(fd, "         linear solver: naive, colamd, scale, row max column max, always;\n");
%!       fputs(fd, "         method: ms, 0.6;\n");
%!       fputs(fd, "         max iterations: 10;\n");
%!       fputs(fd, "         tolerance: tolres;\n");
%!       fputs(fd, "         threads: assembly, 1;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: auto;\n");
%!       fputs(fd, "         nonlinear solver: nox, modified, 10,\n");
%!       fputs(fd, "             keep jacobian matrix,\n");
%!       fputs(fd, "             inner iterations before assembly, 6,\n");
%!       fputs(fd, "             jacobian operator, newton krylov,\n");
%!       fputs(fd, "             solver, line search based,\n");
%!       fputs(fd, "             forcing term, type 2,\n");
%!       fputs(fd, "             direction, newton,\n");
%!       fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!       fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!       fputs(fd, "             linear solver, gmres,\n");
%!       fputs(fd, "             linear solver max iterations, 12,\n");
%!       fputs(fd, "             krylov subspace size, 12;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "     use automatic differentiation;\n");
%!       fputs(fd, "     max iterations: 0;\n");
%!       fputs(fd, "     structural nodes: 1;\n");
%!       fputs(fd, "     joints: 2;\n");
%!       fputs(fd, "     rigid bodies: 1;\n");
%!       fputs(fd, "     gravity;\n");
%!       fputs(fd, "     output meter: closest next, 2 * pi * n0 / abs(omega), forever, 2 * pi / (100 * abs(omega));\n");
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fputs(fd, "         structural: 1, modal, 0., 0., h3, eye, reference, global, h3 * omega2, -h3 * omega1, 0., reference, global, omega1, omega2, omega3, accelerations, yes;\n");
%!       fputs(fd, " end: nodes;\n");
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "         body: 1, 1, m, reference, node, s1, s2, s3, diag, J1, J2, J3;\n");
%!       fputs(fd, "         joint: 1, total pin joint, 1,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position constraint, active, active, active, null,\n");
%!       fputs(fd, "            orientation constraint, inactive, inactive, inactive, null;\n");
%!       fputs(fd, "         joint: 2, modal, 1,\n");
%!       fputs(fd, "            1, from file,\n");
%!       fputs(fd, "            damping from file,\n");
%!       fprintf(fd, "          \"%s.fem\",\n", fname2);
%!       fputs(fd, "            origin node, 1000,\n");
%!       fputs(fd, "            0;\n");
%!       fputs(fd, "            gravity: uniform, 0., 0., -1., 9.81;\n");
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!     end_unwind_protect
%!     MRED = 1;
%!     DRED = 0;
%!     SRED = 1e10;
%!     TRED = zeros(6, 1);
%!     X0 = zeros(3, 1);
%!     URED0 = 0;
%!     DURED_DT0 = 0;
%!     DIAGM = [];
%!     M = 0;
%!     XGC = zeros(3, 1);
%!     JGC = zeros(3, 3);
%!     NODE_LIST = 1000;
%!     mbdyn_pre_write_fem_data([fname2, ".fem"], MRED, DRED, SRED, TRED, X0, URED0, DURED_DT0, DIAGM, M, XGC, JGC, NODE_LIST);
%!     options.output_file = fname2;
%!     if (silent)
%!       options.logfile = [fname2, ".stdout"];
%!     endif
%!     mbdyn_solver_run(fname2, options);
%!     [t2, trajectory2, deformation2, velocity2, acceleration2, node_id2] = mbdyn_post_load_output_struct(fname2);
%! %unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname2);
%!       files = dir([fname2, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%! %end_unwind_protect
%!   fd = -1;
%! %unwind_protect
%!     unwind_protect
%!       [fd, fname3] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_write_fem_data3_XXXXXX"));
%!       if (fd == -1)
%!         error("failed to open temporary file");
%!       endif
%!       fprintf(fd, " set: real omega1 = %e;\n", omega1(i));
%!       fprintf(fd, " set: real omega2 = %e;\n", omega2(i));
%!       fprintf(fd, " set: real omega3 = %e;\n", omega3(i));
%!       fputs(fd, " set: real omega = sqrt(omega1^2 + omega2^2 + omega3^2);\n");
%!       fprintf(fd, " set: real m = %e;\n", m);
%!       fprintf(fd, " set: real s1 = %e;\n", s1);
%!       fprintf(fd, " set: real s2 = %e;\n", s2);
%!       fprintf(fd, " set: real s3 = %e;\n", s3);
%!       fprintf(fd, " set: real J1 = %e;\n", J1);
%!       fprintf(fd, " set: real J2 = %e;\n", J2);
%!       fprintf(fd, " set: real J3 = %e;\n", J3);
%!       fprintf(fd, " set: real h3 = %e;\n", h3);
%!       fprintf(fd, " set: real nsteps = %e;\n", nsteps);
%!       fprintf(fd, " set: real nrev = %e;\n", nrev);
%!       fprintf(fd, " set: real n0 = %e;\n", n0);
%!       fprintf(fd, " set: real tolres = %e;\n", tolres);
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "         problem: initial value;\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "         initial time: 0;\n");
%!       fputs(fd, "         final time: 2. * pi * nrev / abs(omega);\n");
%!       fputs(fd, "         time step: 2. * pi / (nsteps * abs(omega));\n");
%!       fputs(fd, "         linear solver: naive, colamd, scale, row max column max, always;\n");
%!       fputs(fd, "         method: ms, 0.6;\n");
%!       fputs(fd, "         max iterations: 10;\n");
%!       fputs(fd, "         tolerance: tolres;\n");
%!       fputs(fd, "         threads: assembly, 1;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: auto;\n");
%!       fputs(fd, "         nonlinear solver: nox, modified, 10,\n");
%!       fputs(fd, "             keep jacobian matrix,\n");
%!       fputs(fd, "             inner iterations before assembly, 6,\n");
%!       fputs(fd, "             jacobian operator, newton krylov,\n");
%!       fputs(fd, "             solver, line search based,\n");
%!       fputs(fd, "             forcing term, type 2,\n");
%!       fputs(fd, "             direction, newton,\n");
%!       fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!       fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!       fputs(fd, "             linear solver, gmres,\n");
%!       fputs(fd, "             linear solver max iterations, 12,\n");
%!       fputs(fd, "             krylov subspace size, 12;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "     use automatic differentiation;\n");
%!       fputs(fd, "     max iterations: 0;\n");
%!       fputs(fd, "     structural nodes: 1;\n");
%!       fputs(fd, "     joints: 2;\n");
%!       fputs(fd, "     gravity;\n");
%!       fputs(fd, "     output meter: closest next, 2 * pi * n0 / abs(omega), forever, 2 * pi / (100 * abs(omega));\n");
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fputs(fd, "         structural: 1, modal, 0., 0., h3, eye, reference, global, h3 * omega2, -h3 * omega1, 0., reference, global, omega1, omega2, omega3, accelerations, yes;\n");
%!       fputs(fd, " end: nodes;\n");
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "         joint: 1, total pin joint, 1,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position constraint, active, active, active, null,\n");
%!       fputs(fd, "            orientation constraint, inactive, inactive, inactive, null;\n");
%!       fputs(fd, "         joint: 2, modal, 1,\n");
%!       fputs(fd, "            1, from file,\n");
%!       fputs(fd, "            damping from file,\n");
%!       fprintf(fd, "          \"%s.fem\",\n", fname3);
%!       fputs(fd, "            origin node, 1000,\n");
%!       fputs(fd, "            0;\n");
%!       fputs(fd, "            gravity: uniform, 0., 0., -1., 9.81;\n");
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!     end_unwind_protect
%!     MRED = 1;
%!     DRED = 0;
%!     SRED = 1e10;
%!     TRED = zeros(6, 1);
%!     X0 = zeros(3, 1);
%!     URED0 = 0;
%!     DURED_DT0 = 0;
%!     DIAGM = [];
%!     M = m;
%!     XGC = [s1; s2; s3];
%!     JGC = diag([J1, J2, J3]);
%!     NODE_LIST = 1000;
%!     mbdyn_pre_write_fem_data([fname3, ".fem"], MRED, DRED, SRED, TRED, X0, URED0, DURED_DT0, DIAGM, M, XGC, JGC, NODE_LIST);
%!     options.output_file = fname3;
%!     if (silent)
%!       options.logfile = [fname3, ".stdout"];
%!     endif
%!     mbdyn_solver_run(fname3, options);
%!     [t3, trajectory3, deformation3, velocity3, acceleration3, node_id3] = mbdyn_post_load_output_struct(fname3);
%! %unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname3);
%!       files = dir([fname3, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%! %end_unwind_protect
%!   tol = 1e-4;
%!   for j=1:numel(node_id2)
%!     assert_simple(trajectory2{j}, trajectory1{j}, tol * max(max(abs(trajectory1{j}))));
%!     assert_simple(velocity2{j}, velocity1{j}, tol * max(max(abs(velocity1{j}))));
%!     assert_simple(acceleration2{j}, acceleration1{j}, tol * max(max(abs(acceleration1{j}))));
%!     assert_simple(trajectory3{j}, trajectory1{j}, tol * max(max(abs(trajectory1{j}))));
%!     assert_simple(velocity3{j}, velocity1{j}, tol * max(max(abs(velocity1{j}))));
%!     assert_simple(acceleration3{j}, acceleration1{j}, tol * max(max(abs(acceleration1{j}))));
%!   endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
