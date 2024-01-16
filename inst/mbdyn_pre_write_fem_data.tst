%!test
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
%!   unwind_protect
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
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname1);
%!       files = dir([fname1, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
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
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname2);
%!       files = dir([fname2, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
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
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname3);
%!       files = dir([fname3, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
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

%!test
%! ## TEST2
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
%!   unwind_protect
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
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname1);
%!       files = dir([fname1, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
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
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname2);
%!       files = dir([fname2, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
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
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname3);
%!       files = dir([fname3, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
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


%!test
%! ## TEST3
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
%!   unwind_protect
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
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname1);
%!       files = dir([fname1, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
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
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname2);
%!       files = dir([fname2, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
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
%!       fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
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
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname3);
%!       files = dir([fname3, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
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

%!test
%! ## TEST4
%! pkg load mboct-fem-pkg;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   a = 800e-3 / SI_unit_meter;
%!   b = 40e-3 / SI_unit_meter;
%!   c = 10e-3 / SI_unit_meter;
%!   d = 0e-3 / SI_unit_meter;
%!   h = c;
%!   options.interactive = false;
%!   options.plot = true;
%!   options.verbose = false;
%!   options.number_of_beams = int32(40);
%!   options.number_of_threads = mbdyn_solver_num_threads_default();
%!   if (options.plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fprintf(fd, "a = %.16e;\n", a);
%!     fprintf(fd, "b = %.16e;\n", b);
%!     fprintf(fd, "c = %.16e;\n", c);
%!     fprintf(fd, "h = %.16e;\n", h);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(2) = {0,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(3) = {a,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(4) = {a, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Line Loop(1) = {1, 2, 3, 4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,c}{Surface{1}; Layers{Max(1, Round(c / h))}; Recombine;};\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {v1[2]};\n");
%!     fputs(fd, "Physical Surface(2) = {v1[4]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   mesh.material_data.E = 70000e6 / SI_unit_pascal;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 2700 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data.alpha = 0e-5 / (1 / SI_unit_second);
%!   mesh.material_data.beta = 0e-5 / (SI_unit_second);
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   grp_idx_beam = find([[mesh.groups.iso20].id] == 1);
%!   grp_idx_clamp = find([[mesh.groups.quad8].id] == 1);
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_beam).elements) = 1;
%!   cms_opt.number_of_threads = options.number_of_threads;
%!   cms_opt.algorithm = "diag-shift-invert";
%!   cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!   cms_opt.nodes.modal.name = "node_id_modal";
%!   cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!   cms_opt.nodes.interfaces.name = "node_id_interface1";
%!   mesh.nodes(cms_opt.nodes.modal.number, 1:3) = [0, 0, 0];
%!   mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = [a + d, 0, 0];
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, 2, cms_opt.nodes.interfaces.number, "quad8");
%!   cms_opt.refine_max_iter = 30;
%!   cms_opt.pre_scaling = false;
%!   cms_opt.solver = "pardiso";
%!   cms_opt.modes.number = 20;
%!   cms_opt.tolerance_tau = -1;
%!   cms_opt.element.name = "elem_id_modal";
%!   cms_opt.create_binary = true;
%!   cms_opt.use_binary = true;
%!   cms_opt.update_binary = true;
%!   cms_opt.invariants = true;
%!   #cms_opt.enable_KTAU0WP = [true; false(2, 1)];
%!   #cms_opt.enable_KTAU0VP = [true; false(2, 1)];
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, :) = true;
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, :) = true;
%!   [mesh_cms, mat_ass_cms, dof_map_cms, sol_eig_cms, cms_opt, sol_tau_cms] = fem_cms_create2(mesh, load_case_dof, cms_opt);
%!   fem_cms_export(filename, mesh_cms, dof_map_cms, mat_ass_cms, cms_opt);
%!   pert.omega = [1e2; 3e2; 2e2] / (SI_unit_rad / SI_unit_second);
%!   pert.omegadot = [1e5; 1e3; 3e3] / (SI_unit_rad / SI_unit_second^2);
%!   pert.loads = [[1e4; 1e3; 1e2] / (SI_unit_newton);
%!                 [1e2; 1e1; 1e1] / (SI_unit_newton * SI_unit_meter)];
%!   pert.g = [1e4; -1e3; -1e2] / (SI_unit_meter / SI_unit_second^2);
%!   pert.a = [-1e4; 1e3; 1e2] / (SI_unit_meter / SI_unit_second^2);
%!   empty_cell = cell(7, 3, 2);
%!   res = struct("info", empty_cell, ...
%!                "t", empty_cell, ...
%!                "trajectory", empty_cell, ...
%!                "deformation", empty_cell, ...
%!                "velocity", empty_cell, ...
%!                "acceleration", empty_cell, ...
%!                "node_id", empty_cell, ...
%!                "force", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "orientation_description", empty_cell, ...
%!                "drive_id", empty_cell, ...
%!                "drive_value", empty_cell, ...
%!                "modal", empty_cell);
%!   empty_cell = cell(7, 3);
%!   param = struct("omega", empty_cell, ...
%!                  "omegadot", empty_cell, ...
%!                  "F1", empty_cell, ...
%!                  "M1", empty_cell, ...
%!                  "a", empty_cell, ...
%!                  "g", empty_cell, ...
%!                  "t1", empty_cell, ...
%!                  "holonomic", empty_cell);
%!   idx_j = 1:rows(param);
%!   idx_k = 1:columns(param);
%!   for j=idx_j
%!     for k=idx_k
%!       param(j, k).omega = zeros(3, 1);
%!       param(j, k).omegadot = zeros(3, 1);
%!       param(j, k).F1 = zeros(3, 1);
%!       param(j, k).M1 = zeros(3, 1);
%!       param(j, k).a = zeros(3, 1);
%!       param(j, k).g = zeros(3, 1);
%!       param(j, k).holonomic = false;
%!       param(j, k).gamma = zeros(3, 1);
%!       param(j, k).N = 50;
%!       switch (j)
%!         case 1
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).gamma = [20; 45; 30] * pi / 180;
%!         case 2
%!           param(j, k).omega(k) = pert.omega(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           param(j, k).holonomic = true;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 2000;
%!         case 3
%!           param(j, k).omegadot(k) = pert.omegadot(k);
%!           param(j, k).t1 = 1e-2 / SI_unit_second;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 200;
%!         case 4
%!           param(j, k).F1(k) = pert.loads(k);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           switch (k)
%!           case 1
%!             param(j, k).gamma(3) = 5 * pi / 180;
%!           case 2
%!             param(j, k).gamma(3) = 45 * pi / 180;
%!           case 3
%!             param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 5
%!           param(j, k).M1(k) = pert.loads(k + 3);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           param(j, k).gamma(1) = 30 * pi / 180;
%!         case 6
%!           param(j, k).a(k) = pert.a(k);
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).N = 200;
%!           switch (k)
%!           case 1
%!             param(j, k).gamma(3) = 5 * pi / 180;
%!           case 2
%!             param(j, k).gamma(3) = 45 * pi / 180;
%!           case 3
%!             param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 7
%!           param(j, k).g(k) = pert.g(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           switch (k)
%!           case 1
%!             param(j, k).gamma(3) = 5 * pi / 180;
%!           case 2
%!             param(j, k).gamma(3) = 45 * pi / 180;
%!           case 3
%!             param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!       endswitch
%!       for l=1:2
%!         fd = -1;
%!         unwind_protect
%!           filename_mbdyn = sprintf("%s_%d_%d_%d.mbdyn", filename, j, k, l);
%!           [fd, msg] = fopen(filename_mbdyn, "w");
%!           if (fd == -1)
%!             error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!           endif
%!           fprintf(fd, "set: real a = %.16e;\n", a);
%!           fprintf(fd, "set: real b = %.16e;\n", b);
%!           fprintf(fd, "set: real c = %.16e;\n", c);
%!           fprintf(fd, "set: real d = %.16e;\n", d);
%!           for i=1:3
%!             fprintf(fd, "set: real gamma%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).gamma(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGA%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omega(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGAP%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omegadot(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real F1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).F1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real M1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).M1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real a%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).a(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real g%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).g(i));
%!           endfor
%!           fprintf(fd, "set: real t1 = %.16e;\n", param(j, k).t1);
%!           fprintf(fd, "set: integer N = %d;\n", param(j, k).N);
%!           fprintf(fd, "set: integer M = %d;\n", options.number_of_beams);
%!           fputs(fd, "set: integer ref_id_ground = 1;\n");
%!           fputs(fd, "set: integer ref_id_tilt = 2;\n");
%!           fputs(fd, "set: integer joint_id_ground = 1;\n");
%!           fputs(fd, "set: integer force_id1;\n");
%!           fputs(fd, "set: integer torque_id1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_PHI1 = 1;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP1 = 3;\n");
%!           fputs(fd, "set: integer drive_id_PHI2 = 4;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA2 = 5;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP2 = 6;\n");
%!           fputs(fd, "set: integer drive_id_PHIx = 7;\n");
%!           fputs(fd, "set: integer drive_id_PHIy = 8;\n");
%!           fputs(fd, "set: integer drive_id_PHIz = 9;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAx = 10;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAy = 11;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAz = 12;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPx = 13;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPy = 14;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPz = 15;\n");
%!           fputs(fd, "set: integer drive_id_Xx = 16;\n");
%!           fputs(fd, "set: integer drive_id_Xy = 17;\n");
%!           fputs(fd, "set: integer drive_id_Xz = 18;\n");
%!           fputs(fd, "set: integer drive_id_XPx = 19;\n");
%!           fputs(fd, "set: integer drive_id_XPy = 20;\n");
%!           fputs(fd, "set: integer drive_id_XPz = 21;\n");
%!           fputs(fd, "set: integer drive_id_XPPx = 22;\n");
%!           fputs(fd, "set: integer drive_id_XPPy = 23;\n");
%!           fputs(fd, "set: integer drive_id_XPPz = 24;\n");
%!           fputs(fd, "set: integer drive_id_gx = 25;\n");
%!           fputs(fd, "set: integer drive_id_gy = 26;\n");
%!           fputs(fd, "set: integer drive_id_gz = 27;\n");
%!           fputs(fd, "set: integer drive_id_F1x = 28;\n");
%!           fputs(fd, "set: integer drive_id_F1y = 29;\n");
%!           fputs(fd, "set: integer drive_id_F1z = 30;\n");
%!           fputs(fd, "set: integer drive_id_M1x = 31;\n");
%!           fputs(fd, "set: integer drive_id_M1y = 32;\n");
%!           fputs(fd, "set: integer drive_id_M1z = 33;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "set: integer ref_id_modal = 3;\n");
%!               fputs(fd, "set: integer node_id_modal = 1;\n");
%!               fputs(fd, "set: integer ref_id_interface1 = 4;\n");
%!               fputs(fd, "set: integer node_id_interface1 = 2;\n");
%!               fputs(fd, "set: integer elem_id_modal = 2;\n");
%!             case 2
%!               fputs(fd, "set: integer ref_id_beam1 = 3;\n");
%!               fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!               fprintf(fd, "set: real E = %.16e;\n", mesh.material_data.E);
%!               fprintf(fd, "set: real nu = %.16e;\n", mesh.material_data.nu);
%!               fprintf(fd, "set: real rho = %.16e;\n", mesh.material_data.rho);
%!               fprintf(fd, "set: real alpha = %.16e;\n", mesh.material_data.alpha);
%!               fprintf(fd, "set: real beta = %.16e;\n", mesh.material_data.beta);
%!               fputs(fd, "set: real G = E / (2. * (1. + nu));\n");
%!               fputs(fd, "set: real A = b * c;\n");
%!               fputs(fd, "set: real As = 9. / 10. * A;\n");
%!               fputs(fd, "set: real Iy = b * c^3 / 12.;\n");
%!               fputs(fd, "set: real Iz = c * b^3 / 12.;\n");
%!               fputs(fd, "set: real Ip = Iy + Iz;\n");
%!               c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!               w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!               fprintf(fd, "set: real It = %.16e;\n", interp1(w_h, c2, max(c, b) / min(c, b)) * max(c, b) * min(c, b)^3);
%!           endswitch
%!           fputs(fd, "begin: data;\n");
%!           fputs(fd, "        problem: initial value;\n");
%!           fputs(fd, "end: data;\n");
%!           fputs(fd, "begin: initial value;\n");
%!           fputs(fd, "        initial time: 0;\n");
%!           fputs(fd, "        final time: t1;\n");
%!           fputs(fd, "        time step: t1 / N;\n");
%!           fputs(fd, "        max time step: t1 / N;\n");
%!           fputs(fd, "        min time step: t1 / N;\n");
%!           fputs(fd, "        method: ss4, 0.;\n");
%!           fputs(fd, "        tolerance: 1e-4, test, minmax, 1e-6, test,minmax;\n");
%!           fputs(fd, "        max iterations: 1000;\n");
%!           fputs(fd, "        derivatives max iterations: 50;\n");
%!           fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!           fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!           fputs(fd, "        strategy: factor, 0.8, 3, 1.25, 3, 3, 6;\n");
%!           fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!           fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 100;\n");
%!           fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!           fputs(fd, "             keep jacobian matrix,\n");
%!           fputs(fd, "             inner iterations before assembly, 6,\n");
%!           fputs(fd, "             jacobian operator, newton krylov,\n");
%!           fputs(fd, "             solver, line search based,\n");
%!           fputs(fd, "             line search method, backtrack,\n");
%!           fputs(fd, "             recovery step type, constant,\n");
%!           fputs(fd, "             recovery step, 1e-6,\n");
%!           fputs(fd, "             verbose, yes,\n");
%!           fputs(fd, "             forcing term, type 2,\n");
%!           fputs(fd, "             direction, newton,\n");
%!           fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!           fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!           fputs(fd, "             linear solver, gmres,\n");
%!           fputs(fd, "             linear solver max iterations, 300,\n");
%!           fputs(fd, "             minimum step, 1e-12,\n");
%!           fputs(fd, "             krylov subspace size, 300;\n");
%!           fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!           fputs(fd, "        threads: assembly, 1;\n");
%!           fputs(fd, "    eigenanalysis: list, 1, t1,\n");
%!           fputs(fd, "    output matrices, \n");
%!           fprintf(fd, "          parameter, %.16e,\n", 1);
%!           fputs(fd, "    output eigenvectors,\n");
%!           fputs(fd, "        output geometry,\n");
%!           fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1 / SI_unit_second^-1, 100000 / SI_unit_second^-1);
%!           switch (l)
%!           case 1
%!             fputs(fd, "    use lapack, balance, permute, suffix format, \"%02d\";\n");
%!           case 2
%!             fprintf(fd, "    use arpack,%d,%d,0.,suffix format,\"%%02d\";\n", 2 * cms_opt.modes.number, 4 * cms_opt.modes.number + 1);
%!           endswitch
%!           fputs(fd, "end: initial value;\n");
%!           fputs(fd, "begin: control data;\n");
%!           fputs(fd, "        output meter: closest next, 0., forever, t1 / 100.;\n");
%!           switch (l)
%!           case 2
%!             fputs(fd, "        rigid body kinematics: drive,\n");
%!             fputs(fd, "            angular velocity,\n");
%!             fputs(fd, "                   component,\n");
%!             for i=1:3
%!               fprintf(fd, "                postponed, drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!             endfor
%!             fputs(fd, "            acceleration,\n");
%!             fputs(fd, "                   component,\n");
%!             for i=1:3
%!               fprintf(fd, "               postponed, drive_id_XPP%s,\n", {"x", "y", "z"}{i});
%!             endfor
%!             fputs(fd, "            angular acceleration,\n");
%!             fputs(fd, "                   component");
%!             for i=1:3
%!               fprintf(fd, ",\n               postponed, drive_id_OMEGAP%s", {"x","y","z"}{i});
%!             endfor
%!             fputs(fd, ";\n");
%!           endswitch
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       use automatic differentiation;\n");
%!           fputs(fd, "       default output: none, structural nodes, drive callers;\n");
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       output precision: 16;\n");
%!           fputs(fd, "       max iterations: 0;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural nodes:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# interface 1\n");
%!               fputs(fd, "        ;\n");
%!               fputs(fd, "        joints:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# ground\n");
%!               fputs(fd, "        ;\n");
%!             case 2
%!               fputs(fd, "       structural nodes: 2 * M + 1;\n");
%!               fputs(fd, "       rigid bodies: 2 * M + 1;\n");
%!               fputs(fd, "       beams: M;\n");
%!               fputs(fd, "       joints: 1;\n");
%!           endswitch
%!           fputs(fd, "        forces: 2;\n");
%!           fputs(fd, "        gravity;\n");
%!           fputs(fd, "end: control data;\n");
%!           fputs(fd, "drive caller: drive_id_PHI1, string, \"(((pi*Time)/(2*t1)-sin((pi*Time)/t1)/2)*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA1, string, \"sin((pi*Time)/(2*t1))^2\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP1, string, \"(pi*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1)))/t1\";\n");
%!           fputs(fd, "drive caller: drive_id_PHI2, string, \"-(4*sin((pi*Time)/(2*t1))^3*t1^2)/(3*pi^2)\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA2, string, \"-(2*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1))^2*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP2, string, \"-(3*cos((pi*Time)/(2*t1))^2-1)*sin((pi*Time)/(2*t1))\";\n");
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_PHI%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_PHI1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGA1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGAP%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGAP1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_X%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XPP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_g%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, g%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_F1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, F1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_M1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, M1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, "reference: ref_id_ground,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, eye,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, null;\n");
%!           fputs(fd, "reference: ref_id_tilt,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, euler123, gammax, gammay, gammaz,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, null;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "reference: ref_id_modal,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!               fputs(fd, "reference: ref_id_interface1,\n");
%!               fputs(fd, "        reference, ref_id_modal, a + d,  0., 0.,\n");
%!               fputs(fd, "        reference, ref_id_modal, eye,\n");
%!               fputs(fd, "        reference, ref_id_modal, null,\n");
%!               fputs(fd, "        reference, ref_id_modal, null;\n");
%!             case 2
%!               fputs(fd, "reference: ref_id_beam1,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!           endswitch
%!           fputs(fd, "begin: nodes;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural: node_id_modal, modal,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, eye,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, null, accelerations, yes;\n");
%!               fputs(fd, "        structural: node_id_interface1, static,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, eye,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null;\n");
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!                 fprintf(fd, "                reference, ref_id_beam1, 0.5 * a / M * %d, 0., 0.,\n", i - 1);
%!                 fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null, accelerations, yes;\n");
%!               endfor
%!           endswitch
%!           fputs(fd, "end: nodes;\n");
%!           fputs(fd, "begin: elements;\n");
%!           switch (l)
%!           case 1
%!           fputs(fd, "        joint: joint_id_ground, total pin joint,\n");
%!           fprintf(fd, "                %s,\n", {"node_id_modal", "node_id_beam1"}{l});
%!           fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fputs(fd, "               position constraint,\n");
%!           if (~param(j, k).holonomic)
%!             fputs(fd, "                        velocity, velocity, velocity,\n");
%!             fputs(fd, "                        component,\n");
%!             for i=1:3
%!               fprintf(fd, "                      reference, drive_id_XP%s,\n", {"x", "y", "z"}{i});
%!             endfor
%!           else
%!             fputs(fd, "                        active, active, active,\n");
%!             fputs(fd, "                        component,\n");
%!             for i=1:3
%!               fprintf(fd, "                      reference, drive_id_X%s,\n", {"x", "y", "z"}{i});
%!             endfor
%!           endif
%!           fputs(fd, "               orientation constraint,\n");
%!           if (~param(j, k).holonomic)
%!             fputs(fd, "                        angular velocity, angular velocity, angular velocity,\n");
%!             fputs(fd, "                        component");
%!             for i=1:3
%!               fprintf(fd, ",\n                   reference, drive_id_OMEGA%s", {"x","y","z"}{i});
%!             endfor
%!           else
%!             fputs(fd, "                        active, active, active,\n");
%!             fputs(fd, "                        component");
%!             for i=1:3
%!               fprintf(fd, ",\n                   reference, drive_id_PHI%s", {"x","y","z"}{i});
%!             endfor
%!           endif
%!           fputs(fd, ";\n");
%!           case 2
%!             fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!           endswitch
%!           switch (l)
%!             case 1
%!               fprintf(fd, "        include: \"%s.elm\";\n", filename);
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!                 fputs(fd, "               rho * A * a / (2 * M + 1),\n");
%!                 fputs(fd, "               reference, node, null, \n");
%!                 fputs(fd, "               diag,   rho * Ip * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iy * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iz * a / (2 * M + 1),\n");
%!                 fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!               endfor
%!               for i=1:options.number_of_beams
%!                 fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "               linear elastic generic, \n");
%!                 fputs(fd, "               diag, E * A , G * As, G * As, \n");
%!                 fputs(fd, "                     G * It, E * Iy, E * Iz,\n");
%!                 fputs(fd, "               same,\n");
%!                 fputs(fd, "               same;\n");
%!               endfor
%!           endswitch
%!           fprintf(fd, "        force: force_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M"}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_F1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fprintf(fd, "        couple: torque_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M"}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_M1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fputs(fd, "        gravity: uniform, component");
%!           for i=1:3
%!             fprintf(fd, ",\n       reference, drive_id_g%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd,";\n");
%!           fputs(fd, "end: elements;\n");
%!         unwind_protect_cleanup
%!           if (fd ~= -1)
%!             fclose(fd);
%!           endif
%!           fd = -1;
%!         end_unwind_protect
%!         options_mbd.output_file = sprintf("%s_%d_%d_%d_mbd", filename, j, k, l);
%!         if (~options.verbose)
%!           options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!         endif
%!         options_mbd.mbdyn_command = "mbdyn";
%!         options_eig.positive_frequencies = false;
%!         if (options.verbose)
%!           shell(sprintf("cat %s | nl", filename_mbdyn));
%!         endif
%!         res(j, k, l).info = mbdyn_solver_run(filename_mbdyn, options_mbd);
%!         output_file_rel_frame = [options_mbd.output_file, "_rel"];
%!         mbdyn_post_abs_to_rel(1, options_mbd.output_file, output_file_rel_frame, 0);
%!         exts = {".log", ".out"};
%!         for i=1:numel(exts)
%!           [err, msg] = symlink([options_mbd.output_file, exts{i}], [output_file_rel_frame, exts{i}]);
%!           if (err ~= 0)
%!             error("failed to create symlink: %s", msg);
%!           endif
%!         endfor
%!         [res(j, k, l).t, ...
%!          res(j, k, l).trajectory, ...
%!          res(j, k, l).deformation, ...
%!          res(j, k, l).velocity, ...
%!          res(j, k, l).acceleration, ...
%!          res(j, k, l).node_id, ...
%!          res(j, k, l).force, ...
%!          res(j, k, l).force_id, ...
%!          res(j, k, l).force_node_id, ...
%!          res(j, k, l).orientation_description] = mbdyn_post_load_output_struct(output_file_rel_frame);
%!         res(j, k, l).log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!         [res(j, k, l).drive_id, ...
%!         res(j, k, l).drive_value] = mbdyn_post_load_output_drv(options_mbd.output_file, [], numel(res(j, k, l).t));
%!         res(j, k, l).modal = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, 0);
%!       endfor
%!     endfor
%!   endfor
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, :) = true;
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, :) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = options.number_of_threads;
%!   load_case = struct("omega", empty_cell, ...
%!                      "omegadot", empty_cell, ...
%!                      "loads", empty_cell, ...
%!                      "loaded_nodes", empty_cell, ...
%!                      "joints", empty_cell, ...
%!                      "g", empty_cell, ...
%!                      "tau0", empty_cell);
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   for i=1:numel(load_case)
%!     load_case(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!     load_case(i).loads = zeros(1, 6);
%!     load_case(i).omega = zeros(3, 1);
%!     load_case(i).omegadot = zeros(3, 1);
%!     load_case(i).g = zeros(3, 1);
%!     load_case(i).tau0.iso20 = zeros(rows(mesh.elements.iso20), columns(mesh.elements.iso20), 6);
%!   endfor
%!   sol_eig = struct("def", empty_cell, "lambda", empty_cell, "f", empty_cell, "D", empty_cell);
%!   sol_eig_red = struct("lambda_red", empty_cell, "Ured", empty_cell);
%!   for j=idx_j
%!     for k=idx_k
%!       R = euler123_to_rotation_matrix(param(j, k).gamma);
%!       load_case(j, k).omega = R.' * param(j, k).omega;
%!       load_case(j, k).omegadot = R.' * param(j, k).omegadot;
%!       load_case(j, k).loads = [(R.' * param(j, k).F1).', (R.' * param(j, k).M1).'];
%!       load_case(j, k).g = R.' * (param(j, k).g - param(j, k).a);
%!       [mat_ass.M, ...
%!        mat_ass.D, ...
%!        mat_ass.K, ...
%!        mat_ass.KOMEGA, ...
%!        mat_ass.KOMEGA_DOT, ...
%!        mat_ass.DOMEGA, ...
%!        mat_ass.R, ...
%!        mat_ass.mat_info, ...
%!        mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                            dof_map, ...
%!                                            [FEM_MAT_MASS, ...
%!                                             FEM_MAT_DAMPING, ...
%!                                             FEM_MAT_STIFFNESS, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                             FEM_MAT_DAMPING_OMEGA, ...
%!                                             FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case(j, k));
%!       cms_opt.symmetric = false;
%!       sol_statjk = fem_sol_static(mesh, dof_map, mat_ass, cms_opt);
%!       sol_statjk.stress = fem_ass_matrix(mesh, ...
%!                                          dof_map, ...
%!                                          [FEM_VEC_STRESS_CAUCH], ...
%!                                          load_case(j, k), ...
%!                                          sol_statjk);
%!       sol_stat(j, k) = sol_statjk;
%!       load_case(j, k).tau0 = sol_stat(j, k).stress.tau;
%!       mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS_TAU0], ...
%!                                      load_case(j, k));
%!       mat_ass.K += mat_ass.KOMEGA + mat_ass.KOMEGA_DOT + mat_ass.KTAU0;
%!       mat_ass.D += mat_ass.DOMEGA;
%!       sol_eig(j, k) = fem_sol_modal_damped(mesh, ...
%!                                            dof_map, ...
%!                                            mat_ass, ...
%!                                            cms_opt.modes.number, ...
%!                                            cms_opt);
%!       Mred = mat_ass_cms.Mred;
%!       Dred = mat_ass_cms.Dred;
%!       Kred = mat_ass_cms.Kred;
%!       Dred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.DOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA_DOT(dof_map.idx_node, dof_map.idx_node), "Full");
%!       omegaq = [load_case(j, k).omega.^2;
%!                 load_case(j, k).omega(1) * load_case(j, k).omega(2);
%!                 load_case(j, k).omega(2) * load_case(j, k).omega(3);
%!                 load_case(j, k).omega(3) * load_case(j, k).omega(1)];
%!       idx = int32(0);
%!       for i=1:numel(omegaq)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * omegaq(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).omegadot)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).omegadot(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).g)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred -= mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).g(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).loads)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).loads(i);
%!       endfor
%!       [sol_eig_red(j, k).Ured, sol_eig_red(j, k).lambda_red] = fem_sol_eigsd(Kred, Dred, Mred, cms_opt.modes.number, cms_opt);
%!     endfor
%!   endfor
%! tol_abs = [0, 0] / SI_unit_second^-1;
%! tol_rel = [0.3e-2, 3e-2];
%! tol_disp_rel = 3e-2;
%! err_u_modal = err_v_modal = zeros(size(param));
%! printf("deformation/velocity:\n");
%! colors = rainbow(3);
%! width = 1:size(res, 3);
%! linestyle = {"-", "--"};
%! for i=idx_j
%!   for j=idx_k
%!      u_modal = res(i, j, 1).trajectory{end} - res(i, j, 1).trajectory{end}(1, :);
%!      u_beam = res(i, j, 2).trajectory{end} - res(i, j, 2).trajectory{end}(1, :);
%!      v_modal = res(i, j, 1).velocity{end};
%!      v_beam = res(i, j, 2).velocity{end};
%!      if (options.plot)
%!      figure("visible", "off");
%!      hold on;
%!      for k=1:size(res, 3)
%!        for l=1:3
%!          hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l) - res(i, j, k).trajectory{end}(1, l)) * SI_unit_meter);
%!          set(hnd, "color", colors(l, :));
%!          set(hnd, "linewidth", width(k));
%!          set(hnd, "linestyle", linestyle{k});
%!        endfor
%!      endfor
%!      xlabel("t [s]");
%!      ylabel("u [m]");
%!      grid on;
%!      grid minor on;
%!      title(sprintf("linear displacement %d:%d", i, j));
%!      figure("visible", "off");
%!      hold on;
%!      for k=1:size(res, 3)
%!        for l=1:3
%!          hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l + 3) - res(i, j, k).trajectory{end}(1, l + 3)) * 180 / pi);
%!          set(hnd, "color", colors(l, :));
%!          set(hnd, "linewidth", width(k));
%!          set(hnd, "linestyle", linestyle{k});
%!        endfor
%!      endfor
%!      xlabel("t [s]");
%!      ylabel("Phi [deg]");
%!      grid on;
%!      grid minor on;
%!      title(sprintf("angular displacement %d:%d", i, j));
%!      endif
%!      err_u_modal(i, j) = max(max(abs(u_modal - u_beam))) / max(1, max(max(abs(u_beam))));
%!      err_v_modal(i, j) = max(max(abs(v_modal - v_beam))) / max(1, max(max(abs(v_beam))));
%!      printf("%d:%d %.1f%%/%.1f%%\n", i, j, 100 * err_u_modal(i, j), 100 * err_v_modal(i, j));
%!   endfor
%! endfor
%! printf("natural frequencies:\n");
%! MACR = cell(size(param));
%! result_data = struct("f_mbd", cell(size(param)), "f_fem", cell(size(param)));
%! for i=idx_j
%!   for j=idx_k
%!     f_fem = sort(sol_eig(i, j).f(:));
%!     f_fem = f_fem(f_fem > 0);
%!     f_mbd = zeros(rows(f_fem), size(res, 3));
%!     PhiR = zeros(6, rows(f_fem), size(res, 3));
%!     for k=1:size(res, 3)
%!       [f_mbd_k, idx_mbd_k] = sort(res(i, j, k).modal.f(:));
%!       D_mbd_k = res(i, j, k).modal.D(idx_mbd_k);
%!       idx_mbd_k = idx_mbd_k(f_mbd_k > 0);
%!       f_mbd_k = f_mbd_k(f_mbd_k > 0);
%!       idx_mbd_k = idx_mbd_k(1:rows(f_fem));
%!       f_mbd(:, k) = f_mbd_k(1:rows(f_fem));
%!       PhiR(:, :, k) = res(i, j, k).modal.VR(res(i, j, k).modal.idx(end) + (1:6), idx_mbd_k);
%!     endfor
%!     result_data(i, j).f_fem = f_fem;
%!     result_data(i, j).f_mbd = f_mbd;
%!     MACR{i, j} = MACL{i, j} = zeros(rows(f_fem), rows(f_fem));
%!     for k=1:rows(f_fem)
%!       for l=1:rows(f_fem)
%!         MACR{i, j}(k, l) = (PhiR(:, k, 1)' * PhiR(:, k, 2)) * conj(PhiR(:, k, 1)' * PhiR(:, k, 2)) / ((PhiR(:, k, 1)' * PhiR(:, k, 1)) * (PhiR(:, k, 2)' * PhiR(:, k, 2)));
%!       endfor
%!     endfor
%!     printf("%d:%d\n", i, j);
%!     for k=1:rows(f_fem)
%!       printf("%10.2f", f_fem(k) * SI_unit_second^-1);
%!       for l=1:columns(f_mbd)
%!         printf("\t%10.2f", f_mbd(k, l) * SI_unit_second^-1);
%!       endfor
%!       for l=1:columns(f_mbd)
%!         printf("\t%.1f%%", 100 * (f_mbd(k, l) / f_fem(k) - 1));
%!       endfor
%!       printf("\t%.3f", MACR{i, j}(k, k));
%!       fputs(stdout, "\n");
%!     endfor
%!    fputs(stdout, "\n\n");
%!   endfor
%! endfor
%! for i=idx_j
%!   for j=idx_k
%!     for k=1:rows(result_data(i, j).f_fem)
%!       for l=1:columns(result_data(i, j).f_mbd)
%!         assert_simple(result_data(i, j).f_mbd(k, l), result_data(i, j).f_fem(k), tol_abs(l) + tol_rel(l) * abs(result_data(i, j).f_fem(k)));
%!       endfor
%!     endfor
%!   endfor
%! endfor
%! assert_simple(all(all(err_u_modal < tol_disp_rel)));
%! for j=idx_j
%!   for k=idx_k
%!       tol = 2e-2;
%!       [lambda_s] = sortrows([imag(sol_eig(j, k).lambda(:)), real(sol_eig(j, k).lambda(:))],[1,2]);
%!       [lambda_red_s] = sortrows([imag(sol_eig_red(j, k).lambda_red(:)), real(sol_eig_red(j, k).lambda_red(:))],[1,2]);
%!       K = min(20, rows(lambda_s));
%!       lambda_s = 1j * lambda_s(:,1) + lambda_s(:, 2);
%!       lambda_red_s = 1j * lambda_red_s(:, 1) + lambda_red_s(:, 2);
%!       assert_simple(lambda_red_s(1:K), lambda_s(1:K), tol * norm(lambda_s(1:K)));
%!   endfor
%! endfor
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST5
%! ## Oskar Wallrapp, Richard Schwertassek, 1998
%! ## Dynamik flexibler Mehrkoerpersysteme
%! ## chapter 5, table 5.7, page 242
%! ## Natural frequencies of a rotating beam structure
%! pkg load mboct-fem-pkg;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   options.verbose = false;
%!   options.number_of_modes = int32(20);
%!   options.interactive = false;
%!   c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!   w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!   param.num_fem_elem_per_sec = 5;
%!   param.N1 = int32(50);
%!   param.N2 = int32(20);
%!   param.E1 = 7e10 / SI_unit_pascal;
%!   param.E2 = 21e10 / SI_unit_pascal;
%!   param.nu1 = 0.3;
%!   param.nu2 = 0.3;
%!   param.rho1 = 3000 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.rho2 = 7895 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.l1 = 2 / SI_unit_meter;
%!   param.l2 = 0.4 / SI_unit_meter;
%!   param.h1 = 0.009 / SI_unit_meter;
%!   param.h2 = param.h1;
%!   param.w1 = 0.009 / SI_unit_meter;
%!   param.w2 = 0.0095 / SI_unit_meter;
%!   param.OMEGAx = 0 / SI_unit_second^-1;
%!   param.OMEGAy = 0 / SI_unit_second^-1;
%!   param.OMEGAz = 6 / SI_unit_second^-1;
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w1 = %.16e;\n", param.w1);
%!     fprintf(fd, "w2 = %.16e;\n", param.w2);
%!     fprintf(fd, "l1 = %.16e;\n", param.l1);
%!     fprintf(fd, "l2 = %.16e;\n", param.l2);
%!     fprintf(fd, "h1 = %.16e;\n", param.h1);
%!     fprintf(fd, "h = %.16e;\n", min([param.w1, param.w2]) / param.num_fem_elem_per_sec);
%!     fputs(fd, "Point(1) = {0, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(2) = {l1 - 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(3) = {l1 - 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(4) = {l1 + 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(5) = {l1 + 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(6) = {l1 + 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(7) = {l1 - 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(8) = {0, -0.5 * w1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line(9) = {2, 7};\n");
%!     fputs(fd, "Line(10) = {2, 5};\n");
%!     fputs(fd, "Line Loop(1) = {1, 9, 7, 8};\n");
%!     fputs(fd, "Line Loop(2) = {10, 5, 6, 9};\n");
%!     fputs(fd, "Line Loop(3) = {2, 3, 4, 10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Plane Surface(2) = {2};\n");
%!     fputs(fd, "Plane Surface(3) = {3};\n");
%!     fputs(fd, "hw1 = w1 / Round(w1 / h);\n");
%!     fputs(fd, "hw2 = w2 / Round(w2 / h);\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round((l1 - 0.5 * w2) / hw1)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round((l2 - 0.5 * w1) / hw2)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(w2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round((l2 - 0.5 * w1) / hw2)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(5) = Max(1, Round(w1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(6) = Max(1, Round(w2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(7) = Max(1, Round((l1 - 0.5 * w2) / hw1)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(8) = Max(1, Round(w1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(9) = Max(1, Round(w1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(10) = Max(1, Round(w2 / h)) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,h1}{Surface{1,2,3}; Layers{Max(1, Round(h1 / h))}; Recombine;};\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Recombine Surface{2,12};\n");
%!     fputs(fd, "Recombine Surface{3,16};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Transfinite Surface(2) = {};\n");
%!     fputs(fd, "Transfinite Surface(3) = {};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 1) = {7};\n");
%!     fputs(fd, "Physical Volume(\"beam1\", 1) = {1, 2};\n");
%!     fputs(fd, "Physical Volume(\"beam2\", 2) = {3};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "1", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   mesh.material_data = struct("E", cell(1, 2), "nu", cell(1, 2), "rho", cell(1, 2));
%!   mesh.material_data(1).E = param.E1;
%!   mesh.material_data(1).nu = param.nu1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(2).E = param.E2;
%!   mesh.material_data(2).nu = param.nu2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!   grp_idx_beam1 = find([[mesh.groups.iso8].id] == 1);
%!   grp_idx_beam2 = find([[mesh.groups.iso8].id] == 2);
%!   grp_idx_clamp = find([[mesh.groups.iso4].id] == 1);
%!   mesh.materials.iso8(mesh.groups.iso8(grp_idx_beam1).elements) = 1;
%!   mesh.materials.iso8(mesh.groups.iso8(grp_idx_beam2).elements) = 2;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.iso4(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.omega = [param.OMEGAx;
%!                      param.OMEGAy;
%!                      param.OMEGAz];
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = mbdyn_solver_num_threads_default();
%!   [mat_ass.M, ...
%!    mat_ass.Mdiag, ...
%!    mat_ass.K, ...
%!    mat_ass.KOMEGA, ...
%!    mat_ass.DOMEGA, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_MASS_LUMPED, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.tau0 = sol_stat.stress.tau;
%!   mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case);
%!   opt_solver.pre_scaling = true;
%!   opt_solver.refine_max_iter = int32(10);
%!   opt_solver.solver = "pardiso";
%!   opt_solver.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_solver.symmetric = true;
%!   sol_eig(1) = fem_sol_modal(mesh, ...
%!                            dof_map, ...
%!                            mat_ass, ...
%!                            options.number_of_modes, ...
%!                            0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   sol_eig_diag = fem_sol_modal(mesh, ...
%!                                dof_map, ...
%!                                setfield(mat_ass, "M", mat_ass.Mdiag), ...
%!                                options.number_of_modes, ...
%!                                0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0;
%!   mat_ass.D = mat_ass.DOMEGA;
%!   sol_eig(2) = fem_sol_modal_damped(mesh, ...
%!                                     dof_map, ...
%!                                     mat_ass, ...
%!                                     options.number_of_modes, ...
%!                                     opt_solver);
%!   unwind_protect
%!     filename_mbdyn = [filename, ".mbdyn"];
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     mbdyn_pre_write_param_file(fd, param);
%!     fprintf(fd, "set: real t1 = %g;\n", 100 / SI_unit_second);
%!     fputs(fd, "set: integer N = 20;\n");
%!     fputs(fd, "set: real G1 = E1 / (2. * (1. + nu1));\n");
%!     fputs(fd, "set: real G2 = E2 / (2. * (1. + nu2));\n");
%!     fputs(fd, "set: real A1 = w1 * h1;\n");
%!     fputs(fd, "set: real A2 = w2 * h2;\n");
%!     fputs(fd, "set: real As1 = 9. / 10. * A1;\n");
%!     fputs(fd, "set: real As2 = 9. / 10. * A2;\n");
%!     fputs(fd, "set: real Iy1 = w1 * h1^3 / 12.;\n");
%!     fputs(fd, "set: real Iy2 = w2 * h2^3 / 12.;\n");
%!     fputs(fd, "set: real Iz1 = h1 * w1^3 / 12.;\n");
%!     fputs(fd, "set: real Iz2 = h2 * w2^3 / 12.;\n");
%!     fputs(fd, "set: real Ip1 = Iy1 + Iz1;\n");
%!     fputs(fd, "set: real Ip2 = Iy2 + Iz2;\n");
%!     fprintf(fd, "set: real c21 = %.16e;\n", interp1(w_h, c2, max(param.w1, param.h1) / min(param.w1, param.h1)));
%!     fprintf(fd, "set: real c22 = %.16e;\n", interp1(w_h, c2, max(param.w2, param.h2) / min(param.w2, param.h2)));
%!     fputs(fd, "set: real It1 = c21 * h1 * w1^3;\n");
%!     fputs(fd, "set: real It2 = c22 * h2 * w2^3;\n");
%!     fputs(fd, "set: integer ref_id_ground = 1;\n");
%!     fputs(fd, "set: integer ref_id_beam1 = 2;\n");
%!     fputs(fd, "set: integer ref_id_beam2 = 3;\n");
%!     fputs(fd, "set: integer joint_id_ground = 1;\n");
%!     fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer node_id_beam2 = 2 * N1 + node_id_beam1;\n");
%!     fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer body_id_beam2 = 2 * N1 + body_id_beam1 + 1;\n");
%!     fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer beam_id_beam2 = beam_id_beam1 + N1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 2 * t1;\n");
%!     fputs(fd, "        time step: t1 / N;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-5, test,norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             line search method, backtrack,\n");
%!     fputs(fd, "             recovery step type, constant,\n");
%!     fputs(fd, "             recovery step, 1e-6,\n");
%!     fputs(fd, "             verbose, yes,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "           eigenanalysis: list, 2, 0, 2 * t1,\n");
%!     fputs(fd, "           # output matrices, \n");
%!     fprintf(fd, "         parameter, %.16e,\n", 1e-4 / SI_unit_second);
%!     fputs(fd, "           output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1e-3 / SI_unit_second^-1, 200. / SI_unit_second^-1);
%!     fprintf(fd, "           use arpack, %d, %d, 0, suffix format, \"%%02d\";\n", 3 * options.number_of_modes, 10 * options.number_of_modes + 1);
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, t1 / 10.;\n");
%!     fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       structural nodes: 2 * N1 + 1 + 2 * N2;\n");
%!     fputs(fd, "       rigid bodies: 2 * N1 + 1 + 2 * N2 + 1;\n");
%!     fputs(fd, "       beams: N1 + N2;\n");
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_beam1,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, null;\n");
%!     fputs(fd, "reference: ref_id_beam2,\n");
%!     fputs(fd, "        reference, ref_id_beam1, l1,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, 1, 0., 1., 0., 3, 0., 0., 1.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam1, 0.5 * l1 / N1 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null;\n");
%!     endfor
%!     for i=2:(2 * param.N2 + 1)
%!       fprintf(fd, "        structural: node_id_beam2 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam2, 0.5 * l2 / N2 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null;\n");
%!     endfor
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!       fputs(fd, "               rho1 * A1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho1 * Ip1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iy1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iz1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!     endfor
%!     for i=1:(2 * param.N2 + 1)
%!       fprintf(fd, "     body: body_id_beam2 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d,\n", i - 1);
%!       fputs(fd, "               rho2 * A2 * l2 / (2 * N2 + 1), \n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho2 * Ip2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iy2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iz2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam2, eye;\n");
%!     endfor
%!     for i=1:param.N1
%!       fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E1 * A1 , G1 * As1, G1 * As1, \n");
%!       fputs(fd, "                     G1 * It1, E1 * Iy1, E1 * Iz1,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     for i=1:param.N2
%!       fprintf(fd, "        beam3: beam_id_beam2 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E2 * A2 , G2 * As2, G2 * As2, \n");
%!       fputs(fd, "                     G2 * It2, E2 * Iy2, E2 * Iz2,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = [filename, "_mbd"];
%!   if (~options.verbose)
%!     options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
%!   options_mbd.mbdyn_command = "mbdyn";
%!   options_eig.positive_frequencies = false;
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   for i=1:2
%!     res.modal(i) = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, i - 1);
%!   endfor
%!   if (options.interactive)
%!   for j=1:numel(res.modal)
%!     for i=1:numel(res.modal(j).f)
%!       opt_modal.mode_index = i;
%!       opt_modal.scale = 100;
%!       mode_file = [options_mbd.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(options_mbd.output_file, [mode_file, ".mov"], opt_modal, res.modal(j));
%!       [err, msg] = symlink([options_mbd.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = options.interactive;
%!       opt_post.f_runEasyAnim = options.interactive;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!   endfor
%!   endif
%!   idx_mode = [1, 4, 6];
%!   ## table 5.5 (11 elements), table 5.7 (11 elements)
%!   fref =  [0.9614, 1.044;
%!            7.037,  7.478;
%!           16.67,  17.09];
%!   tol = 2e-2;
%!   for i=1:2
%!     fmbd = sort(res.modal(i).f(:)) * SI_unit_second^-1;
%!     fmbd = fmbd(fmbd > 0);
%!     ffem = sort(sol_eig(i).f(:)) * SI_unit_second^-1;
%!     ffem = ffem(ffem >= 0);
%!     Nfem = min(numel(fmbd),numel(ffem));
%!     assert_simple(fmbd(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(ffem(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(fmbd(1:Nfem), ffem(1:Nfem), tol * max(fmbd(1:Nfem)));
%!   endfor
%!   assert_simple(max(abs(sol_eig_diag.f/sol_eig(1).f - 1)) < 0.03);
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST6
%! ## Oskar Wallrapp, Richard Schwertassek, 1998
%! ## Dynamik flexibler Mehrkoerpersysteme
%! ## chapter 5, table 5.7, page 242
%! ## Natural frequencies of a rotating beam structure
%! pkg load mboct-fem-pkg;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   options.verbose = false;
%!   options.number_of_modes = int32(20);
%!   options.interactive = false;
%!   c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!   w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!   param.num_fem_elem_per_sec = 0.5;
%!   param.N1 = int32(50);
%!   param.N2 = int32(20);
%!   param.E1 = 7e10 / SI_unit_pascal;
%!   param.E2 = 21e10 / SI_unit_pascal;
%!   param.nu1 = 0.3;
%!   param.nu2 = 0.3;
%!   param.rho1 = 3000 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.rho2 = 7895 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.l1 = 2 / SI_unit_meter;
%!   param.l2 = 0.4 / SI_unit_meter;
%!   param.h1 = 0.009 / SI_unit_meter;
%!   param.h2 = param.h1;
%!   param.w1 = 0.009 / SI_unit_meter;
%!   param.w2 = 0.0095 / SI_unit_meter;
%!   param.OMEGAx = 0 / SI_unit_second^-1;
%!   param.OMEGAy = 0 / SI_unit_second^-1;
%!   param.OMEGAz = 6 / SI_unit_second^-1;
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w1 = %.16e;\n", param.w1);
%!     fprintf(fd, "w2 = %.16e;\n", param.w2);
%!     fprintf(fd, "l1 = %.16e;\n", param.l1);
%!     fprintf(fd, "l2 = %.16e;\n", param.l2);
%!     fprintf(fd, "h1 = %.16e;\n", param.h1);
%!     fprintf(fd, "h = %.16e;\n", min([param.w1, param.w2, param.h1]) / param.num_fem_elem_per_sec);
%!     fputs(fd, "Point(1) = {0, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(2) = {l1 - 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(3) = {l1 - 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(4) = {l1 + 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(5) = {l1 + 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(6) = {l1 + 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(7) = {l1 - 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(8) = {0, -0.5 * w1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line(9) = {2, 7};\n");
%!     fputs(fd, "Line(10) = {2, 5};\n");
%!     fputs(fd, "Line Loop(1) = {1, 9, 7, 8};\n");
%!     fputs(fd, "Line Loop(2) = {10, 5, 6, 9};\n");
%!     fputs(fd, "Line Loop(3) = {2, 3, 4, 10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Plane Surface(2) = {2};\n");
%!     fputs(fd, "Plane Surface(3) = {3};\n");
%!     fputs(fd, "hw1 = w1 / Round(w1 / h);\n");
%!     fputs(fd, "hw2 = w2 / Round(w2 / h);\n");
%!     fputs(fd, "Transfinite Curve(1) = Round((l1 - 0.5 * w2) / hw1) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Round((l2 - 0.5 * w1) / hw2) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Round((l2 - 0.5 * w1) / hw2) + 1;\n");
%!     fputs(fd, "Transfinite Curve(5) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(6) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(7) = Round((l1 - 0.5 * w2) / hw1) + 1;\n");
%!     fputs(fd, "Transfinite Curve(8) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(9) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(10) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,h1}{Surface{1,2,3}; Layers{Round(h1 / h)}; Recombine;};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Transfinite Surface(2) = {};\n");
%!     fputs(fd, "Transfinite Surface(3) = {};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 1) = {7};\n");
%!     fputs(fd, "Physical Volume(\"beam1\", 1) = {1, 2};\n");
%!     fputs(fd, "Physical Volume(\"beam2\", 2) = {3};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   mesh.material_data = struct("E", cell(1, 2), "nu", cell(1, 2), "rho", cell(1, 2));
%!   mesh.material_data(1).E = param.E1;
%!   mesh.material_data(1).nu = param.nu1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(2).E = param.E2;
%!   mesh.material_data(2).nu = param.nu2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!   grp_idx_beam1 = find([[mesh.groups.penta15].id] == 1);
%!   grp_idx_beam2 = find([[mesh.groups.penta15].id] == 2);
%!   grp_idx_clamp = find([[mesh.groups.quad8].id] == 1);
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_beam1).elements) = 1;
%!   mesh.materials.penta15(mesh.groups.penta15(grp_idx_beam2).elements) = 2;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.omega = [param.OMEGAx;
%!                      param.OMEGAy;
%!                      param.OMEGAz];
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = mbdyn_solver_num_threads_default();
%!   [mat_ass.M, ...
%!    mat_ass.Mdiag, ...
%!    mat_ass.K, ...
%!    mat_ass.KOMEGA, ...
%!    mat_ass.DOMEGA, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_MASS_LUMPED, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.tau0 = sol_stat.stress.tau;
%!   mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case);
%!   opt_solver.pre_scaling = true;
%!   opt_solver.refine_max_iter = int32(10);
%!   opt_solver.solver = "pardiso";
%!   opt_solver.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_solver.symmetric = true;
%!   sol_eig(1) = fem_sol_modal(mesh, ...
%!                            dof_map, ...
%!                            mat_ass, ...
%!                            options.number_of_modes, ...
%!                            0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   sol_eig_diag = fem_sol_modal(mesh, ...
%!                                dof_map, ...
%!                                setfield(mat_ass, "M", mat_ass.Mdiag), ...
%!                                options.number_of_modes, ...
%!                                0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0;
%!   mat_ass.D = mat_ass.DOMEGA;
%!   sol_eig(2) = fem_sol_modal_damped(mesh, ...
%!                                     dof_map, ...
%!                                     mat_ass, ...
%!                                     options.number_of_modes, ...
%!                                     opt_solver);
%!   unwind_protect
%!     filename_mbdyn = [filename, ".mbdyn"];
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     mbdyn_pre_write_param_file(fd, param);
%!     fprintf(fd, "set: real t1 = %g;\n", 100 / SI_unit_second);
%!     fputs(fd, "set: integer N = 20;\n");
%!     fputs(fd, "set: real G1 = E1 / (2. * (1. + nu1));\n");
%!     fputs(fd, "set: real G2 = E2 / (2. * (1. + nu2));\n");
%!     fputs(fd, "set: real A1 = w1 * h1;\n");
%!     fputs(fd, "set: real A2 = w2 * h2;\n");
%!     fputs(fd, "set: real As1 = 9. / 10. * A1;\n");
%!     fputs(fd, "set: real As2 = 9. / 10. * A2;\n");
%!     fputs(fd, "set: real Iy1 = w1 * h1^3 / 12.;\n");
%!     fputs(fd, "set: real Iy2 = w2 * h2^3 / 12.;\n");
%!     fputs(fd, "set: real Iz1 = h1 * w1^3 / 12.;\n");
%!     fputs(fd, "set: real Iz2 = h2 * w2^3 / 12.;\n");
%!     fputs(fd, "set: real Ip1 = Iy1 + Iz1;\n");
%!     fputs(fd, "set: real Ip2 = Iy2 + Iz2;\n");
%!     fprintf(fd, "set: real c21 = %.16e;\n", interp1(w_h, c2, max(param.w1, param.h1) / min(param.w1, param.h1)));
%!     fprintf(fd, "set: real c22 = %.16e;\n", interp1(w_h, c2, max(param.w2, param.h2) / min(param.w2, param.h2)));
%!     fputs(fd, "set: real It1 = c21 * h1 * w1^3;\n");
%!     fputs(fd, "set: real It2 = c22 * h2 * w2^3;\n");
%!     fputs(fd, "set: integer ref_id_ground = 1;\n");
%!     fputs(fd, "set: integer ref_id_beam1 = 2;\n");
%!     fputs(fd, "set: integer ref_id_beam2 = 3;\n");
%!     fputs(fd, "set: integer joint_id_ground = 1;\n");
%!     fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer node_id_beam2 = 2 * N1 + node_id_beam1;\n");
%!     fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer body_id_beam2 = 2 * N1 + body_id_beam1 + 1;\n");
%!     fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer beam_id_beam2 = beam_id_beam1 + N1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 2 * t1;\n");
%!     fputs(fd, "        time step: t1 / N;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-5, test,norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             line search method, backtrack,\n");
%!     fputs(fd, "             recovery step type, constant,\n");
%!     fputs(fd, "             recovery step, 1e-6,\n");
%!     fputs(fd, "             verbose, yes,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "           eigenanalysis: list, 2, 0, 2 * t1,\n");
%!     fputs(fd, "           # output matrices, \n");
%!     fprintf(fd, "         parameter, %.16e,\n", 1e-4 / SI_unit_second);
%!     fputs(fd, "           output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1e-3 / SI_unit_second^-1, 200. / SI_unit_second^-1);
%!     fprintf(fd, "           use arpack, %d, %d, 0, suffix format, \"%%02d\";\n", 3 * options.number_of_modes, 10 * options.number_of_modes + 1);
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, t1 / 10.;\n");
%!     fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       structural nodes: 2 * N1 + 1 + 2 * N2;\n");
%!     fputs(fd, "       rigid bodies: 2 * N1 + 1 + 2 * N2 + 1;\n");
%!     fputs(fd, "       beams: N1 + N2;\n");
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_beam1,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, null;\n");
%!     fputs(fd, "reference: ref_id_beam2,\n");
%!     fputs(fd, "        reference, ref_id_beam1, l1,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, 1, 0., 1., 0., 3, 0., 0., 1.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam1, 0.5 * l1 / N1 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null;\n");
%!     endfor
%!     for i=2:(2 * param.N2 + 1)
%!       fprintf(fd, "        structural: node_id_beam2 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam2, 0.5 * l2 / N2 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null;\n");
%!     endfor
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!       fputs(fd, "               rho1 * A1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho1 * Ip1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iy1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iz1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!     endfor
%!     for i=1:(2 * param.N2 + 1)
%!       fprintf(fd, "     body: body_id_beam2 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d,\n", i - 1);
%!       fputs(fd, "               rho2 * A2 * l2 / (2 * N2 + 1), \n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho2 * Ip2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iy2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iz2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam2, eye;\n");
%!     endfor
%!     for i=1:param.N1
%!       fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E1 * A1 , G1 * As1, G1 * As1, \n");
%!       fputs(fd, "                     G1 * It1, E1 * Iy1, E1 * Iz1,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     for i=1:param.N2
%!       fprintf(fd, "        beam3: beam_id_beam2 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E2 * A2 , G2 * As2, G2 * As2, \n");
%!       fputs(fd, "                     G2 * It2, E2 * Iy2, E2 * Iz2,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = [filename, "_mbd"];
%!   if (~options.verbose)
%!     options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
%!   options_mbd.mbdyn_command = "mbdyn";
%!   options_eig.positive_frequencies = false;
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   for i=1:2
%!     res.modal(i) = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, i - 1);
%!   endfor
%!   if (options.interactive)
%!   for j=1:numel(res.modal)
%!     for i=1:numel(res.modal(j).f)
%!       opt_modal.mode_index = i;
%!       opt_modal.scale = 100;
%!       mode_file = [options_mbd.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(options_mbd.output_file, [mode_file, ".mov"], opt_modal, res.modal(j));
%!       [err, msg] = symlink([options_mbd.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = options.interactive;
%!       opt_post.f_runEasyAnim = options.interactive;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!   endfor
%!   endif
%!   idx_mode = [1, 4, 6];
%!   ## table 5.5 (11 elements), table 5.7 (11 elements)
%!   fref =  [0.9614, 1.044;
%!            7.037,  7.478;
%!           16.67,  17.09];
%!   tol = 2e-2;
%!   for i=1:2
%!     fmbd = sort(res.modal(i).f(:)) * SI_unit_second^-1;
%!     fmbd = fmbd(fmbd > 0);
%!     ffem = sort(sol_eig(i).f(:)) * SI_unit_second^-1;
%!     ffem = ffem(ffem >= 0);
%!     Nfem = min(numel(fmbd),numel(ffem));
%!     assert_simple(fmbd(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(ffem(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(fmbd(1:Nfem), ffem(1:Nfem), tol * max(fmbd(1:Nfem)));
%!   endfor
%!   assert_simple(max(abs(sol_eig_diag.f/sol_eig(1).f - 1)) < 0.025);
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST7
%! ## Oskar Wallrapp, Richard Schwertassek, 1998
%! ## Dynamik flexibler Mehrkoerpersysteme
%! ## chapter 5, table 5.7, page 242
%! ## Natural frequencies of a rotating beam structure
%! pkg load mboct-fem-pkg;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   options.verbose = false;
%!   options.number_of_modes = int32(20);
%!   options.interactive = false;
%!   c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!   w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!   param.num_fem_elem_per_sec = 1;
%!   param.N1 = int32(50);
%!   param.N2 = int32(20);
%!   param.E1 = 7e10 / SI_unit_pascal;
%!   param.E2 = 21e10 / SI_unit_pascal;
%!   param.nu1 = 0.3;
%!   param.nu2 = 0.3;
%!   param.rho1 = 3000 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.rho2 = 7895 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.l1 = 2 / SI_unit_meter;
%!   param.l2 = 0.4 / SI_unit_meter;
%!   param.h1 = 0.009 / SI_unit_meter;
%!   param.h2 = param.h1;
%!   param.w1 = 0.009 / SI_unit_meter;
%!   param.w2 = 0.0095 / SI_unit_meter;
%!   param.OMEGAx = 0 / SI_unit_second^-1;
%!   param.OMEGAy = 0 / SI_unit_second^-1;
%!   param.OMEGAz = 6 / SI_unit_second^-1;
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w1 = %.16e;\n", param.w1);
%!     fprintf(fd, "w2 = %.16e;\n", param.w2);
%!     fprintf(fd, "l1 = %.16e;\n", param.l1);
%!     fprintf(fd, "l2 = %.16e;\n", param.l2);
%!     fprintf(fd, "h1 = %.16e;\n", param.h1);
%!     fprintf(fd, "h = %.16e;\n", min([param.w1, param.w2]) / param.num_fem_elem_per_sec);
%!     fputs(fd, "Point(1) = {0, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(2) = {l1 - 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(3) = {l1 - 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(4) = {l1 + 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(5) = {l1 + 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(6) = {l1 + 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(7) = {l1 - 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(8) = {0, -0.5 * w1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line(9) = {2, 7};\n");
%!     fputs(fd, "Line(10) = {2, 5};\n");
%!     fputs(fd, "Line Loop(1) = {1, 9, 7, 8};\n");
%!     fputs(fd, "Line Loop(2) = {10, 5, 6, 9};\n");
%!     fputs(fd, "Line Loop(3) = {2, 3, 4, 10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Plane Surface(2) = {2};\n");
%!     fputs(fd, "Plane Surface(3) = {3};\n");
%!     fputs(fd, "hw1 = w1 / Round(w1 / h);\n");
%!     fputs(fd, "hw2 = w2 / Round(w2 / h);\n");
%!     fputs(fd, "Transfinite Curve(1) = Round((l1 - 0.5 * w2) / hw1) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Round((l2 - 0.5 * w1) / hw2) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Round((l2 - 0.5 * w1) / hw2) + 1;\n");
%!     fputs(fd, "Transfinite Curve(5) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(6) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(7) = Round((l1 - 0.5 * w2) / hw1) + 1;\n");
%!     fputs(fd, "Transfinite Curve(8) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(9) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(10) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,h1}{Surface{1,2,3}; Layers{Round(h1 / h)}; Recombine;};\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Recombine Surface{2,12};\n");
%!     fputs(fd, "Recombine Surface{3,16};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Transfinite Surface(2) = {};\n");
%!     fputs(fd, "Transfinite Surface(3) = {};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 1) = {7};\n");
%!     fputs(fd, "Physical Volume(\"beam1\", 1) = {1, 2};\n");
%!     fputs(fd, "Physical Volume(\"beam2\", 2) = {3};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   mesh.material_data = struct("E", cell(1, 2), "nu", cell(1, 2), "rho", cell(1, 2));
%!   mesh.material_data(1).E = param.E1;
%!   mesh.material_data(1).nu = param.nu1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(2).E = param.E2;
%!   mesh.material_data(2).nu = param.nu2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   grp_idx_beam1 = find([[mesh.groups.iso20].id] == 1);
%!   grp_idx_beam2 = find([[mesh.groups.iso20].id] == 2);
%!   grp_idx_clamp = find([[mesh.groups.quad8].id] == 1);
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_beam1).elements) = 1;
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_beam2).elements) = 2;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.omega = [param.OMEGAx;
%!                      param.OMEGAy;
%!                      param.OMEGAz];
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = mbdyn_solver_num_threads_default();
%!   [mat_ass.M, ...
%!    mat_ass.Mdiag, ...
%!    mat_ass.K, ...
%!    mat_ass.KOMEGA, ...
%!    mat_ass.DOMEGA, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_MASS_LUMPED, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.tau0 = sol_stat.stress.tau;
%!   mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case);
%!   opt_solver.pre_scaling = true;
%!   opt_solver.refine_max_iter = int32(10);
%!   opt_solver.solver = "pardiso";
%!   opt_solver.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_solver.symmetric = true;
%!   sol_eig(1) = fem_sol_modal(mesh, ...
%!                            dof_map, ...
%!                            mat_ass, ...
%!                            options.number_of_modes, ...
%!                            0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   sol_eig_diag = fem_sol_modal(mesh, ...
%!                                dof_map, ...
%!                                setfield(mat_ass, "M", mat_ass.Mdiag), ...
%!                                options.number_of_modes, ...
%!                                0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0;
%!   mat_ass.D = mat_ass.DOMEGA;
%!   sol_eig(2) = fem_sol_modal_damped(mesh, ...
%!                                     dof_map, ...
%!                                     mat_ass, ...
%!                                     options.number_of_modes, ...
%!                                     opt_solver);
%!   unwind_protect
%!     filename_mbdyn = [filename, ".mbdyn"];
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     mbdyn_pre_write_param_file(fd, param);
%!     fprintf(fd, "set: real t1 = %g;\n", 100 / SI_unit_second);
%!     fputs(fd, "set: integer N = 20;\n");
%!     fputs(fd, "set: real G1 = E1 / (2. * (1. + nu1));\n");
%!     fputs(fd, "set: real G2 = E2 / (2. * (1. + nu2));\n");
%!     fputs(fd, "set: real A1 = w1 * h1;\n");
%!     fputs(fd, "set: real A2 = w2 * h2;\n");
%!     fputs(fd, "set: real As1 = 9. / 10. * A1;\n");
%!     fputs(fd, "set: real As2 = 9. / 10. * A2;\n");
%!     fputs(fd, "set: real Iy1 = w1 * h1^3 / 12.;\n");
%!     fputs(fd, "set: real Iy2 = w2 * h2^3 / 12.;\n");
%!     fputs(fd, "set: real Iz1 = h1 * w1^3 / 12.;\n");
%!     fputs(fd, "set: real Iz2 = h2 * w2^3 / 12.;\n");
%!     fputs(fd, "set: real Ip1 = Iy1 + Iz1;\n");
%!     fputs(fd, "set: real Ip2 = Iy2 + Iz2;\n");
%!     fprintf(fd, "set: real c21 = %.16e;\n", interp1(w_h, c2, max(param.w1, param.h1) / min(param.w1, param.h1)));
%!     fprintf(fd, "set: real c22 = %.16e;\n", interp1(w_h, c2, max(param.w2, param.h2) / min(param.w2, param.h2)));
%!     fputs(fd, "set: real It1 = c21 * h1 * w1^3;\n");
%!     fputs(fd, "set: real It2 = c22 * h2 * w2^3;\n");
%!     fputs(fd, "set: integer ref_id_ground = 1;\n");
%!     fputs(fd, "set: integer ref_id_beam1 = 2;\n");
%!     fputs(fd, "set: integer ref_id_beam2 = 3;\n");
%!     fputs(fd, "set: integer joint_id_ground = 1;\n");
%!     fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer node_id_beam2 = 2 * N1 + node_id_beam1;\n");
%!     fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer body_id_beam2 = 2 * N1 + body_id_beam1 + 1;\n");
%!     fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer beam_id_beam2 = beam_id_beam1 + N1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 2 * t1;\n");
%!     fputs(fd, "        time step: t1 / N;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-5, test,norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             line search method, backtrack,\n");
%!     fputs(fd, "             recovery step type, constant,\n");
%!     fputs(fd, "             recovery step, 1e-6,\n");
%!     fputs(fd, "             verbose, yes,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "           eigenanalysis: list, 2, 0, 2 * t1,\n");
%!     fputs(fd, "           # output matrices, \n");
%!     fprintf(fd, "         parameter, %.16e,\n", 1e-4 / SI_unit_second);
%!     fputs(fd, "           output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1e-3 / SI_unit_second^-1, 200. / SI_unit_second^-1);
%!     fprintf(fd, "           use arpack, %d, %d, 0, suffix format, \"%%02d\";\n", 3 * options.number_of_modes, 10 * options.number_of_modes + 1);
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, t1 / 10.;\n");
%!     fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       structural nodes: 2 * N1 + 1 + 2 * N2;\n");
%!     fputs(fd, "       rigid bodies: 2 * N1 + 1 + 2 * N2 + 1;\n");
%!     fputs(fd, "       beams: N1 + N2;\n");
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_beam1,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, null;\n");
%!     fputs(fd, "reference: ref_id_beam2,\n");
%!     fputs(fd, "        reference, ref_id_beam1, l1,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, 1, 0., 1., 0., 3, 0., 0., 1.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam1, 0.5 * l1 / N1 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null;\n");
%!     endfor
%!     for i=2:(2 * param.N2 + 1)
%!       fprintf(fd, "        structural: node_id_beam2 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam2, 0.5 * l2 / N2 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null;\n");
%!     endfor
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!       fputs(fd, "               rho1 * A1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho1 * Ip1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iy1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iz1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!     endfor
%!     for i=1:(2 * param.N2 + 1)
%!       fprintf(fd, "     body: body_id_beam2 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d,\n", i - 1);
%!       fputs(fd, "               rho2 * A2 * l2 / (2 * N2 + 1), \n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho2 * Ip2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iy2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iz2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam2, eye;\n");
%!     endfor
%!     for i=1:param.N1
%!       fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E1 * A1 , G1 * As1, G1 * As1, \n");
%!       fputs(fd, "                     G1 * It1, E1 * Iy1, E1 * Iz1,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     for i=1:param.N2
%!       fprintf(fd, "        beam3: beam_id_beam2 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E2 * A2 , G2 * As2, G2 * As2, \n");
%!       fputs(fd, "                     G2 * It2, E2 * Iy2, E2 * Iz2,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = [filename, "_mbd"];
%!   if (~options.verbose)
%!     options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
%!   options_mbd.mbdyn_command = "mbdyn";
%!   options_eig.positive_frequencies = false;
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   for i=1:2
%!     res.modal(i) = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, i - 1);
%!   endfor
%!   if (options.interactive)
%!   for j=1:numel(res.modal)
%!     for i=1:numel(res.modal(j).f)
%!       opt_modal.mode_index = i;
%!       opt_modal.scale = 100;
%!       mode_file = [options_mbd.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(options_mbd.output_file, [mode_file, ".mov"], opt_modal, res.modal(j));
%!       [err, msg] = symlink([options_mbd.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = options.interactive;
%!       opt_post.f_runEasyAnim = options.interactive;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!   endfor
%!   endif
%!   idx_mode = [1, 4, 6];
%!   ## table 5.5 (11 elements), table 5.7 (11 elements)
%!   fref =  [0.9614, 1.044;
%!            7.037,  7.478;
%!           16.67,  17.09];
%!   tol = 2e-2;
%!   for i=1:2
%!     fmbd = sort(res.modal(i).f(:)) * SI_unit_second^-1;
%!     fmbd = fmbd(fmbd > 0);
%!     ffem = sort(sol_eig(i).f(:)) * SI_unit_second^-1;
%!     ffem = ffem(ffem >= 0);
%!     Nfem = min(numel(fmbd),numel(ffem));
%!     assert_simple(fmbd(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(ffem(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(fmbd(1:Nfem), ffem(1:Nfem), tol * max(fmbd(1:Nfem)));
%!   endfor
%!   assert_simple(max(abs(sol_eig_diag.f/sol_eig(1).f - 1)) < 0.04);
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST8
%! ## Oskar Wallrapp, Richard Schwertassek, 1998
%! ## Dynamik flexibler Mehrkoerpersysteme
%! ## chapter 5, table 5.7, page 242
%! ## Natural frequencies of a rotating beam structure
%! pkg load mboct-fem-pkg;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   options.verbose = false;
%!   options.number_of_modes = int32(20);
%!   options.interactive = false;
%!   c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!   w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!   param.num_fem_elem_per_sec = 1;
%!   param.N1 = int32(50);
%!   param.N2 = int32(20);
%!   param.E1 = 7e10 / SI_unit_pascal;
%!   param.E2 = 21e10 / SI_unit_pascal;
%!   param.nu1 = 0.3;
%!   param.nu2 = 0.3;
%!   param.rho1 = 3000 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.rho2 = 7895 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.l1 = 2 / SI_unit_meter;
%!   param.l2 = 0.4 / SI_unit_meter;
%!   param.h1 = 0.009 / SI_unit_meter;
%!   param.h2 = param.h1;
%!   param.w1 = 0.009 / SI_unit_meter;
%!   param.w2 = 0.0095 / SI_unit_meter;
%!   param.OMEGAx = 0 / SI_unit_second^-1;
%!   param.OMEGAy = 0 / SI_unit_second^-1;
%!   param.OMEGAz = 6 / SI_unit_second^-1;
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w1 = %.16e;\n", param.w1);
%!     fprintf(fd, "w2 = %.16e;\n", param.w2);
%!     fprintf(fd, "l1 = %.16e;\n", param.l1);
%!     fprintf(fd, "l2 = %.16e;\n", param.l2);
%!     fprintf(fd, "h1 = %.16e;\n", param.h1);
%!     fprintf(fd, "h = %.16e;\n", min([param.w1, param.w2]) / param.num_fem_elem_per_sec);
%!     fputs(fd, "Point(1) = {0, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(2) = {l1 - 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(3) = {l1 - 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(4) = {l1 + 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(5) = {l1 + 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(6) = {l1 + 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(7) = {l1 - 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(8) = {0, -0.5 * w1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line(9) = {2, 7};\n");
%!     fputs(fd, "Line(10) = {2, 5};\n");
%!     fputs(fd, "Line Loop(1) = {1, 9, 7, 8};\n");
%!     fputs(fd, "Line Loop(2) = {10, 5, 6, 9};\n");
%!     fputs(fd, "Line Loop(3) = {2, 3, 4, 10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Plane Surface(2) = {2};\n");
%!     fputs(fd, "Plane Surface(3) = {3};\n");
%!     fputs(fd, "v1 = Extrude{0,0,h1}{Surface{1,2,3};};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 1) = {7};\n");
%!     fputs(fd, "Physical Volume(\"beam1\", 1) = {1, 2};\n");
%!     fputs(fd, "Physical Volume(\"beam2\", 2) = {3};\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{1, 2, 3}; } } = h;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!     fputs(fd, "Mesh.HighOrderThresholdMin=0.9;\n");
%!     fputs(fd, "Mesh.HighOrderThresholdMax=1.1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   opt_mesh.elem_type = {"tet10h", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   mesh.material_data = struct("E", cell(1, 2), "nu", cell(1, 2), "rho", cell(1, 2));
%!   mesh.material_data(1).E = param.E1;
%!   mesh.material_data(1).nu = param.nu1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(2).E = param.E2;
%!   mesh.material_data(2).nu = param.nu2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.materials.tet10h = zeros(rows(mesh.elements.tet10h), 1, "int32");
%!   grp_idx_beam1 = find([[mesh.groups.tet10h].id] == 1);
%!   grp_idx_beam2 = find([[mesh.groups.tet10h].id] == 2);
%!   grp_idx_clamp = find([[mesh.groups.tria6h].id] == 1);
%!   mesh.materials.tet10h(mesh.groups.tet10h(grp_idx_beam1).elements) = 1;
%!   mesh.materials.tet10h(mesh.groups.tet10h(grp_idx_beam2).elements) = 2;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.tria6h(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.omega = [param.OMEGAx;
%!                      param.OMEGAy;
%!                      param.OMEGAz];
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = mbdyn_solver_num_threads_default();
%!   [mat_ass.M, ...
%!    mat_ass.K, ...
%!    mat_ass.KOMEGA, ...
%!    mat_ass.DOMEGA, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.tau0 = sol_stat.stress.tau;
%!   mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case);
%!   opt_solver.pre_scaling = true;
%!   opt_solver.refine_max_iter = int32(10);
%!   opt_solver.solver = "pardiso";
%!   opt_solver.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_solver.symmetric = true;
%!   sol_eig(1) = fem_sol_modal(mesh, ...
%!                              dof_map, ...
%!                              mat_ass, ...
%!                              options.number_of_modes, ...
%!                              0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0;
%!   mat_ass.D = mat_ass.DOMEGA;
%!   sol_eig(2) = fem_sol_modal_damped(mesh, ...
%!                                     dof_map, ...
%!                                     mat_ass, ...
%!                                     options.number_of_modes, ...
%!                                     opt_solver);
%!   unwind_protect
%!     filename_mbdyn = [filename, ".mbdyn"];
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     mbdyn_pre_write_param_file(fd, param);
%!     fprintf(fd, "set: real t1 = %g;\n", 100 / SI_unit_second);
%!     fputs(fd, "set: integer N = 20;\n");
%!     fputs(fd, "set: real G1 = E1 / (2. * (1. + nu1));\n");
%!     fputs(fd, "set: real G2 = E2 / (2. * (1. + nu2));\n");
%!     fputs(fd, "set: real A1 = w1 * h1;\n");
%!     fputs(fd, "set: real A2 = w2 * h2;\n");
%!     fputs(fd, "set: real As1 = 9. / 10. * A1;\n");
%!     fputs(fd, "set: real As2 = 9. / 10. * A2;\n");
%!     fputs(fd, "set: real Iy1 = w1 * h1^3 / 12.;\n");
%!     fputs(fd, "set: real Iy2 = w2 * h2^3 / 12.;\n");
%!     fputs(fd, "set: real Iz1 = h1 * w1^3 / 12.;\n");
%!     fputs(fd, "set: real Iz2 = h2 * w2^3 / 12.;\n");
%!     fputs(fd, "set: real Ip1 = Iy1 + Iz1;\n");
%!     fputs(fd, "set: real Ip2 = Iy2 + Iz2;\n");
%!     fprintf(fd, "set: real c21 = %.16e;\n", interp1(w_h, c2, max(param.w1, param.h1) / min(param.w1, param.h1)));
%!     fprintf(fd, "set: real c22 = %.16e;\n", interp1(w_h, c2, max(param.w2, param.h2) / min(param.w2, param.h2)));
%!     fputs(fd, "set: real It1 = c21 * h1 * w1^3;\n");
%!     fputs(fd, "set: real It2 = c22 * h2 * w2^3;\n");
%!     fputs(fd, "set: integer ref_id_ground = 1;\n");
%!     fputs(fd, "set: integer ref_id_beam1 = 2;\n");
%!     fputs(fd, "set: integer ref_id_beam2 = 3;\n");
%!     fputs(fd, "set: integer joint_id_ground = 1;\n");
%!     fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer node_id_beam2 = 2 * N1 + node_id_beam1;\n");
%!     fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer body_id_beam2 = 2 * N1 + body_id_beam1 + 1;\n");
%!     fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer beam_id_beam2 = beam_id_beam1 + N1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 2 * t1;\n");
%!     fputs(fd, "        time step: t1 / N;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-5, test,norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             line search method, backtrack,\n");
%!     fputs(fd, "             recovery step type, constant,\n");
%!     fputs(fd, "             recovery step, 1e-6,\n");
%!     fputs(fd, "             verbose, yes,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "           eigenanalysis: list, 2, 0, 2 * t1,\n");
%!     fputs(fd, "           # output matrices, \n");
%!     fprintf(fd, "         parameter, %.16e,\n", 1e-4 / SI_unit_second);
%!     fputs(fd, "           output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1e-3 / SI_unit_second^-1, 200. / SI_unit_second^-1);
%!     fprintf(fd, "           use arpack, %d, %d, 0, suffix format, \"%%02d\";\n", 3 * options.number_of_modes, 10 * options.number_of_modes + 1);
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, t1 / 10.;\n");
%!     fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       structural nodes: 2 * N1 + 1 + 2 * N2;\n");
%!     fputs(fd, "       rigid bodies: 2 * N1 + 1 + 2 * N2 + 1;\n");
%!     fputs(fd, "       beams: N1 + N2;\n");
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_beam1,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, null;\n");
%!     fputs(fd, "reference: ref_id_beam2,\n");
%!     fputs(fd, "        reference, ref_id_beam1, l1,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, 1, 0., 1., 0., 3, 0., 0., 1.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam1, 0.5 * l1 / N1 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null;\n");
%!     endfor
%!     for i=2:(2 * param.N2 + 1)
%!       fprintf(fd, "        structural: node_id_beam2 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam2, 0.5 * l2 / N2 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null;\n");
%!     endfor
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!       fputs(fd, "               rho1 * A1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho1 * Ip1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iy1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iz1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!     endfor
%!     for i=1:(2 * param.N2 + 1)
%!       fprintf(fd, "     body: body_id_beam2 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d,\n", i - 1);
%!       fputs(fd, "               rho2 * A2 * l2 / (2 * N2 + 1), \n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho2 * Ip2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iy2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iz2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam2, eye;\n");
%!     endfor
%!     for i=1:param.N1
%!       fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E1 * A1 , G1 * As1, G1 * As1, \n");
%!       fputs(fd, "                     G1 * It1, E1 * Iy1, E1 * Iz1,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     for i=1:param.N2
%!       fprintf(fd, "        beam3: beam_id_beam2 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E2 * A2 , G2 * As2, G2 * As2, \n");
%!       fputs(fd, "                     G2 * It2, E2 * Iy2, E2 * Iz2,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = [filename, "_mbd"];
%!   if (~options.verbose)
%!     options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
%!   options_mbd.mbdyn_command = "mbdyn";
%!   options_eig.positive_frequencies = false;
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   for i=1:2
%!     res.modal(i) = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, i - 1);
%!   endfor
%!   if (options.interactive)
%!   for j=1:numel(res.modal)
%!     for i=1:numel(res.modal(j).f)
%!       opt_modal.mode_index = i;
%!       opt_modal.scale = 100;
%!       mode_file = [options_mbd.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(options_mbd.output_file, [mode_file, ".mov"], opt_modal, res.modal(j));
%!       [err, msg] = symlink([options_mbd.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = options.interactive;
%!       opt_post.f_runEasyAnim = options.interactive;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!   endfor
%!   endif
%!   idx_mode = [1, 4, 6];
%!   ## table 5.5 (11 elements), table 5.7 (11 elements)
%!   fref =  [0.9614, 1.044;
%!            7.037,  7.478;
%!           16.67,  17.09];
%!   tol = 2e-2;
%!   for i=1:2
%!     fmbd = sort(res.modal(i).f(:)) * SI_unit_second^-1;
%!     fmbd = fmbd(fmbd > 0);
%!     ffem = sort(sol_eig(i).f(:)) * SI_unit_second^-1;
%!     ffem = ffem(ffem >= 0);
%!     Nfem = min(numel(fmbd),numel(ffem));
%!     assert_simple(fmbd(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(ffem(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(fmbd(1:Nfem), ffem(1:Nfem), tol * max(fmbd(1:Nfem)));
%!   endfor
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST9
%! ## Oskar Wallrapp, Richard Schwertassek, 1998
%! ## Dynamik flexibler Mehrkoerpersysteme
%! ## chapter 5, table 5.7, page 242
%! ## Natural frequencies of a rotating beam structure
%! pkg load mboct-fem-pkg;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   options.verbose = false;
%!   options.number_of_modes = int32(20);
%!   options.interactive = false;
%!   c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!   w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!   param.num_fem_elem_per_sec = 0.25;
%!   param.N1 = int32(50);
%!   param.N2 = int32(20);
%!   param.E1 = 7e10 / SI_unit_pascal;
%!   param.E2 = 21e10 / SI_unit_pascal;
%!   param.nu1 = 0.3;
%!   param.nu2 = 0.3;
%!   param.rho1 = 3000 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.rho2 = 7895 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.l1 = 2 / SI_unit_meter;
%!   param.l2 = 0.4 / SI_unit_meter;
%!   param.h1 = 0.009 / SI_unit_meter;
%!   param.h2 = param.h1;
%!   param.w1 = 0.009 / SI_unit_meter;
%!   param.w2 = 0.0095 / SI_unit_meter;
%!   param.OMEGAx = 0 / SI_unit_second^-1;
%!   param.OMEGAy = 0 / SI_unit_second^-1;
%!   param.OMEGAz = 6 / SI_unit_second^-1;
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w1 = %.16e;\n", param.w1);
%!     fprintf(fd, "w2 = %.16e;\n", param.w2);
%!     fprintf(fd, "l1 = %.16e;\n", param.l1);
%!     fprintf(fd, "l2 = %.16e;\n", param.l2);
%!     fprintf(fd, "h1 = %.16e;\n", param.h1);
%!     fprintf(fd, "h = %.16e;\n", min([param.w1, param.w2]) / param.num_fem_elem_per_sec);
%!     fputs(fd, "Point(1) = {0, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(2) = {l1 - 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(3) = {l1 - 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(4) = {l1 + 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(5) = {l1 + 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(6) = {l1 + 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(7) = {l1 - 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(8) = {0, -0.5 * w1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line(9) = {2, 7};\n");
%!     fputs(fd, "Line(10) = {2, 5};\n");
%!     fputs(fd, "Line Loop(1) = {1, 9, 7, 8};\n");
%!     fputs(fd, "Line Loop(2) = {10, 5, 6, 9};\n");
%!     fputs(fd, "Line Loop(3) = {2, 3, 4, 10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Plane Surface(2) = {2};\n");
%!     fputs(fd, "Plane Surface(3) = {3};\n");
%!     fputs(fd, "v1 = Extrude{0,0,h1}{Surface{1,2,3};};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 1) = {7};\n");
%!     fputs(fd, "Physical Volume(\"beam1\", 1) = {1, 2};\n");
%!     fputs(fd, "Physical Volume(\"beam2\", 2) = {3};\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{1, 2, 3}; } } = h;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=1;\n");
%!     fputs(fd, "Mesh.HighOrderThresholdMin=0.9;\n");
%!     fputs(fd, "Mesh.HighOrderThresholdMax=1.1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "3", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   opt_mesh.elem_type = {"tet20", "tria10"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_mesh));
%!   mesh.material_data = struct("E", cell(1, 2), "nu", cell(1, 2), "rho", cell(1, 2));
%!   mesh.material_data(1).E = param.E1;
%!   mesh.material_data(1).nu = param.nu1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(2).E = param.E2;
%!   mesh.material_data(2).nu = param.nu2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.materials.tet20 = zeros(rows(mesh.elements.tet20), 1, "int32");
%!   grp_idx_beam1 = find([[mesh.groups.tet20].id] == 1);
%!   grp_idx_beam2 = find([[mesh.groups.tet20].id] == 2);
%!   grp_idx_clamp = find([[mesh.groups.tria10].id] == 1);
%!   mesh.materials.tet20(mesh.groups.tet20(grp_idx_beam1).elements) = 1;
%!   mesh.materials.tet20(mesh.groups.tet20(grp_idx_beam2).elements) = 2;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.tria10(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.omega = [param.OMEGAx;
%!                      param.OMEGAy;
%!                      param.OMEGAz];
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = mbdyn_solver_num_threads_default();
%!   [mat_ass.M, ...
%!    mat_ass.Mdiag, ...
%!    mat_ass.K, ...
%!    mat_ass.KOMEGA, ...
%!    mat_ass.DOMEGA, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_MASS_LUMPED, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.tau0 = sol_stat.stress.tau;
%!   mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case);
%!   opt_solver.pre_scaling = true;
%!   opt_solver.refine_max_iter = int32(10);
%!   opt_solver.solver = "pardiso";
%!   opt_solver.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_solver.symmetric = true;
%!   sol_eig(1) = fem_sol_modal(mesh, ...
%!                            dof_map, ...
%!                            mat_ass, ...
%!                            options.number_of_modes, ...
%!                            0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   sol_eig_diag = fem_sol_modal(mesh, ...
%!                                dof_map, ...
%!                                setfield(mat_ass, "M", mat_ass.Mdiag), ...
%!                                options.number_of_modes, ...
%!                                0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0;
%!   mat_ass.D = mat_ass.DOMEGA;
%!   sol_eig(2) = fem_sol_modal_damped(mesh, ...
%!                                     dof_map, ...
%!                                     mat_ass, ...
%!                                     options.number_of_modes, ...
%!                                     opt_solver);
%!   unwind_protect
%!     filename_mbdyn = [filename, ".mbdyn"];
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     mbdyn_pre_write_param_file(fd, param);
%!     fprintf(fd, "set: real t1 = %g;\n", 100 / SI_unit_second);
%!     fputs(fd, "set: integer N = 20;\n");
%!     fputs(fd, "set: real G1 = E1 / (2. * (1. + nu1));\n");
%!     fputs(fd, "set: real G2 = E2 / (2. * (1. + nu2));\n");
%!     fputs(fd, "set: real A1 = w1 * h1;\n");
%!     fputs(fd, "set: real A2 = w2 * h2;\n");
%!     fputs(fd, "set: real As1 = 9. / 10. * A1;\n");
%!     fputs(fd, "set: real As2 = 9. / 10. * A2;\n");
%!     fputs(fd, "set: real Iy1 = w1 * h1^3 / 12.;\n");
%!     fputs(fd, "set: real Iy2 = w2 * h2^3 / 12.;\n");
%!     fputs(fd, "set: real Iz1 = h1 * w1^3 / 12.;\n");
%!     fputs(fd, "set: real Iz2 = h2 * w2^3 / 12.;\n");
%!     fputs(fd, "set: real Ip1 = Iy1 + Iz1;\n");
%!     fputs(fd, "set: real Ip2 = Iy2 + Iz2;\n");
%!     fprintf(fd, "set: real c21 = %.16e;\n", interp1(w_h, c2, max(param.w1, param.h1) / min(param.w1, param.h1)));
%!     fprintf(fd, "set: real c22 = %.16e;\n", interp1(w_h, c2, max(param.w2, param.h2) / min(param.w2, param.h2)));
%!     fputs(fd, "set: real It1 = c21 * h1 * w1^3;\n");
%!     fputs(fd, "set: real It2 = c22 * h2 * w2^3;\n");
%!     fputs(fd, "set: integer ref_id_ground = 1;\n");
%!     fputs(fd, "set: integer ref_id_beam1 = 2;\n");
%!     fputs(fd, "set: integer ref_id_beam2 = 3;\n");
%!     fputs(fd, "set: integer joint_id_ground = 1;\n");
%!     fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer node_id_beam2 = 2 * N1 + node_id_beam1;\n");
%!     fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer body_id_beam2 = 2 * N1 + body_id_beam1 + 1;\n");
%!     fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer beam_id_beam2 = beam_id_beam1 + N1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 2 * t1;\n");
%!     fputs(fd, "        time step: t1 / N;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-5, test,norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             line search method, backtrack,\n");
%!     fputs(fd, "             recovery step type, constant,\n");
%!     fputs(fd, "             recovery step, 1e-6,\n");
%!     fputs(fd, "             verbose, yes,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "           eigenanalysis: list, 2, 0, 2 * t1,\n");
%!     fputs(fd, "           # output matrices, \n");
%!     fprintf(fd, "         parameter, %.16e,\n", 1e-4 / SI_unit_second);
%!     fputs(fd, "           output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1e-3 / SI_unit_second^-1, 200. / SI_unit_second^-1);
%!     fprintf(fd, "           use arpack, %d, %d, 0, suffix format, \"%%02d\";\n", 3 * options.number_of_modes, 10 * options.number_of_modes + 1);
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, t1 / 10.;\n");
%!     fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       structural nodes: 2 * N1 + 1 + 2 * N2;\n");
%!     fputs(fd, "       rigid bodies: 2 * N1 + 1 + 2 * N2 + 1;\n");
%!     fputs(fd, "       beams: N1 + N2;\n");
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_beam1,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, null;\n");
%!     fputs(fd, "reference: ref_id_beam2,\n");
%!     fputs(fd, "        reference, ref_id_beam1, l1,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, 1, 0., 1., 0., 3, 0., 0., 1.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam1, 0.5 * l1 / N1 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null;\n");
%!     endfor
%!     for i=2:(2 * param.N2 + 1)
%!       fprintf(fd, "        structural: node_id_beam2 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam2, 0.5 * l2 / N2 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null;\n");
%!     endfor
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!       fputs(fd, "               rho1 * A1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho1 * Ip1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iy1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iz1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!     endfor
%!     for i=1:(2 * param.N2 + 1)
%!       fprintf(fd, "     body: body_id_beam2 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d,\n", i - 1);
%!       fputs(fd, "               rho2 * A2 * l2 / (2 * N2 + 1), \n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho2 * Ip2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iy2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iz2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam2, eye;\n");
%!     endfor
%!     for i=1:param.N1
%!       fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E1 * A1 , G1 * As1, G1 * As1, \n");
%!       fputs(fd, "                     G1 * It1, E1 * Iy1, E1 * Iz1,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     for i=1:param.N2
%!       fprintf(fd, "        beam3: beam_id_beam2 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E2 * A2 , G2 * As2, G2 * As2, \n");
%!       fputs(fd, "                     G2 * It2, E2 * Iy2, E2 * Iz2,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = [filename, "_mbd"];
%!   if (~options.verbose)
%!     options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
%!   options_mbd.mbdyn_command = "mbdyn";
%!   options_eig.positive_frequencies = false;
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   for i=1:2
%!     res.modal(i) = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, i - 1);
%!   endfor
%!   if (options.interactive)
%!   for j=1:numel(res.modal)
%!     for i=1:numel(res.modal(j).f)
%!       opt_modal.mode_index = i;
%!       opt_modal.scale = 100;
%!       mode_file = [options_mbd.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(options_mbd.output_file, [mode_file, ".mov"], opt_modal, res.modal(j));
%!       [err, msg] = symlink([options_mbd.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = options.interactive;
%!       opt_post.f_runEasyAnim = options.interactive;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!   endfor
%!   endif
%!   idx_mode = [1, 4, 6];
%!   ## table 5.5 (11 elements), table 5.7 (11 elements)
%!   fref =  [0.9614, 1.044;
%!            7.037,  7.478;
%!           16.67,  17.09];
%!   tol = 2e-2;
%!   for i=1:2
%!     fmbd = sort(res.modal(i).f(:)) * SI_unit_second^-1;
%!     fmbd = fmbd(fmbd > 0);
%!     ffem = sort(sol_eig(i).f(:)) * SI_unit_second^-1;
%!     ffem = ffem(ffem >= 0);
%!     Nfem = min(numel(fmbd),numel(ffem));
%!     assert_simple(fmbd(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(ffem(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(fmbd(1:Nfem), ffem(1:Nfem), tol * max(fmbd(1:Nfem)));
%!   endfor
%!   assert_simple(max(abs(sol_eig_diag.f / sol_eig(1).f - 1)) < 0.03);
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST10
%! pkg load mboct-fem-pkg;
%! printf("fem_cms_create2: test7\n");
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   a = 800e-3 / SI_unit_meter;
%!   b = 40e-3 / SI_unit_meter;
%!   c = 10e-3 / SI_unit_meter;
%!   d = 0e-3 / SI_unit_meter;
%!   h = 2 * c;
%!   options.interactive = false;
%!   options.plot = true;
%!   options.verbose = false;
%!   options.number_of_beams = int32(40);
%!   options.number_of_threads = mbdyn_solver_num_threads_default();
%!   if (options.plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fprintf(fd, "a = %.16e;\n", a);
%!     fprintf(fd, "b = %.16e;\n", b);
%!     fprintf(fd, "c = %.16e;\n", c);
%!     fprintf(fd, "h = %.16e;\n", h);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(2) = {0,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(3) = {a,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(4) = {a, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Line Loop(1) = {1, 2, 3, 4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,c}{Surface{1}; Layers{Max(1, Round(c / h))}; Recombine;};\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {v1[2]};\n");
%!     fputs(fd, "Physical Surface(2) = {v1[4]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   mesh.material_data.E = 70000e6 / SI_unit_pascal;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 2700 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data.alpha = 0e-5 / (1 / SI_unit_second);
%!   mesh.material_data.beta = 0e-5 / (SI_unit_second);
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   grp_idx_beam = find([[mesh.groups.iso20].id] == 1);
%!   grp_idx_clamp = find([[mesh.groups.quad8].id] == 1);
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_beam).elements) = 1;
%!   cms_opt.number_of_threads = options.number_of_threads;
%!   cms_opt.algorithm = "diag-shift-invert";
%!   cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!   cms_opt.nodes.modal.name = "node_id_modal";
%!   cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!   cms_opt.nodes.interfaces.name = "node_id_interface1";
%!   mesh.nodes(cms_opt.nodes.modal.number, 1:3) = [0, 0, 0];
%!   mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = [a + d, 0, 0];
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, [2], [cms_opt.nodes.interfaces.number], "quad8");
%!   cms_opt.refine_max_iter = 30;
%!   cms_opt.pre_scaling = false;
%!   cms_opt.solver = "umfpack"; ## Need to test at least one case with an unsymmetric solver
%!   cms_opt.modes.number = 20;
%!   cms_opt.tolerance_tau = -1;
%!   cms_opt.element.name = "elem_id_modal";
%!   cms_opt.create_binary = true;
%!   cms_opt.use_binary = true;
%!   cms_opt.update_binary = true;
%!   cms_opt.invariants = true;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, :) = true; ## Avoid singular matrix
%!   [mesh_cms, mat_ass_cms, dof_map_cms, sol_eig_cms, cms_opt, sol_tau_cms] = fem_cms_create2(mesh, load_case_dof, cms_opt);
%!   fem_cms_export(filename, mesh_cms, dof_map_cms, mat_ass_cms, cms_opt);
%!   pert.omega = [1e2; 3e2; 2e2] / (SI_unit_rad / SI_unit_second);
%!   pert.omegadot = [1e5; 1e3; 3e3] / (SI_unit_rad / SI_unit_second^2);
%!   pert.loads = [[1e4; 1e3; 1e2] / (SI_unit_newton);
%!                 [1e2; 1e1; 1e1] / (SI_unit_newton * SI_unit_meter)];
%!   pert.g = [1e4; -1e3; -1e2] / (SI_unit_meter / SI_unit_second^2);
%!   pert.a = [-1e4; 1e3; 1e2] / (SI_unit_meter / SI_unit_second^2);
%!   empty_cell = cell(7, 3, 2);
%!   res = struct("info", empty_cell, ...
%!                "t", empty_cell, ...
%!                "trajectory", empty_cell, ...
%!                "deformation", empty_cell, ...
%!                "velocity", empty_cell, ...
%!                "acceleration", empty_cell, ...
%!                "node_id", empty_cell, ...
%!                "force", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "orientation_description", empty_cell, ...
%!                "drive_id", empty_cell, ...
%!                "drive_value", empty_cell, ...
%!                "modal", empty_cell);
%!   empty_cell = cell(7, 3);
%!   param = struct("omega", empty_cell, ...
%!                  "omegadot", empty_cell, ...
%!                  "F1", empty_cell, ...
%!                  "M1", empty_cell, ...
%!                  "a", empty_cell, ...
%!                  "g", empty_cell, ...
%!                  "t1", empty_cell, ...
%!                  "holonomic", empty_cell);
%!   idx_j = 1:rows(param);
%!   idx_k = 1:columns(param);
%!   for j=idx_j
%!     for k=idx_k
%!       param(j, k).omega = zeros(3, 1);
%!       param(j, k).omegadot = zeros(3, 1);
%!       param(j, k).F1 = zeros(3, 1);
%!       param(j, k).M1 = zeros(3, 1);
%!       param(j, k).a = zeros(3, 1);
%!       param(j, k).g = zeros(3, 1);
%!       param(j, k).holonomic = false;
%!       param(j, k).gamma = zeros(3, 1);
%!       param(j, k).N = 50;
%!       switch (j)
%!         case 1
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).gamma = [20; 45; 30] * pi / 180;
%!         case 2
%!           param(j, k).omega(k) = pert.omega(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           param(j, k).holonomic = true;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 2000;
%!         case 3
%!           param(j, k).omegadot(k) = pert.omegadot(k);
%!           param(j, k).t1 = 1e-2 / SI_unit_second;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 200;
%!         case 4
%!           param(j, k).F1(k) = pert.loads(k);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           switch (k)
%!             case 1
%!               param(j, k).gamma(3) = 5 * pi / 180;
%!             case 2
%!               param(j, k).gamma(3) = 45 * pi / 180;
%!             case 3
%!               param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 5
%!           param(j, k).M1(k) = pert.loads(k + 3);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           param(j, k).gamma(1) = 30 * pi / 180;
%!         case 6
%!           param(j, k).a(k) = pert.a(k);
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).N = 200;
%!           switch (k)
%!             case 1
%!               param(j, k).gamma(3) = 5 * pi / 180;
%!             case 2
%!               param(j, k).gamma(3) = 45 * pi / 180;
%!             case 3
%!               param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 7
%!           param(j, k).g(k) = pert.g(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           switch (k)
%!             case 1
%!               param(j, k).gamma(3) = 5 * pi / 180;
%!             case 2
%!               param(j, k).gamma(3) = 45 * pi / 180;
%!             case 3
%!               param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!       endswitch
%!       for l=1:3
%!         switch (l)
%!           case 3
%!             nodes_file = sprintf("%s_%d_%d_%d.nod", filename, j, k, l);
%!             csl_file = sprintf("%s_%d_%d_%d.csl", filename, j, k, l);
%!             elem_file = sprintf("%s_%d_%d_%d.elm", filename, j, k, l);
%!             opt_mbd_mesh = struct();
%!             opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_clamp";
%!             opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!             opt_mbd_mesh.struct_nodes.type(cms_opt.nodes.interfaces.number) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!             opt_mbd_mesh.struct_nodes.type(cms_opt.nodes.modal.number) = MBDYN_NODE_TYPE_STATIC_STRUCT_DISP;
%!             opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!             opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!             load_case_empty = struct();
%!             opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!         endswitch
%!         fd = -1;
%!         unwind_protect
%!           filename_mbdyn = sprintf("%s_%d_%d_%d.mbdyn", filename, j, k, l);
%!           [fd, msg] = fopen(filename_mbdyn, "w");
%!           if (fd == -1)
%!             error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!           endif
%!           fprintf(fd, "set: real a = %.16e;\n", a);
%!           fprintf(fd, "set: real b = %.16e;\n", b);
%!           fprintf(fd, "set: real c = %.16e;\n", c);
%!           fprintf(fd, "set: real d = %.16e;\n", d);
%!           for i=1:3
%!             fprintf(fd, "set: real gamma%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).gamma(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGA%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omega(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGAP%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omegadot(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real F1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).F1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real M1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).M1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real a%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).a(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real g%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).g(i));
%!           endfor
%!           fprintf(fd, "set: real t1 = %.16e;\n", param(j, k).t1);
%!           fprintf(fd, "set: integer N = %d;\n", param(j, k).N);
%!           fprintf(fd, "set: integer M = %d;\n", options.number_of_beams);
%!           fputs(fd, "set: integer ref_id_ground = 1;\n");
%!           fputs(fd, "set: integer ref_id_tilt = 2;\n");
%!           fputs(fd, "set: integer joint_id_ground = 1;\n");
%!           fputs(fd, "set: integer force_id1;\n");
%!           fputs(fd, "set: integer torque_id1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_PHI1 = 1;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP1 = 3;\n");
%!           fputs(fd, "set: integer drive_id_PHI2 = 4;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA2 = 5;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP2 = 6;\n");
%!           fputs(fd, "set: integer drive_id_PHIx = 7;\n");
%!           fputs(fd, "set: integer drive_id_PHIy = 8;\n");
%!           fputs(fd, "set: integer drive_id_PHIz = 9;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAx = 10;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAy = 11;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAz = 12;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPx = 13;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPy = 14;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPz = 15;\n");
%!           fputs(fd, "set: integer drive_id_Xx = 16;\n");
%!           fputs(fd, "set: integer drive_id_Xy = 17;\n");
%!           fputs(fd, "set: integer drive_id_Xz = 18;\n");
%!           fputs(fd, "set: integer drive_id_XPx = 19;\n");
%!           fputs(fd, "set: integer drive_id_XPy = 20;\n");
%!           fputs(fd, "set: integer drive_id_XPz = 21;\n");
%!           fputs(fd, "set: integer drive_id_XPPx = 22;\n");
%!           fputs(fd, "set: integer drive_id_XPPy = 23;\n");
%!           fputs(fd, "set: integer drive_id_XPPz = 24;\n");
%!           fputs(fd, "set: integer drive_id_gx = 25;\n");
%!           fputs(fd, "set: integer drive_id_gy = 26;\n");
%!           fputs(fd, "set: integer drive_id_gz = 27;\n");
%!           fputs(fd, "set: integer drive_id_F1x = 28;\n");
%!           fputs(fd, "set: integer drive_id_F1y = 29;\n");
%!           fputs(fd, "set: integer drive_id_F1z = 30;\n");
%!           fputs(fd, "set: integer drive_id_M1x = 31;\n");
%!           fputs(fd, "set: integer drive_id_M1y = 32;\n");
%!           fputs(fd, "set: integer drive_id_M1z = 33;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "set: integer ref_id_modal = 3;\n");
%!               fputs(fd, "set: integer node_id_modal = 1;\n");
%!               fputs(fd, "set: integer ref_id_interface1 = 4;\n");
%!               fputs(fd, "set: integer node_id_interface1 = 2;\n");
%!               fputs(fd, "set: integer elem_id_modal = 2;\n");
%!             case 2
%!               fputs(fd, "set: integer ref_id_beam1 = 3;\n");
%!               fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!               fprintf(fd, "set: real E = %.16e;\n", mesh.material_data.E);
%!               fprintf(fd, "set: real nu = %.16e;\n", mesh.material_data.nu);
%!               fprintf(fd, "set: real rho = %.16e;\n", mesh.material_data.rho);
%!               fprintf(fd, "set: real alpha = %.16e;\n", mesh.material_data.alpha);
%!               fprintf(fd, "set: real beta = %.16e;\n", mesh.material_data.beta);
%!               fputs(fd, "set: real G = E / (2. * (1. + nu));\n");
%!               fputs(fd, "set: real A = b * c;\n");
%!               fputs(fd, "set: real As = 9. / 10. * A;\n");
%!               fputs(fd, "set: real Iy = b * c^3 / 12.;\n");
%!               fputs(fd, "set: real Iz = c * b^3 / 12.;\n");
%!               fputs(fd, "set: real Ip = Iy + Iz;\n");
%!               c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!               w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!               fprintf(fd, "set: real It = %.16e;\n", interp1(w_h, c2, max(c, b) / min(c, b)) * max(c, b) * min(c, b)^3);
%!             case 3
%!               fputs(fd, "set: integer ref_id_clamp = 3;\n");
%!               fprintf(fd, "set: integer node_id_interface1 = %d;\n", cms_opt.nodes.interfaces.number);
%!           endswitch
%!           fputs(fd, "begin: data;\n");
%!           fputs(fd, "        problem: initial value;\n");
%!           fputs(fd, "end: data;\n");
%!           fputs(fd, "begin: initial value;\n");
%!           fputs(fd, "        initial time: 0;\n");
%!           fputs(fd, "        final time: t1;\n");
%!           fputs(fd, "        time step: t1 / N;\n");
%!           fputs(fd, "        max time step: t1 / N;\n");
%!           fputs(fd, "        min time step: t1 / N;\n");
%!           fputs(fd, "        method: ss4, 0.;\n");
%!           fputs(fd, "        tolerance: 1e-4, test, minmax, 1e-6, test,minmax;\n");
%!           fputs(fd, "        max iterations: 1000;\n");
%!           fputs(fd, "        derivatives max iterations: 50;\n");
%!           fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!           fputs(fd, "        derivatives tolerance: 1e-6, 1;\n");
%!           fputs(fd, "        strategy: factor, 0.8, 3, 1.25, 3, 3, 6;\n");
%!           fputs(fd, "        output: iterations, solver condition number, stat, yes, cpu time;\n");
%!           fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!           fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!           fputs(fd, "             keep jacobian matrix,\n");
%!           fputs(fd, "             inner iterations before assembly, 6,\n");
%!           fputs(fd, "             jacobian operator, newton krylov,\n");
%!           fputs(fd, "             solver, line search based,\n");
%!           fputs(fd, "             line search method, backtrack,\n");
%!           fputs(fd, "             recovery step type, constant,\n");
%!           fputs(fd, "             recovery step, 1e-6,\n");
%!           fputs(fd, "             verbose, yes,\n");
%!           fputs(fd, "             forcing term, type 2,\n");
%!           fputs(fd, "             direction, newton,\n");
%!           fputs(fd, "             weighted rms absolute tolerance, 0,\n");
%!           fputs(fd, "             weighted rms relative tolerance, 0,\n");
%!           fputs(fd, "             linear solver, gmres,\n");
%!           fputs(fd, "             linear solver max iterations, 30,\n");
%!           fputs(fd, "             minimum step, 1e-12,\n");
%!           fputs(fd, "             krylov subspace size, 30;\n");
%!           fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!           switch (l)
%!           case 3
%!             fprintf(fd, "        threads: assembly, %d;\n", cms_opt.number_of_threads);
%!           otherwise
%!             fputs(fd, "        threads: assembly, 1;\n");
%!           endswitch
%!           fputs(fd, "    eigenanalysis: list, 1, t1,\n");
%!           fputs(fd, "    output eigenvectors,\n");
%!           fputs(fd, "    # output matrices, \n");
%!           fputs(fd, "        output geometry,\n");
%!           fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1. / SI_unit_second^-1, 10000 / SI_unit_second^-1);
%!           switch (l)
%!             case 1
%!               fputs(fd, "    use lapack, balance, permute, suffix format, \"%02d\";\n");
%!             case {2, 3}
%!               fprintf(fd, "    use arpack,%d,%d,0.,suffix format,\"%%02d\";\n", cms_opt.modes.number, 2 * cms_opt.modes.number + 1);
%!           endswitch
%!           fputs(fd, "end: initial value;\n");
%!           fputs(fd, "begin: control data;\n");
%!           fputs(fd, "        output meter: closest next, 0., forever, t1 / 20.;\n");
%!           switch (l)
%!             case {2, 3}
%!               fputs(fd, "        rigid body kinematics: drive,\n");
%!               fputs(fd, "            angular velocity,\n");
%!               fputs(fd, "                   component,\n");
%!               for i=1:3
%!                 fprintf(fd, "                postponed, drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!               endfor
%!               fputs(fd, "            acceleration,\n");
%!               fputs(fd, "                   component,\n");
%!               for i=1:3
%!                 fprintf(fd, "               postponed, drive_id_XPP%s,\n", {"x", "y", "z"}{i});
%!               endfor
%!               fputs(fd, "            angular acceleration,\n");
%!               fputs(fd, "                   component");
%!               for i=1:3
%!                 fprintf(fd, ",\n               postponed, drive_id_OMEGAP%s", {"x","y","z"}{i});
%!               endfor
%!               fputs(fd, ";\n");
%!           endswitch
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       use automatic differentiation;\n");
%!           fputs(fd, "       default output: none, structural nodes, drive callers;\n");
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       output precision: 16;\n");
%!           fputs(fd, "       max iterations: 0;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural nodes:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# interface 1\n");
%!               fputs(fd, "        ;\n");
%!               fputs(fd, "        joints:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# ground\n");
%!               fputs(fd, "        ;\n");
%!               fputs(fd, "        forces: 2;\n");
%!             case 2
%!               fputs(fd, "       structural nodes: 2 * M + 1;\n");
%!               fputs(fd, "       rigid bodies: 2 * M + 1;\n");
%!               fputs(fd, "       beams: M;\n");
%!               fputs(fd, "       joints: 1;\n");
%!               fputs(fd, "       forces: 2;\n");
%!             case 3
%!               fprintf(fd, "     structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!               fprintf(fd, "     solids: %d;\n", opt_mbd_mesh.solids.number);
%!               fprintf(fd, "     genels: %d;\n", opt_mbd_mesh.genels.number);
%!               fprintf(fd, "     joints: %d;\n", opt_mbd_mesh.joints.number);
%!               fprintf(fd, "     forces: %d;\n", opt_mbd_mesh.forces.number + 2);
%!           endswitch
%!           fputs(fd, "        gravity;\n");
%!           fputs(fd, "end: control data;\n");
%!           fputs(fd, "drive caller: drive_id_PHI1, string, \"(((pi*Time)/(2*t1)-sin((pi*Time)/t1)/2)*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA1, string, \"sin((pi*Time)/(2*t1))^2\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP1, string, \"(pi*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1)))/t1\";\n");
%!           fputs(fd, "drive caller: drive_id_PHI2, string, \"-(4*sin((pi*Time)/(2*t1))^3*t1^2)/(3*pi^2)\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA2, string, \"-(2*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1))^2*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP2, string, \"-(3*cos((pi*Time)/(2*t1))^2-1)*sin((pi*Time)/(2*t1))\";\n");
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_PHI%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_PHI1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGA1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGAP%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGAP1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_X%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XPP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_g%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, g%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_F1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, F1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_M1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, M1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, "reference: ref_id_ground,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, eye,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, null;\n");
%!           fputs(fd, "reference: ref_id_tilt,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, euler123, gammax, gammay, gammaz,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, null;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "reference: ref_id_modal,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!               fputs(fd, "reference: ref_id_interface1,\n");
%!               fputs(fd, "        reference, ref_id_modal, a + d,  0., 0.,\n");
%!               fputs(fd, "        reference, ref_id_modal, eye,\n");
%!               fputs(fd, "        reference, ref_id_modal, null,\n");
%!               fputs(fd, "        reference, ref_id_modal, null;\n");
%!             case 2
%!               fputs(fd, "reference: ref_id_beam1,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!             case 3
%!               fputs(fd, "reference: ref_id_clamp,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!               fprintf(fd, "include: \"%s\";\n", csl_file);
%!           endswitch
%!           fputs(fd, "begin: nodes;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural: node_id_modal, modal,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, eye,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, null, accelerations, yes;\n");
%!               fputs(fd, "        structural: node_id_interface1, static,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, eye,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null;\n");
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!                 fprintf(fd, "                reference, ref_id_beam1, 0.5 * a / M * %d, 0., 0.,\n", i - 1);
%!                 fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null, accelerations, yes;\n");
%!               endfor
%!             case 3
%!               fprintf(fd, "include: \"%s\";\n", nodes_file);
%!           endswitch
%!           fputs(fd, "end: nodes;\n");
%!           fputs(fd, "begin: elements;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        joint: joint_id_ground, total pin joint,\n");
%!               fprintf(fd, "                %s,\n", {"node_id_modal", "node_id_beam1"}{l});
%!               fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fputs(fd, "               position constraint,\n");
%!               if (~param(j, k).holonomic)
%!                 fputs(fd, "                        velocity, velocity, velocity,\n");
%!                 fputs(fd, "                        component,\n");
%!                 for i=1:3
%!                   fprintf(fd, "                      reference, drive_id_XP%s,\n", {"x", "y", "z"}{i});
%!                 endfor
%!               else
%!                 fputs(fd, "                        active, active, active,\n");
%!                 fputs(fd, "                        component,\n");
%!                 for i=1:3
%!                   fprintf(fd, "                      reference, drive_id_X%s,\n", {"x", "y", "z"}{i});
%!                 endfor
%!               endif
%!               fputs(fd, "               orientation constraint,\n");
%!               if (~param(j, k).holonomic)
%!                 fputs(fd, "                        angular velocity, angular velocity, angular velocity,\n");
%!                 fputs(fd, "                        component");
%!                 for i=1:3
%!                   fprintf(fd, ",\n                   reference, drive_id_OMEGA%s", {"x","y","z"}{i});
%!                 endfor
%!               else
%!                 fputs(fd, "                        active, active, active,\n");
%!                 fputs(fd, "                        component");
%!                 for i=1:3
%!                   fprintf(fd, ",\n                   reference, drive_id_PHI%s", {"x","y","z"}{i});
%!                 endfor
%!               endif
%!               fputs(fd, ";\n");
%!             case 2
%!               fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!           endswitch
%!           switch (l)
%!             case 1
%!               fprintf(fd, "        include: \"%s.elm\";\n", filename);
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!                 fputs(fd, "               rho * A * a / (2 * M + 1),\n");
%!                 fputs(fd, "               reference, node, null, \n");
%!                 fputs(fd, "               diag,   rho * Ip * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iy * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iz * a / (2 * M + 1),\n");
%!                 fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!               endfor
%!               for i=1:options.number_of_beams
%!                 fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "               linear elastic generic, \n");
%!                 fputs(fd, "               diag, E * A , G * As, G * As, \n");
%!                 fputs(fd, "                     G * It, E * Iy, E * Iz,\n");
%!                 fputs(fd, "               same,\n");
%!                 fputs(fd, "               same;\n");
%!               endfor
%!             case 3
%!               fprintf(fd, "include: \"%s\";\n", elem_file);
%!           endswitch
%!           fprintf(fd, "        force: force_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M", sprintf("%d", cms_opt.nodes.interfaces.number)}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_F1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fprintf(fd, "        couple: torque_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M", sprintf("%d", cms_opt.nodes.interfaces.number)}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_M1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fputs(fd, "        gravity: uniform, component");
%!           for i=1:3
%!             fprintf(fd, ",\n       reference, drive_id_g%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd,";\n");
%!           fputs(fd, "end: elements;\n");
%!         unwind_protect_cleanup
%!           if (fd ~= -1)
%!             fclose(fd);
%!           endif
%!           fd = -1;
%!         end_unwind_protect
%!         options_mbd.output_file = sprintf("%s_%d_%d_%d_mbd", filename, j, k, l);
%!         if (~options.verbose)
%!           options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!         endif
%!         options_mbd.mbdyn_command = "mbdyn -C";
%!         options_eig.positive_frequencies = false;
%!         if (options.verbose)
%!           shell(sprintf("cat %s | nl", filename_mbdyn));
%!           switch (l)
%!           case 3
%!             shell(sprintf("cat %s | nl", csl_file));
%!             shell(sprintf("cat %s | nl", nodes_file));
%!             shell(sprintf("cat %s | nl", elem_file));
%!           endswitch
%!         endif
%!         res(j, k, l).info = mbdyn_solver_run(filename_mbdyn, options_mbd);
%!         output_file_rel_frame = [options_mbd.output_file, "_rel"];
%!         mbdyn_post_abs_to_rel(1, options_mbd.output_file, output_file_rel_frame, 0);
%!         exts = {".log", ".out"};
%!         for i=1:numel(exts)
%!           [err, msg] = symlink([options_mbd.output_file, exts{i}], [output_file_rel_frame, exts{i}]);
%!           if (err ~= 0)
%!             error("failed to create symlink: %s", msg);
%!           endif
%!         endfor
%!         [res(j, k, l).t, ...
%!          res(j, k, l).trajectory, ...
%!          res(j, k, l).deformation, ...
%!          res(j, k, l).velocity, ...
%!          res(j, k, l).acceleration, ...
%!          res(j, k, l).node_id, ...
%!          res(j, k, l).force, ...
%!          res(j, k, l).force_id, ...
%!          res(j, k, l).force_node_id, ...
%!          res(j, k, l).orientation_description] = mbdyn_post_load_output_struct(output_file_rel_frame);
%!         res(j, k, l).log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!         [res(j, k, l).drive_id, ...
%!          res(j, k, l).drive_value] = mbdyn_post_load_output_drv(options_mbd.output_file, [], numel(res(j, k, l).t));
%!         res(j, k, l).modal = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, 0);
%!       endfor
%!     endfor
%!   endfor
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, :) = true; ## Avoid singular matrix
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = options.number_of_threads;
%!   load_case = struct("omega", empty_cell, ...
%!                      "omegadot", empty_cell, ...
%!                      "loads", empty_cell, ...
%!                      "loaded_nodes", empty_cell, ...
%!                      "joints", empty_cell, ...
%!                      "g", empty_cell, ...
%!                      "tau0", empty_cell);
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   for i=1:numel(load_case)
%!     load_case(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!     load_case(i).loads = zeros(1, 6);
%!     load_case(i).omega = zeros(3, 1);
%!     load_case(i).omegadot = zeros(3, 1);
%!     load_case(i).g = zeros(3, 1);
%!     load_case(i).tau0.iso20 = zeros(rows(mesh.elements.iso20), columns(mesh.elements.iso20), 6);
%!   endfor
%!   sol_eig = struct("def", empty_cell, "lambda", empty_cell, "f", empty_cell, "D", empty_cell);
%!   sol_eig_red = struct("lambda_red", empty_cell, "Ured", empty_cell);
%!   for j=idx_j
%!     for k=idx_k
%!       R = euler123_to_rotation_matrix(param(j, k).gamma);
%!       load_case(j, k).omega = R.' * param(j, k).omega;
%!       load_case(j, k).omegadot = R.' * param(j, k).omegadot;
%!       load_case(j, k).loads = [(R.' * param(j, k).F1).', (R.' * param(j, k).M1).'];
%!       load_case(j, k).g = R.' * (param(j, k).g - param(j, k).a);
%!       [mat_ass.M, ...
%!        mat_ass.D, ...
%!        mat_ass.K, ...
%!        mat_ass.KOMEGA, ...
%!        mat_ass.KOMEGA_DOT, ...
%!        mat_ass.DOMEGA, ...
%!        mat_ass.R, ...
%!        mat_ass.mat_info, ...
%!        mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                            dof_map, ...
%!                                            [FEM_MAT_MASS, ...
%!                                             FEM_MAT_DAMPING, ...
%!                                             FEM_MAT_STIFFNESS, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                             FEM_MAT_DAMPING_OMEGA, ...
%!                                             FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case(j, k));
%!       cms_opt.symmetric = false;
%!       sol_statjk = fem_sol_static(mesh, dof_map, mat_ass, cms_opt);
%!       sol_statjk.stress = fem_ass_matrix(mesh, ...
%!                                          dof_map, ...
%!                                          [FEM_VEC_STRESS_CAUCH], ...
%!                                          load_case(j, k), ...
%!                                          sol_statjk);
%!       sol_stat(j, k) = sol_statjk;
%!       load_case(j, k).tau0 = sol_stat(j, k).stress.tau;
%!       mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS_TAU0], ...
%!                                      load_case(j, k));
%!       mat_ass.K += mat_ass.KOMEGA + mat_ass.KOMEGA_DOT + mat_ass.KTAU0;
%!       mat_ass.D += mat_ass.DOMEGA;
%!       sol_eig(j, k) = fem_sol_modal_damped(mesh, ...
%!                                            dof_map, ...
%!                                            mat_ass, ...
%!                                            cms_opt.modes.number, ...
%!                                            cms_opt);
%!       Mred = mat_ass_cms.Mred;
%!       Dred = mat_ass_cms.Dred;
%!       Kred = mat_ass_cms.Kred;
%!       Dred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.DOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA_DOT(dof_map.idx_node, dof_map.idx_node), "Full");
%!       omegaq = [load_case(j, k).omega.^2;
%!                 load_case(j, k).omega(1) * load_case(j, k).omega(2);
%!                 load_case(j, k).omega(2) * load_case(j, k).omega(3);
%!                 load_case(j, k).omega(3) * load_case(j, k).omega(1)];
%!       idx = int32(0);
%!       for i=1:numel(omegaq)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * omegaq(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).omegadot)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).omegadot(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).g)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred -= mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).g(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).loads)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).loads(i);
%!       endfor
%!       [sol_eig_red(j, k).Ured, sol_eig_red(j, k).lambda_red] = fem_sol_eigsd(Kred, Dred, Mred, cms_opt.modes.number, cms_opt);
%!     endfor
%!   endfor
%!   tol_abs = [0, 0, 0] / SI_unit_second^-1;
%!   tol_rel = [0.3e-2, 4.5e-2, 3e-2];
%!   tol_disp_rel = 4e-2;
%!   err_u_modal = err_v_modal = zeros(size(param));
%!   printf("deformation/velocity:\n");
%!   colors = rainbow(3);
%!   width = 1:size(res, 3);
%!   linestyle = {"-", "--", "-."};
%!   for i=idx_j
%!     for j=idx_k
%!       u_modal = res(i, j, 1).trajectory{end} - res(i, j, 1).trajectory{end}(1, :);
%!       u_beam = res(i, j, 2).trajectory{end} - res(i, j, 2).trajectory{end}(1, :);
%!       v_modal = res(i, j, 1).velocity{end};
%!       v_beam = res(i, j, 2).velocity{end};
%!       if (options.plot)
%!         figure("visible", "off");
%!         hold on;
%!         for k=1:size(res, 3)
%!           for l=1:3
%!             hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l) - res(i, j, k).trajectory{end}(1, l)) * SI_unit_meter);
%!             set(hnd, "color", colors(l, :));
%!             set(hnd, "linewidth", width(k));
%!             set(hnd, "linestyle", linestyle{k});
%!           endfor
%!         endfor
%!         xlabel("t [s]");
%!         ylabel("u [m]");
%!         grid on;
%!         grid minor on;
%!         title(sprintf("linear displacement %d:%d", i, j));
%!         figure("visible", "off");
%!         hold on;
%!         for k=1:size(res, 3)
%!           for l=1:3
%!             hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l + 3) - res(i, j, k).trajectory{end}(1, l + 3)) * 180 / pi);
%!             set(hnd, "color", colors(l, :));
%!             set(hnd, "linewidth", width(k));
%!             set(hnd, "linestyle", linestyle{k});
%!           endfor
%!         endfor
%!         xlabel("t [s]");
%!         ylabel("Phi [deg]");
%!         grid on;
%!         grid minor on;
%!         title(sprintf("angular displacement %d:%d", i, j));
%!       endif
%!       err_u_modal(i, j) = max(max(abs(u_modal - u_beam))) / max(1, max(max(abs(u_beam))));
%!       err_v_modal(i, j) = max(max(abs(v_modal - v_beam))) / max(1, max(max(abs(v_beam))));
%!       printf("%d:%d %.1f%%/%.1f%%\n", i, j, 100 * err_u_modal(i, j), 100 * err_v_modal(i, j));
%!     endfor
%!   endfor
%!   printf("natural frequencies:\n");
%!   MACR = cell(size(param));
%!   result_data = struct("f_mbd", cell(size(param)), "f_fem", cell(size(param)));
%!   for i=idx_j
%!     for j=idx_k
%!       f_fem = sort(sol_eig(i, j).f(:));
%!       f_fem = f_fem(f_fem > 0);
%!       f_mbd = zeros(rows(f_fem), size(res, 3));
%!       PhiR = zeros(6, rows(f_fem), size(res, 3));
%!       for k=1:size(res, 3)
%!         [f_mbd_k, idx_mbd_k] = sort(res(i, j, k).modal.f(:));
%!         D_mbd_k = res(i, j, k).modal.D(idx_mbd_k);
%!         idx_mbd_k = idx_mbd_k(f_mbd_k > 0);
%!         f_mbd_k = f_mbd_k(f_mbd_k > 0);
%!         idx_mbd_k = idx_mbd_k(1:rows(f_fem));
%!         f_mbd(:, k) = f_mbd_k(1:rows(f_fem));
%!         PhiR(:, :, k) = res(i, j, k).modal.VR(res(i, j, k).modal.idx(end) + (1:6), idx_mbd_k);
%!       endfor
%!       result_data(i, j).f_fem = f_fem;
%!       result_data(i, j).f_mbd = f_mbd;
%!       MACR{i, j} = MACL{i, j} = zeros(rows(f_fem), rows(f_fem));
%!       for k=1:rows(f_fem)
%!         for l=1:rows(f_fem)
%!           MACR{i, j}(k, l) = (PhiR(:, k, 1)' * PhiR(:, k, 2)) * conj(PhiR(:, k, 1)' * PhiR(:, k, 2)) / ((PhiR(:, k, 1)' * PhiR(:, k, 1)) * (PhiR(:, k, 2)' * PhiR(:, k, 2)));
%!         endfor
%!       endfor
%!       printf("%d:%d\n", i, j);
%!       for k=1:rows(f_fem)
%!         printf("%10.2f", f_fem(k) * SI_unit_second^-1);
%!         for l=1:columns(f_mbd)
%!           printf("\t%10.2f", f_mbd(k, l) * SI_unit_second^-1);
%!         endfor
%!         for l=1:columns(f_mbd)
%!           printf("\t%.1f%%", 100 * (f_mbd(k, l) / f_fem(k) - 1));
%!         endfor
%!         printf("\t%.3f", MACR{i, j}(k, k));
%!         fputs(stdout, "\n");
%!       endfor
%!       fputs(stdout, "\n\n");
%!     endfor
%!   endfor
%!   for i=idx_j
%!     for j=idx_k
%!       for k=1:rows(result_data(i, j).f_fem)
%!         for l=1:columns(result_data(i, j).f_mbd)
%!           assert_simple(result_data(i, j).f_mbd(k, l), result_data(i, j).f_fem(k), tol_abs(l) + tol_rel(l) * abs(result_data(i, j).f_fem(k)));
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%!   assert_simple(all(all(err_u_modal < tol_disp_rel)));
%!   for j=idx_j
%!     for k=idx_k
%!       tol = 2e-2;
%!       [lambda_s] = sortrows([imag(sol_eig(j, k).lambda(:)), real(sol_eig(j, k).lambda(:))],[1,2]);
%!       [lambda_red_s] = sortrows([imag(sol_eig_red(j, k).lambda_red(:)), real(sol_eig_red(j, k).lambda_red(:))],[1,2]);
%!       K = min(20, rows(lambda_s));
%!       lambda_s = 1j * lambda_s(:,1) + lambda_s(:, 2);
%!       lambda_red_s = 1j * lambda_red_s(:, 1) + lambda_red_s(:, 2);
%!       assert_simple(lambda_red_s(1:K), lambda_s(1:K), tol * norm(lambda_s(1:K)));
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST11
%! pkg load mboct-fem-pkg;
%! printf("fem_cms_create2: test8\n");
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1;
%!   SI_unit_kilogram = 1e3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   geometry.l = 4e-3 / SI_unit_meter;
%!   geometry.w = 6e-3 / SI_unit_meter;
%!   geometry.h = 8e-3 / SI_unit_meter;
%!   h = [geometry.l; geometry.w; geometry.h];
%!   t1 = 1;
%!   dt = t1 / 40;
%!   model = "static";
%!   method = "implicit euler";
%!   options.verbose = false;
%!   material.E = 210000e6 / SI_unit_pascal;
%!   material.ET = 2100e6 / SI_unit_pascal;
%!   material.sigmayv = 235000e6 / SI_unit_pascal;
%!   material.nu = 0.3;
%!   material.G = material.E / (2 * (1 + material.nu));
%!   material.kappa = material.E / (3 * (1 - 2 * material.nu));
%!   material.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   material.beta = 1 / SI_unit_second^-1;
%!   material.delta = 0.25;
%!   material.theta = 0.5;
%!   material.tau = 0.5 / SI_unit_second;
%!   elem_types = {"tet10h", ...
%!                 "tet10upc", ...
%!                 "iso8", ...
%!                 "iso8upc", ...
%!                 "iso20", ...
%!                 "iso20r", ...
%!                 "iso20upc", ...
%!                 "iso20upcr", ...
%!                 "iso27", ...
%!                 "penta15", ...
%!                 "penta15upc", ...
%!                };
%!   mat_types = {"linear elastic generic", ...
%!                "hookean linear elastic isotropic", ...
%!                "neo hookean elastic", ...
%!                "bilinear isotropic hardening", ...
%!                "mooney rivlin elastic", ...
%!                "linear viscoelastic generic", ...
%!                "hookean linear viscoelastic isotropic", ...
%!                "neo hookean viscoelastic", ...
%!                "linear viscoelastic maxwell1", ...
%!                "linear viscoelastic maxwelln", ...
%!               };
%!   f_transfinite_mesh = [true, ...
%!                         false, ...
%!                        ];
%!   boundary_cond = { #"symmetry", ...
%!                     "three point", ...
%!                     #"two surfaces one line", ...
%!                   };
%!   load_type = {"traction", "pressure", "prestrain"};
%!   sigma = material.sigmayv * diag([1.5, 0.3, 0.9, 0.2, 0.2, 0.2]);
%!   epsilon0 = [1; 2; 3; 0.4; 0.5; 0.6];
%!   for idx_sigma=1:columns(sigma)
%!     sigmav = sqrt(sum(sigma(1:3, idx_sigma).^2) - (sigma(1, idx_sigma) * sigma(2, idx_sigma) + sigma(2, idx_sigma) * sigma(3, idx_sigma) + sigma(1, idx_sigma) * sigma(3, idx_sigma)) + 3 * sum(sigma(4:6, idx_sigma).^2));
%!     for idx_load_type=1:numel(load_type)
%!       switch (load_type{idx_load_type})
%!         case "pressure"
%!           if (norm(sigma(1:3, idx_sigma)) == 0)
%!             continue; ## Uniform pressure cannot cause shear stress
%!           endif
%!         case "prestrain"
%!           if (idx_sigma > columns(epsilon0))
%!             continue;
%!           endif
%!       endswitch
%!       for idx_boundary_cond=1:numel(boundary_cond)
%!         switch (boundary_cond{idx_boundary_cond})
%!           case "symmetry"
%!             switch (load_type{idx_load_type})
%!             case "prestrain"
%!               sym_cond = norm(epsilon0(4:6, idx_sigma));
%!             otherwise
%!               sym_cond = norm(sigma(4:6, idx_sigma));
%!             endswitch
%!             if (sym_cond)
%!               continue; ## Deformation cannot be symmetric because of shear stress
%!             endif
%!         endswitch
%!         for idx_transfinite=1:numel(f_transfinite_mesh)
%!           for idx_mat_type=1:numel(mat_types)
%!             material.type = mat_types{idx_mat_type};
%!             for idx_elem_type=1:numel(elem_types)
%!               t2 = 0;
%!               elem_type = elem_types{idx_elem_type};
%!               switch (elem_type)
%!               case {"iso8upc", "iso20upc", "iso20upcr", "penta15upc", "tet10upc"}
%!                 switch (material.type)
%!                 case {"hookean linear elastic isotropic", "mooney rivlin elastic", "bilinear isotropic hardening"}
%!                 otherwise
%!                   ## incompressible version of constitutive law not implemented yet
%!                   continue;
%!                 endswitch
%!               endswitch
%!               switch (load_type{idx_load_type})
%!               case "prestrain"
%!                 switch (mat_types{idx_mat_type})
%!                 case {"linear elastic generic", "hookean linear elastic isotropic", "neo hookean elastic", "mooney rivlin elastic", "linear viscoelastic generic", "hookean linear viscoelastic isotropic", "neo hookean viscoelastic", "bilinear isotropic hardening"}
%!                 otherwise
%!                   ## prestrain is not implemented yet
%!                   continue;
%!                 endswitch
%!                 switch (mat_types{idx_mat_type})
%!                 case {"linear viscoelastic generic", "hookean linear viscoelastic isotropic", "neo hookean viscoelastic"}
%!                   t2 = 1000;
%!                 endswitch
%!                 switch (mat_types{idx_mat_type})
%!                 case "neo hookean viscoelastic"
%!                   continue; ## FIXME: test not passed yet
%!                 endswitch
%!               endswitch
%!               switch (elem_type)
%!                 case {"iso20upcr", "iso20r"}
%!                   if (f_transfinite_mesh(idx_transfinite))
%!                     elem_factor_h = [0.5; 2; 2]; ## avoid hourglass instability
%!                   else
%!                     elem_factor_h = [0.5; 0.5; 0.5];
%!                   endif
%!                 otherwise
%!                   elem_factor_h = [2; 2; 2];
%!               endswitch
%!               switch (mat_types{idx_mat_type})
%!                 case {"neo hookean elastic", "neo hookean viscoelastic", "mooney rivlin elastic", "linear viscoelastic generic", "hookean linear viscoelastic isotropic"}
%!                   switch(boundary_cond{idx_boundary_cond})
%!                   case {"three point", "two surfaces one line"}
%!                     switch (idx_sigma)
%!                     case {4, 5, 6}
%!                       ## shear deformation with those materials and elements not passed yet because the Jacobian may become singular
%!                       switch (elem_type)
%!                       case {"tet10h", "tet10upc", "penta15", "penta15upc"}
%!                         continue;
%!                       otherwise
%!                         if (~f_transfinite_mesh(idx_transfinite))
%!                           continue;
%!                         endif
%!                       endswitch
%!                     endswitch
%!                   endswitch
%!               endswitch
%!               file_prefix = sprintf("%s_%d_%d_%d_%d_%d_%d", filename, idx_sigma, idx_load_type, idx_boundary_cond, idx_transfinite, idx_mat_type, idx_elem_type);
%!               geo_file = [file_prefix, "_gmsh.geo"];
%!               mesh_file = [file_prefix, "_gmsh.msh"];
%!               nodes_file = [file_prefix, "_mbd.nod"];
%!               elem_file = [file_prefix, "_mbd.elm"];
%!               set_file = [file_prefix, "_mbd.set"];
%!               csl_file = [file_prefix, "_mbd.csl"];
%!               control_file = [file_prefix, "_mbd.con"];
%!               initial_value_file = [file_prefix, "_mbd.inv"];
%!               input_file = [file_prefix, "_mbd_inp.mbdyn"];
%!               output_file = [file_prefix, "_mbd_out"];
%!               opt_mbd.output_file = output_file;
%!               if (~options.verbose)
%!                 opt_mbd.logfile = [opt_mbd.output_file, ".stdout"];
%!               endif
%!               opt_mbd.mbdyn_command = "mbdyn -C";
%!               opt_mbd.f_run_mbdyn = true;
%!               switch (elem_type)
%!                 case {"iso8", "iso8upc"}
%!                   mesh_order = 1;
%!                   elem_type_solid = {elem_type};
%!                   elem_type_surf = {"iso4"};
%!                 case "iso20"
%!                   mesh_order = 2;
%!                   elem_type_solid = {"iso20", "penta15"};
%!                   elem_type_surf = {"quad8", "tria6h"};
%!                 case "iso20upc"
%!                   mesh_order = 2;
%!                   elem_type_solid = {elem_type};
%!                   elem_type_surf = {"quad8"};
%!                 case "iso20upcr"
%!                   mesh_order = 2;
%!                   elem_type_solid = {elem_type};
%!                   elem_type_surf = {"quad8r"};
%!                 case "iso27"
%!                   mesh_order = 2;
%!                   elem_type_solid = {"iso27"};
%!                   elem_type_surf = {"quad9"};
%!                 case "iso20r"
%!                   mesh_order = 2;
%!                   elem_type_solid = {"iso20r", "penta15"};
%!                   elem_type_surf = {"quad8r", "tria6h"};
%!                 case {"penta15", "penta15upc"}
%!                   mesh_order = 2;
%!                   elem_type_solid = {elem_type};
%!                   elem_type_surf = {"quad8", "tria6h"};
%!                 case {"tet10h", "tet10upc"}
%!                   mesh_order = 2;
%!                   elem_type_solid = {elem_type};
%!                   elem_type_surf = {"tria6h"};
%!                 otherwise
%!                   error("unknown element type \"%s\"", elem_type);
%!               endswitch
%!               fd = -1;
%!               unwind_protect
%!                 [fd, msg] = fopen(geo_file, "w");
%!                 if (fd == -1)
%!                   error("failed to open file \"%s.geo\"", geo_file);
%!                 endif
%!                 fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!                 fprintf(fd, "a=%g;\n", geometry.l);
%!                 fprintf(fd, "b=%g;\n", geometry.w);
%!                 fprintf(fd, "c=%g;\n", geometry.h);
%!                 for i=1:3
%!                   fprintf(fd, "h%s = %g;\n", {"x","y","z"}{i}, h(i) * elem_factor_h(i));
%!                 endfor
%!                 switch (elem_type)
%!                   case {"iso20", "iso20upc", "iso20upcr", "iso20r", "penta15", "penta15upc"}
%!                     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!                   otherwise
%!                     fputs(fd, "Mesh.SecondOrderIncomplete=0;\n");
%!                 endswitch
%!                 fprintf(fd, "Mesh.ElementOrder = %d;\n", mesh_order);
%!                 fputs(fd, "Point(1) = { 0, -0.5 * b, -0.5 * c};\n");
%!                 fputs(fd, "Point(2) = { a, -0.5 * b, -0.5 * c};\n");
%!                 fputs(fd, "Point(3) = { a,  0.5 * b, -0.5 * c};\n");
%!                 fputs(fd, "Point(4) = { 0,  0.5 * b, -0.5 * c};\n");
%!                 fputs(fd, "Line(1) = {4,3};\n");
%!                 fputs(fd, "Line(2) = {3,2};\n");
%!                 fputs(fd, "Line(3) = {2,1};\n");
%!                 fputs(fd, "Line(4) = {1,4};\n");
%!                 switch (elem_type)
%!                 case {"iso20upcr", "iso20r"}
%!                   num_layers = 2; ## because of hourglass instability
%!                 otherwise
%!                   num_layers = 1;
%!                 endswitch
%!                 fprintf(fd, "num_layers = %d;\n", num_layers);
%!                 if (f_transfinite_mesh(idx_transfinite))
%!                   fprintf(fd, "Transfinite Curve(1) = Max(num_layers + 1, Round(a / hx) + 1);\n");
%!                   fprintf(fd, "Transfinite Curve(2) = Max(num_layers + 1, Round(b / hy) + 1);\n");
%!                   fprintf(fd, "Transfinite Curve(3) = Max(num_layers + 1, Round(a / hx) + 1);\n");
%!                   fprintf(fd, "Transfinite Curve(4) = Max(num_layers + 1, Round(b / hy) + 1);\n");
%!                 endif
%!                 fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!                 fputs(fd, "Plane Surface(6) = {5};\n");
%!                 if (f_transfinite_mesh(idx_transfinite))
%!                   fputs(fd, "Transfinite Surface(6) = {2,3,4,1};\n");
%!                 endif
%!                 fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!                 switch (elem_type)
%!                   case {"tet10h", "tet10upc"}
%!                     fputs(fd, "  Surface{6};\n");
%!                   otherwise
%!                     fprintf(fd, "  Surface{6}; Layers{Max(num_layers, Round(c/hz))}; Recombine;\n");
%!                 endswitch
%!                 fputs(fd, "};\n");
%!                 f_unstruct_mesh_size = false;
%!                 switch (elem_type)
%!                   case {"iso8", "iso8upc", "iso20", "iso20r", "iso20upc", "iso20upcr", "iso27"}
%!                     f_unstruct_mesh_size = ~f_transfinite_mesh(idx_transfinite);
%!                     fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!                   otherwise
%!                     f_unstruct_mesh_size = true;
%!                 endswitch
%!                 if (f_unstruct_mesh_size)
%!                   fprintf(fd, "MeshSize{PointsOf{Volume{tmp[1]};}} = %.16e;\n", 2 * mean(h .* elem_factor_h));
%!                 endif
%!                 switch (elem_type)
%!                   case {"tet10h", "tet10upc"}
%!                     if (~f_transfinite_mesh(idx_transfinite))
%!                       fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!                       fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!                     endif
%!                 endswitch
%!                 fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!                 fputs(fd, "Physical Surface(\"clamp\",2) = {tmp[4]};\n");
%!                 fputs(fd, "Physical Surface(\"load+x\",3) = {tmp[2]};\n");
%!                 fputs(fd, "Physical Surface(\"load+y\",6) = {tmp[5]};\n");
%!                 fputs(fd, "Physical Surface(\"load+z\",7) = {tmp[0]};\n");
%!                 fputs(fd, "Physical Surface(\"load-x\",8) = {tmp[4]};\n");
%!                 fputs(fd, "Physical Surface(\"load-y\",9) = {tmp[3]};\n");
%!                 fputs(fd, "Physical Surface(\"load-z\",10) = {6};\n");
%!                 fputs(fd, "Physical Surface(\"symmetry-xy\",4) = {6};\n");
%!                 fputs(fd, "Physical Surface(\"symmetry-xz\",5) = {tmp[3]};\n");
%!               unwind_protect_cleanup
%!                 if (fd ~= -1)
%!                   fclose(fd);
%!                 endif
%!               end_unwind_protect
%!               pid = spawn("gmsh", {"-format", "msh2", "-3", geo_file});
%!               status = spawn_wait(pid);
%!               if (status ~= 0)
%!                 warning("gmsh failed with status %d", status);
%!               endif
%!               opt_msh.elem_type = {elem_type_solid{:}, elem_type_surf{:}};
%!               mesh = fem_pre_mesh_reorder(fem_pre_mesh_import(mesh_file, "gmsh", opt_msh));
%!               opt_mbd_mesh = struct();
%!               switch (model)
%!                 case "dynamic"
%!                   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!                 case "static"
%!                   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!               endswitch
%!               grp_idx_volume = zeros(1, numel(elem_type_solid), "int32");
%!               for i=1:numel(elem_type_solid)
%!                 if (~isfield(mesh.groups, elem_type_solid{i}))
%!                   continue;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_solid{i}).id] == 1);
%!                 if (isempty(idx))
%!                   continue;
%!                 endif
%!                 grp_idx_volume(i) = idx;
%!               endfor
%!               grp_idx_load_px = grp_idx_load_py = grp_idx_load_pz = grp_idx_clamp = grp_idx_symmetry_xy = grp_idx_symmetry_xz = zeros(1, numel(elem_type_surf), "int32");
%!               grp_idx_load_mx = grp_idx_load_my = grp_idx_load_mz = zeros(1, numel(elem_type_surf), "int32");
%!               for i=1:numel(elem_type_surf)
%!                 if (~isfield(mesh.groups, elem_type_surf{i}))
%!                   continue;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 2);
%!                 if (~isempty(idx))
%!                   grp_idx_clamp(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 3);
%!                 if (~isempty(idx))
%!                   grp_idx_load_px(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 4);
%!                 if (~isempty(idx))
%!                   grp_idx_symmetry_xy(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 5);
%!                 if (~isempty(idx))
%!                   grp_idx_symmetry_xz(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 6);
%!                 if (~isempty(idx))
%!                   grp_idx_load_py(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 7);
%!                 if (~isempty(idx))
%!                   grp_idx_load_pz(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 8);
%!                 if (~isempty(idx))
%!                   grp_idx_load_mx(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 9);
%!                 if (~isempty(idx))
%!                   grp_idx_load_my(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 10);
%!                 if (~isempty(idx))
%!                   grp_idx_load_mz(i) = idx;
%!                 endif
%!               endfor
%!               load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!               grp_idx_p1 = find((mesh.nodes(:, 1) == 0) & (mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 3) == -0.5 * geometry.h));
%!               grp_idx_p2 = find((mesh.nodes(:, 1) == geometry.l) & (mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 3) == -0.5 * geometry.h));
%!               grp_idx_p3 = find((mesh.nodes(:, 1) == geometry.l) & (mesh.nodes(:, 2) == 0.5 * geometry.w) & (mesh.nodes(:, 3) == -0.5 * geometry.h));
%!               grp_idx_p4 = find((mesh.nodes(:, 1) == 0) & (mesh.nodes(:, 2) == 0.5 * geometry.w) & (mesh.nodes(:, 3) == -0.5 * geometry.h));
%!               grp_idx_p5 = find((mesh.nodes(:, 1) == 0) & (mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 3) == 0.5 * geometry.h));
%!               grp_idx_p6 = find((mesh.nodes(:, 1) == geometry.l) & (mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 3) == 0.5 * geometry.h));
%!               grp_idx_p7 = find((mesh.nodes(:, 1) == geometry.l) & (mesh.nodes(:, 2) == 0.5 * geometry.w) & (mesh.nodes(:, 3) == 0.5 * geometry.h));
%!               grp_idx_p8 = find((mesh.nodes(:, 1) == 0) & (mesh.nodes(:, 2) == 0.5 * geometry.w) & (mesh.nodes(:, 3) == 0.5 * geometry.h));
%!               switch (boundary_cond{idx_boundary_cond})
%!                 case "symmetry"
%!                   for i=1:numel(elem_type_surf)
%!                     if (grp_idx_clamp(i))
%!                       load_case_dof.locked_dof(getfield(mesh.groups, elem_type_surf{i})(grp_idx_clamp(i)).nodes, 1) = true;
%!                     endif
%!                     if (grp_idx_symmetry_xz(i))
%!                       load_case_dof.locked_dof(getfield(mesh.groups, elem_type_surf{i})(grp_idx_symmetry_xz(i)).nodes, 2) = true;
%!                     endif
%!                     if (grp_idx_symmetry_xy(i))
%!                       load_case_dof.locked_dof(getfield(mesh.groups, elem_type_surf{i})(grp_idx_symmetry_xy(i)).nodes, 3) = true;
%!                     endif
%!                   endfor
%!                 case "three point"
%!                   load_case_dof.locked_dof(grp_idx_p1, 1:3) = true;
%!                   load_case_dof.locked_dof(grp_idx_p2, 2:3) = true;
%!                   load_case_dof.locked_dof(grp_idx_p4, 3) = true;
%!                 case "two surfaces one line"
%!                     switch (idx_sigma)
%!                     case {1, 2, 3, 4}
%!                       load_case_dof.locked_dof(mesh.nodes(:, 2) == -0.5 * geometry.w, 2) = true;
%!                       load_case_dof.locked_dof(mesh.nodes(:, 3) == -0.5 * geometry.h, 3) = true;
%!                       load_case_dof.locked_dof((mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 1) == 0), 1) = true;
%!                     case 5
%!                       load_case_dof.locked_dof(mesh.nodes(:, 3) == -0.5 * geometry.h, 3) = true;
%!                       load_case_dof.locked_dof(mesh.nodes(:, 1) == 0, 1) = true;
%!                       load_case_dof.locked_dof((mesh.nodes(:, 3) == -0.5 * geometry.h) & (mesh.nodes(:, 2) == -0.5 * geometry.w), 2) = true;
%!                     case 6
%!                       load_case_dof.locked_dof(mesh.nodes(:, 3) == -0.5 * geometry.h, 3) = true;
%!                       load_case_dof.locked_dof(mesh.nodes(:, 2) == -0.5 * geometry.w, 2) = true;
%!                       load_case_dof.locked_dof((mesh.nodes(:, 3) == -0.5 * geometry.h) & (mesh.nodes(:, 1) == 0), 1) = true;
%!                     otherwise
%!                       error("invalid boundary condition");
%!                     endswitch
%!                 otherwise
%!                   error("unkown boundary condition");
%!               endswitch
%!               mesh.material_data = material;
%!               switch (mesh.material_data.type)
%!               case "linear viscoelastic maxwelln"
%!                 mesh.material_data.tau = repmat(mesh.material_data.tau, 1, 2);
%!                 mesh.material_data.theta = repmat(mesh.material_data.theta / 2, 1, 2);
%!               endswitch
%!               mesh.materials = struct();
%!               for i=1:numel(elem_type_solid)
%!                 if (~isfield(mesh.elements, elem_type_solid{i}))
%!                   continue;
%!                 endif
%!                 elem_mat = zeros(rows(getfield(mesh.elements, elem_type_solid{i})), 1, "int32");
%!                 elem_mat(getfield(mesh.groups, elem_type_solid{i})(grp_idx_volume(i)).elements) = 1;
%!                 mesh.materials = setfield(mesh.materials, elem_type_solid{i}, elem_mat);
%!               endfor
%!               opt_mbd_mesh.forces.time_function = "time";
%!               opt_mbd_mesh.surface_loads.time_function = opt_mbd_mesh.forces.time_function;
%!               load_case.pressure = struct();
%!               load_case.traction = struct();
%!               load_case.traction_abs = struct();
%!               for i=1:numel(grp_idx_load_px)
%!                 if (~isfield(mesh.elements, elem_type_surf{i}))
%!                   continue;
%!                 endif
%!                 if (grp_idx_load_px(i) > 0)
%!                   elem_px = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_px(i)).elements;
%!                 else
%!                   elem_px = [];
%!                 endif
%!                 if (grp_idx_load_py(i) > 0)
%!                   elem_py = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_py(i)).elements;
%!                 else
%!                   elem_py = [];
%!                 endif
%!                 if (grp_idx_load_pz(i) > 0)
%!                   elem_pz = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_pz(i)).elements;
%!                 else
%!                   elem_pz = [];
%!                 endif
%!                 if (grp_idx_load_mx(i) > 0)
%!                   elem_mx = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_mx(i)).elements;
%!                 else
%!                   elem_mx = [];
%!                 endif
%!                 if (grp_idx_load_my(i) > 0)
%!                   elem_my = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_my(i)).elements;
%!                 else
%!                   elem_my = [];
%!                 endif
%!                 if (grp_idx_load_mz(i) > 0)
%!                   elem_mz = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_mz(i)).elements;
%!                 else
%!                   elem_mz = [];
%!                 endif
%!                 elem_nodes = [getfield(mesh.elements, elem_type_surf{i})(elem_px, :);
%!                               getfield(mesh.elements, elem_type_surf{i})(elem_py, :);
%!                               getfield(mesh.elements, elem_type_surf{i})(elem_pz, :);
%!                               getfield(mesh.elements, elem_type_surf{i})(elem_mx, :);
%!                               getfield(mesh.elements, elem_type_surf{i})(elem_my, :);
%!                               getfield(mesh.elements, elem_type_surf{i})(elem_mz, :)];
%!                 switch (load_type{idx_load_type})
%!                   case "pressure"
%!                     elem_press = [repmat(-sigma(1, idx_sigma), numel(elem_px), columns(elem_nodes));
%!                                   repmat(-sigma(2, idx_sigma), numel(elem_py), columns(elem_nodes));
%!                                   repmat(-sigma(3, idx_sigma), numel(elem_pz), columns(elem_nodes));
%!                                   repmat(-sigma(1, idx_sigma), numel(elem_mx), columns(elem_nodes));
%!                                   repmat(-sigma(2, idx_sigma), numel(elem_my), columns(elem_nodes));
%!                                   repmat(-sigma(3, idx_sigma), numel(elem_mz), columns(elem_nodes))];
%!                     load_case.pressure = setfield(load_case.pressure, ...
%!                                                   elem_type_surf{i}, ...
%!                                                   struct("elements", elem_nodes, ...
%!                                                          "p", elem_press));
%!                   case {"traction", "traction_abs"}
%!                     elem_trac = zeros(rows(elem_nodes), columns(elem_nodes), 3);
%!                     ioffset = 0;               elem_trac(ioffset + (1:numel(elem_px)), :, 1) += sigma(1, idx_sigma);
%!                     ioffset += numel(elem_px); elem_trac(ioffset + (1:numel(elem_py)), :, 2) += sigma(2, idx_sigma);
%!                     ioffset += numel(elem_py); elem_trac(ioffset + (1:numel(elem_pz)), :, 3) += sigma(3, idx_sigma);
%!                     ioffset += numel(elem_pz); elem_trac(ioffset + (1:numel(elem_mx)), :, 1) -= sigma(1, idx_sigma);
%!                     ioffset += numel(elem_mx); elem_trac(ioffset + (1:numel(elem_my)), :, 2) -= sigma(2, idx_sigma);
%!                     ioffset += numel(elem_my); elem_trac(ioffset + (1:numel(elem_mz)), :, 3) -= sigma(3, idx_sigma);
%!                     ioffset += numel(elem_mz);
%!                     ioffset = 0;               elem_trac(ioffset + (1:numel(elem_px)), :, 3) += sigma(6, idx_sigma);
%!                     ioffset += numel(elem_px); elem_trac(ioffset + (1:numel(elem_py)), :, 2) += 0;
%!                     ioffset += numel(elem_py); elem_trac(ioffset + (1:numel(elem_pz)), :, 1) += sigma(6, idx_sigma);
%!                     ioffset += numel(elem_pz); elem_trac(ioffset + (1:numel(elem_mx)), :, 3) -= sigma(6, idx_sigma);
%!                     ioffset += numel(elem_mx); elem_trac(ioffset + (1:numel(elem_my)), :, 2) -= 0;
%!                     ioffset += numel(elem_my); elem_trac(ioffset + (1:numel(elem_mz)), :, 1) -= sigma(6, idx_sigma);
%!                     ioffset += numel(elem_mz);
%!                     ioffset = 0;               elem_trac(ioffset + (1:numel(elem_px)), :, 2) += sigma(4, idx_sigma);
%!                     ioffset += numel(elem_px); elem_trac(ioffset + (1:numel(elem_py)), :, 1) += sigma(4, idx_sigma);
%!                     ioffset += numel(elem_py); elem_trac(ioffset + (1:numel(elem_pz)), :, 1) += 0;
%!                     ioffset += numel(elem_pz); elem_trac(ioffset + (1:numel(elem_mx)), :, 2) -= sigma(4, idx_sigma);
%!                     ioffset += numel(elem_mx); elem_trac(ioffset + (1:numel(elem_my)), :, 1) -= sigma(4, idx_sigma);
%!                     ioffset += numel(elem_my); elem_trac(ioffset + (1:numel(elem_mz)), :, 1) -= 0;
%!                     ioffset += numel(elem_mz);
%!                     ioffset = 0;               elem_trac(ioffset + (1:numel(elem_px)), :, 3) += 0;
%!                     ioffset += numel(elem_px); elem_trac(ioffset + (1:numel(elem_py)), :, 3) += sigma(5, idx_sigma);
%!                     ioffset += numel(elem_py); elem_trac(ioffset + (1:numel(elem_pz)), :, 2) += sigma(5, idx_sigma);
%!                     ioffset += numel(elem_pz); elem_trac(ioffset + (1:numel(elem_mx)), :, 2) -= 0;
%!                     ioffset += numel(elem_mx); elem_trac(ioffset + (1:numel(elem_my)), :, 3) -= sigma(5, idx_sigma);
%!                     ioffset += numel(elem_my); elem_trac(ioffset + (1:numel(elem_mz)), :, 2) -= sigma(5, idx_sigma);
%!                     ioffset += numel(elem_mz);
%!                     load_case = setfield(load_case, ...
%!                                          load_type{idx_load_type}, ...
%!                                          setfield(getfield(load_case, load_type{idx_load_type}), ...
%!                                                   elem_type_surf{i}, ...
%!                                                   struct("elements", elem_nodes, ...
%!                                                          "f", elem_trac)));
%!                 case "prestrain"
%!                     mesh.material_data.extra_data = ", prestrain, reference, tpl_drive_id_epsilon0";
%!                 endswitch
%!               endfor
%!               opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!               opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!               opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!               idx_joint = int32(0);
%!               unwind_protect
%!                 [fd, msg] = fopen(set_file, "wt");
%!                 if (fd == -1)
%!                   error("failed to open file \"%s\": %s", set_file, msg);
%!                 endif
%!                 fprintf(fd, "set: integer number_of_nodes = %d;\n", opt_mbd_mesh.struct_nodes.number);
%!                 fprintf(fd, "set: integer number_of_nodes_hydraulic = %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!                 fprintf(fd, "set: integer number_of_solids = %d;\n", opt_mbd_mesh.solids.number);
%!                 fprintf(fd, "set: integer number_of_genels = %d;\n", opt_mbd_mesh.genels.number);
%!                 fprintf(fd, "set: integer number_of_forces = %d;\n", opt_mbd_mesh.forces.number);
%!                 fprintf(fd, "set: integer number_of_joints = %d;\n", idx_joint);
%!                 fprintf(fd, "set: integer number_of_surface_loads = %d;\n", opt_mbd_mesh.surface_loads.number);
%!                 fprintf(fd, "set: real t1 = %.16e;\n", t1);
%!                 fprintf(fd, "set: real t2 = %.16e;\n", t2);
%!                 fprintf(fd, "set: real dt = %.16e;\n", dt);
%!                 switch (load_type{idx_load_type})
%!                 case "prestrain"
%!                   fprintf(fd, "set: integer tpl_drive_id_epsilon0 = 2001;\n");
%!                   for i=1:6
%!                     fprintf(fd, "set: real epsilon0_%d = %.16e;\n", i, epsilon0(i, idx_sigma));
%!                   endfor
%!                 endswitch
%!               unwind_protect_cleanup
%!                 if (fd ~= -1)
%!                   fclose(fd);
%!                 endif
%!                 fd = -1;
%!               end_unwind_protect
%!               unwind_protect
%!                 [fd, msg] = fopen(control_file, "wt");
%!                 if (fd == -1)
%!                   error("failed to open file \"%s\": %s", control_file, msg);
%!                 endif
%!                 switch (model)
%!                   case "static"
%!                     fprintf(fd, "model: %s;\n", model);
%!                   case "dynamic"
%!                     fprintf(fd, "# model: %s;\n", model);
%!                 endswitch
%!               unwind_protect_cleanup
%!                 if (fd ~= -1)
%!                   fclose(fd);
%!                 endif
%!                 fd = -1;
%!               end_unwind_protect
%!               unwind_protect
%!                 [fd, msg] = fopen(initial_value_file, "wt");
%!                 if (fd == -1)
%!                   error("failed to open file \"%s\": %s", initial_value_file, msg);
%!                 endif
%!                 fprintf(fd, "method: %s;\n", method);
%!               unwind_protect_cleanup
%!                 if (fd ~= -1)
%!                   fclose(fd);
%!                 endif
%!                 fd = -1;
%!               end_unwind_protect
%!               fd = -1;
%!               unwind_protect
%!                 [fd, msg] = fopen(input_file, "wt");
%!                 if (fd == -1)
%!                   error("failed to open file \"%s\": %s", input_file, msg);
%!                 endif
%!                 fprintf(fd, "include: \"%s\";\n", set_file);
%!                 fprintf(fd, "begin: data;\n");
%!                 fprintf(fd, "        problem: initial value; # the default\n");
%!                 fprintf(fd, "end: data;\n");
%!                 fprintf(fd, "begin: initial value;\n");
%!                 fprintf(fd, "        initial time: 0;\n");
%!                 fprintf(fd, "        final time: t1 + t2;\n");
%!                 fprintf(fd, "        time step: dt;\n");
%!                 if (t2 > 0)
%!                   fprintf(fd, "        strategy: change, piecewise linear, 4, 0., dt, t1, dt, t1 + 0.1 * (t2 - t1), (t2 - t1) / t1 * dt, t2, (t2 - t1) / t1 * dt;\n");
%!                 endif
%!                 fprintf(fd, "        max iterations: 100;\n");
%!                 tolerance_type = "sepnorm";
%!                 switch (elem_type)
%!                 case {"tet10upc", "iso8upc", "iso20upc", "iso20upcr", "penta15upc"}
%!                   ## FIXME: test "sepnorm" does not work well with u/p-c elements; use norm instead
%!                   tolerance_type = "norm";
%!                 endswitch
%!                 switch (load_type{idx_load_type})
%!                 case "prestrain"
%!                   ## because there is no load
%!                   tolerance_type = "norm";
%!                 endswitch
%!                 switch (tolerance_type)
%!                 case "norm"
%!                   fprintf(fd, "        tolerance: 1e-6, test, norm, 1e-12, test, norm;\n");
%!                 case "sepnorm"
%!                   fprintf(fd, "        tolerance: 1e-12, test, sepnorm, 1e-12, test, norm;\n");
%!                 endswitch
%!                 fprintf(fd, "        output: messages;\n");
%!                 fprintf(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!                 fprintf(fd, "        nonlinear solver: nox,\n");
%!                 fprintf(fd, "                          modified, 30,\n");
%!                 fprintf(fd, "                          keep jacobian matrix,\n");
%!                 fprintf(fd, "                          use preconditioner as solver, no,\n");
%!                 fprintf(fd, "                          linesearch method, backtrack,\n");
%!                 fprintf(fd, "                          direction, newton,\n");
%!                 fprintf(fd, "                          jacobian operator, newton krylov,\n");
%!                 fprintf(fd, "                          forcing term, constant,\n");
%!                 fprintf(fd, "                          linear solver tolerance, 1e-12,\n");
%!                 fprintf(fd, "                          inner iterations before assembly, 15,\n");
%!                 fprintf(fd, "                          linear solver max iterations, 300,\n");
%!                 fprintf(fd, "                          krylov subspace size, 300,\n");
%!                 fprintf(fd, "                          minimum step, 1e-12,\n");
%!                 fprintf(fd, "                          recovery step type, constant,\n");
%!                 fprintf(fd, "                          recovery step, 1e-12,\n");
%!                 fprintf(fd, "                          verbose, 3,\n");
%!                 fprintf(fd, "                          print convergence info, no;\n");
%!                 fprintf(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!                 fprintf(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!                 fprintf(fd, "        derivatives tolerance: 1e-5, 1e-5;\n");
%!                 fprintf(fd, "        derivatives max iterations: 10;\n");
%!                 fprintf(fd, "        threads: assembly, 1;\n");
%!                 fprintf(fd, "        threads: solver, 1;\n");
%!                 fprintf(fd, "        output: cpu time;\n");
%!                 fprintf(fd, "        include: \"%s\";\n", initial_value_file);
%!                 fprintf(fd, "end: initial value;\n");
%!                 fprintf(fd, "begin: control data;\n");
%!                 fprintf(fd, "       output meter: closest next, t1 + t2, forever, dt;\n");
%!                 fprintf(fd, "       skip initial joint assembly;\n");
%!                 fprintf(fd, "       output precision: 16;\n");
%!                 fprintf(fd, "       include: \"%s\";\n", control_file);
%!                 fprintf(fd, "       default output: all, solids, accelerations;\n");
%!                 fprintf(fd, "       structural nodes: number_of_nodes;\n");
%!                 fprintf(fd, "       hydraulic nodes: number_of_nodes_hydraulic;\n");
%!                 fprintf(fd, "       solids: number_of_solids;\n");
%!                 fprintf(fd, "       genels: number_of_genels;\n");
%!                 fprintf(fd, "       forces: number_of_forces;\n");
%!                 fprintf(fd, "       surface loads: number_of_surface_loads;\n");
%!                 fprintf(fd, "       joints: number_of_joints;\n");
%!                 fprintf(fd, "       use automatic differentiation;\n");
%!                 fprintf(fd, "end: control data;\n");
%!                 switch (load_type{idx_load_type})
%!                 case "prestrain"
%!                    fprintf(fd, "template drive caller: tpl_drive_id_epsilon0, 6,\n");
%!                    fprintf(fd, " green lagrange strain, single,\n");
%!                    for i=1:6
%!                      fprintf(fd, "  epsilon0_%d,\n", i);
%!                    endfor
%!                    fprintf(fd, "  min, 2, const, 1., mult, time, const, 1. / t1;\n");
%!                 endswitch
%!                 fprintf(fd, "include: \"%s\";\n", csl_file);
%!                 fprintf(fd, "begin: nodes;\n");
%!                 fprintf(fd, "       include: \"%s\";\n", nodes_file);
%!                 fprintf(fd, "end: nodes;\n");
%!                 fprintf(fd, "begin: elements;\n");
%!                 fprintf(fd, "       include: \"%s\";\n", elem_file);
%!                 fprintf(fd, "end: elements;\n");
%!               unwind_protect_cleanup
%!                 if (fd ~= -1)
%!                   fclose(fd);
%!                 endif
%!                 fd = -1;
%!               end_unwind_protect
%!               if (options.verbose)
%!                 shell(sprintf("cat \"%s\" | nl", set_file));
%!                 shell(sprintf("cat \"%s\" | nl", input_file));
%!                 shell(sprintf("cat \"%s\" | nl", nodes_file));
%!                 shell(sprintf("cat \"%s\" | nl", csl_file));
%!                 shell(sprintf("cat \"%s\" | nl", elem_file));
%!               endif
%!               fprintf(stderr, "element type: %s\n", elem_type);
%!               info = mbdyn_solver_run(input_file, opt_mbd);
%!               [mesh_sol, sol] = mbdyn_post_load_output_sol(output_file);
%!               [genel_id, genel_data] = mbdyn_post_load_output([output_file, ".gen"], 1, [], numel(sol.t), 1);
%!               switch (load_type{idx_load_type})
%!               case {"pressure", "traction", "traction_abs"}
%!                 [surfl_id, surfl_data] = mbdyn_post_load_output([output_file, ".prl"], 3, [], numel(sol.t), 3);
%!               endswitch
%!               F11 = (sol.def(grp_idx_p2, 1, :)(:) - sol.def(grp_idx_p1, 1, :)(:)) ./ (mesh.nodes(grp_idx_p2, 1, :)(:) - mesh.nodes(grp_idx_p1, 1, end)(:)) + 1;
%!               F22 = (sol.def(grp_idx_p4, 2, :)(:) - sol.def(grp_idx_p1, 2, :)(:)) ./ (mesh.nodes(grp_idx_p4, 2, :)(:) - mesh.nodes(grp_idx_p1, 2, end)(:)) + 1;
%!               F33 = (sol.def(grp_idx_p5, 3, :)(:) - sol.def(grp_idx_p1, 3, :)(:)) ./ (mesh.nodes(grp_idx_p5, 3, :)(:) - mesh.nodes(grp_idx_p1, 3, end)(:)) + 1;
%!               F12 = (sol.def(grp_idx_p4, 1, :)(:) - sol.def(grp_idx_p1, 1, :)(:)) ./ (mesh.nodes(grp_idx_p4, 2, :)(:) - mesh.nodes(grp_idx_p1, 2, end)(:));
%!               F31 = (sol.def(grp_idx_p2, 3, :)(:) - sol.def(grp_idx_p1, 3, :)(:)) ./ (mesh.nodes(grp_idx_p2, 1, :)(:) - mesh.nodes(grp_idx_p1, 1, end)(:));
%!               F32 = (sol.def(grp_idx_p4, 3, :)(:) - sol.def(grp_idx_p1, 3, :)(:)) ./ (mesh.nodes(grp_idx_p4, 2, :)(:) - mesh.nodes(grp_idx_p1, 2, end)(:));
%!               F21 = (sol.def(grp_idx_p2, 2, :)(:) - sol.def(grp_idx_p1, 2, :)(:)) ./ (mesh.nodes(grp_idx_p2, 1, :)(:) - mesh.nodes(grp_idx_p1, 1, end)(:));
%!               F13 = (sol.def(grp_idx_p5, 1, :)(:) - sol.def(grp_idx_p1, 1, :)(:)) ./ (mesh.nodes(grp_idx_p5, 3, :)(:) - mesh.nodes(grp_idx_p1, 3, end)(:));
%!               F23 = (sol.def(grp_idx_p5, 2, :)(:) - sol.def(grp_idx_p1, 2, :)(:)) ./ (mesh.nodes(grp_idx_p5, 3, :)(:) - mesh.nodes(grp_idx_p1, 3, end)(:));
%!               F = zeros(3, 3, numel(sol.t));
%!               for i=1:numel(sol.t)
%!                 F(1, 1, :) = F11;
%!                 F(1, 2, :) = F12;
%!                 F(1, 3, :) = F13;
%!                 F(2, 1, :) = F21;
%!                 F(2, 2, :) = F22;
%!                 F(2, 3, :) = F23;
%!                 F(3, 1, :) = F31;
%!                 F(3, 2, :) = F32;
%!                 F(3, 3, :) = F33;
%!               endfor
%!               G = C = zeros(size(F));
%!               for i=1:numel(sol.t)
%!                 G(:, :, i) = 0.5 * (F(:, :, i).' * F(:, :, i) - eye(3));
%!                 C(:, :, i) = F(:, :, i).' * F(:, :, i);
%!               endfor
%!               Epsilon = epsilon = zeros(6, size(G, 3));
%!               for i=1:3
%!                 epsilon(i, :) = sqrt(1 + 2 * G(i, i, :)(:).') - 1;
%!                 Epsilon(i, :) = G(i, i, :)(:).';
%!               endfor
%!               epsilon(4, :) = 2 * G(1, 2, :)(:).' ./ ((1 + epsilon(1, :)) .* (1 + epsilon(2, :)));
%!               epsilon(5, :) = 2 * G(2, 3, :)(:).' ./ ((1 + epsilon(2, :)) .* (1 + epsilon(3, :)));
%!               epsilon(6, :) = 2 * G(3, 1, :)(:).' ./ ((1 + epsilon(3, :)) .* (1 + epsilon(1, :)));
%!               Epsilon(4, :) = 2 * G(1, 2, :)(:).';
%!               Epsilon(5, :) = 2 * G(2, 3, :)(:).';
%!               Epsilon(6, :) = 2 * G(3, 1, :)(:).';
%!               tau_ref = sigma(:, idx_sigma);
%!               sin_gamma = zeros(3, 1);
%!               sin_gamma_cnt = 0;
%!               for i=1:numel(elem_type_solid)
%!                 if (~isfield(sol.strain.epsilon, elem_type_solid{i}))
%!                   continue;
%!                 endif
%!                 sin_gamma += mean(mean(getfield(sol.strain.epsilon, elem_type_solid{i})(:, :, 4:6, end), 1), 2)(:);
%!                 ++sin_gamma_cnt;
%!               endfor
%!               sin_gamma /= sin_gamma_cnt;
%!               cos_gamma = sqrt(1 - sin_gamma.^2);
%!               tan_gamma = sin_gamma ./ cos_gamma;
%!               tau_ref(1) += 2 * sigma(6, idx_sigma) * tan_gamma(3) + 2 * sigma(4, idx_sigma) * tan_gamma(1);
%!               tau_ref(2) += 2 * sigma(5, idx_sigma) * tan_gamma(2);
%!               tol_epsilon = 1e-8;
%!               tol_sigma = 1e-10;
%!               tol = 1e-9;
%!               tol_F = 1e-8;
%!               switch (load_type{idx_load_type})
%!               case "prestrain"
%!                 for j=1:numel(elem_type_solid)
%!                   if (~isfield(sol.strain.epsilon, elem_type_solid{j}))
%!                     continue;
%!                   endif
%!                   epsilon_res = getfield(sol.strain.epsilon, elem_type_solid{j});
%!                   sigma_res = getfield(sol.stress.tau, elem_type_solid{j});
%!                   assert_simple(max(max(max(max(abs(sigma_res))))) < tol_sigma * material.E);
%!                   for i=1:numel(sol.t)
%!                     for k=1:size(epsilon_res, 1)
%!                       for l=1:size(epsilon_res, 2)
%!                         assert_simple(epsilon_res(k, l, :, i)(:), epsilon0(:, idx_sigma), tol_epsilon * norm(epsilon0(:, idx_sigma)));
%!                       endfor
%!                     endfor
%!                   endfor
%!                 endfor
%!                 continue;
%!               endswitch
%!               Fsurfl_sum = zeros(3, 1);
%!               for i=1:numel(surfl_data)
%!                 Fsurfl_sum += surfl_data{i}(end, :)(:);
%!               endfor
%!               for i=1:numel(elem_type_solid)
%!                 if (~isfield(mesh.elements, elem_type_solid{i}))
%!                   continue;
%!                 endif
%!                 tau_res = getfield(sol.stress.tau, elem_type_solid{i});
%!                 for j=1:size(tau_res, 1)
%!                   for k=1:size(tau_res, 2)
%!                     assert_simple(tau_res(j, k, :, end)(:), tau_ref, tol * norm(tau_ref));
%!                   endfor
%!                 endfor
%!               endfor
%!               Fref = max(abs([geometry.w * geometry.h * tau_ref(1), ...
%!                               geometry.l * geometry.h * tau_ref(2), ...
%!                               geometry.l * geometry.w * tau_ref(3)]));
%!               for i=1:numel(genel_data)
%!                 assert_simple(all(all(abs(genel_data{i}) < tol_F * Fref)));
%!               endfor
%!               assert_simple(norm(Fsurfl_sum) < tol_F * Fref);
%!               for j=1:numel(elem_type_solid)
%!                 if (~isfield(sol.strain.epsilon, elem_type_solid{j}))
%!                   continue;
%!                 endif
%!                 epsilon_res = getfield(sol.strain.epsilon, elem_type_solid{j});
%!                 for i=1:numel(sol.t)
%!                   for k=1:size(epsilon_res, 1)
%!                     for l=1:size(epsilon_res, 2)
%!                       assert_simple(epsilon_res(k, l, :, i)(:), epsilon(:, i), tol_epsilon);
%!                     endfor
%!                   endfor
%!                 endfor
%!               endfor
%!               S_res = zeros(3, 3, numel(sol.t));
%!               switch (mesh.material_data.type)
%!                 case {"linear elastic generic", "hookean linear elastic isotropic", "bilinear isotropic hardening"}
%!                   switch (mesh.material_data.type)
%!                     case {"bilinear isotropic hardening"}
%!                       if (sigmav > mesh.material_data.sigmayv)
%!                         continue;
%!                       endif
%!                   endswitch
%!                   H = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!                   sigma_res = H * Epsilon;
%!                   for i=1:3
%!                     S_res(i, i, :) = sigma_res(i, :);
%!                   endfor
%!                   S_res(1, 2, :) = S_res(2, 1, :) = sigma_res(4, :);
%!                   S_res(2, 3, :) = S_res(3, 2, :) = sigma_res(5, :);
%!                   S_res(3, 1, :) = S_res(1, 3, :) = sigma_res(6, :);
%!                 case {"neo hookean elastic", "mooney rivlin elastic"}
%!                   mu = mesh.material_data.E / (2 * (1 + mesh.material_data.nu));
%!                   lambda = mesh.material_data.E * mesh.material_data.nu / ((1 + mesh.material_data.nu ) * (1 - 2 * mesh.material_data.nu));
%!                   for i=1:numel(sol.t)
%!                     IC = trace(C(:, :, i));
%!                     IIC = 1/2 * (trace(C(:, :, i))^2 - trace(C(:, :, i)^2));
%!                     IIIC = det(C(:, :, i));
%!                     invC = inv(C(:, :, i));
%!                     for k=1:3
%!                       for l=1:3
%!                         switch (mesh.material_data.type)
%!                         case "neo hookean elastic"
%!                           S_res(k, l, i) = mu * (k == l) + (lambda * (IIIC - sqrt(IIIC)) - mu) * invC(k, l);
%!                         case {"mooney rivlin elastic"}
%!                           C1 = mesh.material_data.G / (2 * (1 + mesh.material_data.delta));
%!                           C2 = mesh.material_data.delta * C1;
%!                           S_res(k, l, i) = 2 * (C1 * IIIC^(-1/3) * (k == l) + C2 * IIIC^(-2/3) * (IC * (k == l) - C(k, l, i)) + (1/2 * mesh.material_data.kappa * (IIIC - sqrt(IIIC)) - 1/3 * C1 * IC * IIIC^(-1/3) - 2/3 * C2 * IIC * IIIC^(-2/3)) * invC(k, l));
%!                         endswitch
%!                       endfor
%!                     endfor
%!                   endfor
%!               endswitch
%!               switch (mesh.material_data.type)
%!               case {"linear viscoelastic generic", "hookean linear viscoelastic isotropic", "neo hookean viscoelastic", "linear viscoelastic maxwelln", "linear viscoelastic maxwell1"}
%!                 ## TODO: viscoelastic case is not handled yet
%!               otherwise
%!                 tau_res = zeros(3, 3, numel(sol.t));
%!                 for i=1:numel(sol.t)
%!                   tau_res(:, :, i) = F(:, :, i) * S_res(:, :, i) * F(:, :, i).' / det(F(:, :, i));
%!                 endfor
%!                 Tau_res = zeros(6, numel(sol.t));
%!                 for i=1:3
%!                   Tau_res(i, :) = tau_res(i, i, :);
%!                 endfor
%!                 Tau_res(4, :) = tau_res(1, 2, :);
%!                 Tau_res(5, :) = tau_res(2, 3, :);
%!                 Tau_res(6, :) = tau_res(3, 1, :);
%!                 assert_simple(Tau_res(:, end), tau_ref, tol * norm(tau_ref));
%!               endswitch
%!             endfor
%!           endfor
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST12
%! pkg load mboct-fem-pkg;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   a = 800e-3 / SI_unit_meter;
%!   b = 40e-3 / SI_unit_meter;
%!   c = 10e-3 / SI_unit_meter;
%!   d = 0e-3 / SI_unit_meter;
%!   h = c;
%!   options.interactive = false;
%!   options.plot = true;
%!   options.verbose = false;
%!   options.number_of_beams = int32(40);
%!   options.number_of_threads = mbdyn_solver_num_threads_default();
%!   if (options.plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fprintf(fd, "a = %.16e;\n", a);
%!     fprintf(fd, "b = %.16e;\n", b);
%!     fprintf(fd, "c = %.16e;\n", c);
%!     fprintf(fd, "h = %.16e;\n", h);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(2) = {0,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(3) = {a,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(4) = {a, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Line Loop(1) = {1, 2, 3, 4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,c}{Surface{1}; Layers{Max(1, Round(c / h))}; Recombine;};\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {v1[2]};\n");
%!     fputs(fd, "Physical Surface(2) = {v1[4]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   mesh.material_data.E = 70000e6 / SI_unit_pascal;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 2700 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data.alpha = 0e-5 / (1 / SI_unit_second);
%!   mesh.material_data.beta = 0e-5 / (SI_unit_second);
%!   mesh.materials.iso27 = zeros(rows(mesh.elements.iso27), 1, "int32");
%!   grp_idx_beam = find([[mesh.groups.iso27].id] == 1);
%!   grp_idx_clamp = find([[mesh.groups.quad9].id] == 1);
%!   mesh.materials.iso27(mesh.groups.iso27(grp_idx_beam).elements) = 1;
%!   cms_opt.number_of_threads = options.number_of_threads;
%!   cms_opt.algorithm = "diag-shift-invert";
%!   cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!   cms_opt.nodes.modal.name = "node_id_modal";
%!   cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!   cms_opt.nodes.interfaces.name = "node_id_interface1";
%!   mesh.nodes(cms_opt.nodes.modal.number, 1:3) = [0, 0, 0];
%!   mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = [a + d, 0, 0];
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, 2, cms_opt.nodes.interfaces.number, "quad9");
%!   cms_opt.refine_max_iter = 30;
%!   cms_opt.pre_scaling = false;
%!   cms_opt.solver = "pardiso";
%!   cms_opt.modes.number = 20;
%!   cms_opt.tolerance_tau = -1;
%!   cms_opt.element.name = "elem_id_modal";
%!   cms_opt.create_binary = true;
%!   cms_opt.use_binary = true;
%!   cms_opt.update_binary = true;
%!   cms_opt.invariants = true;
%!   #cms_opt.enable_KTAU0WP = [true; false(2, 1)];
%!   #cms_opt.enable_KTAU0VP = [true; false(2, 1)];
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, :) = true;
%!   load_case_dof.locked_dof(mesh.groups.quad9(grp_idx_clamp).nodes, :) = true;
%!   [mesh_cms, mat_ass_cms, dof_map_cms, sol_eig_cms, cms_opt, sol_tau_cms] = fem_cms_create2(mesh, load_case_dof, cms_opt);
%!   fem_cms_export(filename, mesh_cms, dof_map_cms, mat_ass_cms, cms_opt);
%!   pert.omega = [1e2; 3e2; 2e2] / (SI_unit_rad / SI_unit_second);
%!   pert.omegadot = [1e5; 1e3; 3e3] / (SI_unit_rad / SI_unit_second^2);
%!   pert.loads = [[1e4; 1e3; 1e2] / (SI_unit_newton);
%!                 [1e2; 1e1; 1e1] / (SI_unit_newton * SI_unit_meter)];
%!   pert.g = [1e4; -1e3; -1e2] / (SI_unit_meter / SI_unit_second^2);
%!   pert.a = [-1e4; 1e3; 1e2] / (SI_unit_meter / SI_unit_second^2);
%!   empty_cell = cell(7, 3, 2);
%!   res = struct("info", empty_cell, ...
%!                "t", empty_cell, ...
%!                "trajectory", empty_cell, ...
%!                "deformation", empty_cell, ...
%!                "velocity", empty_cell, ...
%!                "acceleration", empty_cell, ...
%!                "node_id", empty_cell, ...
%!                "force", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "orientation_description", empty_cell, ...
%!                "drive_id", empty_cell, ...
%!                "drive_value", empty_cell, ...
%!                "modal", empty_cell);
%!   empty_cell = cell(7, 3);
%!   param = struct("omega", empty_cell, ...
%!                  "omegadot", empty_cell, ...
%!                  "F1", empty_cell, ...
%!                  "M1", empty_cell, ...
%!                  "a", empty_cell, ...
%!                  "g", empty_cell, ...
%!                  "t1", empty_cell, ...
%!                  "holonomic", empty_cell);
%!   idx_j = 1:rows(param);
%!   idx_k = 1:columns(param);
%!   for j=idx_j
%!     for k=idx_k
%!       param(j, k).omega = zeros(3, 1);
%!       param(j, k).omegadot = zeros(3, 1);
%!       param(j, k).F1 = zeros(3, 1);
%!       param(j, k).M1 = zeros(3, 1);
%!       param(j, k).a = zeros(3, 1);
%!       param(j, k).g = zeros(3, 1);
%!       param(j, k).holonomic = false;
%!       param(j, k).gamma = zeros(3, 1);
%!       param(j, k).N = 50;
%!       switch (j)
%!         case 1
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).gamma = [20; 45; 30] * pi / 180;
%!         case 2
%!           param(j, k).omega(k) = pert.omega(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           param(j, k).holonomic = true;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 2000;
%!         case 3
%!           param(j, k).omegadot(k) = pert.omegadot(k);
%!           param(j, k).t1 = 1e-2 / SI_unit_second;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 200;
%!         case 4
%!           param(j, k).F1(k) = pert.loads(k);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           switch (k)
%!           case 1
%!             param(j, k).gamma(3) = 5 * pi / 180;
%!           case 2
%!             param(j, k).gamma(3) = 45 * pi / 180;
%!           case 3
%!             param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 5
%!           param(j, k).M1(k) = pert.loads(k + 3);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           param(j, k).gamma(1) = 30 * pi / 180;
%!         case 6
%!           param(j, k).a(k) = pert.a(k);
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).N = 200;
%!           switch (k)
%!           case 1
%!             param(j, k).gamma(3) = 5 * pi / 180;
%!           case 2
%!             param(j, k).gamma(3) = 45 * pi / 180;
%!           case 3
%!             param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 7
%!           param(j, k).g(k) = pert.g(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           switch (k)
%!           case 1
%!             param(j, k).gamma(3) = 5 * pi / 180;
%!           case 2
%!             param(j, k).gamma(3) = 45 * pi / 180;
%!           case 3
%!             param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!       endswitch
%!       for l=1:2
%!         fd = -1;
%!         unwind_protect
%!           filename_mbdyn = sprintf("%s_%d_%d_%d.mbdyn", filename, j, k, l);
%!           [fd, msg] = fopen(filename_mbdyn, "w");
%!           if (fd == -1)
%!             error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!           endif
%!           fprintf(fd, "set: real a = %.16e;\n", a);
%!           fprintf(fd, "set: real b = %.16e;\n", b);
%!           fprintf(fd, "set: real c = %.16e;\n", c);
%!           fprintf(fd, "set: real d = %.16e;\n", d);
%!           for i=1:3
%!             fprintf(fd, "set: real gamma%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).gamma(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGA%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omega(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGAP%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omegadot(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real F1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).F1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real M1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).M1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real a%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).a(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real g%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).g(i));
%!           endfor
%!           fprintf(fd, "set: real t1 = %.16e;\n", param(j, k).t1);
%!           fprintf(fd, "set: integer N = %d;\n", param(j, k).N);
%!           fprintf(fd, "set: integer M = %d;\n", options.number_of_beams);
%!           fputs(fd, "set: integer ref_id_ground = 1;\n");
%!           fputs(fd, "set: integer ref_id_tilt = 2;\n");
%!           fputs(fd, "set: integer joint_id_ground = 1;\n");
%!           fputs(fd, "set: integer force_id1;\n");
%!           fputs(fd, "set: integer torque_id1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_PHI1 = 1;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP1 = 3;\n");
%!           fputs(fd, "set: integer drive_id_PHI2 = 4;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA2 = 5;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP2 = 6;\n");
%!           fputs(fd, "set: integer drive_id_PHIx = 7;\n");
%!           fputs(fd, "set: integer drive_id_PHIy = 8;\n");
%!           fputs(fd, "set: integer drive_id_PHIz = 9;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAx = 10;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAy = 11;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAz = 12;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPx = 13;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPy = 14;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPz = 15;\n");
%!           fputs(fd, "set: integer drive_id_Xx = 16;\n");
%!           fputs(fd, "set: integer drive_id_Xy = 17;\n");
%!           fputs(fd, "set: integer drive_id_Xz = 18;\n");
%!           fputs(fd, "set: integer drive_id_XPx = 19;\n");
%!           fputs(fd, "set: integer drive_id_XPy = 20;\n");
%!           fputs(fd, "set: integer drive_id_XPz = 21;\n");
%!           fputs(fd, "set: integer drive_id_XPPx = 22;\n");
%!           fputs(fd, "set: integer drive_id_XPPy = 23;\n");
%!           fputs(fd, "set: integer drive_id_XPPz = 24;\n");
%!           fputs(fd, "set: integer drive_id_gx = 25;\n");
%!           fputs(fd, "set: integer drive_id_gy = 26;\n");
%!           fputs(fd, "set: integer drive_id_gz = 27;\n");
%!           fputs(fd, "set: integer drive_id_F1x = 28;\n");
%!           fputs(fd, "set: integer drive_id_F1y = 29;\n");
%!           fputs(fd, "set: integer drive_id_F1z = 30;\n");
%!           fputs(fd, "set: integer drive_id_M1x = 31;\n");
%!           fputs(fd, "set: integer drive_id_M1y = 32;\n");
%!           fputs(fd, "set: integer drive_id_M1z = 33;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "set: integer ref_id_modal = 3;\n");
%!               fputs(fd, "set: integer node_id_modal = 1;\n");
%!               fputs(fd, "set: integer ref_id_interface1 = 4;\n");
%!               fputs(fd, "set: integer node_id_interface1 = 2;\n");
%!               fputs(fd, "set: integer elem_id_modal = 2;\n");
%!             case 2
%!               fputs(fd, "set: integer ref_id_beam1 = 3;\n");
%!               fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!               fprintf(fd, "set: real E = %.16e;\n", mesh.material_data.E);
%!               fprintf(fd, "set: real nu = %.16e;\n", mesh.material_data.nu);
%!               fprintf(fd, "set: real rho = %.16e;\n", mesh.material_data.rho);
%!               fprintf(fd, "set: real alpha = %.16e;\n", mesh.material_data.alpha);
%!               fprintf(fd, "set: real beta = %.16e;\n", mesh.material_data.beta);
%!               fputs(fd, "set: real G = E / (2. * (1. + nu));\n");
%!               fputs(fd, "set: real A = b * c;\n");
%!               fputs(fd, "set: real As = 9. / 10. * A;\n");
%!               fputs(fd, "set: real Iy = b * c^3 / 12.;\n");
%!               fputs(fd, "set: real Iz = c * b^3 / 12.;\n");
%!               fputs(fd, "set: real Ip = Iy + Iz;\n");
%!               c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!               w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!               fprintf(fd, "set: real It = %.16e;\n", interp1(w_h, c2, max(c, b) / min(c, b)) * max(c, b) * min(c, b)^3);
%!           endswitch
%!           fputs(fd, "begin: data;\n");
%!           fputs(fd, "        problem: initial value;\n");
%!           fputs(fd, "end: data;\n");
%!           fputs(fd, "begin: initial value;\n");
%!           fputs(fd, "        initial time: 0;\n");
%!           fputs(fd, "        final time: t1;\n");
%!           fputs(fd, "        time step: t1 / N;\n");
%!           fputs(fd, "        max time step: t1 / N;\n");
%!           fputs(fd, "        min time step: t1 / N;\n");
%!           fputs(fd, "        method: ss4, 0.;\n");
%!           fputs(fd, "        tolerance: 1e-4, test, minmax, 1e-6, test,minmax;\n");
%!           fputs(fd, "        max iterations: 1000;\n");
%!           fputs(fd, "        derivatives max iterations: 50;\n");
%!           fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!           fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!           fputs(fd, "        strategy: factor, 0.8, 3, 1.25, 3, 3, 6;\n");
%!           fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!           fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 100;\n");
%!           fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!           fputs(fd, "             keep jacobian matrix,\n");
%!           fputs(fd, "             inner iterations before assembly, 6,\n");
%!           fputs(fd, "             jacobian operator, newton krylov,\n");
%!           fputs(fd, "             solver, line search based,\n");
%!           fputs(fd, "             line search method, backtrack,\n");
%!           fputs(fd, "             recovery step type, constant,\n");
%!           fputs(fd, "             recovery step, 1e-6,\n");
%!           fputs(fd, "             verbose, yes,\n");
%!           fputs(fd, "             forcing term, type 2,\n");
%!           fputs(fd, "             direction, newton,\n");
%!           fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!           fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!           fputs(fd, "             linear solver, gmres,\n");
%!           fputs(fd, "             linear solver max iterations, 300,\n");
%!           fputs(fd, "             minimum step, 1e-12,\n");
%!           fputs(fd, "             krylov subspace size, 300;\n");
%!           fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!           fputs(fd, "        threads: assembly, 1;\n");
%!           fputs(fd, "    eigenanalysis: list, 1, t1,\n");
%!           fputs(fd, "    output matrices, \n");
%!           fprintf(fd, "          parameter, %.16e,\n", 1);
%!           fputs(fd, "    output eigenvectors,\n");
%!           fputs(fd, "        output geometry,\n");
%!           fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1 / SI_unit_second^-1, 100000 / SI_unit_second^-1);
%!           switch (l)
%!           case 1
%!             fputs(fd, "    use lapack, balance, permute, suffix format, \"%02d\";\n");
%!           case 2
%!             fprintf(fd, "    use arpack,%d,%d,0.,suffix format,\"%%02d\";\n", 2 * cms_opt.modes.number, 4 * cms_opt.modes.number + 1);
%!           endswitch
%!           fputs(fd, "end: initial value;\n");
%!           fputs(fd, "begin: control data;\n");
%!           fputs(fd, "        output meter: closest next, 0., forever, t1 / 100.;\n");
%!           switch (l)
%!           case 2
%!             fputs(fd, "        rigid body kinematics: drive,\n");
%!             fputs(fd, "            angular velocity,\n");
%!             fputs(fd, "                   component,\n");
%!             for i=1:3
%!               fprintf(fd, "                postponed, drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!             endfor
%!             fputs(fd, "            acceleration,\n");
%!             fputs(fd, "                   component,\n");
%!             for i=1:3
%!               fprintf(fd, "               postponed, drive_id_XPP%s,\n", {"x", "y", "z"}{i});
%!             endfor
%!             fputs(fd, "            angular acceleration,\n");
%!             fputs(fd, "                   component");
%!             for i=1:3
%!               fprintf(fd, ",\n               postponed, drive_id_OMEGAP%s", {"x","y","z"}{i});
%!             endfor
%!             fputs(fd, ";\n");
%!           endswitch
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       use automatic differentiation;\n");
%!           fputs(fd, "       default output: none, structural nodes, drive callers;\n");
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       output precision: 16;\n");
%!           fputs(fd, "       max iterations: 0;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural nodes:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# interface 1\n");
%!               fputs(fd, "        ;\n");
%!               fputs(fd, "        joints:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# ground\n");
%!               fputs(fd, "        ;\n");
%!             case 2
%!               fputs(fd, "       structural nodes: 2 * M + 1;\n");
%!               fputs(fd, "       rigid bodies: 2 * M + 1;\n");
%!               fputs(fd, "       beams: M;\n");
%!               fputs(fd, "       joints: 1;\n");
%!           endswitch
%!           fputs(fd, "        forces: 2;\n");
%!           fputs(fd, "        gravity;\n");
%!           fputs(fd, "end: control data;\n");
%!           fputs(fd, "drive caller: drive_id_PHI1, string, \"(((pi*Time)/(2*t1)-sin((pi*Time)/t1)/2)*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA1, string, \"sin((pi*Time)/(2*t1))^2\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP1, string, \"(pi*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1)))/t1\";\n");
%!           fputs(fd, "drive caller: drive_id_PHI2, string, \"-(4*sin((pi*Time)/(2*t1))^3*t1^2)/(3*pi^2)\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA2, string, \"-(2*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1))^2*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP2, string, \"-(3*cos((pi*Time)/(2*t1))^2-1)*sin((pi*Time)/(2*t1))\";\n");
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_PHI%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_PHI1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGA1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGAP%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGAP1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_X%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XPP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_g%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, g%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_F1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, F1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_M1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, M1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, "reference: ref_id_ground,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, eye,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, null;\n");
%!           fputs(fd, "reference: ref_id_tilt,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, euler123, gammax, gammay, gammaz,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, null;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "reference: ref_id_modal,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!               fputs(fd, "reference: ref_id_interface1,\n");
%!               fputs(fd, "        reference, ref_id_modal, a + d,  0., 0.,\n");
%!               fputs(fd, "        reference, ref_id_modal, eye,\n");
%!               fputs(fd, "        reference, ref_id_modal, null,\n");
%!               fputs(fd, "        reference, ref_id_modal, null;\n");
%!             case 2
%!               fputs(fd, "reference: ref_id_beam1,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!           endswitch
%!           fputs(fd, "begin: nodes;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural: node_id_modal, modal,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, eye,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, null, accelerations, yes;\n");
%!               fputs(fd, "        structural: node_id_interface1, static,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, eye,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null;\n");
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!                 fprintf(fd, "                reference, ref_id_beam1, 0.5 * a / M * %d, 0., 0.,\n", i - 1);
%!                 fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null, accelerations, yes;\n");
%!               endfor
%!           endswitch
%!           fputs(fd, "end: nodes;\n");
%!           fputs(fd, "begin: elements;\n");
%!           switch (l)
%!           case 1
%!           fputs(fd, "        joint: joint_id_ground, total pin joint,\n");
%!           fprintf(fd, "                %s,\n", {"node_id_modal", "node_id_beam1"}{l});
%!           fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!           fputs(fd, "               position constraint,\n");
%!           if (~param(j, k).holonomic)
%!             fputs(fd, "                        velocity, velocity, velocity,\n");
%!             fputs(fd, "                        component,\n");
%!             for i=1:3
%!               fprintf(fd, "                      reference, drive_id_XP%s,\n", {"x", "y", "z"}{i});
%!             endfor
%!           else
%!             fputs(fd, "                        active, active, active,\n");
%!             fputs(fd, "                        component,\n");
%!             for i=1:3
%!               fprintf(fd, "                      reference, drive_id_X%s,\n", {"x", "y", "z"}{i});
%!             endfor
%!           endif
%!           fputs(fd, "               orientation constraint,\n");
%!           if (~param(j, k).holonomic)
%!             fputs(fd, "                        angular velocity, angular velocity, angular velocity,\n");
%!             fputs(fd, "                        component");
%!             for i=1:3
%!               fprintf(fd, ",\n                   reference, drive_id_OMEGA%s", {"x","y","z"}{i});
%!             endfor
%!           else
%!             fputs(fd, "                        active, active, active,\n");
%!             fputs(fd, "                        component");
%!             for i=1:3
%!               fprintf(fd, ",\n                   reference, drive_id_PHI%s", {"x","y","z"}{i});
%!             endfor
%!           endif
%!           fputs(fd, ";\n");
%!           case 2
%!             fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!           endswitch
%!           switch (l)
%!             case 1
%!               fprintf(fd, "        include: \"%s.elm\";\n", filename);
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!                 fputs(fd, "               rho * A * a / (2 * M + 1),\n");
%!                 fputs(fd, "               reference, node, null, \n");
%!                 fputs(fd, "               diag,   rho * Ip * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iy * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iz * a / (2 * M + 1),\n");
%!                 fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!               endfor
%!               for i=1:options.number_of_beams
%!                 fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "               linear elastic generic, \n");
%!                 fputs(fd, "               diag, E * A , G * As, G * As, \n");
%!                 fputs(fd, "                     G * It, E * Iy, E * Iz,\n");
%!                 fputs(fd, "               same,\n");
%!                 fputs(fd, "               same;\n");
%!               endfor
%!           endswitch
%!           fprintf(fd, "        force: force_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M"}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_F1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fprintf(fd, "        couple: torque_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M"}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_M1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fputs(fd, "        gravity: uniform, component");
%!           for i=1:3
%!             fprintf(fd, ",\n       reference, drive_id_g%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd,";\n");
%!           fputs(fd, "end: elements;\n");
%!         unwind_protect_cleanup
%!           if (fd ~= -1)
%!             fclose(fd);
%!           endif
%!           fd = -1;
%!         end_unwind_protect
%!         options_mbd.output_file = sprintf("%s_%d_%d_%d_mbd", filename, j, k, l);
%!         if (~options.verbose)
%!           options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!         endif
%!         options_mbd.mbdyn_command = "mbdyn";
%!         options_eig.positive_frequencies = false;
%!         if (options.verbose)
%!           shell(sprintf("cat %s | nl", filename_mbdyn));
%!         endif
%!         res(j, k, l).info = mbdyn_solver_run(filename_mbdyn, options_mbd);
%!         output_file_rel_frame = [options_mbd.output_file, "_rel"];
%!         mbdyn_post_abs_to_rel(1, options_mbd.output_file, output_file_rel_frame, 0);
%!         exts = {".log", ".out"};
%!         for i=1:numel(exts)
%!           [err, msg] = symlink([options_mbd.output_file, exts{i}], [output_file_rel_frame, exts{i}]);
%!           if (err ~= 0)
%!             error("failed to create symlink: %s", msg);
%!           endif
%!         endfor
%!         [res(j, k, l).t, ...
%!          res(j, k, l).trajectory, ...
%!          res(j, k, l).deformation, ...
%!          res(j, k, l).velocity, ...
%!          res(j, k, l).acceleration, ...
%!          res(j, k, l).node_id, ...
%!          res(j, k, l).force, ...
%!          res(j, k, l).force_id, ...
%!          res(j, k, l).force_node_id, ...
%!          res(j, k, l).orientation_description] = mbdyn_post_load_output_struct(output_file_rel_frame);
%!         res(j, k, l).log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!         [res(j, k, l).drive_id, ...
%!         res(j, k, l).drive_value] = mbdyn_post_load_output_drv(options_mbd.output_file, [], numel(res(j, k, l).t));
%!         res(j, k, l).modal = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, 0);
%!       endfor
%!     endfor
%!   endfor
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, :) = true;
%!   load_case_dof.locked_dof(mesh.groups.quad9(grp_idx_clamp).nodes, :) = true;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = options.number_of_threads;
%!   load_case = struct("omega", empty_cell, ...
%!                      "omegadot", empty_cell, ...
%!                      "loads", empty_cell, ...
%!                      "loaded_nodes", empty_cell, ...
%!                      "joints", empty_cell, ...
%!                      "g", empty_cell, ...
%!                      "tau0", empty_cell);
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   for i=1:numel(load_case)
%!     load_case(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!     load_case(i).loads = zeros(1, 6);
%!     load_case(i).omega = zeros(3, 1);
%!     load_case(i).omegadot = zeros(3, 1);
%!     load_case(i).g = zeros(3, 1);
%!     load_case(i).tau0.iso27 = zeros(rows(mesh.elements.iso27), columns(mesh.elements.iso27), 6);
%!   endfor
%!   sol_eig = struct("def", empty_cell, "lambda", empty_cell, "f", empty_cell, "D", empty_cell);
%!   sol_eig_red = struct("lambda_red", empty_cell, "Ured", empty_cell);
%!   for j=idx_j
%!     for k=idx_k
%!       R = euler123_to_rotation_matrix(param(j, k).gamma);
%!       load_case(j, k).omega = R.' * param(j, k).omega;
%!       load_case(j, k).omegadot = R.' * param(j, k).omegadot;
%!       load_case(j, k).loads = [(R.' * param(j, k).F1).', (R.' * param(j, k).M1).'];
%!       load_case(j, k).g = R.' * (param(j, k).g - param(j, k).a);
%!       [mat_ass.M, ...
%!        mat_ass.D, ...
%!        mat_ass.K, ...
%!        mat_ass.KOMEGA, ...
%!        mat_ass.KOMEGA_DOT, ...
%!        mat_ass.DOMEGA, ...
%!        mat_ass.R, ...
%!        mat_ass.mat_info, ...
%!        mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                            dof_map, ...
%!                                            [FEM_MAT_MASS, ...
%!                                             FEM_MAT_DAMPING, ...
%!                                             FEM_MAT_STIFFNESS, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                             FEM_MAT_DAMPING_OMEGA, ...
%!                                             FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case(j, k));
%!       cms_opt.symmetric = false;
%!       sol_statjk = fem_sol_static(mesh, dof_map, mat_ass, cms_opt);
%!       sol_statjk.stress = fem_ass_matrix(mesh, ...
%!                                          dof_map, ...
%!                                          [FEM_VEC_STRESS_CAUCH], ...
%!                                          load_case(j, k), ...
%!                                          sol_statjk);
%!       sol_stat(j, k) = sol_statjk;
%!       load_case(j, k).tau0 = sol_stat(j, k).stress.tau;
%!       mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS_TAU0], ...
%!                                      load_case(j, k));
%!       mat_ass.K += mat_ass.KOMEGA + mat_ass.KOMEGA_DOT + mat_ass.KTAU0;
%!       mat_ass.D += mat_ass.DOMEGA;
%!       sol_eig(j, k) = fem_sol_modal_damped(mesh, ...
%!                                            dof_map, ...
%!                                            mat_ass, ...
%!                                            cms_opt.modes.number, ...
%!                                            cms_opt);
%!       Mred = mat_ass_cms.Mred;
%!       Dred = mat_ass_cms.Dred;
%!       Kred = mat_ass_cms.Kred;
%!       Dred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.DOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA_DOT(dof_map.idx_node, dof_map.idx_node), "Full");
%!       omegaq = [load_case(j, k).omega.^2;
%!                 load_case(j, k).omega(1) * load_case(j, k).omega(2);
%!                 load_case(j, k).omega(2) * load_case(j, k).omega(3);
%!                 load_case(j, k).omega(3) * load_case(j, k).omega(1)];
%!       idx = int32(0);
%!       for i=1:numel(omegaq)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * omegaq(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).omegadot)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).omegadot(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).g)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred -= mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).g(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).loads)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).loads(i);
%!       endfor
%!       [sol_eig_red(j, k).Ured, sol_eig_red(j, k).lambda_red] = fem_sol_eigsd(Kred, Dred, Mred, cms_opt.modes.number, cms_opt);
%!     endfor
%!   endfor
%! tol_abs = [0, 0] / SI_unit_second^-1;
%! tol_rel = [0.3e-2, 3e-2];
%! tol_disp_rel = 3e-2;
%! err_u_modal = err_v_modal = zeros(size(param));
%! printf("deformation/velocity:\n");
%! colors = rainbow(3);
%! width = 1:size(res, 3);
%! linestyle = {"-", "--"};
%! for i=idx_j
%!   for j=idx_k
%!      u_modal = res(i, j, 1).trajectory{end} - res(i, j, 1).trajectory{end}(1, :);
%!      u_beam = res(i, j, 2).trajectory{end} - res(i, j, 2).trajectory{end}(1, :);
%!      v_modal = res(i, j, 1).velocity{end};
%!      v_beam = res(i, j, 2).velocity{end};
%!      if (options.plot)
%!      figure("visible", "off");
%!      hold on;
%!      for k=1:size(res, 3)
%!        for l=1:3
%!          hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l) - res(i, j, k).trajectory{end}(1, l)) * SI_unit_meter);
%!          set(hnd, "color", colors(l, :));
%!          set(hnd, "linewidth", width(k));
%!          set(hnd, "linestyle", linestyle{k});
%!        endfor
%!      endfor
%!      xlabel("t [s]");
%!      ylabel("u [m]");
%!      grid on;
%!      grid minor on;
%!      title(sprintf("linear displacement %d:%d", i, j));
%!      figure("visible", "off");
%!      hold on;
%!      for k=1:size(res, 3)
%!        for l=1:3
%!          hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l + 3) - res(i, j, k).trajectory{end}(1, l + 3)) * 180 / pi);
%!          set(hnd, "color", colors(l, :));
%!          set(hnd, "linewidth", width(k));
%!          set(hnd, "linestyle", linestyle{k});
%!        endfor
%!      endfor
%!      xlabel("t [s]");
%!      ylabel("Phi [deg]");
%!      grid on;
%!      grid minor on;
%!      title(sprintf("angular displacement %d:%d", i, j));
%!      endif
%!      err_u_modal(i, j) = max(max(abs(u_modal - u_beam))) / max(1, max(max(abs(u_beam))));
%!      err_v_modal(i, j) = max(max(abs(v_modal - v_beam))) / max(1, max(max(abs(v_beam))));
%!      printf("%d:%d %.1f%%/%.1f%%\n", i, j, 100 * err_u_modal(i, j), 100 * err_v_modal(i, j));
%!   endfor
%! endfor
%! printf("natural frequencies:\n");
%! MACR = cell(size(param));
%! result_data = struct("f_mbd", cell(size(param)), "f_fem", cell(size(param)));
%! for i=idx_j
%!   for j=idx_k
%!     f_fem = sort(sol_eig(i, j).f(:));
%!     f_fem = f_fem(f_fem > 0);
%!     f_mbd = zeros(rows(f_fem), size(res, 3));
%!     PhiR = zeros(6, rows(f_fem), size(res, 3));
%!     for k=1:size(res, 3)
%!       [f_mbd_k, idx_mbd_k] = sort(res(i, j, k).modal.f(:));
%!       D_mbd_k = res(i, j, k).modal.D(idx_mbd_k);
%!       idx_mbd_k = idx_mbd_k(f_mbd_k > 0);
%!       f_mbd_k = f_mbd_k(f_mbd_k > 0);
%!       idx_mbd_k = idx_mbd_k(1:rows(f_fem));
%!       f_mbd(:, k) = f_mbd_k(1:rows(f_fem));
%!       PhiR(:, :, k) = res(i, j, k).modal.VR(res(i, j, k).modal.idx(end) + (1:6), idx_mbd_k);
%!     endfor
%!     result_data(i, j).f_fem = f_fem;
%!     result_data(i, j).f_mbd = f_mbd;
%!     MACR{i, j} = MACL{i, j} = zeros(rows(f_fem), rows(f_fem));
%!     for k=1:rows(f_fem)
%!       for l=1:rows(f_fem)
%!         MACR{i, j}(k, l) = (PhiR(:, k, 1)' * PhiR(:, k, 2)) * conj(PhiR(:, k, 1)' * PhiR(:, k, 2)) / ((PhiR(:, k, 1)' * PhiR(:, k, 1)) * (PhiR(:, k, 2)' * PhiR(:, k, 2)));
%!       endfor
%!     endfor
%!     printf("%d:%d\n", i, j);
%!     for k=1:rows(f_fem)
%!       printf("%10.2f", f_fem(k) * SI_unit_second^-1);
%!       for l=1:columns(f_mbd)
%!         printf("\t%10.2f", f_mbd(k, l) * SI_unit_second^-1);
%!       endfor
%!       for l=1:columns(f_mbd)
%!         printf("\t%.1f%%", 100 * (f_mbd(k, l) / f_fem(k) - 1));
%!       endfor
%!       printf("\t%.3f", MACR{i, j}(k, k));
%!       fputs(stdout, "\n");
%!     endfor
%!    fputs(stdout, "\n\n");
%!   endfor
%! endfor
%! for i=idx_j
%!   for j=idx_k
%!     for k=1:rows(result_data(i, j).f_fem)
%!       for l=1:columns(result_data(i, j).f_mbd)
%!         assert_simple(result_data(i, j).f_mbd(k, l), result_data(i, j).f_fem(k), tol_abs(l) + tol_rel(l) * abs(result_data(i, j).f_fem(k)));
%!       endfor
%!     endfor
%!   endfor
%! endfor
%! assert_simple(all(all(err_u_modal < tol_disp_rel)));
%! for j=idx_j
%!   for k=idx_k
%!       tol = 2e-2;
%!       [lambda_s] = sortrows([imag(sol_eig(j, k).lambda(:)), real(sol_eig(j, k).lambda(:))],[1,2]);
%!       [lambda_red_s] = sortrows([imag(sol_eig_red(j, k).lambda_red(:)), real(sol_eig_red(j, k).lambda_red(:))],[1,2]);
%!       K = min(20, rows(lambda_s));
%!       lambda_s = 1j * lambda_s(:,1) + lambda_s(:, 2);
%!       lambda_red_s = 1j * lambda_red_s(:, 1) + lambda_red_s(:, 2);
%!       assert_simple(lambda_red_s(1:K), lambda_s(1:K), tol * norm(lambda_s(1:K)));
%!   endfor
%! endfor
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST13
%! pkg load mboct-fem-pkg;
%! ## Oskar Wallrapp, Richard Schwertassek, 1998
%! ## Dynamik flexibler Mehrkoerpersysteme
%! ## chapter 5, table 5.7, page 242
%! ## Natural frequencies of a rotating beam structure
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   options.verbose = false;
%!   options.number_of_modes = int32(20);
%!   options.interactive = false;
%!   c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!   w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!   param.num_fem_elem_per_sec = 1;
%!   param.N1 = int32(50);
%!   param.N2 = int32(20);
%!   param.E1 = 7e10 / SI_unit_pascal;
%!   param.E2 = 21e10 / SI_unit_pascal;
%!   param.nu1 = 0.3;
%!   param.nu2 = 0.3;
%!   param.rho1 = 3000 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.rho2 = 7895 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.l1 = 2 / SI_unit_meter;
%!   param.l2 = 0.4 / SI_unit_meter;
%!   param.h1 = 0.009 / SI_unit_meter;
%!   param.h2 = param.h1;
%!   param.w1 = 0.009 / SI_unit_meter;
%!   param.w2 = 0.0095 / SI_unit_meter;
%!   param.OMEGAx = 0 / SI_unit_second^-1;
%!   param.OMEGAy = 0 / SI_unit_second^-1;
%!   param.OMEGAz = 6 / SI_unit_second^-1;
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "w1 = %.16e;\n", param.w1);
%!     fprintf(fd, "w2 = %.16e;\n", param.w2);
%!     fprintf(fd, "l1 = %.16e;\n", param.l1);
%!     fprintf(fd, "l2 = %.16e;\n", param.l2);
%!     fprintf(fd, "h1 = %.16e;\n", param.h1);
%!     fprintf(fd, "h = %.16e;\n", min([param.w1, param.w2]) / param.num_fem_elem_per_sec);
%!     fputs(fd, "Point(1) = {0, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(2) = {l1 - 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(3) = {l1 - 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(4) = {l1 + 0.5 * w2, l2, 0};\n");
%!     fputs(fd, "Point(5) = {l1 + 0.5 * w2, 0.5 * w1, 0};\n");
%!     fputs(fd, "Point(6) = {l1 + 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(7) = {l1 - 0.5 * w2, -0.5 * w1, 0};\n");
%!     fputs(fd, "Point(8) = {0, -0.5 * w1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 1};\n");
%!     fputs(fd, "Line(9) = {2, 7};\n");
%!     fputs(fd, "Line(10) = {2, 5};\n");
%!     fputs(fd, "Line Loop(1) = {1, 9, 7, 8};\n");
%!     fputs(fd, "Line Loop(2) = {10, 5, 6, 9};\n");
%!     fputs(fd, "Line Loop(3) = {2, 3, 4, 10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Plane Surface(2) = {2};\n");
%!     fputs(fd, "Plane Surface(3) = {3};\n");
%!     fputs(fd, "hw1 = w1 / Round(w1 / h);\n");
%!     fputs(fd, "hw2 = w2 / Round(w2 / h);\n");
%!     fputs(fd, "Transfinite Curve(1) = Round((l1 - 0.5 * w2) / hw1) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Round((l2 - 0.5 * w1) / hw2) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Round((l2 - 0.5 * w1) / hw2) + 1;\n");
%!     fputs(fd, "Transfinite Curve(5) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(6) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(7) = Round((l1 - 0.5 * w2) / hw1) + 1;\n");
%!     fputs(fd, "Transfinite Curve(8) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(9) = Round(w1 / h) + 1;\n");
%!     fputs(fd, "Transfinite Curve(10) = Round(w2 / h) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,h1}{Surface{1,2,3}; Layers{Round(h1 / h)}; Recombine;};\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Recombine Surface{2,12};\n");
%!     fputs(fd, "Recombine Surface{3,16};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Transfinite Surface(2) = {};\n");
%!     fputs(fd, "Transfinite Surface(3) = {};\n");
%!     fputs(fd, "Physical Surface(\"clamp\", 1) = {7};\n");
%!     fputs(fd, "Physical Volume(\"beam1\", 1) = {1, 2};\n");
%!     fputs(fd, "Physical Volume(\"beam2\", 2) = {3};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=0;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   mesh.material_data = struct("E", cell(1, 2), "nu", cell(1, 2), "rho", cell(1, 2));
%!   mesh.material_data(1).E = param.E1;
%!   mesh.material_data(1).nu = param.nu1;
%!   mesh.material_data(1).rho = param.rho1;
%!   mesh.material_data(2).E = param.E2;
%!   mesh.material_data(2).nu = param.nu2;
%!   mesh.material_data(2).rho = param.rho2;
%!   mesh.materials.iso27 = zeros(rows(mesh.elements.iso27), 1, "int32");
%!   grp_idx_beam1 = find([[mesh.groups.iso27].id] == 1);
%!   grp_idx_beam2 = find([[mesh.groups.iso27].id] == 2);
%!   grp_idx_clamp = find([[mesh.groups.quad9].id] == 1);
%!   mesh.materials.iso27(mesh.groups.iso27(grp_idx_beam1).elements) = 1;
%!   mesh.materials.iso27(mesh.groups.iso27(grp_idx_beam2).elements) = 2;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad9(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case.omega = [param.OMEGAx;
%!                      param.OMEGAy;
%!                      param.OMEGAz];
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = mbdyn_solver_num_threads_default();
%!   [mat_ass.M, ...
%!    mat_ass.Mdiag, ...
%!    mat_ass.K, ...
%!    mat_ass.KOMEGA, ...
%!    mat_ass.DOMEGA, ...
%!    mat_ass.R] = fem_ass_matrix(mesh, ...
%!                                dof_map, ...
%!                                [FEM_MAT_MASS, ...
%!                                 FEM_MAT_MASS_LUMPED, ...
%!                                 FEM_MAT_STIFFNESS, ...
%!                                 FEM_MAT_STIFFNESS_OMEGA, ...
%!                                 FEM_MAT_DAMPING_OMEGA, ...
%!                                 FEM_VEC_LOAD_CONSISTENT], ...
%!                                load_case);
%!   sol_stat = fem_sol_static(mesh, dof_map, mat_ass);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case, ...
%!                                    sol_stat);
%!   load_case.tau0 = sol_stat.stress.tau;
%!   mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_STIFFNESS_TAU0], ...
%!                                  load_case);
%!   opt_solver.pre_scaling = true;
%!   opt_solver.refine_max_iter = int32(10);
%!   opt_solver.solver = "pardiso";
%!   opt_solver.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_solver.symmetric = true;
%!   sol_eig(1) = fem_sol_modal(mesh, ...
%!                            dof_map, ...
%!                            mat_ass, ...
%!                            options.number_of_modes, ...
%!                            0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   sol_eig_diag = fem_sol_modal(mesh, ...
%!                                dof_map, ...
%!                                setfield(mat_ass, "M", mat_ass.Mdiag), ...
%!                                options.number_of_modes, ...
%!                                0, 0, "shift-invert", opt_solver.solver, opt_solver.number_of_threads);
%!   mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0;
%!   mat_ass.D = mat_ass.DOMEGA;
%!   sol_eig(2) = fem_sol_modal_damped(mesh, ...
%!                                     dof_map, ...
%!                                     mat_ass, ...
%!                                     options.number_of_modes, ...
%!                                     opt_solver);
%!   unwind_protect
%!     filename_mbdyn = [filename, ".mbdyn"];
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     mbdyn_pre_write_param_file(fd, param);
%!     fprintf(fd, "set: real t1 = %g;\n", 100 / SI_unit_second);
%!     fputs(fd, "set: integer N = 20;\n");
%!     fputs(fd, "set: real G1 = E1 / (2. * (1. + nu1));\n");
%!     fputs(fd, "set: real G2 = E2 / (2. * (1. + nu2));\n");
%!     fputs(fd, "set: real A1 = w1 * h1;\n");
%!     fputs(fd, "set: real A2 = w2 * h2;\n");
%!     fputs(fd, "set: real As1 = 9. / 10. * A1;\n");
%!     fputs(fd, "set: real As2 = 9. / 10. * A2;\n");
%!     fputs(fd, "set: real Iy1 = w1 * h1^3 / 12.;\n");
%!     fputs(fd, "set: real Iy2 = w2 * h2^3 / 12.;\n");
%!     fputs(fd, "set: real Iz1 = h1 * w1^3 / 12.;\n");
%!     fputs(fd, "set: real Iz2 = h2 * w2^3 / 12.;\n");
%!     fputs(fd, "set: real Ip1 = Iy1 + Iz1;\n");
%!     fputs(fd, "set: real Ip2 = Iy2 + Iz2;\n");
%!     fprintf(fd, "set: real c21 = %.16e;\n", interp1(w_h, c2, max(param.w1, param.h1) / min(param.w1, param.h1)));
%!     fprintf(fd, "set: real c22 = %.16e;\n", interp1(w_h, c2, max(param.w2, param.h2) / min(param.w2, param.h2)));
%!     fputs(fd, "set: real It1 = c21 * h1 * w1^3;\n");
%!     fputs(fd, "set: real It2 = c22 * h2 * w2^3;\n");
%!     fputs(fd, "set: integer ref_id_ground = 1;\n");
%!     fputs(fd, "set: integer ref_id_beam1 = 2;\n");
%!     fputs(fd, "set: integer ref_id_beam2 = 3;\n");
%!     fputs(fd, "set: integer joint_id_ground = 1;\n");
%!     fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer node_id_beam2 = 2 * N1 + node_id_beam1;\n");
%!     fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer body_id_beam2 = 2 * N1 + body_id_beam1 + 1;\n");
%!     fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!     fputs(fd, "set: integer beam_id_beam2 = beam_id_beam1 + N1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 2 * t1;\n");
%!     fputs(fd, "        time step: t1 / N;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-5, test,norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             line search method, backtrack,\n");
%!     fputs(fd, "             recovery step type, constant,\n");
%!     fputs(fd, "             recovery step, 1e-6,\n");
%!     fputs(fd, "             verbose, yes,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "           eigenanalysis: list, 2, 0, 2 * t1,\n");
%!     fputs(fd, "           # output matrices, \n");
%!     fprintf(fd, "         parameter, %.16e,\n", 1e-4 / SI_unit_second);
%!     fputs(fd, "           output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1e-3 / SI_unit_second^-1, 200. / SI_unit_second^-1);
%!     fprintf(fd, "           use arpack, %d, %d, 0, suffix format, \"%%02d\";\n", 3 * options.number_of_modes, 10 * options.number_of_modes + 1);
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, t1 / 10.;\n");
%!     fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!     fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default output: none, structural nodes;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       structural nodes: 2 * N1 + 1 + 2 * N2;\n");
%!     fputs(fd, "       rigid bodies: 2 * N1 + 1 + 2 * N2 + 1;\n");
%!     fputs(fd, "       beams: N1 + N2;\n");
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_beam1,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, eye,\n");
%!     fputs(fd, "        reference, ref_id_ground, null,\n");
%!     fputs(fd, "        reference, ref_id_ground, null;\n");
%!     fputs(fd, "reference: ref_id_beam2,\n");
%!     fputs(fd, "        reference, ref_id_beam1, l1,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, 1, 0., 1., 0., 3, 0., 0., 1.,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null,\n");
%!     fputs(fd, "        reference, ref_id_beam1, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam1, 0.5 * l1 / N1 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null,\n");
%!       fputs(fd, "                reference, ref_id_beam1, null;\n");
%!     endfor
%!     for i=2:(2 * param.N2 + 1)
%!       fprintf(fd, "        structural: node_id_beam2 + %d, dynamic,\n", i - 1);
%!       fprintf(fd, "                reference, ref_id_beam2, 0.5 * l2 / N2 * %d, 0., 0.,\n", i - 1);
%!       fputs(fd, "                reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null,\n");
%!       fputs(fd, "                reference, ref_id_beam2, null;\n");
%!     endfor
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!     for i=1:(2 * param.N1 + 1)
%!       fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!       fputs(fd, "               rho1 * A1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho1 * Ip1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iy1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "                       rho1 * Iz1 * l1 / (2 * N1 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!     endfor
%!     for i=1:(2 * param.N2 + 1)
%!       fprintf(fd, "     body: body_id_beam2 + %d, \n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d,\n", i - 1);
%!       fputs(fd, "               rho2 * A2 * l2 / (2 * N2 + 1), \n");
%!       fputs(fd, "               reference, node, null, \n");
%!       fputs(fd, "               diag,   rho2 * Ip2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iy2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "                       rho2 * Iz2 * l2 / (2 * N2 + 1),\n");
%!       fputs(fd, "               orientation, reference, ref_id_beam2, eye;\n");
%!     endfor
%!     for i=1:param.N1
%!       fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E1 * A1 , G1 * As1, G1 * As1, \n");
%!       fputs(fd, "                     G1 * It1, E1 * Iy1, E1 * Iz1,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     for i=1:param.N2
%!       fprintf(fd, "        beam3: beam_id_beam2 + %d,\n", i - 1);
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fprintf(fd, "             node_id_beam2 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!       fputs(fd, "               orientation, reference, node, eye,\n");
%!       fputs(fd, "               reference, ref_id_beam2, eye,\n");
%!       fputs(fd, "               linear elastic generic, \n");
%!       fputs(fd, "               diag, E2 * A2 , G2 * As2, G2 * As2, \n");
%!       fputs(fd, "                     G2 * It2, E2 * Iy2, E2 * Iz2,\n");
%!       fputs(fd, "               same,\n");
%!       fputs(fd, "               same;\n");
%!     endfor
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = [filename, "_mbd"];
%!   if (~options.verbose)
%!     options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
%!   options_mbd.mbdyn_command = "mbdyn";
%!   options_eig.positive_frequencies = false;
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   res.log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   for i=1:2
%!     res.modal(i) = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, i - 1);
%!   endfor
%!   if (options.interactive)
%!   for j=1:numel(res.modal)
%!     for i=1:numel(res.modal(j).f)
%!       opt_modal.mode_index = i;
%!       opt_modal.scale = 100;
%!       mode_file = [options_mbd.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(options_mbd.output_file, [mode_file, ".mov"], opt_modal, res.modal(j));
%!       [err, msg] = symlink([options_mbd.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = options.interactive;
%!       opt_post.f_runEasyAnim = options.interactive;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!   endfor
%!   endif
%!   idx_mode = [1, 4, 6];
%!   ## table 5.5 (11 elements), table 5.7 (11 elements)
%!   fref =  [0.9614, 1.044;
%!            7.037,  7.478;
%!           16.67,  17.09];
%!   tol = 2e-2;
%!   for i=1:2
%!     fmbd = sort(res.modal(i).f(:)) * SI_unit_second^-1;
%!     fmbd = fmbd(fmbd > 0);
%!     ffem = sort(sol_eig(i).f(:)) * SI_unit_second^-1;
%!     ffem = ffem(ffem >= 0);
%!     Nfem = min(numel(fmbd),numel(ffem));
%!     assert_simple(fmbd(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(ffem(idx_mode), fref(:, i), tol * max(fref(:,i)));
%!     assert_simple(fmbd(1:Nfem), ffem(1:Nfem), tol * max(fmbd(1:Nfem)));
%!   endfor
%!   assert_simple(max(abs(sol_eig_diag.f/sol_eig(1).f - 1)) < 0.04);
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST14
%! pkg load mboct-fem-pkg;
%! printf("fem_cms_create2: test7\n");
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   a = 800e-3 / SI_unit_meter;
%!   b = 40e-3 / SI_unit_meter;
%!   c = 10e-3 / SI_unit_meter;
%!   d = 0e-3 / SI_unit_meter;
%!   h = 2 * c;
%!   options.interactive = false;
%!   options.plot = true;
%!   options.verbose = false;
%!   options.number_of_beams = int32(40);
%!   options.number_of_threads = mbdyn_solver_num_threads_default();
%!   if (options.plot)
%!     close all;
%!   endif
%!   fd = -1;
%!   filename_geo = [filename, "_gmsh.geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(filename_geo, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", filename_geo);
%!     endif
%!     fprintf(fd, "a = %.16e;\n", a);
%!     fprintf(fd, "b = %.16e;\n", b);
%!     fprintf(fd, "c = %.16e;\n", c);
%!     fprintf(fd, "h = %.16e;\n", h);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(2) = {0,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(3) = {a,  0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Point(4) = {a, -0.5 * b, -0.5 * c};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Line Loop(1) = {1, 2, 3, 4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(a / h)) + 1;\n");
%!     fputs(fd, "v1 = Extrude{0,0,c}{Surface{1}; Layers{Max(1, Round(c / h))}; Recombine;};\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Transfinite Surface(1) = {};\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {v1[2]};\n");
%!     fputs(fd, "Physical Surface(2) = {v1[4]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.interactive)
%!     pid = spawn("gmsh", {filename_geo});
%!     status = spawn_wait(pid);
%!   endif
%!   pid = spawn("gmsh", {"-format", "msh2", ...
%!                        "-3", ...
%!                        "-order", "2", ...
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   mesh.material_data.E = 70000e6 / SI_unit_pascal;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = 2700 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data.alpha = 0e-5 / (1 / SI_unit_second);
%!   mesh.material_data.beta = 0e-5 / (SI_unit_second);
%!   mesh.materials.iso27 = zeros(rows(mesh.elements.iso27), 1, "int32");
%!   grp_idx_beam = find([[mesh.groups.iso27].id] == 1);
%!   grp_idx_clamp = find([[mesh.groups.quad9].id] == 1);
%!   mesh.materials.iso27(mesh.groups.iso27(grp_idx_beam).elements) = 1;
%!   cms_opt.number_of_threads = options.number_of_threads;
%!   cms_opt.algorithm = "diag-shift-invert";
%!   cms_opt.nodes.modal.number = rows(mesh.nodes) + 2;
%!   cms_opt.nodes.modal.name = "node_id_modal";
%!   cms_opt.nodes.interfaces.number = rows(mesh.nodes) + 1;
%!   cms_opt.nodes.interfaces.name = "node_id_interface1";
%!   mesh.nodes(cms_opt.nodes.modal.number, 1:3) = [0, 0, 0];
%!   mesh.nodes(cms_opt.nodes.interfaces.number, 1:3) = [a + d, 0, 0];
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, [2], [cms_opt.nodes.interfaces.number], "quad9");
%!   cms_opt.refine_max_iter = 30;
%!   cms_opt.pre_scaling = false;
%!   cms_opt.solver = "pardiso";
%!   cms_opt.modes.number = 20;
%!   cms_opt.tolerance_tau = -1;
%!   cms_opt.element.name = "elem_id_modal";
%!   cms_opt.create_binary = true;
%!   cms_opt.use_binary = true;
%!   cms_opt.update_binary = true;
%!   cms_opt.invariants = true;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad9(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, :) = true; ## Avoid singular matrix
%!   [mesh_cms, mat_ass_cms, dof_map_cms, sol_eig_cms, cms_opt, sol_tau_cms] = fem_cms_create2(mesh, load_case_dof, cms_opt);
%!   fem_cms_export(filename, mesh_cms, dof_map_cms, mat_ass_cms, cms_opt);
%!   pert.omega = [1e2; 3e2; 2e2] / (SI_unit_rad / SI_unit_second);
%!   pert.omegadot = [1e5; 1e3; 3e3] / (SI_unit_rad / SI_unit_second^2);
%!   pert.loads = [[1e4; 1e3; 1e2] / (SI_unit_newton);
%!                 [1e2; 1e1; 1e1] / (SI_unit_newton * SI_unit_meter)];
%!   pert.g = [1e4; -1e3; -1e2] / (SI_unit_meter / SI_unit_second^2);
%!   pert.a = [-1e4; 1e3; 1e2] / (SI_unit_meter / SI_unit_second^2);
%!   empty_cell = cell(7, 3, 2);
%!   res = struct("info", empty_cell, ...
%!                "t", empty_cell, ...
%!                "trajectory", empty_cell, ...
%!                "deformation", empty_cell, ...
%!                "velocity", empty_cell, ...
%!                "acceleration", empty_cell, ...
%!                "node_id", empty_cell, ...
%!                "force", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "force_node_id", empty_cell, ...
%!                "orientation_description", empty_cell, ...
%!                "drive_id", empty_cell, ...
%!                "drive_value", empty_cell, ...
%!                "modal", empty_cell);
%!   empty_cell = cell(7, 3);
%!   param = struct("omega", empty_cell, ...
%!                  "omegadot", empty_cell, ...
%!                  "F1", empty_cell, ...
%!                  "M1", empty_cell, ...
%!                  "a", empty_cell, ...
%!                  "g", empty_cell, ...
%!                  "t1", empty_cell, ...
%!                  "holonomic", empty_cell);
%!   idx_j = 1:rows(param);
%!   idx_k = 1:columns(param);
%!   for j=idx_j
%!     for k=idx_k
%!       param(j, k).omega = zeros(3, 1);
%!       param(j, k).omegadot = zeros(3, 1);
%!       param(j, k).F1 = zeros(3, 1);
%!       param(j, k).M1 = zeros(3, 1);
%!       param(j, k).a = zeros(3, 1);
%!       param(j, k).g = zeros(3, 1);
%!       param(j, k).holonomic = false;
%!       param(j, k).gamma = zeros(3, 1);
%!       param(j, k).N = 50;
%!       switch (j)
%!         case 1
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).gamma = [20; 45; 30] * pi / 180;
%!         case 2
%!           param(j, k).omega(k) = pert.omega(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           param(j, k).holonomic = true;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 2000;
%!         case 3
%!           param(j, k).omegadot(k) = pert.omegadot(k);
%!           param(j, k).t1 = 1e-2 / SI_unit_second;
%!           param(j, k).gamma(1) = 45 * pi / 180;
%!           param(j, k).N = 200;
%!         case 4
%!           param(j, k).F1(k) = pert.loads(k);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           switch (k)
%!             case 1
%!               param(j, k).gamma(3) = 5 * pi / 180;
%!             case 2
%!               param(j, k).gamma(3) = 45 * pi / 180;
%!             case 3
%!               param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 5
%!           param(j, k).M1(k) = pert.loads(k + 3);
%!           param(j, k).t1 = 10 / SI_unit_second;
%!           param(j, k).gamma(1) = 30 * pi / 180;
%!         case 6
%!           param(j, k).a(k) = pert.a(k);
%!           param(j, k).t1 = 1e-3 / SI_unit_second;
%!           param(j, k).N = 200;
%!           switch (k)
%!             case 1
%!               param(j, k).gamma(3) = 5 * pi / 180;
%!             case 2
%!               param(j, k).gamma(3) = 45 * pi / 180;
%!             case 3
%!               param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!         case 7
%!           param(j, k).g(k) = pert.g(k);
%!           param(j, k).t1 = 1 / SI_unit_second;
%!           switch (k)
%!             case 1
%!               param(j, k).gamma(3) = 5 * pi / 180;
%!             case 2
%!               param(j, k).gamma(3) = 45 * pi / 180;
%!             case 3
%!               param(j, k).gamma(1) = 80 * pi / 180;
%!           endswitch
%!       endswitch
%!       for l=1:3
%!         switch (l)
%!           case 3
%!             nodes_file = sprintf("%s_%d_%d_%d.nod", filename, j, k, l);
%!             csl_file = sprintf("%s_%d_%d_%d.csl", filename, j, k, l);
%!             elem_file = sprintf("%s_%d_%d_%d.elm", filename, j, k, l);
%!             opt_mbd_mesh = struct();
%!             opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_clamp";
%!             opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!             opt_mbd_mesh.struct_nodes.type(cms_opt.nodes.interfaces.number) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!             opt_mbd_mesh.struct_nodes.type(cms_opt.nodes.modal.number) = MBDYN_NODE_TYPE_STATIC_STRUCT_DISP;
%!             opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!             opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!             load_case_empty = struct();
%!             opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!         endswitch
%!         fd = -1;
%!         unwind_protect
%!           filename_mbdyn = sprintf("%s_%d_%d_%d.mbdyn", filename, j, k, l);
%!           [fd, msg] = fopen(filename_mbdyn, "w");
%!           if (fd == -1)
%!             error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!           endif
%!           fprintf(fd, "set: real a = %.16e;\n", a);
%!           fprintf(fd, "set: real b = %.16e;\n", b);
%!           fprintf(fd, "set: real c = %.16e;\n", c);
%!           fprintf(fd, "set: real d = %.16e;\n", d);
%!           for i=1:3
%!             fprintf(fd, "set: real gamma%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).gamma(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGA%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omega(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real OMEGAP%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).omegadot(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real F1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).F1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real M1%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).M1(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real a%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).a(i));
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "set: real g%s = %.16e;\n", {"x","y","z"}{i}, param(j, k).g(i));
%!           endfor
%!           fprintf(fd, "set: real t1 = %.16e;\n", param(j, k).t1);
%!           fprintf(fd, "set: integer N = %d;\n", param(j, k).N);
%!           fprintf(fd, "set: integer M = %d;\n", options.number_of_beams);
%!           fputs(fd, "set: integer ref_id_ground = 1;\n");
%!           fputs(fd, "set: integer ref_id_tilt = 2;\n");
%!           fputs(fd, "set: integer joint_id_ground = 1;\n");
%!           fputs(fd, "set: integer force_id1;\n");
%!           fputs(fd, "set: integer torque_id1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_PHI1 = 1;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA1 = 2;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP1 = 3;\n");
%!           fputs(fd, "set: integer drive_id_PHI2 = 4;\n");
%!           fputs(fd, "set: integer drive_id_OMEGA2 = 5;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAP2 = 6;\n");
%!           fputs(fd, "set: integer drive_id_PHIx = 7;\n");
%!           fputs(fd, "set: integer drive_id_PHIy = 8;\n");
%!           fputs(fd, "set: integer drive_id_PHIz = 9;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAx = 10;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAy = 11;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAz = 12;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPx = 13;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPy = 14;\n");
%!           fputs(fd, "set: integer drive_id_OMEGAPz = 15;\n");
%!           fputs(fd, "set: integer drive_id_Xx = 16;\n");
%!           fputs(fd, "set: integer drive_id_Xy = 17;\n");
%!           fputs(fd, "set: integer drive_id_Xz = 18;\n");
%!           fputs(fd, "set: integer drive_id_XPx = 19;\n");
%!           fputs(fd, "set: integer drive_id_XPy = 20;\n");
%!           fputs(fd, "set: integer drive_id_XPz = 21;\n");
%!           fputs(fd, "set: integer drive_id_XPPx = 22;\n");
%!           fputs(fd, "set: integer drive_id_XPPy = 23;\n");
%!           fputs(fd, "set: integer drive_id_XPPz = 24;\n");
%!           fputs(fd, "set: integer drive_id_gx = 25;\n");
%!           fputs(fd, "set: integer drive_id_gy = 26;\n");
%!           fputs(fd, "set: integer drive_id_gz = 27;\n");
%!           fputs(fd, "set: integer drive_id_F1x = 28;\n");
%!           fputs(fd, "set: integer drive_id_F1y = 29;\n");
%!           fputs(fd, "set: integer drive_id_F1z = 30;\n");
%!           fputs(fd, "set: integer drive_id_M1x = 31;\n");
%!           fputs(fd, "set: integer drive_id_M1y = 32;\n");
%!           fputs(fd, "set: integer drive_id_M1z = 33;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "set: integer ref_id_modal = 3;\n");
%!               fputs(fd, "set: integer node_id_modal = 1;\n");
%!               fputs(fd, "set: integer ref_id_interface1 = 4;\n");
%!               fputs(fd, "set: integer node_id_interface1 = 2;\n");
%!               fputs(fd, "set: integer elem_id_modal = 2;\n");
%!             case 2
%!               fputs(fd, "set: integer ref_id_beam1 = 3;\n");
%!               fputs(fd, "set: integer node_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer body_id_beam1 = 1;\n");
%!               fputs(fd, "set: integer beam_id_beam1 = 1;\n");
%!               fprintf(fd, "set: real E = %.16e;\n", mesh.material_data.E);
%!               fprintf(fd, "set: real nu = %.16e;\n", mesh.material_data.nu);
%!               fprintf(fd, "set: real rho = %.16e;\n", mesh.material_data.rho);
%!               fprintf(fd, "set: real alpha = %.16e;\n", mesh.material_data.alpha);
%!               fprintf(fd, "set: real beta = %.16e;\n", mesh.material_data.beta);
%!               fputs(fd, "set: real G = E / (2. * (1. + nu));\n");
%!               fputs(fd, "set: real A = b * c;\n");
%!               fputs(fd, "set: real As = 9. / 10. * A;\n");
%!               fputs(fd, "set: real Iy = b * c^3 / 12.;\n");
%!               fputs(fd, "set: real Iz = c * b^3 / 12.;\n");
%!               fputs(fd, "set: real Ip = Iy + Iz;\n");
%!               c2  = [0.141, 0.166, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.312, 0.33];
%!               w_h = [    1,   1.2,   1.5,     2,   2.5,     3,     4,     5,    10,  inf];
%!               fprintf(fd, "set: real It = %.16e;\n", interp1(w_h, c2, max(c, b) / min(c, b)) * max(c, b) * min(c, b)^3);
%!             case 3
%!               fputs(fd, "set: integer ref_id_clamp = 3;\n");
%!               fprintf(fd, "set: integer node_id_interface1 = %d;\n", cms_opt.nodes.interfaces.number);
%!           endswitch
%!           fputs(fd, "begin: data;\n");
%!           fputs(fd, "        problem: initial value;\n");
%!           fputs(fd, "end: data;\n");
%!           fputs(fd, "begin: initial value;\n");
%!           fputs(fd, "        initial time: 0;\n");
%!           fputs(fd, "        final time: t1;\n");
%!           fputs(fd, "        time step: t1 / N;\n");
%!           fputs(fd, "        max time step: t1 / N;\n");
%!           fputs(fd, "        min time step: t1 / N;\n");
%!           fputs(fd, "        method: ss4, 0.;\n");
%!           fputs(fd, "        tolerance: 1e-4, test, minmax, 1e-6, test,minmax;\n");
%!           fputs(fd, "        max iterations: 1000;\n");
%!           fputs(fd, "        derivatives max iterations: 50;\n");
%!           fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!           fputs(fd, "        derivatives tolerance: 1e-6, 1;\n");
%!           fputs(fd, "        strategy: factor, 0.8, 3, 1.25, 3, 3, 6;\n");
%!           fputs(fd, "        output: iterations, solver condition number, stat, yes, cpu time;\n");
%!           fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!           fputs(fd, "        nonlinear solver: nox, modified, 100,\n");
%!           fputs(fd, "             keep jacobian matrix,\n");
%!           fputs(fd, "             inner iterations before assembly, 6,\n");
%!           fputs(fd, "             jacobian operator, newton krylov,\n");
%!           fputs(fd, "             solver, line search based,\n");
%!           fputs(fd, "             line search method, backtrack,\n");
%!           fputs(fd, "             recovery step type, constant,\n");
%!           fputs(fd, "             recovery step, 1e-6,\n");
%!           fputs(fd, "             verbose, yes,\n");
%!           fputs(fd, "             forcing term, type 2,\n");
%!           fputs(fd, "             direction, newton,\n");
%!           fputs(fd, "             weighted rms absolute tolerance, 0,\n");
%!           fputs(fd, "             weighted rms relative tolerance, 0,\n");
%!           fputs(fd, "             linear solver, gmres,\n");
%!           fputs(fd, "             linear solver max iterations, 30,\n");
%!           fputs(fd, "             minimum step, 1e-12,\n");
%!           fputs(fd, "             krylov subspace size, 30;\n");
%!           fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!           switch (l)
%!           case 3
%!             fprintf(fd, "        threads: assembly, %d;\n", cms_opt.number_of_threads);
%!           otherwise
%!             fputs(fd, "        threads: assembly, 1;\n");
%!           endswitch
%!           fputs(fd, "    eigenanalysis: list, 1, t1,\n");
%!           fputs(fd, "    output eigenvectors,\n");
%!           fputs(fd, "    # output matrices, \n");
%!           fputs(fd, "        output geometry,\n");
%!           fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", 1. / SI_unit_second^-1, 10000 / SI_unit_second^-1);
%!           switch (l)
%!             case 1
%!               fputs(fd, "    use lapack, balance, permute, suffix format, \"%02d\";\n");
%!             case {2, 3}
%!               fprintf(fd, "    use arpack,%d,%d,0.,suffix format,\"%%02d\";\n", cms_opt.modes.number, 2 * cms_opt.modes.number + 1);
%!           endswitch
%!           fputs(fd, "end: initial value;\n");
%!           fputs(fd, "begin: control data;\n");
%!           fputs(fd, "        output meter: closest next, 0., forever, t1 / 20.;\n");
%!           switch (l)
%!             case {2, 3}
%!               fputs(fd, "        rigid body kinematics: drive,\n");
%!               fputs(fd, "            angular velocity,\n");
%!               fputs(fd, "                   component,\n");
%!               for i=1:3
%!                 fprintf(fd, "                postponed, drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!               endfor
%!               fputs(fd, "            acceleration,\n");
%!               fputs(fd, "                   component,\n");
%!               for i=1:3
%!                 fprintf(fd, "               postponed, drive_id_XPP%s,\n", {"x", "y", "z"}{i});
%!               endfor
%!               fputs(fd, "            angular acceleration,\n");
%!               fputs(fd, "                   component");
%!               for i=1:3
%!                 fprintf(fd, ",\n               postponed, drive_id_OMEGAP%s", {"x","y","z"}{i});
%!               endfor
%!               fputs(fd, ";\n");
%!           endswitch
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       use automatic differentiation;\n");
%!           fputs(fd, "       default output: none, structural nodes, drive callers;\n");
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       output precision: 16;\n");
%!           fputs(fd, "       max iterations: 0;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural nodes:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# interface 1\n");
%!               fputs(fd, "        ;\n");
%!               fputs(fd, "        joints:\n");
%!               fputs(fd, "                +1		# modal\n");
%!               fputs(fd, "                +1		# ground\n");
%!               fputs(fd, "        ;\n");
%!               fputs(fd, "        forces: 2;\n");
%!             case 2
%!               fputs(fd, "       structural nodes: 2 * M + 1;\n");
%!               fputs(fd, "       rigid bodies: 2 * M + 1;\n");
%!               fputs(fd, "       beams: M;\n");
%!               fputs(fd, "       joints: 1;\n");
%!               fputs(fd, "       forces: 2;\n");
%!             case 3
%!               fprintf(fd, "     structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!               fprintf(fd, "     solids: %d;\n", opt_mbd_mesh.solids.number);
%!               fprintf(fd, "     genels: %d;\n", opt_mbd_mesh.genels.number);
%!               fprintf(fd, "     joints: %d;\n", opt_mbd_mesh.joints.number);
%!               fprintf(fd, "     forces: %d;\n", opt_mbd_mesh.forces.number + 2);
%!           endswitch
%!           fputs(fd, "        gravity;\n");
%!           fputs(fd, "end: control data;\n");
%!           fputs(fd, "drive caller: drive_id_PHI1, string, \"(((pi*Time)/(2*t1)-sin((pi*Time)/t1)/2)*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA1, string, \"sin((pi*Time)/(2*t1))^2\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP1, string, \"(pi*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1)))/t1\";\n");
%!           fputs(fd, "drive caller: drive_id_PHI2, string, \"-(4*sin((pi*Time)/(2*t1))^3*t1^2)/(3*pi^2)\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGA2, string, \"-(2*cos((pi*Time)/(2*t1))*sin((pi*Time)/(2*t1))^2*t1)/pi\";\n");
%!           fputs(fd, "drive caller: drive_id_OMEGAP2, string, \"-(3*cos((pi*Time)/(2*t1))^2-1)*sin((pi*Time)/(2*t1))\";\n");
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_PHI%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_PHI1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGA%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGA1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_OMEGAP%s,\n", {"x","y","z"}{i});
%!             fputs(fd, "  array, 2,\n");
%!             fprintf(fd, "     mult, const, OMEGA%s, reference, drive_id_OMEGAP1,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, OMEGAP%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_X%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_PHI2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGA2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_XPP%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, a%s, reference, drive_id_OMEGAP2;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_g%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, g%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_F1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, F1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           for i=1:3
%!             fprintf(fd, "drive caller: drive_id_M1%s,\n", {"x","y","z"}{i});
%!             fprintf(fd, "     mult, const, M1%s, reference, drive_id_OMEGA1;\n", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, "reference: ref_id_ground,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, eye,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, null;\n");
%!           fputs(fd, "reference: ref_id_tilt,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, euler123, gammax, gammay, gammaz,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, null;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "reference: ref_id_modal,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!               fputs(fd, "reference: ref_id_interface1,\n");
%!               fputs(fd, "        reference, ref_id_modal, a + d,  0., 0.,\n");
%!               fputs(fd, "        reference, ref_id_modal, eye,\n");
%!               fputs(fd, "        reference, ref_id_modal, null,\n");
%!               fputs(fd, "        reference, ref_id_modal, null;\n");
%!             case 2
%!               fputs(fd, "reference: ref_id_beam1,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!             case 3
%!               fputs(fd, "reference: ref_id_clamp,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, eye,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null,\n");
%!               fputs(fd, "        reference, ref_id_tilt, null;\n");
%!               fprintf(fd, "include: \"%s\";\n", csl_file);
%!           endswitch
%!           fputs(fd, "begin: nodes;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        structural: node_id_modal, modal,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, eye,\n");
%!               fputs(fd, "                reference, ref_id_modal, null,\n");
%!               fputs(fd, "                reference, ref_id_modal, null, accelerations, yes;\n");
%!               fputs(fd, "        structural: node_id_interface1, static,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, eye,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null;\n");
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "        structural: node_id_beam1 + %d, dynamic,\n", i - 1);
%!                 fprintf(fd, "                reference, ref_id_beam1, 0.5 * a / M * %d, 0., 0.,\n", i - 1);
%!                 fputs(fd, "                reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null,\n");
%!                 fputs(fd, "                reference, ref_id_beam1, null, accelerations, yes;\n");
%!               endfor
%!             case 3
%!               fprintf(fd, "include: \"%s\";\n", nodes_file);
%!           endswitch
%!           fputs(fd, "end: nodes;\n");
%!           fputs(fd, "begin: elements;\n");
%!           switch (l)
%!             case 1
%!               fputs(fd, "        joint: joint_id_ground, total pin joint,\n");
%!               fprintf(fd, "                %s,\n", {"node_id_modal", "node_id_beam1"}{l});
%!               fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        position, reference, %s, null,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        position orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fprintf(fd, "                        rotation orientation, reference, %s, eye,\n", {"ref_id_ground", "ref_id_ground"}{l});
%!               fputs(fd, "               position constraint,\n");
%!               if (~param(j, k).holonomic)
%!                 fputs(fd, "                        velocity, velocity, velocity,\n");
%!                 fputs(fd, "                        component,\n");
%!                 for i=1:3
%!                   fprintf(fd, "                      reference, drive_id_XP%s,\n", {"x", "y", "z"}{i});
%!                 endfor
%!               else
%!                 fputs(fd, "                        active, active, active,\n");
%!                 fputs(fd, "                        component,\n");
%!                 for i=1:3
%!                   fprintf(fd, "                      reference, drive_id_X%s,\n", {"x", "y", "z"}{i});
%!                 endfor
%!               endif
%!               fputs(fd, "               orientation constraint,\n");
%!               if (~param(j, k).holonomic)
%!                 fputs(fd, "                        angular velocity, angular velocity, angular velocity,\n");
%!                 fputs(fd, "                        component");
%!                 for i=1:3
%!                   fprintf(fd, ",\n                   reference, drive_id_OMEGA%s", {"x","y","z"}{i});
%!                 endfor
%!               else
%!                 fputs(fd, "                        active, active, active,\n");
%!                 fputs(fd, "                        component");
%!                 for i=1:3
%!                   fprintf(fd, ",\n                   reference, drive_id_PHI%s", {"x","y","z"}{i});
%!                 endfor
%!               endif
%!               fputs(fd, ";\n");
%!             case 2
%!               fputs(fd, "joint: joint_id_ground, clamp, node_id_beam1, node, node;\n");
%!           endswitch
%!           switch (l)
%!             case 1
%!               fprintf(fd, "        include: \"%s.elm\";\n", filename);
%!             case 2
%!               for i=1:(2 * options.number_of_beams + 1)
%!                 fprintf(fd, "     body: body_id_beam1 + %d, \n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d,\n", i - 1);
%!                 fputs(fd, "               rho * A * a / (2 * M + 1),\n");
%!                 fputs(fd, "               reference, node, null, \n");
%!                 fputs(fd, "               diag,   rho * Ip * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iy * a / (2 * M + 1),\n");
%!                 fputs(fd, "                       rho * Iz * a / (2 * M + 1),\n");
%!                 fputs(fd, "               orientation, reference, ref_id_beam1, eye;\n");
%!               endfor
%!               for i=1:options.number_of_beams
%!                 fprintf(fd, "        beam3: beam_id_beam1 + %d,\n", i - 1);
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1));
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 1);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fprintf(fd, "             node_id_beam1 + %d, position, reference, node, null,\n", 2 * (i - 1) + 2);
%!                 fputs(fd, "               orientation, reference, node, eye,\n");
%!                 fputs(fd, "               reference, ref_id_beam1, eye,\n");
%!                 fputs(fd, "               linear elastic generic, \n");
%!                 fputs(fd, "               diag, E * A , G * As, G * As, \n");
%!                 fputs(fd, "                     G * It, E * Iy, E * Iz,\n");
%!                 fputs(fd, "               same,\n");
%!                 fputs(fd, "               same;\n");
%!               endfor
%!             case 3
%!               fprintf(fd, "include: \"%s\";\n", elem_file);
%!           endswitch
%!           fprintf(fd, "        force: force_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M", sprintf("%d", cms_opt.nodes.interfaces.number)}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_F1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fprintf(fd, "        couple: torque_id1, absolute, %s,\n", {"node_id_interface1", "node_id_beam1 + 2 * M", sprintf("%d", cms_opt.nodes.interfaces.number)}{l});
%!           fputs(fd, "               position, reference, node, null,\n");
%!           fputs(fd, "                  component");
%!           for i=1:3
%!             fprintf(fd, ",\n             reference, drive_id_M1%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd, ";\n");
%!           fputs(fd, "        gravity: uniform, component");
%!           for i=1:3
%!             fprintf(fd, ",\n       reference, drive_id_g%s", {"x","y","z"}{i});
%!           endfor
%!           fputs(fd,";\n");
%!           fputs(fd, "end: elements;\n");
%!         unwind_protect_cleanup
%!           if (fd ~= -1)
%!             fclose(fd);
%!           endif
%!           fd = -1;
%!         end_unwind_protect
%!         options_mbd.output_file = sprintf("%s_%d_%d_%d_mbd", filename, j, k, l);
%!         if (~options.verbose)
%!           options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!         endif
%!         options_mbd.mbdyn_command = "mbdyn -C";
%!         options_eig.positive_frequencies = false;
%!         if (options.verbose)
%!           shell(sprintf("cat %s | nl", filename_mbdyn));
%!           switch (l)
%!           case 3
%!             shell(sprintf("cat %s | nl", csl_file));
%!             shell(sprintf("cat %s | nl", nodes_file));
%!             shell(sprintf("cat %s | nl", elem_file));
%!           endswitch
%!         endif
%!         res(j, k, l).info = mbdyn_solver_run(filename_mbdyn, options_mbd);
%!         output_file_rel_frame = [options_mbd.output_file, "_rel"];
%!         mbdyn_post_abs_to_rel(1, options_mbd.output_file, output_file_rel_frame, 0);
%!         exts = {".log", ".out"};
%!         for i=1:numel(exts)
%!           [err, msg] = symlink([options_mbd.output_file, exts{i}], [output_file_rel_frame, exts{i}]);
%!           if (err ~= 0)
%!             error("failed to create symlink: %s", msg);
%!           endif
%!         endfor
%!         [res(j, k, l).t, ...
%!          res(j, k, l).trajectory, ...
%!          res(j, k, l).deformation, ...
%!          res(j, k, l).velocity, ...
%!          res(j, k, l).acceleration, ...
%!          res(j, k, l).node_id, ...
%!          res(j, k, l).force, ...
%!          res(j, k, l).force_id, ...
%!          res(j, k, l).force_node_id, ...
%!          res(j, k, l).orientation_description] = mbdyn_post_load_output_struct(output_file_rel_frame);
%!         res(j, k, l).log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!         [res(j, k, l).drive_id, ...
%!          res(j, k, l).drive_value] = mbdyn_post_load_output_drv(options_mbd.output_file, [], numel(res(j, k, l).t));
%!         res(j, k, l).modal = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig, 0);
%!       endfor
%!     endfor
%!   endfor
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(mesh.groups.quad9(grp_idx_clamp).nodes, 1:3) = true;
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, :) = true; ## Avoid singular matrix
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = options.number_of_threads;
%!   load_case = struct("omega", empty_cell, ...
%!                      "omegadot", empty_cell, ...
%!                      "loads", empty_cell, ...
%!                      "loaded_nodes", empty_cell, ...
%!                      "joints", empty_cell, ...
%!                      "g", empty_cell, ...
%!                      "tau0", empty_cell);
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   for i=1:numel(load_case)
%!     load_case(i).loaded_nodes = cms_opt.nodes.interfaces.number;
%!     load_case(i).loads = zeros(1, 6);
%!     load_case(i).omega = zeros(3, 1);
%!     load_case(i).omegadot = zeros(3, 1);
%!     load_case(i).g = zeros(3, 1);
%!     load_case(i).tau0.iso27 = zeros(rows(mesh.elements.iso27), columns(mesh.elements.iso27), 6);
%!   endfor
%!   sol_eig = struct("def", empty_cell, "lambda", empty_cell, "f", empty_cell, "D", empty_cell);
%!   sol_eig_red = struct("lambda_red", empty_cell, "Ured", empty_cell);
%!   for j=idx_j
%!     for k=idx_k
%!       R = euler123_to_rotation_matrix(param(j, k).gamma);
%!       load_case(j, k).omega = R.' * param(j, k).omega;
%!       load_case(j, k).omegadot = R.' * param(j, k).omegadot;
%!       load_case(j, k).loads = [(R.' * param(j, k).F1).', (R.' * param(j, k).M1).'];
%!       load_case(j, k).g = R.' * (param(j, k).g - param(j, k).a);
%!       [mat_ass.M, ...
%!        mat_ass.D, ...
%!        mat_ass.K, ...
%!        mat_ass.KOMEGA, ...
%!        mat_ass.KOMEGA_DOT, ...
%!        mat_ass.DOMEGA, ...
%!        mat_ass.R, ...
%!        mat_ass.mat_info, ...
%!        mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                            dof_map, ...
%!                                            [FEM_MAT_MASS, ...
%!                                             FEM_MAT_DAMPING, ...
%!                                             FEM_MAT_STIFFNESS, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA, ...
%!                                             FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                             FEM_MAT_DAMPING_OMEGA, ...
%!                                             FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case(j, k));
%!       cms_opt.symmetric = false;
%!       sol_statjk = fem_sol_static(mesh, dof_map, mat_ass, cms_opt);
%!       sol_statjk.stress = fem_ass_matrix(mesh, ...
%!                                          dof_map, ...
%!                                          [FEM_VEC_STRESS_CAUCH], ...
%!                                          load_case(j, k), ...
%!                                          sol_statjk);
%!       sol_stat(j, k) = sol_statjk;
%!       load_case(j, k).tau0 = sol_stat(j, k).stress.tau;
%!       mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                      dof_map, ...
%!                                      [FEM_MAT_STIFFNESS_TAU0], ...
%!                                      load_case(j, k));
%!       mat_ass.K += mat_ass.KOMEGA + mat_ass.KOMEGA_DOT + mat_ass.KTAU0;
%!       mat_ass.D += mat_ass.DOMEGA;
%!       sol_eig(j, k) = fem_sol_modal_damped(mesh, ...
%!                                            dof_map, ...
%!                                            mat_ass, ...
%!                                            cms_opt.modes.number, ...
%!                                            cms_opt);
%!       Mred = mat_ass_cms.Mred;
%!       Dred = mat_ass_cms.Dred;
%!       Kred = mat_ass_cms.Kred;
%!       Dred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.DOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA(dof_map.idx_node, dof_map.idx_node), "Full");
%!       Kred += fem_cms_matrix_trans(mat_ass_cms.Tred, mat_ass.KOMEGA_DOT(dof_map.idx_node, dof_map.idx_node), "Full");
%!       omegaq = [load_case(j, k).omega.^2;
%!                 load_case(j, k).omega(1) * load_case(j, k).omega(2);
%!                 load_case(j, k).omega(2) * load_case(j, k).omega(3);
%!                 load_case(j, k).omega(3) * load_case(j, k).omega(1)];
%!       idx = int32(0);
%!       for i=1:numel(omegaq)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * omegaq(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).omegadot)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).omegadot(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).g)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred -= mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).g(i);
%!       endfor
%!       for i=1:numel(load_case(j, k).loads)
%!         l = find(cms_opt.index_KTAU0red == ++idx);
%!         if (isempty(l))
%!           continue;
%!         endif
%!         Kred += mat_ass_cms.KTAU0red(:, :, l) * load_case(j, k).loads(i);
%!       endfor
%!       [sol_eig_red(j, k).Ured, sol_eig_red(j, k).lambda_red] = fem_sol_eigsd(Kred, Dred, Mred, cms_opt.modes.number, cms_opt);
%!     endfor
%!   endfor
%!   tol_abs = [0, 0, 0] / SI_unit_second^-1;
%!   tol_rel = [0.3e-2, 4.5e-2, 3e-2];
%!   tol_disp_rel = 4e-2;
%!   err_u_modal = err_v_modal = zeros(size(param));
%!   printf("deformation/velocity:\n");
%!   colors = rainbow(3);
%!   width = 1:size(res, 3);
%!   linestyle = {"-", "--", "-."};
%!   for i=idx_j
%!     for j=idx_k
%!       u_modal = res(i, j, 1).trajectory{end} - res(i, j, 1).trajectory{end}(1, :);
%!       u_beam = res(i, j, 2).trajectory{end} - res(i, j, 2).trajectory{end}(1, :);
%!       v_modal = res(i, j, 1).velocity{end};
%!       v_beam = res(i, j, 2).velocity{end};
%!       if (options.plot)
%!         figure("visible", "off");
%!         hold on;
%!         for k=1:size(res, 3)
%!           for l=1:3
%!             hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l) - res(i, j, k).trajectory{end}(1, l)) * SI_unit_meter);
%!             set(hnd, "color", colors(l, :));
%!             set(hnd, "linewidth", width(k));
%!             set(hnd, "linestyle", linestyle{k});
%!           endfor
%!         endfor
%!         xlabel("t [s]");
%!         ylabel("u [m]");
%!         grid on;
%!         grid minor on;
%!         title(sprintf("linear displacement %d:%d", i, j));
%!         figure("visible", "off");
%!         hold on;
%!         for k=1:size(res, 3)
%!           for l=1:3
%!             hnd = plot(res(i, j, k).t * SI_unit_second, (res(i, j, k).trajectory{end}(:, l + 3) - res(i, j, k).trajectory{end}(1, l + 3)) * 180 / pi);
%!             set(hnd, "color", colors(l, :));
%!             set(hnd, "linewidth", width(k));
%!             set(hnd, "linestyle", linestyle{k});
%!           endfor
%!         endfor
%!         xlabel("t [s]");
%!         ylabel("Phi [deg]");
%!         grid on;
%!         grid minor on;
%!         title(sprintf("angular displacement %d:%d", i, j));
%!       endif
%!       err_u_modal(i, j) = max(max(abs(u_modal - u_beam))) / max(1, max(max(abs(u_beam))));
%!       err_v_modal(i, j) = max(max(abs(v_modal - v_beam))) / max(1, max(max(abs(v_beam))));
%!       printf("%d:%d %.1f%%/%.1f%%\n", i, j, 100 * err_u_modal(i, j), 100 * err_v_modal(i, j));
%!     endfor
%!   endfor
%!   printf("natural frequencies:\n");
%!   MACR = cell(size(param));
%!   result_data = struct("f_mbd", cell(size(param)), "f_fem", cell(size(param)));
%!   for i=idx_j
%!     for j=idx_k
%!       f_fem = sort(sol_eig(i, j).f(:));
%!       f_fem = f_fem(f_fem > 0);
%!       f_mbd = zeros(rows(f_fem), size(res, 3));
%!       PhiR = zeros(6, rows(f_fem), size(res, 3));
%!       for k=1:size(res, 3)
%!         [f_mbd_k, idx_mbd_k] = sort(res(i, j, k).modal.f(:));
%!         D_mbd_k = res(i, j, k).modal.D(idx_mbd_k);
%!         idx_mbd_k = idx_mbd_k(f_mbd_k > 0);
%!         f_mbd_k = f_mbd_k(f_mbd_k > 0);
%!         idx_mbd_k = idx_mbd_k(1:rows(f_fem));
%!         f_mbd(:, k) = f_mbd_k(1:rows(f_fem));
%!         PhiR(:, :, k) = res(i, j, k).modal.VR(res(i, j, k).modal.idx(end) + (1:6), idx_mbd_k);
%!       endfor
%!       result_data(i, j).f_fem = f_fem;
%!       result_data(i, j).f_mbd = f_mbd;
%!       MACR{i, j} = MACL{i, j} = zeros(rows(f_fem), rows(f_fem));
%!       for k=1:rows(f_fem)
%!         for l=1:rows(f_fem)
%!           MACR{i, j}(k, l) = (PhiR(:, k, 1)' * PhiR(:, k, 2)) * conj(PhiR(:, k, 1)' * PhiR(:, k, 2)) / ((PhiR(:, k, 1)' * PhiR(:, k, 1)) * (PhiR(:, k, 2)' * PhiR(:, k, 2)));
%!         endfor
%!       endfor
%!       printf("%d:%d\n", i, j);
%!       for k=1:rows(f_fem)
%!         printf("%10.2f", f_fem(k) * SI_unit_second^-1);
%!         for l=1:columns(f_mbd)
%!           printf("\t%10.2f", f_mbd(k, l) * SI_unit_second^-1);
%!         endfor
%!         for l=1:columns(f_mbd)
%!           printf("\t%.1f%%", 100 * (f_mbd(k, l) / f_fem(k) - 1));
%!         endfor
%!         printf("\t%.3f", MACR{i, j}(k, k));
%!         fputs(stdout, "\n");
%!       endfor
%!       fputs(stdout, "\n\n");
%!     endfor
%!   endfor
%!   for i=idx_j
%!     for j=idx_k
%!       for k=1:rows(result_data(i, j).f_fem)
%!         for l=1:columns(result_data(i, j).f_mbd)
%!           assert_simple(result_data(i, j).f_mbd(k, l), result_data(i, j).f_fem(k), tol_abs(l) + tol_rel(l) * abs(result_data(i, j).f_fem(k)));
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%!   assert_simple(all(all(err_u_modal < tol_disp_rel)));
%!   for j=idx_j
%!     for k=idx_k
%!       tol = 2e-2;
%!       [lambda_s] = sortrows([imag(sol_eig(j, k).lambda(:)), real(sol_eig(j, k).lambda(:))],[1,2]);
%!       [lambda_red_s] = sortrows([imag(sol_eig_red(j, k).lambda_red(:)), real(sol_eig_red(j, k).lambda_red(:))],[1,2]);
%!       K = min(20, rows(lambda_s));
%!       lambda_s = 1j * lambda_s(:,1) + lambda_s(:, 2);
%!       lambda_red_s = 1j * lambda_red_s(:, 1) + lambda_red_s(:, 2);
%!       assert_simple(lambda_red_s(1:K), lambda_s(1:K), tol * norm(lambda_s(1:K)));
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 15
%! pkg load mboct-fem-pkg;
%! close all;
%! function R2j = norm_ref_frame(R2)
%!   e1 = R2(:, 1);
%!   e3 = [0; 0; 1];
%!   e2 = cross(e3, e1);
%!   e1 = cross(e2, e3);
%!   R2j = [e1 / norm(e1), e2 / norm(e2), e3 / norm(e3)];
%! endfunction
%! function [omega, r] = rotordynamics_test_case(param, options, SI_unit)
%!   cms_opt.algorithm = "shift-invert";
%!   cms_opt.refine_max_iter = int32(3);
%!   cms_opt.element.name = "elem_id_rotor";
%!   cms_opt.nodes.modal.name = "node_id_rotor";
%!   cms_opt.nodes.interfaces(1).name = "node_id_bearing1";
%!   cms_opt.nodes.interfaces(2).name = "node_id_bearing2";
%!   cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!   cms_opt.verbose = options.verbose;
%!   enable_filenames = [options.f_enable_modal, options.f_enable_beam, options.f_enable_solid];
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     mbdyn_filename_suffix = {"cms", "beam", "solid"};
%!     for i=1:numel(mbdyn_filename_suffix)
%!       mbdyn_filenames{i} = [filename, mbdyn_filename_suffix{i}, ".mbdyn"];
%!     endfor
%!     fd = -1;
%!     if (options.f_enable_modal)
%!       unwind_protect
%!         fd = fopen(mbdyn_filenames{1}, "w");
%!         if (fd == -1)
%!           error("failed to open file \"%s\"", mbdyn_filenames{1});
%!         endif
%!         fputs(fd, "include: \"${MBDYN_ROTOR_DYN_CMS_PARAM_FILE}\";\n");
%!         fputs(fd, "set: integer ref_id_ground = 1001;\n");
%!         fputs(fd, "set: integer ref_id_shaft = 1002;\n");
%!         fputs(fd, "set: integer ref_id_rotor = 1003;\n");
%!         fputs(fd, "set: integer ref_id_bearing1 = 1004;\n");
%!         fputs(fd, "set: integer ref_id_bearing2 = 1005;\n");
%!         fputs(fd, "set: integer node_id_rotor = 2002;\n");
%!         fputs(fd, "set: integer node_id_bearing1 = 2003;\n");
%!         fputs(fd, "set: integer node_id_bearing2 = 2004;\n");
%!         fputs(fd, "set: integer body_id_unbalance = 3000;\n");
%!         fputs(fd, "set: integer elem_id_rotor = 3005;\n");
%!         fputs(fd, "set: integer joint_id_bearing1 = 3006;\n");
%!         fputs(fd, "set: integer joint_id_bearing2 = 3007;\n");
%!         fputs(fd, "set: integer joint_id_drive = 3008;\n");
%!         fputs(fd, "set: integer elem_id_inertia = 3010;\n");
%!         fputs(fd, "set: integer drive_id_rotor_speed = 5001;\n");
%!         fputs(fd, "set: integer drive_id_time_step = 5002;\n");
%!         fputs(fd, "set: real initial_time = 0.;\n");
%!         fputs(fd, "set: real final_time = 2. * pi * n / (abs(omega1 + omega0) / 2.);\n");
%!         fputs(fd, "begin: data;\n");
%!         fputs(fd, "        problem: initial value;\n");
%!         fputs(fd, "end: data;\n");
%!         fputs(fd, "begin: initial value;\n");
%!         fputs(fd, "        initial time: initial_time;\n");
%!         fputs(fd, "        final time: final_time;\n");
%!         fputs(fd, "        time step: 2. * pi / (360. * abs(omega0));\n");
%!         fputs(fd, "        strategy: change, postponed, drive_id_time_step;\n");
%!         fputs(fd, "        method: ms, 0.6;\n");
%!         fputs(fd, "        tolerance: 1e-6, test, norm, 1e-6, test, norm;\n");
%!         fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!         fputs(fd, "        max iterations: 100;\n");
%!         fputs(fd, "        derivatives max iterations: 10;\n");
%!         fputs(fd, "        derivatives coefficient: auto;\n");
%!         fputs(fd, "        output: iterations;\n");
%!         fputs(fd, "        linear solver: naive, colamd, scale, row max column max, always, max iterations, 100;\n");
%!         fputs(fd, "        nonlinear solver: nox, modified, 10,\n");
%!         fputs(fd, "             keep jacobian matrix,\n");
%!         fputs(fd, "             inner iterations before assembly, 6,\n");
%!         fputs(fd, "             jacobian operator, newton krylov,\n");
%!         fputs(fd, "             solver, line search based,\n");
%!         fputs(fd, "             forcing term, type 2,\n");
%!         fputs(fd, "             direction, newton,\n");
%!         fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!         fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!         fputs(fd, "             linear solver, gmres,\n");
%!         fputs(fd, "             linear solver max iterations, 12,\n");
%!         fputs(fd, "             krylov subspace size, 12;\n");
%!         fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!         fputs(fd, "        threads: assembly, 1;\n");
%!         fputs(fd, "end: initial value;\n");
%!         fputs(fd, "begin: control data;\n");
%!         fputs(fd, "       default output: none, structural nodes;\n");
%!         fputs(fd, "       default orientation: euler123;\n");
%!         fputs(fd, "       output precision: 16;\n");
%!         fputs(fd, "       max iterations: 0;\n");
%!         fputs(fd, "        structural nodes:\n");
%!         fputs(fd, "                +1           # modal\n");
%!         fputs(fd, "                +1           # interface 1\n");
%!         fputs(fd, "                +1              # interface 2\n");
%!         fputs(fd, "        ;\n");
%!         fputs(fd, "        joints:\n");
%!         fputs(fd, "                +1           # modal\n");
%!         fputs(fd, "                +1           # bearing1\n");
%!         fputs(fd, "                +1           # bearing2\n");
%!         fputs(fd, "                +1              # drive\n");
%!         fputs(fd, "        ;\n");
%!         fputs(fd, "        rigid bodies: 1;\n");
%!         fputs(fd, "end: control data;\n");
%!         fputs(fd, "reference: ref_id_ground,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, eye,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, null;\n");
%!         fputs(fd, "reference: ref_id_shaft,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         fputs(fd, "        reference, ref_id_ground, eye,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         fputs(fd, "        reference, ref_id_ground, 0., 0., omega0;\n");
%!         fputs(fd, "reference: ref_id_rotor,\n");
%!         fputs(fd, "        reference, ref_id_shaft, 0., 0., o,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing1,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   -0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing2,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "begin: nodes;\n");
%!         fputs(fd, "        structural: node_id_rotor, modal,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null;\n");
%!         fputs(fd, "        structural: node_id_bearing1, static,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null;\n");
%!         fputs(fd, "        structural: node_id_bearing2, static,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null;\n");
%!         fputs(fd, "end: nodes;\n");
%!         fputs(fd, "begin: elements;\n");
%!         fputs(fd, "       drive caller: drive_id_rotor_speed, string, \"(omega1 - omega0) / (final_time - initial_time) * Time + omega0\", output, yes;\n");
%!         fputs(fd, "       drive caller: drive_id_time_step, string, \"2. * pi / (max(1., 360. * abs(model::drive(drive_id_rotor_speed, Time))))\";\n");
%!         fputs(fd, "       body: body_id_unbalance,\n");
%!         fputs(fd, "             node_id_rotor,\n");
%!         fputs(fd, "                dm,\n");
%!         fputs(fd, "                reference, ref_id_rotor, dr, 0., 0.,\n");
%!         fputs(fd, "                diag, 0., 0., 0.;\n");
%!         fputs(fd, "        joint: joint_id_bearing1, total pin joint,\n");
%!         fputs(fd, "                node_id_bearing1,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, active,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "        joint: joint_id_bearing2, total pin joint,\n");
%!         fputs(fd, "               node_id_bearing2,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "               position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, inactive,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "     joint: joint_id_drive, angular velocity,\n");
%!         fputs(fd, "             # node label\n");
%!         fputs(fd, "             node_id_rotor, \n");
%!         fputs(fd, "             # direction\n");
%!         fputs(fd, "             0.,0.,1.,\n");
%!         fputs(fd, "             # angular velocity\n");
%!         fputs(fd, "             reference, drive_id_rotor_speed;\n");
%!         fputs(fd, "        include: \"${MBDYN_ROTOR_DYN_CMS_ELEM_FILE}\";\n");
%!         fputs(fd, "        inertia: elem_id_inertia,\n");
%!         fputs(fd, "                 position, reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                 orientation, reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                 body, all, joint, all, loadable, all,\n");
%!         fputs(fd, "                 output, both;\n");
%!         fputs(fd, "end: elements;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!     endif
%!     fd = -1;
%!     if (options.f_enable_beam)
%!       unwind_protect
%!         fd = fopen(mbdyn_filenames{2}, "w");
%!         if (fd == -1)
%!           error("failed to open file \"%s\"", mbdyn_filenames{2});
%!         endif
%!         fputs(fd, "include: \"${MBDYN_ROTOR_DYN_CMS_PARAM_FILE}\";\n");
%!         fputs(fd, "set: integer ref_id_ground = 1001;\n");
%!         fputs(fd, "set: integer ref_id_shaft = 1002;\n");
%!         fputs(fd, "set: integer ref_id_rotor = 1003;\n");
%!         fputs(fd, "set: integer ref_id_bearing1 = 1004;\n");
%!         fputs(fd, "set: integer ref_id_bearing2 = 1005;\n");
%!         fputs(fd, "set: integer ref_id_shaft_section = 1006;\n");
%!         fputs(fd, "set: integer ref_id_bearing1_center = 1007;\n");
%!         fputs(fd, "set: integer ref_id_bearing2_center = 1008;\n");
%!         fputs(fd, "set: integer node_id_rotor = 2001;\n");
%!         fputs(fd, "set: integer node_id_bearing1 = 2002;\n");
%!         fputs(fd, "set: integer node_id_bearing1_center = 2003;\n");
%!         fputs(fd, "set: integer node_id_bearing2 = 2004;\n");
%!         fputs(fd, "set: integer node_id_bearing2_center = 2005;\n");
%!         fputs(fd, "set: integer body_id_rotor = 3001;\n");
%!         fputs(fd, "set: integer body_id_bearing1 = 3002;\n");
%!         fputs(fd, "set: integer body_id_bearing1_center = 3003;\n");
%!         fputs(fd, "set: integer body_id_bearing2 = 3004;\n");
%!         fputs(fd, "set: integer body_id_bearing2_center = 3005;\n");
%!         fputs(fd, "set: integer beam_id_shaft1 = 4001;\n");
%!         fputs(fd, "set: integer beam_id_shaft2 = 4002;\n");
%!         fputs(fd, "set: integer joint_id_bearing1 = 3006;\n");
%!         fputs(fd, "set: integer joint_id_bearing2 = 3007;\n");
%!         fputs(fd, "set: integer joint_id_drive = 3008;\n");
%!         fputs(fd, "set: integer elem_id_inertia = 4001;\n");
%!         fputs(fd, "set: integer drive_id_rotor_speed = 5001;\n");
%!         fputs(fd, "set: integer drive_id_time_step = 5002;\n");
%!         fputs(fd, "set: real l1 = 0.5 * l + o - 0.5 * w;\n");
%!         fputs(fd, "set: real l2 = 0.5 * l - o - 0.5 * w;\n");
%!         fputs(fd, "set: real rotor_m = rho * D^2 * pi / 4. * w;\n");
%!         fputs(fd, "set: real rotor_Jx = rotor_m * ((0.5 * D)^2 + w^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real rotor_Jy = rotor_Jx;\n");
%!         fputs(fd, "set: real rotor_Jz = rotor_m * (0.5 * D)^2 / 2.;\n");
%!         fputs(fd, "set: real shaft_rotor1_m = rho * d^2 * pi / 4. * (l1 / 4.);\n");
%!         fputs(fd, "set: real shaft_rotor1_Jx = shaft_rotor1_m * ((0.5 * d)^2 + (l1 / 4.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_rotor1_Jy = shaft_rotor1_Jx;\n");
%!         fputs(fd, "set: real shaft_rotor1_Jz = shaft_rotor1_m * (0.5 * d)^2 / 2.;\n");
%!         fputs(fd, "set: real shaft_rotor2_m = rho * d^2 * pi / 4. * (l2 / 4.);\n");
%!         fputs(fd, "set: real shaft_rotor2_Jx = shaft_rotor2_m * ((0.5 * d)^2 + (l2 / 4.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_rotor2_Jy = shaft_rotor2_Jx;\n");
%!         fputs(fd, "set: real shaft_rotor2_Jz = shaft_rotor2_m * (0.5 * d)^2 / 2.;\n");
%!         fputs(fd, "set: real shaft_bearing1_m = rho * d^2 * pi / 4. * (l1 / 4.);\n");
%!         fputs(fd, "set: real shaft_bearing1_Jx = shaft_bearing1_m * ((0.5 * d)^2 + (l1 / 4.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_bearing1_Jy = shaft_bearing1_Jx;\n");
%!         fputs(fd, "set: real shaft_bearing1_Jz = shaft_bearing1_m * (0.5 * d)^2 / 2.;\n");
%!         fputs(fd, "set: real shaft_bearing1_center_m = rho * d^2 * pi / 4. * (l1 / 2.);\n");
%!         fputs(fd, "set: real shaft_bearing1_center_Jx = shaft_bearing1_center_m * ((0.5 * d)^2 + (l1 / 2.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_bearing1_center_Jy = shaft_bearing1_center_Jx;\n");
%!         fputs(fd, "set: real shaft_bearing1_center_Jz = shaft_bearing1_center_m * ((0.5 * d)^2) / 2.;\n");
%!         fputs(fd, "set: real shaft_bearing2_m = rho * d^2 * pi / 4. * (l2 / 4.);\n");
%!         fputs(fd, "set: real shaft_bearing2_Jx = shaft_bearing2_m * ((0.5 * d)^2 + (l2 / 4.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_bearing2_Jy = shaft_bearing2_Jx;\n");
%!         fputs(fd, "set: real shaft_bearing2_Jz = shaft_bearing2_m * (0.5 * d)^2 / 2.;\n");
%!         fputs(fd, "set: real shaft_bearing2_center_m = rho * d^2 * pi / 4. * (l2 / 2.);\n");
%!         fputs(fd, "set: real shaft_bearing2_center_Jx = shaft_bearing2_center_m * ((0.5 * d)^2 + (l2 / 2.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_bearing2_center_Jy = shaft_bearing2_center_Jx;\n");
%!         fputs(fd, "set: real shaft_bearing2_center_Jz = shaft_bearing2_center_m * ((0.5 * d)^2) / 2.;\n");
%!         fputs(fd, "set: real shaft_A = d^2 * pi / 4;\n");
%!         fputs(fd, "set: real shaft_As = 9. / 10. * shaft_A;\n");
%!         fputs(fd, "set: real shaft_Iy = d^4 * pi / 64.;\n");
%!         fputs(fd, "set: real shaft_Iz = shaft_Iy;\n");
%!         fputs(fd, "set: real shaft_Ip = shaft_Iy + shaft_Iz;\n");
%!         fputs(fd, "set: real shaft_It = shaft_Ip;\n");
%!         fputs(fd, "set: real shaft_E = E;\n");
%!         fputs(fd, "set: real shaft_nu = nu;\n");
%!         fputs(fd, "set: real shaft_G = shaft_E / (2 * (1 + shaft_nu));\n");
%!         fputs(fd, "set: real shaft_rho = rho;\n");
%!         fputs(fd, "set: real shaft_damping_ratio = beta;\n");
%!         fputs(fd, "set: real initial_time = 0.;\n");
%!         fputs(fd, "set: real final_time = 2. * pi * n / (abs(omega1 + omega0) / 2.);\n");
%!         fputs(fd, "begin: data;\n");
%!         fputs(fd, "        problem: initial value;\n");
%!         fputs(fd, "end: data;\n");
%!         fputs(fd, "begin: initial value;\n");
%!         fputs(fd, "        initial time: initial_time;\n");
%!         fputs(fd, "        final time: final_time;\n");
%!         fputs(fd, "        time step: 2. * pi / (360. * abs(omega0));\n");
%!         fputs(fd, "        strategy: change, postponed, drive_id_time_step;\n");
%!         fputs(fd, "        method: ms, 0.6;\n");
%!         fputs(fd, "        tolerance: 1e-6, test, norm, 1e-6, test, norm;\n");
%!         fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!         fputs(fd, "        max iterations: 100;\n");
%!         fputs(fd, "        derivatives max iterations: 10;\n");
%!         fputs(fd, "        derivatives coefficient: auto;\n");
%!         fputs(fd, "        output: iterations;\n");
%!         fputs(fd, "        linear solver: naive, colamd, scale, row max column max, always, max iterations, 100;\n");
%!         fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!         fputs(fd, "        threads: assembly, 1;\n");
%!         fputs(fd, "        nonlinear solver: nox, modified, 10,\n");
%!         fputs(fd, "             keep jacobian matrix,\n");
%!         fputs(fd, "             inner iterations before assembly, 6,\n");
%!         fputs(fd, "             jacobian operator, newton krylov,\n");
%!         fputs(fd, "             solver, line search based,\n");
%!         fputs(fd, "             forcing term, type 2,\n");
%!         fputs(fd, "             direction, newton,\n");
%!         fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!         fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!         fputs(fd, "             linear solver, gmres,\n");
%!         fputs(fd, "             linear solver max iterations, 12,\n");
%!         fputs(fd, "             krylov subspace size, 12;\n");
%!         fputs(fd, "end: initial value;\n");
%!         fputs(fd, "begin: control data;\n");
%!         fputs(fd, "       use automatic differentiation;\n");
%!         fputs(fd, "       default output: none, structural nodes;\n");
%!         fputs(fd, "       default orientation: euler123;\n");
%!         fputs(fd, "       output precision: 16;\n");
%!         fputs(fd, "       max iterations: 0;\n");
%!         fputs(fd, "        structural nodes: 5;\n");
%!         fputs(fd, "        joints: 3;\n");
%!         fputs(fd, "        beams: 2;\n");
%!         fputs(fd, "        rigid bodies: 5;\n");
%!         if (options.f_rbk)
%!           fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!           fputs(fd, "        component,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           array, 2,\n");
%!           fputs(fd, "             mult, time, const, (omega1 - omega0) / (final_time - initial_time),\n");
%!           fputs(fd, "             const, omega0,\n");
%!           fputs(fd, "        angular acceleration,\n");
%!           fputs(fd, "        component,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           const, (omega1 - omega0) / (final_time - initial_time);\n");
%!         endif
%!         fputs(fd, "end: control data;\n");
%!         fputs(fd, "reference: ref_id_ground,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, eye,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, null;\n");
%!         fputs(fd, "reference: ref_id_shaft,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         fputs(fd, "        reference, ref_id_ground, eye,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         if (options.f_rbk)
%!           fputs(fd, "        reference, ref_id_ground, null;\n");
%!         else
%!           fputs(fd, "        reference, ref_id_ground, 0., 0., omega0;\n");
%!         endif
%!         fputs(fd, "reference: ref_id_rotor,\n");
%!         fputs(fd, "        reference, ref_id_shaft, 0., 0., o,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_shaft_section,\n");
%!         fputs(fd, "           reference, ref_id_shaft, null,\n");
%!         fputs(fd, "           reference, ref_id_shaft, 1, 0., 0., 1.,\n");
%!         fputs(fd, "                                    2, 0., 1., 0.,\n");
%!         fputs(fd, "           reference, ref_id_shaft, null,\n");
%!         fputs(fd, "           reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing1,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   -0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing1_center,\n");
%!         fputs(fd, "        reference, ref_id_bearing1,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.5 * l1,\n");
%!         fputs(fd, "        reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "        reference, ref_id_rotor, null,\n");
%!         fputs(fd, "        reference, ref_id_rotor, null;\n");
%!         fputs(fd, "reference: ref_id_bearing2,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing2_center,\n");
%!         fputs(fd, "        reference, ref_id_bearing2,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   -0.5 * l2,\n");
%!         fputs(fd, "        reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "        reference, ref_id_rotor, null,\n");
%!         fputs(fd, "        reference, ref_id_rotor, null;\n");
%!         fputs(fd, "begin: nodes;\n");
%!         fputs(fd, "        structural: node_id_rotor, dynamic,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null;\n");
%!         fputs(fd, "        structural: node_id_bearing1, dynamic,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null;\n");
%!         fputs(fd, "        structural: node_id_bearing1_center, dynamic,\n");
%!         fputs(fd, "                reference, ref_id_bearing1_center, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1_center, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing1_center, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1_center, null;\n");
%!         fputs(fd, "        structural: node_id_bearing2, dynamic,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null;\n");
%!         fputs(fd, "        structural: node_id_bearing2_center, dynamic,\n");
%!         fputs(fd, "                reference, ref_id_bearing2_center, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2_center, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing2_center, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2_center, null;\n");
%!         fputs(fd, "end: nodes;\n");
%!         fputs(fd, "begin: elements;\n");
%!         fputs(fd, "       drive caller: drive_id_rotor_speed, string, \"(omega1 - omega0) / (final_time - initial_time) * Time + omega0\", output, yes;\n");
%!         fputs(fd, "       drive caller: drive_id_time_step, string, \"2. * pi / (max(1., 360. * abs(model::drive(drive_id_rotor_speed, Time))))\";\n");
%!         fputs(fd, "       body: body_id_rotor,\n");
%!         fputs(fd, "             node_id_rotor,\n");
%!         fputs(fd, "             condense, 4,\n");
%!         fputs(fd, "                rotor_m,\n");
%!         fputs(fd, "                  reference, node, null,\n");
%!         fputs(fd, "                  diag, rotor_Jx, rotor_Jy, rotor_Jz,\n");
%!         fputs(fd, "                shaft_rotor1_m,\n");
%!         fputs(fd, "                  reference, node, 0., 0., -w / 2. - l1 / 8.,\n");
%!         fputs(fd, "                  diag, shaft_rotor1_Jx, shaft_rotor1_Jy, shaft_rotor1_Jz,\n");
%!         fputs(fd, "                shaft_rotor2_m,\n");
%!         fputs(fd, "                  reference, node, 0., 0., w / 2. + l2 / 8.,\n");
%!         fputs(fd, "                  diag, shaft_rotor2_Jx, shaft_rotor2_Jy, shaft_rotor2_Jz,                \n");
%!         fputs(fd, "                dm,\n");
%!         fputs(fd, "                  reference, node, dr, 0., 0.,\n");
%!         fputs(fd, "                  diag, 0., 0., 0.;\n");
%!         fputs(fd, "        body: body_id_bearing1,\n");
%!         fputs(fd, "              node_id_bearing1,\n");
%!         fputs(fd, "              shaft_bearing1_m,\n");
%!         fputs(fd, "              reference, node, 0., 0., l1 / 8.,\n");
%!         fputs(fd, "              diag, shaft_bearing1_Jx, shaft_bearing1_Jy, shaft_bearing1_Jz;\n");
%!         fputs(fd, "        body: body_id_bearing1_center,\n");
%!         fputs(fd, "              node_id_bearing1_center,\n");
%!         fputs(fd, "              shaft_bearing1_center_m,\n");
%!         fputs(fd, "              reference, node, null,\n");
%!         fputs(fd, "              diag, shaft_bearing1_center_Jx, shaft_bearing1_center_Jy, shaft_bearing1_center_Jz;\n");
%!         fputs(fd, "        body: body_id_bearing2,\n");
%!         fputs(fd, "              node_id_bearing2,\n");
%!         fputs(fd, "              shaft_bearing2_m,\n");
%!         fputs(fd, "              reference, node, 0., 0., -l2 / 8.,\n");
%!         fputs(fd, "              diag, shaft_bearing2_Jx, shaft_bearing2_Jy, shaft_bearing2_Jz;\n");
%!         fputs(fd, "        body: body_id_bearing2_center,\n");
%!         fputs(fd, "              node_id_bearing2_center,\n");
%!         fputs(fd, "              shaft_bearing2_center_m,\n");
%!         fputs(fd, "              reference, node, null,\n");
%!         fputs(fd, "              diag, shaft_bearing2_center_Jx, shaft_bearing2_center_Jy, shaft_bearing2_center_Jz;\n");
%!         fputs(fd, "    beam3: beam_id_shaft1,\n");
%!         fputs(fd, "                # node 1\n");
%!         fputs(fd, "                node_id_bearing1, position, reference, node, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # node 2\n");
%!         fputs(fd, "                node_id_bearing1_center, position, reference, node, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # node 3,\n");
%!         fputs(fd, "                node_id_rotor, position, reference, node, 0., 0., -0.5 * w,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # orientation matrix section I\n");
%!         fputs(fd, "                reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # constitutive law section I\n");
%!         fputs(fd, "                linear viscoelastic generic,\n");
%!         fputs(fd, "                diag, shaft_E * shaft_A , shaft_G * shaft_As, shaft_G * shaft_As,\n");
%!         fputs(fd, "                      shaft_G * shaft_It, shaft_E * shaft_Iy, shaft_E * shaft_Iz,\n");
%!         fputs(fd, "                proportional, shaft_damping_ratio,\n");
%!         fputs(fd, "                # orientation matrix section II\n");
%!         fputs(fd, "                same,\n");
%!         fputs(fd, "                # constitutive law section II\n");
%!         fputs(fd, "                same;\n");
%!         fputs(fd, "    beam3: beam_id_shaft2,\n");
%!         fputs(fd, "                # node 1\n");
%!         fputs(fd, "                node_id_rotor, position, reference, node, 0., 0., 0.5 * w,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,                \n");
%!         fputs(fd, "                # node 2\n");
%!         fputs(fd, "                node_id_bearing2_center, position, reference, node, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # node 3,\n");
%!         fputs(fd, "                node_id_bearing2, position, reference, node, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # orientation matrix section I\n");
%!         fputs(fd, "                reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # constitutive law section I\n");
%!         fputs(fd, "                linear viscoelastic generic,\n");
%!         fputs(fd, "                diag, shaft_E * shaft_A , shaft_G * shaft_As, shaft_G * shaft_As,\n");
%!         fputs(fd, "                      shaft_G * shaft_It, shaft_E * shaft_Iy, shaft_E * shaft_Iz,\n");
%!         fputs(fd, "                proportional, shaft_damping_ratio,\n");
%!         fputs(fd, "                # orientation matrix section II\n");
%!         fputs(fd, "                same,\n");
%!         fputs(fd, "                # constitutive law section II\n");
%!         fputs(fd, "                same;\n");
%!         fputs(fd, "        joint: joint_id_bearing1, total pin joint,\n");
%!         fputs(fd, "                node_id_bearing1,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, active,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "        joint: joint_id_bearing2, total pin joint,\n");
%!         fputs(fd, "               node_id_bearing2,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "               position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, inactive,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "     joint: joint_id_drive, angular velocity,\n");
%!         fputs(fd, "             # node label\n");
%!         fputs(fd, "             node_id_rotor, \n");
%!         fputs(fd, "             # direction\n");
%!         fputs(fd, "             0.,0.,1.,\n");
%!         fputs(fd, "             # angular velocity\n");
%!         if (options.f_rbk)
%!           fputs(fd, "           null;\n");
%!         else
%!           fputs(fd, "           reference, drive_id_rotor_speed;\n");
%!         endif
%!         fputs(fd, "        inertia: elem_id_inertia,\n");
%!         fputs(fd, "                 position, reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                 orientation, reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                 body, all, joint, all, loadable, all,\n");
%!         fputs(fd, "                 output, both;\n");
%!         fputs(fd, "end: elements;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!     endif
%!     fd = -1;
%!     if (options.f_enable_modal || options.f_enable_solid)
%!       geometry_file = [filename, ".geo"];
%!       unwind_protect
%!         [fd, msg] = fopen(geometry_file, "w");
%!         if (fd == -1)
%!           error("failed to open file \"%s\"", geometry_file);
%!         endif
%!         fn = fieldnames(param);
%!         for i=1:length(fn)
%!           fprintf(fd, "%s = %g;\n", fn{i}, getfield(param, fn{i}));
%!         endfor
%!         fputs(fd, "SetFactory(\"Built-in\");\n");
%!         fputs(fd, "Point(1) = {0, 0, -0.5 * l};\n");
%!         fputs(fd, "Point(2) = {0.5 * d, 0, -0.5 * l};\n");
%!         fputs(fd, "Point(3) = {0.5 * d, 0, -0.5 * w + o};\n");
%!         fputs(fd, "Point(4) = {0.5 * D, 0, -0.5 * w + o};\n");
%!         fputs(fd, "Point(5) = {0.5 * D, 0, 0.5 * w + o};\n");
%!         fputs(fd, "Point(6) = {0.5 * d, 0, 0.5 * w + o};\n");
%!         fputs(fd, "Point(7) = {0.5 * d, 0, 0.5 * l};\n");
%!         fputs(fd, "Point(8) = {0, 0, 0.5 * l};\n");
%!         fputs(fd, "Point(9) = {0, 0, 0.5 * w + o};\n");
%!         fputs(fd, "Point(10) = {0, 0, -0.5 * w + o};\n");
%!         fputs(fd, "Line(1) = {1, 2};\n");
%!         fputs(fd, "Line(2) = {2, 3};\n");
%!         fputs(fd, "Line(3) = {3, 4};\n");
%!         fputs(fd, "Line(4) = {4, 5};\n");
%!         fputs(fd, "Line(5) = {5, 6};\n");
%!         fputs(fd, "Line(6) = {6, 7};\n");
%!         fputs(fd, "Line(7) = {7, 8};\n");
%!         fputs(fd, "Line(8) = {8, 9};\n");
%!         fputs(fd, "Line(9) = {9, 10};\n");
%!         fputs(fd, "Line(10) = {10, 1};\n");
%!         fputs(fd, "Line(11) = {10,3};\n");
%!         fputs(fd, "Line(12) = {3,6};\n");
%!         fputs(fd, "Line(13) = {6,9};\n");
%!         fputs(fd, "Transfinite Curve(1) = Max(1, Round(0.5 * d / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(2) = Max(1, Round((0.5 * l - 0.5 * w + o) / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(3) = Max(1, Round((0.5 * (D - d)) / h2)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(4) = Max(1, Round(w / h2)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(5) = Max(1, Round((0.5 * (D - d)) / h2)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(6) = Max(1, Round((0.5 * l - 0.5 * w - o) / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(7) = Max(1, Round(0.5 * d / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(8) = Max(1, Round((0.5 * l - 0.5 * w - o) / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(9) = Max(1, Round(w / h2)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(10) = Max(1, Round((0.5 * l - 0.5 * w + o) / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(11) = Max(1, Round(0.5 * d / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(12) = Max(1, Round(w / h2)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(13) = Max(1, Round(0.5 * d / h1)) + 1;\n");
%!         fputs(fd, "Line Loop(1) = {1, 2, -11, 10};\n");
%!         fputs(fd, "Line Loop(2) = {11, 12, 13, 9};\n");
%!         fputs(fd, "Line Loop(3) = {3, 4, 5, -12};\n");
%!         fputs(fd, "Line Loop(4) = {-13, 6, 7, 8};\n");
%!         fputs(fd, "Plane Surface(1) = {1};\n");
%!         fputs(fd, "Plane Surface(2) = {2};\n");
%!         fputs(fd, "Plane Surface(3) = {3};\n");
%!         fputs(fd, "Plane Surface(4) = {4};\n");
%!         fputs(fd, "Transfinite Surface(1) = {1, 2, 3, 10};\n");
%!         fputs(fd, "Transfinite Surface(2) = {10, 3, 6, 9};\n");
%!         fputs(fd, "Transfinite Surface(3) = {3, 4, 5, 6};\n");
%!         fputs(fd, "Transfinite Surface(4) = {9, 6, 7, 8};\n");
%!         fputs(fd, "vol11[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{1}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol21[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol11[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol31[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol21[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol41[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol31[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "Recombine Surface{1, vol11[0]};\n");
%!         fputs(fd, "Recombine Surface{vol11[0], vol21[0]};\n");
%!         fputs(fd, "Recombine Surface{vol21[0], vol31[0]};\n");
%!         fputs(fd, "Recombine Surface{vol31[0], vol41[0]};\n");
%!         fputs(fd, "vol12[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{2}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol22[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol12[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol32[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol22[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol42[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol32[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "Recombine Surface{2, vol12[0]};\n");
%!         fputs(fd, "Recombine Surface{vol12[0], vol22[0]};\n");
%!         fputs(fd, "Recombine Surface{vol22[0], vol32[0]};\n");
%!         fputs(fd, "Recombine Surface{vol32[0], vol42[0]};\n");
%!         fputs(fd, "vol13[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{3}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol23[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol13[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol33[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol23[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol43[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol33[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "Recombine Surface{3, vol13[0]};\n");
%!         fputs(fd, "Recombine Surface{vol13[0], vol23[0]};\n");
%!         fputs(fd, "Recombine Surface{vol23[0], vol33[0]};\n");
%!         fputs(fd, "Recombine Surface{vol33[0], vol43[0]};\n");
%!         fputs(fd, "vol14[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{4}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol24[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol14[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol34[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol24[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol44[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol34[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "Recombine Surface{4, vol14[0]};\n");
%!         fputs(fd, "Recombine Surface{vol14[0], vol24[0]};\n");
%!         fputs(fd, "Recombine Surface{vol24[0], vol34[0]};\n");
%!         fputs(fd, "Recombine Surface{vol34[0], vol44[0]};\n");
%!         fputs(fd, "Coherence;\n");
%!         switch (options.elem_type)
%!           case "iso8"
%!             fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!             fputs(fd, "Mesh.ElementOrder = 1;");
%!           case "iso20"
%!             fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!             fputs(fd, "Mesh.ElementOrder = 2;");
%!         endswitch
%!         fputs(fd, "Mesh 3;\n");
%!         fputs(fd, "Coherence Mesh;\n");
%!         fputs(fd, "Physical Surface(\"bearing1\", 1) = {21, 38, 55, 72};\n");
%!         fputs(fd, "Physical Surface(\"bearing2\", 2) = {260, 243, 277, 294};\n");
%!         fputs(fd, "Physical Surface(\"rotor\", 3) = {202, 180, 224, 158};\n");
%!         fputs(fd, "Physical Volume(\"shaft\", 1) = {1, 2, 4, 3, 14, 13, 15, 16, 5, 6, 8, 7};\n");
%!         fputs(fd, "Physical Volume(\"disc\", 2) = {9, 10, 12, 11};\n");
%!         fputs(fd, "Mesh.Format = 1;\n");
%!         fprintf(fd, "Save \"%s.msh\";\n", filename);
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!       ##spawn_wait(spawn("gmsh", {geometry_file})); return;
%!       fprintf(stderr, "meshing ...\n");
%!       pid = spawn("gmsh", {"-format", "msh2", "-0",geometry_file});
%!       status = spawn_wait(pid);
%!       if (status ~= 0)
%!         warning("gmsh failed with status %d", status);
%!       endif
%!       fprintf(stderr, "loading mesh ...\n");
%!       switch (options.elem_type)
%!         case "iso20"
%!           opt_mesh.elem_type = {"iso20", "penta15", "quad8", "tria6h"};
%!         case "iso8"
%!           opt_mesh.elem_type = {"iso8", "penta6", "iso4", "tria3"};
%!       endswitch
%!       mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh",opt_mesh));
%!       cms_opt.modes.number = int32(options.number_of_modes);
%!       cms_opt.nodes.modal.number = int32(rows(mesh.nodes) + 1);
%!       cms_opt.nodes.interfaces(1).number = int32(rows(mesh.nodes) + 2);
%!       cms_opt.nodes.interfaces(2).number = int32(rows(mesh.nodes) + 3);
%!       cms_opt.invariants = true;
%!       mesh.nodes(cms_opt.nodes.modal.number, :) = [0, 0, param.o, 0, 0, 0];
%!       mesh.nodes([cms_opt.nodes.interfaces.number], :) = [0, 0, -0.5 * param.l, 0, 0, 0;
%!                                                           0, 0,  0.5 * param.l, 0, 0, 0];
%!       switch (options.elem_type)
%!         case "iso20"
%!           mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 1, cms_opt.nodes.interfaces(1).number, "tria6h");
%!           mesh.elements.rbe3(3) = fem_pre_mesh_rbe3_from_surf(mesh, 2, cms_opt.nodes.interfaces(2).number, "tria6h");
%!           mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 3, cms_opt.nodes.modal.number, "quad8");
%!         case "iso8"
%!           mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 1, cms_opt.nodes.interfaces(1).number, "iso4");
%!           mesh.elements.rbe3(3) = fem_pre_mesh_rbe3_from_surf(mesh, 2, cms_opt.nodes.interfaces(2).number, "iso4");
%!           mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 3, cms_opt.nodes.modal.number, "iso4");
%!       endswitch
%!       load_case_dof.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
%!       switch (options.elem_type)
%!         case "iso20"
%!           mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!           mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!           mesh.materials.iso20([mesh.groups.iso20(find([mesh.groups.iso20.id] == 1)).elements]) = 1;
%!           mesh.materials.iso20([mesh.groups.iso20(find([mesh.groups.iso20.id] == 2)).elements]) = 2;
%!           mesh.materials.penta15([mesh.groups.penta15(find([mesh.groups.penta15.id] == 1)).elements]) = 1;
%!           mesh.materials.penta15([mesh.groups.penta15(find([mesh.groups.penta15.id] == 2)).elements]) = 2;
%!         case "iso8"
%!           mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!           mesh.materials.iso8([mesh.groups.iso8(find([mesh.groups.iso8.id] == 1)).elements]) = 1;
%!           mesh.materials.iso8([mesh.groups.iso8(find([mesh.groups.iso8.id] == 2)).elements]) = 2;
%!       endswitch
%!       mesh.material_data(1).rho = param.rho;
%!       mesh.material_data(1).E = param.E;
%!       mesh.material_data(1).nu = param.nu;
%!       mesh.material_data(2).rho = param.rho;
%!       mesh.material_data(2).E = 100 * param.E;
%!       mesh.material_data(2).nu = param.nu;
%!     endif
%!     if (options.f_enable_modal)
%!       fprintf(stderr, "building cms element ...\n");
%!       [mesh_cms, mat_ass, dof_map, sol_eig, cms_opt] = fem_cms_create2(mesh, load_case_dof, cms_opt);
%!       mesh_post_pro_file = sprintf("%s_post.msh", filename);
%!       fem_post_mesh_export(mesh_post_pro_file, mesh_cms);
%!       for j=1:size(sol_eig.def, 3)
%!         eig_post_pro_file_mode{j} = sprintf("%s_eig_def_%03d.msh", filename, j);
%!         fem_post_sol_step_export(eig_post_pro_file_mode{j}, sol_eig, j, j, sol_eig.f(j), options.scale_def / max(norm(sol_eig.def(:, 1:3, j), "rows")));
%!       endfor
%!       eig_post_pro_file = sprintf("%s_modes_post.geo", filename);
%!       fd = -1;
%!       unwind_protect
%!         [fd, msg] = fopen(eig_post_pro_file, "w");
%!         if (fd == -1)
%!           error("failed to open file \"%s\"", eig_post_pro_file);
%!         endif
%!         fprintf(fd, "Merge \"%s\";\n", mesh_post_pro_file);
%!         for j=1:numel(eig_post_pro_file_mode)
%!           fprintf(fd, "Merge \"%s\";\n", eig_post_pro_file_mode{j});
%!         endfor
%!         fputs(fd, "View.Type = 1;\n");
%!         fputs(fd, "View.VectorType = 5;\n");
%!         fputs(fd, "View.Visible = 1;\n");
%!         fputs(fd, "View.DisplacementFactor = 1;\n");
%!         fputs(fd, "View.ShowTime = 6;\n");
%!         fputs(fd, "View.ShowElement = 1;\n");
%!         fputs(fd, "View.IntervalsType = 3;\n");
%!         fputs(fd, "View.NbIso = 20;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!       mat_ass.Dred = param.alpha * mat_ass.Mred + param.beta * mat_ass.Kred;
%!       fem_cms_export([filename, "_cms"], mesh_cms, dof_map, mat_ass, cms_opt);
%!     endif
%!     if (options.f_enable_solid)
%!       opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!       opt_mbd_mesh.struct_nodes.type([cms_opt.nodes.modal.number]) = MBDYN_NODE_TYPE_DYNAMIC_STRUCT;
%!       opt_mbd_mesh.struct_nodes.type([cms_opt.nodes.interfaces.number]) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!       opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_shaft";
%!       solid_csl_file = [filename, "_solid.csl"];
%!       solid_nodes_file = [filename, "_solid.nod"];
%!       solid_elem_file = [filename, "_solid.elm"];
%!       load_case = struct();
%!       assert_simple(mesh.nodes(cms_opt.nodes.modal.number, 1:3), [0, 0, param.o]);
%!       assert_simple(mesh.nodes(cms_opt.nodes.interfaces(1).number, 1:3), [0, 0, -0.5 * param.l]);
%!       assert_simple(mesh.nodes(cms_opt.nodes.interfaces(2).number, 1:3), [0, 0, 0.5 * param.l]);
%!       opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, solid_csl_file, opt_mbd_mesh);
%!       opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, solid_nodes_file, opt_mbd_mesh);
%!       opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, solid_elem_file, opt_mbd_mesh);
%!       fd = -1;
%!       unwind_protect
%!         fd = fopen(mbdyn_filenames{3}, "w");
%!         if (fd == -1)
%!           error("failed to open file \"%s\"", mbdyn_filenames{2});
%!         endif
%!         fputs(fd, "include: \"${MBDYN_ROTOR_DYN_CMS_PARAM_FILE}\";\n");
%!         fputs(fd, "set: integer ref_id_ground = 1001;\n");
%!         fputs(fd, "set: integer ref_id_shaft = 1002;\n");
%!         fputs(fd, "set: integer ref_id_rotor = 1003;\n");
%!         fputs(fd, "set: integer ref_id_bearing1 = 1004;\n");
%!         fputs(fd, "set: integer ref_id_bearing2 = 1005;\n");
%!         fprintf(fd, "set: integer node_id_rotor = %d;\n", cms_opt.nodes.modal.number);
%!         fprintf(fd, "set: integer node_id_bearing1 = %d;\n", cms_opt.nodes.interfaces(1).number);
%!         fprintf(fd, "set: integer node_id_bearing2 = %d;\n", cms_opt.nodes.interfaces(2).number);
%!         fputs(fd, "set: integer body_id_unbalance = 2001;\n");
%!         fputs(fd, "set: integer joint_id_bearing1 = 3006;\n");
%!         fputs(fd, "set: integer joint_id_bearing2 = 3007;\n");
%!         fputs(fd, "set: integer joint_id_drive = 3008;\n");
%!         fputs(fd, "set: integer elem_id_inertia = 4001;\n");
%!         fputs(fd, "set: integer drive_id_rotor_speed = 5001;\n");
%!         fputs(fd, "set: integer drive_id_time_step = 5002;\n");
%!         fputs(fd, "set: real initial_time = 0.;\n");
%!         fputs(fd, "set: real final_time = 2. * pi * n / (abs(omega1 + omega0) / 2.);\n");
%!         fputs(fd, "begin: data;\n");
%!         fputs(fd, "        problem: initial value;\n");
%!         fputs(fd, "end: data;\n");
%!         fputs(fd, "begin: initial value;\n");
%!         fputs(fd, "        initial time: initial_time;\n");
%!         fputs(fd, "        final time: final_time;\n");
%!         fputs(fd, "        time step: 2. * pi / (360. * abs(omega0));\n");
%!         fputs(fd, "        strategy: change, postponed, drive_id_time_step;\n");
%!         fputs(fd, "        method: ms, 0.6;\n");
%!         fputs(fd, "        tolerance: 1e-5, test, sepnorm, 1e-5, test, norm;\n");
%!         fputs(fd, "        derivatives tolerance: 1e-3, 1e-3;\n");
%!         fputs(fd, "        max iterations: 100;\n");
%!         fputs(fd, "        derivatives max iterations: 10;\n");
%!         fputs(fd, "        derivatives coefficient: auto;\n");
%!         fputs(fd, "        output: iterations, cpu time, solver condition number, stat, yes;\n");
%!         fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!         fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!         fprintf(fd, "      threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!         fprintf(fd, "      threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!         fputs(fd, "        nonlinear solver: nox, modified, 25,\n");
%!         fputs(fd, "             keep jacobian matrix,\n");
%!         fputs(fd, "             inner iterations before assembly, 12,\n");
%!         fputs(fd, "             use preconditioner as solver, no,\n");
%!         fputs(fd, "             jacobian operator, newton krylov,\n");
%!         fputs(fd, "             solver, line search based,\n");
%!         fputs(fd, "             forcing term, type 2,\n");
%!         fputs(fd, "             forcing term min tolerance, 1e-10,\n");
%!         fputs(fd, "             forcing term max tolerance, 1e-8,\n");
%!         fputs(fd, "             direction, newton,\n");
%!         fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!         fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!         fputs(fd, "             linear solver, gmres,\n");
%!         fputs(fd, "             linear solver max iterations, 24,\n");
%!         fputs(fd, "             krylov subspace size, 24;\n");
%!         fputs(fd, "end: initial value;\n");
%!         fputs(fd, "begin: control data;\n");
%!         fputs(fd, "       use automatic differentiation;\n");
%!         fputs(fd, "       default output: none, structural nodes;\n");
%!         fputs(fd, "       default orientation: euler123;\n");
%!         fputs(fd, "       output precision: 16;\n");
%!         fputs(fd, "       max iterations: 0;\n");
%!         fprintf(fd, "        structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!         fprintf(fd, "        solids: %d;\n", opt_mbd_mesh.solids.number);
%!         fprintf(fd, "        joints: 3 + %d;\n", opt_mbd_mesh.joints.number);
%!         fputs(fd,   "        rigid bodies: 1;\n");
%!         if (options.f_rbk)
%!           fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!           fputs(fd, "        component,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           array, 2,\n");
%!           fputs(fd, "             mult, time, const, (omega1 - omega0) / (final_time - initial_time),\n");
%!           fputs(fd, "             const, omega0,\n");
%!           fputs(fd, "        angular acceleration,\n");
%!           fputs(fd, "        component,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           const, (omega1 - omega0) / (final_time - initial_time);\n");
%!         endif
%!         fputs(fd, "end: control data;\n");
%!         fputs(fd, "reference: ref_id_ground,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, eye,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, null;\n");
%!         fputs(fd, "reference: ref_id_shaft,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         fputs(fd, "        reference, ref_id_ground, eye,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         if (options.f_rbk)
%!           fputs(fd, "        reference, ref_id_ground, null;\n");
%!         else
%!           fputs(fd, "        reference, ref_id_ground, 0., 0., omega0;\n");
%!         endif
%!         fputs(fd, "reference: ref_id_rotor,\n");
%!         fputs(fd, "        reference, ref_id_shaft, 0., 0., o,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing1,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   -0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing2,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "begin: nodes;\n");
%!         fprintf(fd, "include: \"%s\";\n", solid_nodes_file);
%!         fputs(fd, "end: nodes;\n");
%!         fprintf(fd, "include: \"%s\";\n", solid_csl_file);
%!         fputs(fd, "begin: elements;\n");
%!         fputs(fd, "       drive caller: drive_id_rotor_speed, string, \"(omega1 - omega0) / (final_time - initial_time) * Time + omega0\", output, yes;\n");
%!         fputs(fd, "       drive caller: drive_id_time_step, string, \"2. * pi / (max(1., 360. * abs(model::drive(drive_id_rotor_speed, Time))))\";\n");
%!         fprintf(fd, "      include: \"%s\";\n", solid_elem_file);
%!         fputs(fd, "       body: body_id_unbalance,\n");
%!         fputs(fd, "             node_id_rotor,\n");
%!         fputs(fd, "                dm,\n");
%!         fputs(fd, "                  reference, ref_id_rotor, dr, 0., 0.,\n");
%!         fputs(fd, "                  diag, 0., 0., 0.;\n");
%!         fputs(fd, "        joint: joint_id_bearing1, total pin joint,\n");
%!         fputs(fd, "                node_id_bearing1,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, active,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "        joint: joint_id_bearing2, total pin joint,\n");
%!         fputs(fd, "               node_id_bearing2,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "               position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, inactive,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "     joint: joint_id_drive, angular velocity,\n");
%!         fputs(fd, "             # node label\n");
%!         fputs(fd, "             node_id_rotor, \n");
%!         fputs(fd, "             # direction\n");
%!         fputs(fd, "             0.,0.,1.,\n");
%!         fputs(fd, "             # angular velocity\n");
%!         if (~options.f_rbk)
%!           fputs(fd, "           reference, drive_id_rotor_speed;\n");
%!         else
%!           fputs(fd, "           null;\n");
%!         endif
%!         fputs(fd, "        inertia: elem_id_inertia,\n");
%!         fputs(fd, "                 position, reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                 orientation, reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                 body, all, solid, all, loadable, all,\n");
%!         fputs(fd, "                 output, both;\n");
%!         fputs(fd, "end: elements;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!     endif
%!     for i=1:numel(mbdyn_filenames)
%!       if (~enable_filenames(i))
%!         continue;
%!       endif
%!       options_mbd(i).output_file = sprintf("%s_%d", filename, i);
%!       options_mbd(i).mbdyn_command = options.mbdyn_command;
%!       if (~options.verbose)
%!         options_mbd(i).logfile = sprintf("%s_%d.stdout", filename, i);
%!       endif
%!       options_mbd(i).f_run_mbdyn2easyanim = false;
%!       param_file = sprintf("%s_%d.set", filename, i);
%!       putenv("MBDYN_ROTOR_DYN_CMS_ELEM_FILE", [filename, "_cms.elm"]);
%!       putenv("MBDYN_ROTOR_DYN_CMS_PARAM_FILE", param_file);
%!       mbdyn_pre_write_param_file(param_file, param);
%!       mbdyn_solver_run(mbdyn_filenames{i}, options_mbd(i));
%!       res(i).log_dat = mbdyn_post_load_log(options_mbd(i).output_file);
%!       [res(i).t, res(i).trajectory, res(i).deformation, res(i).velocity, res(i).acceleration, res(i).node_id] = mbdyn_post_load_output_struct(options_mbd(i).output_file);
%!       res(i).log_dat.vars = mbdyn_post_id_to_index(res(i), res(i).log_dat.vars);
%!       [res(i).drive_id, res(i).drive_data] = mbdyn_post_load_output_drv(options_mbd(i).output_file);
%!     endfor
%!     r = omega = cell(1, numel(mbdyn_filenames));
%!     figure("visible", "off");
%!     hold on;
%!     for i=1:numel(mbdyn_filenames)
%!       if (~enable_filenames(i))
%!         continue;
%!       endif
%!       omega{i} = res(i).drive_data{find(res(i).drive_id == res(i).log_dat.vars.drive_id_rotor_speed)};
%!       r{i} = norm(res(i).trajectory{res(i).log_dat.vars.node_idx_rotor}(:, 1:2), "rows");
%!       plot(omega{i} * 30 / pi * (1 / SI_unit.second), 1e3 * r{i} * SI_unit.meter, sprintf("-;%s;%d", printable_title(mbdyn_filename_suffix{i}), i));
%!     endfor
%!     xlabel("n [rpm]");
%!     ylabel("r [mm]");
%!     grid on;
%!     grid minor on;
%!     title("resonance curve center of disk versus speed - magnitude");
%!     figure_list();
%!     tol = 0.5e-2;
%!     for i=1:numel(mbdyn_filenames)
%!       if (~enable_filenames(i))
%!         continue;
%!       endif
%!       for j=1:numel(mbdyn_filenames)
%!         if (~enable_filenames(j))
%!           continue;
%!         endif
%!         assert_simple(interp1(omega{i}, r{i}, omega{j}, "pchip", "extrap"), r{j}, tol * max(abs(r{j})));
%!       endfor
%!     endfor
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   end_unwind_protect
%! endfunction
%!
%! ## Define the unit system
%! SI_unit.meter = 1e-3;
%! SI_unit.second = 1e-3;
%! SI_unit.kilogram = 1e-3;
%! SI_unit.newton = SI_unit.kilogram * SI_unit.meter / SI_unit.second^2;
%! SI_unit.pascal = SI_unit.newton / SI_unit.meter^2;
%! param.alpha = 0 / (1 / SI_unit.second);
%! param.beta = 1e-7 / (SI_unit.second);
%! param.h1 = 5e-3 / SI_unit.meter; ## fine mesh size
%! param.h2 = 5e-3 / SI_unit.meter; ## fine mesh size
%! param.l = 350e-3 / SI_unit.meter; ## bearing distance
%! param.d = 10e-3 / SI_unit.meter; ## shaft diameter
%! param.D = 150e-3 / SI_unit.meter; ##disk diameter
%! param.w = 15e-3 / SI_unit.meter; ## disk width
%! param.o = 75e-3 / SI_unit.meter; ## disk offset
%! param.ecg = 1e-3 * param.D;
%! param.dm = 1e-6 / SI_unit.kilogram;
%! param.E = 210000e6 / SI_unit.pascal;
%! param.nu = 0.3;
%! param.rho = 7850 / (SI_unit.kilogram / SI_unit.meter^3);
%! m1 = param.D^2 * pi / 4 * param.w * param.rho;
%! param.dr = param.ecg * (m1 + param.dm) / param.dm;
%! param.omega0 = 1000 * pi / 30 / (1 / SI_unit.second);
%! param.omega1 = 10000 * pi / 30 / (1 / SI_unit.second);
%! param.n = 10;
%! options.number_of_modes = int32(10);
%! options.scale_def = 10e-3;
%! options.geo_tol = sqrt(eps);
%! options.code.use_package = false;
%! options.mbdyn_command = "mbdyn";
%! options.f_run_mbdyn = [true, true];
%! options.verbose = false;
%! options.elem_type = "iso20";
%! options.f_rbk = false;
%! options.f_enable_beam = true;
%! options.f_enable_modal = true;
%! options.f_enable_solid = false; ## disabled because of long execution time
%! [omega1, r1] = rotordynamics_test_case(param, options, SI_unit);
%! options.f_rbk = true;
%! options.f_enable_modal = false;
%! [omega2, r2] = rotordynamics_test_case(param, options, SI_unit);
%! assert_simple(r2{2}, r1{2}, 1e-3 * max(abs(r1{2})));
%! options.elem_type = "iso20";
%! options.f_enable_beam = false; ## disabled because of coarse mesh and linear elements
%! options.f_enable_solid = true;
%! options.f_enable_modal = true;
%! param.h1 = 10e-3 / SI_unit.meter; ## coarse mesh size
%! param.h2 = 100e-3 / SI_unit.meter; ## coarse mesh size
%! [omega3, r3] = rotordynamics_test_case(param, options, SI_unit);
%! options.f_rbk = false;
%! options.f_enable_modal = false;
%! [omega4, r4] = rotordynamics_test_case(param, options, SI_unit);
%! assert_simple(r4{3}, r3{3}, 1e-3 * max(abs(r3{3})));

%!test
%! ## TEST 16
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E12 = 80000e6 / SI_unit_pascal;
%! param.E3 = 1000 * param.E12;
%! param.nu = 0.3;
%! param.rho = 1700 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.d2 = 15e-3 / SI_unit_meter;
%! param.d3 = 150e-3 / SI_unit_meter;
%! param.l1 = 120e-3 / SI_unit_meter;
%! param.l2 = 100e-3 / SI_unit_meter;
%! param.l3 = 15e-3 / SI_unit_meter;
%! param.h = 1.5e-3 / SI_unit_meter;
%! param.h3 = 10e-3 / SI_unit_meter;
%! OMEGA = 2 * pi * [0, 50,   0,   0;
%!                   0,  0, 300,   0;
%!                   0,  0,   0,   0] / (1 / SI_unit_second);
%! OMEGADOT = [0,  0,  0,  2e5;
%!             0,  0,  0,    0;
%!             0,  0,  0,    0];
%! idx = 1:columns(OMEGA);
%! OMEGA = OMEGA(:, idx);
%! OMEGADOT = OMEGADOT(:, idx);
%! opt_solver.pre_scaling = true;
%! opt_solver.refine_max_iter = int32(100);
%! opt_solver.solver = "pardiso";
%! opt_solver.number_of_threads = mbdyn_solver_num_threads_default();
%! opt_solver.symmetric = false; ## FEM_MAT_STIFFNESS_OMEGA_DOT makes it unsymmetric
%! I1 = param.d1^4 * pi / 64;
%! I2 = param.d2^4 * pi / 64;

%! alpha = d11 = param.l1 * param.l2^2 / (3 * param.E12 * I1) + (param.l2^3 - param.l3^3) / (3 * param.E12 * I2);
%! delta = gamma = d12 = param.l1 * param.l2 / (3 * param.E12 * I1) + (param.l2^2 - param.l3^2) / (2 * param.E12 * I2);
%! beta = d22 = param.l1 / (3 * param.E12 * I1) + (param.l2 - param.l3) / (param.E12 * I2);
%! m = param.rho * pi / 4 * param.d3^2 * 2 * param.l3;
%! Ja = m * (3 * param.d3^2 + 4 * (2 * param.l3)^2) / 48;
%! Jp = m * param.d3^2 / 8;
%! lambda0_ref = sqrt((alpha * m + beta * Ja) / (2 * m * Ja * (alpha * beta - gamma^2)) * (1 + [-1, 1] * sqrt(1 - (4 * m * Ja * (alpha * beta - gamma^2))/(alpha * m + beta * Ja)^2)));
%! lambda_ref = zeros(4, columns(OMEGA));
%! for i=1:columns(OMEGA)
%!   lambda_ref(:, i) = roots([m * Ja * (alpha * beta - gamma^2);
%!                             -m * Jp * OMEGA(1, i) * (alpha * beta - gamma^2);
%!                             -(alpha * m + beta * Ja);
%!                             beta * Jp * OMEGA(1, i);
%!                             1]);
%! endfor
%! lambda_ref = real(lambda_ref);
%! options.number_of_modes = int32(10);
%! fref = zeros(5, columns(OMEGA));
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   modal = struct("mbdyn", cell(1, columns(OMEGA)));
%!   for i=1:columns(OMEGA)
%!     unwind_protect
%!       [fd, msg] = fopen(mbdyn_file, "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s\"", mbdyn_file);
%!       endif
%!       fprintf(fd, " set: real E = %.16e;\n", param.E12);
%!       fprintf(fd, " set: real nu = %.16e;\n", param.nu);
%!       fprintf(fd, " set: real rho = %.16e;\n", param.rho);
%!       fprintf(fd, " set: real d1 = %.16e;\n", param.d1);
%!       fprintf(fd, " set: real d2 = %.16e;\n", param.d2);
%!       fprintf(fd, " set: real d3 = %.16e;\n", param.d3);
%!       fprintf(fd, " set: real l1 = %.16e;\n", param.l1);
%!       fprintf(fd, " set: real l2 = %.16e;\n", param.l2);
%!       fprintf(fd, " set: real l3 = %.16e;\n", param.l3);
%!       fputs(fd, " set: real G = E / (2. * (1 + nu));\n");
%!       fputs(fd, " set: real A1 = d1^2 * pi / 4.;\n");
%!       fputs(fd, " set: real A2 = d2^2 * pi / 4.;\n");
%!       fputs(fd, " set: real As1 = 9./10. * A1;\n");
%!       fputs(fd, " set: real As2 = 9./10. * A2;\n");
%!       fputs(fd, " set: real Iy1 = d1^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iy2 = d2^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iz1 = Iy1;\n");
%!       fputs(fd, " set: real Iz2 = Iy2;\n");
%!       fputs(fd, " set: real It1 = Iy1 + Iz1;\n");
%!       fputs(fd, " set: real It2 = Iy2 + Iz2;\n");
%!       fputs(fd, " set: real m = rho * pi / 4 * d3^2 * 2 * l3;\n");
%!       fputs(fd, " set: real Ja = m * (3 * d3^2 + 4 * (2 * l3)^2) / 48;\n");
%!       fputs(fd, " set: real Jp = m * d3^2 / 8;\n");
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGA%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGA(j, i));
%!       endfor
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGADOT%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGADOT(j, i));
%!       endfor
%!       fputs(fd, " set: real t1 = 1;\n");
%!       fputs(fd, " set: real N = 20000;\n");
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "    problem: initial value; # the default\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "    initial time: 0;\n");
%!       fputs(fd, "    final time: t1;\n");
%!       fputs(fd, "    time step: t1 / N;\n");
%!       fputs(fd, "         threads: disable;\n");
%!       fputs(fd, "    max iterations: 10000;\n");
%!       fputs(fd, "    tolerance: 1.e-6;\n");
%!       fputs(fd, "    linear solver: umfpack, scale, row max column max, always, max iterations, 3;\n");
%!       fputs(fd, "    method: implicit euler;\n");
%!       fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!       fputs(fd, "         output: iterations;\n");
%!       fputs(fd, "    eigenanalysis: t1,\n");
%!       fputs(fd, "    output matrices, \n");
%!       fputs(fd, "    parameter, 1e-3,\n");
%!       fputs(fd, "    output eigenvectors,\n");
%!       fputs(fd, "         output geometry,\n");
%!       fputs(fd, "         lower frequency limit, 0.0001, upper frequency limit, 1000,\n");
%!       fputs(fd, "    use lapack,balance,permute,suffix format, \"%02d\";\n");
%!       fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "    output meter: closest next, 0., forever, t1 / 100;\n");
%!       fputs(fd, "        use automatic differentiation;\n");
%!       fputs(fd, "        rigid body kinematics: const, angular velocity, OMEGAx, OMEGAy, OMEGAz;\n");
%!       fputs(fd, "        rigid body kinematics: const, angular acceleration, OMEGADOTx, OMEGADOTy, OMEGADOTz;\n");
%!       fputs(fd, "    structural nodes: 9;\n");
%!       fputs(fd, "    rigid bodies: 1;\n");
%!       fputs(fd, "    beams: 4;\n");
%!       fputs(fd, "    joints: 2;\n");
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fputs(fd, "    structural: 1, static, \n");
%!       fputs(fd, "            reference, global, null, \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 2, static, \n");
%!       fputs(fd, "            reference, global, 0.25 * l1, 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 3, static, \n");
%!       fputs(fd, "            reference, global, 0.5 * l1, 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 4, static, \n");
%!       fputs(fd, "            reference, global, 0.75 * l1, 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 5, static, \n");
%!       fputs(fd, "            reference, global, l1, 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 6, static, \n");
%!       fputs(fd, "            reference, global, l1 + 0.25 * (l2 - l3), 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 7, static, \n");
%!       fputs(fd, "            reference, global, l1 + 0.5 * (l2 - l3), 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 8, static, \n");
%!       fputs(fd, "            reference, global, l1 + 0.75 * (l2 - l3), 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null;\n");
%!       fputs(fd, "    structural: 9, dynamic, \n");
%!       fputs(fd, "            reference, global, l1 + l2, 0., 0., \n");
%!       fputs(fd, "            reference, global, eye, \n");
%!       fputs(fd, "            reference, global, null,\n");
%!       fputs(fd, "            reference, global, null; \n");
%!       fputs(fd, " end: nodes;\n");
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "    body: 1, \n");
%!       fputs(fd, "            9,\n");
%!       fputs(fd, "            m, \n");
%!       fputs(fd, "            null, \n");
%!       fputs(fd, "            diag,   Jp, \n");
%!       fputs(fd, "                    Ja, \n");
%!       fputs(fd, "                    Ja,\n");
%!       fputs(fd, "            orientation, reference, global, eye;\n");
%!       fputs(fd, "        beam3: 1,\n");
%!       fputs(fd, "            1, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            2, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            3, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            reference, global, eye,\n");
%!       fputs(fd, "            linear elastic generic, \n");
%!       fputs(fd, "            diag, E * A1 , G * As1, G * As1, \n");
%!       fputs(fd, "                  G * It1, E * Iy1, E * Iz1,\n");
%!       fputs(fd, "            same,\n");
%!       fputs(fd, "            same;\n");
%!       fputs(fd, "        beam3: 2,\n");
%!       fputs(fd, "            3, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            4, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            5, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            reference, global, eye,\n");
%!       fputs(fd, "            linear elastic generic, \n");
%!       fputs(fd, "            diag, E * A1 , G * As1, G * As1, \n");
%!       fputs(fd, "                  G * It1, E * Iy1, E * Iz1,\n");
%!       fputs(fd, "            same,\n");
%!       fputs(fd, "            same;\n");
%!       fputs(fd, "        beam3: 3,\n");
%!       fputs(fd, "            5, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            6, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            7, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            reference, global, eye,\n");
%!       fputs(fd, "            linear elastic generic, \n");
%!       fputs(fd, "            diag, E * A2 , G * As2, G * As2, \n");
%!       fputs(fd, "                  G * It2, E * Iy2, E * Iz2,\n");
%!       fputs(fd, "            same,\n");
%!       fputs(fd, "            same;\n");
%!       fputs(fd, "       beam3: 4,\n");
%!       fputs(fd, "            7, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            8, position, reference, node, null,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            9, position, reference, global, l1 + l2 - l3, 0., 0.,\n");
%!       fputs(fd, "            orientation, reference, node, eye,\n");
%!       fputs(fd, "            reference, global, eye,\n");
%!       fputs(fd, "            linear elastic generic, \n");
%!       fputs(fd, "            diag, E * A2 , G * As2, G * As2, \n");
%!       fputs(fd, "                  G * It2, E * Iy2, E * Iz2,\n");
%!       fputs(fd, "            same,\n");
%!       fputs(fd, "            same;\n");
%!       fputs(fd, "    joint: 1, total pin joint,\n");
%!       fputs(fd, "                    1,\n");
%!       fputs(fd, "                            position,		reference, global, null,\n");
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            position,		reference, global, null,\n");
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fputs(fd, "    joint: 2, total pin joint,\n");
%!       fputs(fd, "                    5,\n");
%!       fputs(fd, "                            position,		reference, global, l1, 0., 0.,\n");
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            position,		reference, global, l1, 0., 0.,\n");
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    inactive, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!       fd = -1;
%!     end_unwind_protect
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!     opt_mbdyn.output_file = [filename, "_mbdyn"];
%!     info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!     modal(i).mbdyn = mbdyn_post_load_output_eig(opt_mbdyn.output_file);
%!     for j=1:numel(modal(i).mbdyn.f)
%!       opt_modal.mode_index = j;
%!       opt_modal.scale = 100;
%!       mode_file = [opt_mbdyn.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(opt_mbdyn.output_file, [mode_file, ".mov"], opt_modal, modal(i).mbdyn);
%!       [err, msg] = symlink([opt_mbdyn.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = true;
%!       opt_post.f_runEasyAnim = false;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!     fref(:, i) = modal(i).mbdyn.f(:);
%!   endfor
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Geometry.Tolerance = 1e-6;\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "s2[] = Extrude {l2 - l3, 0, 0}{ Line{s1[0]}; Layers{Ceil((l2 - l3) / h)}; Recombine; };\n");
%!     fputs(fd, "s3[] = Extrude {2 * l3, 0, 0}{ Line{s2[0]}; Layers{Ceil(2 * l3 / h3)}; Recombine; };\n");
%!     fputs(fd, "s4[] = Extrude {0, 0.5 * (d2 - d1), 0}{ Line{s2[3],s3[3]}; Layers{Ceil(0.5 * (d2 - d1) / h)}; Recombine; };\n");
%!     fputs(fd, "s6[] = Extrude {0, 0.5 * (d3 - d2), 0}{ Line{s4[4]}; Layers{Ceil(0.5 * (d3 - d1) / h3)}; Recombine; };\n");
%!     fputs(fd, "se0[] = {s1[1], s2[1], s3[1], s4[1], s4[5], s6[1]};\n");
%!     fputs(fd, "v1[] = {};\n");
%!     fputs(fd, "se1[] = {};\n");
%!     fputs(fd, "For i In {0:#se0[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se0[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v1[i] = vtmp[1];\n");
%!     fputs(fd, "  se1[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v2[] = {};\n");
%!     fputs(fd, "se2[] = {};\n");
%!     fputs(fd, "For i In {0:#se1[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se1[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v2[i] = vtmp[1];\n");
%!     fputs(fd, "  se2[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v3[] = {};\n");
%!     fputs(fd, "se3[] = {};\n");
%!     fputs(fd, "For i In {0:#se2[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se2[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v3[i] = vtmp[1];\n");
%!     fputs(fd, "  se3[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v4[] = {};\n");
%!     fputs(fd, "se4[] = {};\n");
%!     fputs(fd, "For i In {0:#se3[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se3[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v4[i] = vtmp[1];\n");
%!     fputs(fd, "  se4[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v[] = {v1[], v2[], v3[], v4[]};\n");
%!     fputs(fd, "Physical Volume(1) = {19, 1, 13, 7};\n");
%!     fputs(fd, "Physical Volume(2) = {22, 4, 16, 10, 14, 8, 20, 2};\n");
%!     fputs(fd, "Physical Volume(3) = {18, 12, 24, 6, 23, 5, 17, 11, 21, 3, 15, 9};\n");
%!     fputs(fd, "Physical Surface(1) = {89, 8, 62, 35};\n");
%!     fputs(fd, "Physical Surface(2) = {19, 100, 73, 46};\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   e1 = [1; 0.5; 0.4];
%!   e2 = [0; 1; 0];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   R = [e1, e2, e3];
%!   R *= diag(1 ./ norm(R, "cols"));
%!   T = [R.', zeros(3, 3);
%!        zeros(3, 3), R.'];
%!   mesh.nodes *= T;
%!   mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, "iso4");
%!   mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, "iso4");
%!   mesh.elements.joints(2).nodes = node_idx_bearing2;
%!   mesh.elements.joints(2).C = [0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0] * T;
%!   mesh.elements.joints(1).nodes = node_idx_bearing1;
%!   mesh.elements.joints(1).C = [1, 0, 0, 0, 0, 0;
%!                                0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0;
%!                                0, 0, 0, 1, 0, 0] * T;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
%!   mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!   for i=1:3
%!     mesh.materials.iso8([mesh.groups.iso8(find([[mesh.groups.iso8.id] == i])).elements]) = i;
%!   endfor
%!   mesh.material_data(1).rho = 0;
%!   mesh.material_data(1).E = param.E12;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(2).rho = 0;
%!   mesh.material_data(2).E = param.E12;
%!   mesh.material_data(2).nu = param.nu;
%!   mesh.material_data(3).rho = param.rho;
%!   mesh.material_data(3).E = param.E3;
%!   mesh.material_data(3).nu = param.nu;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = opt_solver.number_of_threads;
%!   load_case = fem_pre_load_case_create_empty(columns(OMEGA));
%!   empty_cell = cell(1, columns(OMEGA));
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   sol_eig = struct("lambda", empty_cell, "f", empty_cell, "def", empty_cell, "D", empty_cell);
%!   mat_red = struct("Mred", empty_cell, "Kred", empty_cell, "Dred", empty_cell);
%!   for i=1:columns(OMEGA)
%!     fprintf(stderr, "n=%.0frpm\n", norm(OMEGA(:, i)) * 30 / pi / SI_unit_second);
%!     load_case(i).omega = R * OMEGA(:, i);
%!     load_case(i).omegadot = R * OMEGADOT(:, i);
%!     [mat_ass.M, ...
%!      mat_ass.K, ...
%!      mat_ass.KOMEGA, ...
%!      mat_ass.KOMEGA_DOT, ...
%!      mat_ass.DOMEGA, ...
%!      mat_ass.R, ...
%!      mat_ass.dm, ...
%!      mat_ass.S, ...
%!      mat_ass.J] = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_MASS, ...
%!                                   FEM_MAT_STIFFNESS, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                   FEM_MAT_DAMPING_OMEGA, ...
%!                                   FEM_VEC_LOAD_CONSISTENT, ...
%!                                   FEM_SCA_TOT_MASS, ...
%!                                   FEM_VEC_INERTIA_M1, ...
%!                                   FEM_MAT_INERTIA_J], ...
%!                                  load_case(i));
%!     Tred = zeros(rows(mesh.nodes) * 3, 6);
%!     for j=1:rows(mesh.nodes)
%!       Tred((j - 1) * 3 + (1:3), 1:3) = eye(3);
%!       Tred((j - 1) * 3 + (1:3), 4:6) = -skew(mesh.nodes(j, 1:3));
%!     endfor
%!     idx = dof_map.ndof(:, 1:3).'(:);
%!     mat_red(i).Mred = Tred.' * mat_ass.M(idx, idx) * Tred;
%!     mat_red(i).Kred = Tred.' * mat_ass.K(idx, idx) * Tred;
%!     mat_red(i).Dred = Tred.' * mat_ass.DOMEGA(idx, idx) * Tred;
%!     sol_stat(i).def = fem_sol_static(mesh, dof_map, mat_ass).def;
%!     sol_stat(i).stress = fem_ass_matrix(mesh, ...
%!                                         dof_map, ...
%!                                         [FEM_VEC_STRESS_CAUCH], ...
%!                                         load_case(i), ...
%!                                         sol_stat(i));
%!     load_case_pre_stress.tau0 = sol_stat(i).stress.tau;
%!     mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS_TAU0], ...
%!                                    load_case_pre_stress);
%!     mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0 + mat_ass.KOMEGA_DOT;
%!     mat_ass.D = mat_ass.DOMEGA;
%!     sol_eig(i) = fem_sol_modal_damped(mesh, ...
%!                                       dof_map, ...
%!                                       mat_ass, ...
%!                                       options.number_of_modes, ...
%!                                       opt_solver);
%!   endfor
%!   f = zeros(options.number_of_modes, numel(sol_eig));
%!   for i=1:columns(f)
%!     f(:, i) = sort(sol_eig(i).f(:));
%!   endfor
%!   tol = 2e-2;
%!   assert_simple(f(floor(end/2+1):end,:), fref, tol * max(max(abs(fref))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 17
%! pkg load mboct-fem-pkg;
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1e-3;
%! SI_unit_kilogram = 1e-3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E12 = 80000e6 / SI_unit_pascal;
%! param.E3 =  210000e6 / SI_unit_pascal;
%! param.nu = 0.3;
%! param.rho = 1700 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.d2 = 15e-3 / SI_unit_meter;
%! param.d3 = 150e-3 / SI_unit_meter;
%! param.l1 = 120e-3 / SI_unit_meter;
%! param.l2 = 100e-3 / SI_unit_meter;
%! param.l3 = 15e-3 / SI_unit_meter;
%! param.h = 10e-3 / SI_unit_meter;
%! param.h3 = 10e-3 / SI_unit_meter;
%! OMEGA = 2 * pi * [0, 50,   0,   0;
%!                   0,  0, 150,   0;
%!                   0,  0,   0,   0] / (1 / SI_unit_second);
%! OMEGADOT = [0,  0,  0,  2e5;
%!             0,  0,  0,    0;
%!             0,  0,  0,    0] / (1 / SI_unit_second^2);
%! idx = 1:columns(OMEGA);
%! OMEGA = OMEGA(:, idx);
%! OMEGADOT = OMEGADOT(:, idx);
%! options.post_proc_modes = false;
%! options.verbose = false;
%! opt_solver.pre_scaling = true;
%! opt_solver.refine_max_iter = int32(100);
%! opt_solver.solver = "pardiso";
%! opt_solver.number_of_threads = mbdyn_solver_num_threads_default();
%! opt_solver.symmetric = false; ## FEM_MAT_STIFFNESS_OMEGA_DOT makes it unsymmetric
%! I1 = param.d1^4 * pi / 64;
%! I2 = param.d2^4 * pi / 64;
%! alpha = d11 = param.l1 * param.l2^2 / (3 * param.E12 * I1) + (param.l2^3 - param.l3^3) / (3 * param.E12 * I2);
%! delta = gamma = d12 = param.l1 * param.l2 / (3 * param.E12 * I1) + (param.l2^2 - param.l3^2) / (2 * param.E12 * I2);
%! beta = d22 = param.l1 / (3 * param.E12 * I1) + (param.l2 - param.l3) / (param.E12 * I2);
%! m = param.rho * pi / 4 * param.d3^2 * 2 * param.l3;
%! Ja = m * (3 * param.d3^2 + 4 * (2 * param.l3)^2) / 48;
%! Jp = m * param.d3^2 / 8;
%! lambda0_ref = sqrt((alpha * m + beta * Ja) / (2 * m * Ja * (alpha * beta - gamma^2)) * (1 + [-1, 1] * sqrt(1 - (4 * m * Ja * (alpha * beta - gamma^2))/(alpha * m + beta * Ja)^2)));
%! lambda_ref = zeros(4, columns(OMEGA));
%! for i=1:columns(OMEGA)
%!   lambda_ref(:, i) = roots([m * Ja * (alpha * beta - gamma^2);
%!                             -m * Jp * OMEGA(1, i) * (alpha * beta - gamma^2);
%!                             -(alpha * m + beta * Ja);
%!                             beta * Jp * OMEGA(1, i);
%!                             1]);
%! endfor
%! lambda_ref = real(lambda_ref);
%! options.number_of_modes = int32(10);
%! fref = zeros(5, columns(OMEGA));
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   modal = struct("mbdyn", cell(1, columns(OMEGA)));
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Geometry.Tolerance = 1e-6;\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "s2[] = Extrude {l2 - l3, 0, 0}{ Line{s1[0]}; Layers{Ceil((l2 - l3) / h)}; Recombine; };\n");
%!     fputs(fd, "s3[] = Extrude {2 * l3, 0, 0}{ Line{s2[0]}; Layers{Ceil(2 * l3 / h3)}; Recombine; };\n");
%!     fputs(fd, "s4[] = Extrude {0, 0.5 * (d2 - d1), 0}{ Line{s2[3],s3[3]}; Layers{Ceil(0.5 * (d2 - d1) / h)}; Recombine; };\n");
%!     fputs(fd, "s6[] = Extrude {0, 0.5 * (d3 - d2), 0}{ Line{s4[4]}; Layers{Ceil(0.5 * (d3 - d1) / h3)}; Recombine; };\n");
%!     fputs(fd, "se0[] = {s1[1], s2[1], s3[1], s4[1], s4[5], s6[1]};\n");
%!     fputs(fd, "v1[] = {};\n");
%!     fputs(fd, "se1[] = {};\n");
%!     fputs(fd, "For i In {0:#se0[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se0[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v1[i] = vtmp[1];\n");
%!     fputs(fd, "  se1[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v2[] = {};\n");
%!     fputs(fd, "se2[] = {};\n");
%!     fputs(fd, "For i In {0:#se1[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se1[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v2[i] = vtmp[1];\n");
%!     fputs(fd, "  se2[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v3[] = {};\n");
%!     fputs(fd, "se3[] = {};\n");
%!     fputs(fd, "For i In {0:#se2[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se2[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v3[i] = vtmp[1];\n");
%!     fputs(fd, "  se3[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v4[] = {};\n");
%!     fputs(fd, "se4[] = {};\n");
%!     fputs(fd, "For i In {0:#se3[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se3[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v4[i] = vtmp[1];\n");
%!     fputs(fd, "  se4[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v[] = {v1[], v2[], v3[], v4[]};\n");
%!     fputs(fd, "Physical Volume(1) = {19, 1, 13, 7};\n");
%!     fputs(fd, "Physical Volume(2) = {22, 4, 16, 10, 14, 8, 20, 2};\n");
%!     fputs(fd, "Physical Volume(3) = {18, 12, 24, 6, 23, 5, 17, 11, 21, 3, 15, 9};\n");
%!     fputs(fd, "Physical Surface(1) = {89, 8, 62, 35};\n");
%!     fputs(fd, "Physical Surface(2) = {19, 100, 73, 46};\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   e1 = [1; 0.5; 0.4];
%!   e2 = [0; 1; 0];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   R = [e1, e2, e3];
%!   R *= diag(1 ./ norm(R, "cols"));
%!   R = eye(3);
%!   T = [R.', zeros(3, 3);
%!        zeros(3, 3), R.'];
%!   mesh.nodes *= T;
%!   mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, "iso4");
%!   mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, "iso4");
%!   mesh.elements.joints(2).nodes = node_idx_bearing2;
%!   mesh.elements.joints(2).C = [0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0] * T;
%!   mesh.elements.joints(1).nodes = node_idx_bearing1;
%!   mesh.elements.joints(1).C = [1, 0, 0, 0, 0, 0;
%!                                0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0;
%!                                0, 0, 0, 1, 0, 0] * T;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
%!   mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!   for i=1:3
%!     mesh.materials.iso8([mesh.groups.iso8(find([[mesh.groups.iso8.id] == i])).elements]) = i;
%!   endfor
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E12;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(2).rho = param.rho;
%!   mesh.material_data(2).E = param.E12;
%!   mesh.material_data(2).nu = param.nu;
%!   mesh.material_data(3).rho = param.rho;
%!   mesh.material_data(3).E = param.E3;
%!   mesh.material_data(3).nu = param.nu;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = opt_solver.number_of_threads;
%!   load_case = fem_pre_load_case_create_empty(columns(OMEGA));
%!   empty_cell = cell(1, columns(OMEGA));
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   sol_eig = struct("lambda", empty_cell, "f", empty_cell, "def", empty_cell, "D", empty_cell);
%!   mat_red = struct("Mred", empty_cell, "Kred", empty_cell, "Dred", empty_cell);
%!   for i=1:columns(OMEGA)
%!     fprintf(stderr, "n=%.0frpm\n", norm(OMEGA(:, i)) * 30 / pi / SI_unit_second);
%!     load_case(i).omega = R * OMEGA(:, i);
%!     load_case(i).omegadot = R * OMEGADOT(:, i);
%!     [mat_ass.M, ...
%!      mat_ass.K, ...
%!      mat_ass.KOMEGA, ...
%!      mat_ass.KOMEGA_DOT, ...
%!      mat_ass.DOMEGA, ...
%!      mat_ass.R, ...
%!      mat_ass.dm, ...
%!      mat_ass.S, ...
%!      mat_ass.J] = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_MASS, ...
%!                                   FEM_MAT_STIFFNESS, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                   FEM_MAT_DAMPING_OMEGA, ...
%!                                   FEM_VEC_LOAD_CONSISTENT, ...
%!                                   FEM_SCA_TOT_MASS, ...
%!                                   FEM_VEC_INERTIA_M1, ...
%!                                   FEM_MAT_INERTIA_J], ...
%!                                  load_case(i));
%!     Tred = zeros(rows(mesh.nodes) * 3, 6);
%!     for j=1:rows(mesh.nodes)
%!       Tred((j - 1) * 3 + (1:3), 1:3) = eye(3);
%!       Tred((j - 1) * 3 + (1:3), 4:6) = -skew(mesh.nodes(j, 1:3));
%!     endfor
%!     idx = dof_map.ndof(:, 1:3).'(:);
%!     mat_red(i).Mred = Tred.' * mat_ass.M(idx, idx) * Tred;
%!     mat_red(i).Kred = Tred.' * mat_ass.K(idx, idx) * Tred;
%!     mat_red(i).Dred = Tred.' * mat_ass.DOMEGA(idx, idx) * Tred;
%!     sol_stat(i).def = fem_sol_static(mesh, dof_map, mat_ass).def;
%!     sol_stat(i).stress = fem_ass_matrix(mesh, ...
%!                                         dof_map, ...
%!                                         [FEM_VEC_STRESS_CAUCH], ...
%!                                         load_case(i), ...
%!                                         sol_stat(i));
%!     load_case_pre_stress.tau0 = sol_stat(i).stress.tau;
%!     mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS_TAU0], ...
%!                                    load_case_pre_stress);
%!     mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0 + mat_ass.KOMEGA_DOT;
%!     mat_ass.D = mat_ass.DOMEGA;
%!     sol_eig(i) = fem_sol_modal_damped(mesh, ...
%!                                       dof_map, ...
%!                                       mat_ass, ...
%!                                       options.number_of_modes, ...
%!                                       opt_solver);
%!   endfor
%!   for i=1:columns(OMEGA)
%!     opt_mbd_mesh = struct();
%!     opt_mbd_mesh.joints.number = 2;
%!     opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!     opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!     opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!     opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!     load_case_empty = struct();
%!     opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!     unwind_protect
%!       [fd, msg] = fopen(mbdyn_file, "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s\"", mbdyn_file);
%!       endif
%!       fprintf(fd, " set: real E = %.16e;\n", param.E12);
%!       fprintf(fd, " set: real nu = %.16e;\n", param.nu);
%!       fprintf(fd, " set: real rho = %.16e;\n", param.rho);
%!       fprintf(fd, " set: real d1 = %.16e;\n", param.d1);
%!       fprintf(fd, " set: real d2 = %.16e;\n", param.d2);
%!       fprintf(fd, " set: real d3 = %.16e;\n", param.d3);
%!       fprintf(fd, " set: real l1 = %.16e;\n", param.l1);
%!       fprintf(fd, " set: real l2 = %.16e;\n", param.l2);
%!       fprintf(fd, " set: real l3 = %.16e;\n", param.l3);
%!       fputs(fd, " set: real G = E / (2. * (1 + nu));\n");
%!       fputs(fd, " set: real A1 = d1^2 * pi / 4.;\n");
%!       fputs(fd, " set: real A2 = d2^2 * pi / 4.;\n");
%!       fputs(fd, " set: real As1 = 9./10. * A1;\n");
%!       fputs(fd, " set: real As2 = 9./10. * A2;\n");
%!       fputs(fd, " set: real Iy1 = d1^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iy2 = d2^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iz1 = Iy1;\n");
%!       fputs(fd, " set: real Iz2 = Iy2;\n");
%!       fputs(fd, " set: real It1 = Iy1 + Iz1;\n");
%!       fputs(fd, " set: real It2 = Iy2 + Iz2;\n");
%!       fputs(fd, " set: real m = rho * pi / 4 * d3^2 * 2 * l3;\n");
%!       fputs(fd, " set: real Ja = m * (3 * d3^2 + 4 * (2 * l3)^2) / 48;\n");
%!       fputs(fd, " set: real Jp = m * d3^2 / 8;\n");
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGA%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGA(j, i));
%!       endfor
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGADOT%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGADOT(j, i));
%!       endfor
%!       fprintf(fd, " set: real t1 = %g;\n", 1000 / SI_unit_second);
%!       fputs(fd, " set: real N = 1000;\n");
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "    problem: initial value; # the default\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "    initial time: 0;\n");
%!       fputs(fd, "    final time: t1;\n");
%!       fputs(fd, "    time step: t1 / N;\n");
%!       fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!       fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!       fputs(fd, "    max iterations: 10000;\n");
%!       fputs(fd, "    tolerance: 1.e-3, 1e-3;\n");
%!       fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!       fputs(fd, "    method: bdf;\n");
%!       fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!       fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!       fputs(fd, "    eigenanalysis: t1,\n");
%!       fputs(fd, "    # output matrices, \n");
%!       fputs(fd, "    # parameter, 1e-3, ## use default estimate\n");
%!       fputs(fd, "    output eigenvectors,\n");
%!       fputs(fd, "         # output geometry,\n");
%!       fprintf(fd, "         lower frequency limit, %e,\n", 0.01 / (SI_unit_second^-1));
%!       fprintf(fd, "         upper frequency limit, %e,\n", 1000 / (SI_unit_second^-1));
%!       fprintf(fd, "    use arpack,%d,%d,0.,suffix format, \"%%02d\";\n", 2 * options.number_of_modes, options.number_of_modes * 20);
%!       fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-8, forcing term max tolerance, 1e-3;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!       fputs(fd, "        use automatic differentiation;\n");
%!       fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!       fputs(fd, "        component,\n");
%!       fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "        angular acceleration,\n");
%!       fputs(fd, "        component,\n");
%!       fputs(fd, "           string, \"OMEGADOTx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGADOTy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGADOTz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!       fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!       fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!       fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", nodes_file);
%!       fputs(fd, " end: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", csl_file);
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "    joint: 1, total pin joint,\n");
%!       fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fputs(fd, "    joint: 2, total pin joint,\n");
%!       fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    inactive, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fprintf(fd, "include: \"%s\";\n", elem_file);
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!       fd = -1;
%!     end_unwind_protect
%!     if (~options.verbose)
%!       opt_mbdyn.logfile = [filename, ".stdout"];
%!     endif
%!     opt_mbdyn.output_file = [filename, "_mbdyn"];
%!     info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!     [mesh_sol(i), sol(i)] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!     modal(i).mbdyn = mbdyn_post_load_output_eig(opt_mbdyn.output_file);
%!     if (options.post_proc_modes)
%!     for j=1:numel(modal(i).mbdyn.f)
%!       opt_modal.mode_index = j;
%!       opt_modal.scale = 100;
%!       mode_file = [opt_mbdyn.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(opt_mbdyn.output_file, [mode_file, ".mov"], opt_modal, modal(i).mbdyn);
%!       [err, msg] = symlink([opt_mbdyn.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = true;
%!       opt_post.f_runEasyAnim = false;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!     endif
%!     fref(:, i) = modal(i).mbdyn.f(1:rows(fref));
%!   endfor
%!   f = zeros(options.number_of_modes, numel(sol_eig));
%!   for i=1:columns(f)
%!     f(:, i) = sort(sol_eig(i).f(:));
%!   endfor
%!   tol = 2e-2;
%!   assert_simple(f(floor(end/2+1):end,:), fref, tol * max(max(abs(fref))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 18
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1e-3;
%! SI_unit_kilogram = 1e-3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E12 = 80000e6 / SI_unit_pascal;
%! param.E3 =  210000e6 / SI_unit_pascal;
%! param.nu = 0.3;
%! param.rho = 1700 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.d2 = 15e-3 / SI_unit_meter;
%! param.d3 = 150e-3 / SI_unit_meter;
%! param.l1 = 120e-3 / SI_unit_meter;
%! param.l2 = 100e-3 / SI_unit_meter;
%! param.l3 = 15e-3 / SI_unit_meter;
%! param.h = 40e-3 / SI_unit_meter;
%! param.h3 = 10e-3 / SI_unit_meter;
%! OMEGA = 2 * pi * [0, 50,   0,   0;
%!                   0,  0, 150,   0;
%!                   0,  0,   0,   0] / (1 / SI_unit_second);
%! OMEGADOT = [0,  0,  0,  2e5;
%!             0,  0,  0,    0;
%!             0,  0,  0,    0] / (1 / SI_unit_second^2);
%! idx = 1:columns(OMEGA);
%! OMEGA = OMEGA(:, idx);
%! OMEGADOT = OMEGADOT(:, idx);
%! options.post_proc_modes = false;
%! options.verbose = false;
%! opt_solver.pre_scaling = true;
%! opt_solver.refine_max_iter = int32(100);
%! opt_solver.solver = "pardiso";
%! opt_solver.number_of_threads = mbdyn_solver_num_threads_default();
%! opt_solver.symmetric = false; ## FEM_MAT_STIFFNESS_OMEGA_DOT makes it unsymmetric
%! I1 = param.d1^4 * pi / 64;
%! I2 = param.d2^4 * pi / 64;
%! alpha = d11 = param.l1 * param.l2^2 / (3 * param.E12 * I1) + (param.l2^3 - param.l3^3) / (3 * param.E12 * I2);
%! delta = gamma = d12 = param.l1 * param.l2 / (3 * param.E12 * I1) + (param.l2^2 - param.l3^2) / (2 * param.E12 * I2);
%! beta = d22 = param.l1 / (3 * param.E12 * I1) + (param.l2 - param.l3) / (param.E12 * I2);
%! m = param.rho * pi / 4 * param.d3^2 * 2 * param.l3;
%! Ja = m * (3 * param.d3^2 + 4 * (2 * param.l3)^2) / 48;
%! Jp = m * param.d3^2 / 8;
%! lambda0_ref = sqrt((alpha * m + beta * Ja) / (2 * m * Ja * (alpha * beta - gamma^2)) * (1 + [-1, 1] * sqrt(1 - (4 * m * Ja * (alpha * beta - gamma^2))/(alpha * m + beta * Ja)^2)));
%! lambda_ref = zeros(4, columns(OMEGA));
%! for i=1:columns(OMEGA)
%!   lambda_ref(:, i) = roots([m * Ja * (alpha * beta - gamma^2);
%!                             -m * Jp * OMEGA(1, i) * (alpha * beta - gamma^2);
%!                             -(alpha * m + beta * Ja);
%!                             beta * Jp * OMEGA(1, i);
%!                             1]);
%! endfor
%! lambda_ref = real(lambda_ref);
%! options.number_of_modes = int32(10);
%! fref = zeros(5, columns(OMEGA));
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   modal = struct("mbdyn", cell(1, columns(OMEGA)));
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Geometry.Tolerance = 1e-6;\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "s2[] = Extrude {l2 - l3, 0, 0}{ Line{s1[0]}; Layers{Ceil((l2 - l3) / h)}; Recombine; };\n");
%!     fputs(fd, "s3[] = Extrude {2 * l3, 0, 0}{ Line{s2[0]}; Layers{Ceil(2 * l3 / h3)}; Recombine; };\n");
%!     fputs(fd, "s4[] = Extrude {0, 0.5 * (d2 - d1), 0}{ Line{s2[3],s3[3]}; Layers{Ceil(0.5 * (d2 - d1) / h)}; Recombine; };\n");
%!     fputs(fd, "s6[] = Extrude {0, 0.5 * (d3 - d2), 0}{ Line{s4[4]}; Layers{Ceil(0.5 * (d3 - d1) / h3)}; Recombine; };\n");
%!     fputs(fd, "se0[] = {s1[1], s2[1], s3[1], s4[1], s4[5], s6[1]};\n");
%!     fputs(fd, "v1[] = {};\n");
%!     fputs(fd, "se1[] = {};\n");
%!     fputs(fd, "For i In {0:#se0[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se0[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v1[i] = vtmp[1];\n");
%!     fputs(fd, "  se1[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v2[] = {};\n");
%!     fputs(fd, "se2[] = {};\n");
%!     fputs(fd, "For i In {0:#se1[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se1[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v2[i] = vtmp[1];\n");
%!     fputs(fd, "  se2[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v3[] = {};\n");
%!     fputs(fd, "se3[] = {};\n");
%!     fputs(fd, "For i In {0:#se2[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se2[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v3[i] = vtmp[1];\n");
%!     fputs(fd, "  se3[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v4[] = {};\n");
%!     fputs(fd, "se4[] = {};\n");
%!     fputs(fd, "For i In {0:#se3[] - 1}\n");
%!     fputs(fd, "  vtmp[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi/2}{ Surface{se3[i]}; Layers{Ceil(Pi/2 * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "  v4[i] = vtmp[1];\n");
%!     fputs(fd, "  se4[i] = vtmp[0];\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "v[] = {v1[], v2[], v3[], v4[]};\n");
%!     fputs(fd, "Physical Volume(1) = {19, 1, 13, 7};\n");
%!     fputs(fd, "Physical Volume(2) = {22, 4, 16, 10, 14, 8, 20, 2};\n");
%!     fputs(fd, "Physical Volume(3) = {18, 12, 24, 6, 23, 5, 17, 11, 21, 3, 15, 9};\n");
%!     fputs(fd, "Physical Surface(1) = {89, 8, 62, 35};\n");
%!     fputs(fd, "Physical Surface(2) = {19, 100, 73, 46};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"quad8", "iso20", "penta15", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   e1 = [1; 0.5; 0.4];
%!   e2 = [0; 1; 0];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   R = [e1, e2, e3];
%!   R *= diag(1 ./ norm(R, "cols"));
%!   R = eye(3);
%!   T = [R.', zeros(3, 3);
%!        zeros(3, 3), R.'];
%!   mesh.nodes *= T;
%!   mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, "quad8");
%!   mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, "tria6h");
%!   mesh.elements.joints(2).nodes = node_idx_bearing2;
%!   mesh.elements.joints(2).C = [0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0] * T;
%!   mesh.elements.joints(1).nodes = node_idx_bearing1;
%!   mesh.elements.joints(1).C = [1, 0, 0, 0, 0, 0;
%!                                0, 1, 0, 0, 0, 0;
%!                                0, 0, 1, 0, 0, 0;
%!                                0, 0, 0, 1, 0, 0] * T;
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!   for i=1:3
%!     mesh.materials.iso20([mesh.groups.iso20(find([[mesh.groups.iso20.id] == i])).elements]) = i;
%!     mesh.materials.penta15([mesh.groups.penta15(find([[mesh.groups.penta15.id] == i])).elements]) = i;
%!   endfor
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E12;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(2).rho = param.rho;
%!   mesh.material_data(2).E = param.E12;
%!   mesh.material_data(2).nu = param.nu;
%!   mesh.material_data(3).rho = param.rho;
%!   mesh.material_data(3).E = param.E3;
%!   mesh.material_data(3).nu = param.nu;
%!   dof_map = fem_ass_dof_map(mesh, load_case_dof);
%!   dof_map.parallel.threads_ass = opt_solver.number_of_threads;
%!   load_case = fem_pre_load_case_create_empty(columns(OMEGA));
%!   empty_cell = cell(1, columns(OMEGA));
%!   sol_stat = struct("def", empty_cell, "stress", empty_cell);
%!   sol_eig = struct("lambda", empty_cell, "f", empty_cell, "def", empty_cell, "D", empty_cell);
%!   mat_red = struct("Mred", empty_cell, "Kred", empty_cell, "Dred", empty_cell);
%!   for i=1:columns(OMEGA)
%!     fprintf(stderr, "n=%.0frpm\n", norm(OMEGA(:, i)) * 30 / pi / SI_unit_second);
%!     load_case(i).omega = R * OMEGA(:, i);
%!     load_case(i).omegadot = R * OMEGADOT(:, i);
%!     [mat_ass.M, ...
%!      mat_ass.K, ...
%!      mat_ass.KOMEGA, ...
%!      mat_ass.KOMEGA_DOT, ...
%!      mat_ass.DOMEGA, ...
%!      mat_ass.R, ...
%!      mat_ass.dm, ...
%!      mat_ass.S, ...
%!      mat_ass.J] = fem_ass_matrix(mesh, ...
%!                                  dof_map, ...
%!                                  [FEM_MAT_MASS, ...
%!                                   FEM_MAT_STIFFNESS, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA, ...
%!                                   FEM_MAT_STIFFNESS_OMEGA_DOT, ...
%!                                   FEM_MAT_DAMPING_OMEGA, ...
%!                                   FEM_VEC_LOAD_CONSISTENT, ...
%!                                   FEM_SCA_TOT_MASS, ...
%!                                   FEM_VEC_INERTIA_M1, ...
%!                                   FEM_MAT_INERTIA_J], ...
%!                                  load_case(i));
%!     Tred = zeros(rows(mesh.nodes) * 3, 6);
%!     for j=1:rows(mesh.nodes)
%!       Tred((j - 1) * 3 + (1:3), 1:3) = eye(3);
%!       Tred((j - 1) * 3 + (1:3), 4:6) = -skew(mesh.nodes(j, 1:3));
%!     endfor
%!     idx = dof_map.ndof(:, 1:3).'(:);
%!     mat_red(i).Mred = Tred.' * mat_ass.M(idx, idx) * Tred;
%!     mat_red(i).Kred = Tred.' * mat_ass.K(idx, idx) * Tred;
%!     mat_red(i).Dred = Tred.' * mat_ass.DOMEGA(idx, idx) * Tred;
%!     sol_stat(i).def = fem_sol_static(mesh, dof_map, mat_ass).def;
%!     sol_stat(i).stress = fem_ass_matrix(mesh, ...
%!                                         dof_map, ...
%!                                         [FEM_VEC_STRESS_CAUCH], ...
%!                                         load_case(i), ...
%!                                         sol_stat(i));
%!     load_case_pre_stress.tau0 = sol_stat(i).stress.tau;
%!     mat_ass.KTAU0 = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_MAT_STIFFNESS_TAU0], ...
%!                                    load_case_pre_stress);
%!     mat_ass.K += mat_ass.KOMEGA + mat_ass.KTAU0 + mat_ass.KOMEGA_DOT;
%!     mat_ass.D = mat_ass.DOMEGA;
%!     sol_eig(i) = fem_sol_modal_damped(mesh, ...
%!                                       dof_map, ...
%!                                       mat_ass, ...
%!                                       options.number_of_modes, ...
%!                                       opt_solver);
%!   endfor
%!   for i=1:columns(OMEGA)
%!     opt_mbd_mesh = struct();
%!     opt_mbd_mesh.joints.number = 2;
%!     opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!     opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!     opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!     opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!     load_case_empty = struct();
%!     opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!     unwind_protect
%!       [fd, msg] = fopen(mbdyn_file, "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s\"", mbdyn_file);
%!       endif
%!       fprintf(fd, " set: real E = %.16e;\n", param.E12);
%!       fprintf(fd, " set: real nu = %.16e;\n", param.nu);
%!       fprintf(fd, " set: real rho = %.16e;\n", param.rho);
%!       fprintf(fd, " set: real d1 = %.16e;\n", param.d1);
%!       fprintf(fd, " set: real d2 = %.16e;\n", param.d2);
%!       fprintf(fd, " set: real d3 = %.16e;\n", param.d3);
%!       fprintf(fd, " set: real l1 = %.16e;\n", param.l1);
%!       fprintf(fd, " set: real l2 = %.16e;\n", param.l2);
%!       fprintf(fd, " set: real l3 = %.16e;\n", param.l3);
%!       fputs(fd, " set: real G = E / (2. * (1 + nu));\n");
%!       fputs(fd, " set: real A1 = d1^2 * pi / 4.;\n");
%!       fputs(fd, " set: real A2 = d2^2 * pi / 4.;\n");
%!       fputs(fd, " set: real As1 = 9./10. * A1;\n");
%!       fputs(fd, " set: real As2 = 9./10. * A2;\n");
%!       fputs(fd, " set: real Iy1 = d1^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iy2 = d2^4 * pi / 64.;\n");
%!       fputs(fd, " set: real Iz1 = Iy1;\n");
%!       fputs(fd, " set: real Iz2 = Iy2;\n");
%!       fputs(fd, " set: real It1 = Iy1 + Iz1;\n");
%!       fputs(fd, " set: real It2 = Iy2 + Iz2;\n");
%!       fputs(fd, " set: real m = rho * pi / 4 * d3^2 * 2 * l3;\n");
%!       fputs(fd, " set: real Ja = m * (3 * d3^2 + 4 * (2 * l3)^2) / 48;\n");
%!       fputs(fd, " set: real Jp = m * d3^2 / 8;\n");
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGA%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGA(j, i));
%!       endfor
%!       for j=1:3
%!         fprintf(fd, " set: real OMEGADOT%s = %.16e;\n", {"x", "y", "z"}{j}, OMEGADOT(j, i));
%!       endfor
%!       fprintf(fd, " set: real t1 = %g;\n", 1000 / SI_unit_second);
%!       fputs(fd, " set: real N = 1000;\n");
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "    problem: initial value; # the default\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "    initial time: 0;\n");
%!       fputs(fd, "    final time: t1;\n");
%!       fputs(fd, "    time step: t1 / N;\n");
%!       fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!       fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!       fputs(fd, "    max iterations: 10000;\n");
%!       fputs(fd, "    tolerance: 1.e-3, 1e-3;\n");
%!       fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!       fputs(fd, "    method: bdf;\n");
%!       fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!       fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!       fputs(fd, "    eigenanalysis: t1,\n");
%!       fputs(fd, "    # output matrices, \n");
%!       fputs(fd, "    # parameter, 1e-3, ## use default estimate\n");
%!       fputs(fd, "    output eigenvectors,\n");
%!       fputs(fd, "         # output geometry,\n");
%!       fprintf(fd, "         lower frequency limit, %e,\n", 0.01 / (SI_unit_second^-1));
%!       fprintf(fd, "         upper frequency limit, %e,\n", 1000 / (SI_unit_second^-1));
%!       fprintf(fd, "    use arpack,%d,%d,0.,suffix format, \"%%02d\";\n", 2 * options.number_of_modes, options.number_of_modes * 10);
%!       fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-8, forcing term max tolerance, 1e-3, inner iterations before assembly, 30;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!       fputs(fd, "        use automatic differentiation;\n");
%!       fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!       fputs(fd, "        component,\n");
%!       fputs(fd, "           string, \"OMEGAx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGAy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGAz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "        angular acceleration,\n");
%!       fputs(fd, "        component,\n");
%!       fputs(fd, "           string, \"OMEGADOTx * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGADOTy * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\",\n");
%!       fputs(fd, "           string, \"OMEGADOTz * (1. - cos(pi/2 * (Time / t1 * (Time <= t1) + (Time > t1)))^2)\";\n");
%!       fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!       fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!       fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", nodes_file);
%!       fputs(fd, " end: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", csl_file);
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "    joint: 1, total pin joint,\n");
%!       fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fputs(fd, "    joint: 2, total pin joint,\n");
%!       fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fprintf(fd, "                            position,		reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!       fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!       fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!       fputs(fd, "                    position constraint, \n");
%!       fputs(fd, "                                    inactive, \n");
%!       fputs(fd, "                                    active, \n");
%!       fputs(fd, "                                    active,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                    orientation constraint,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                                    inactive,\n");
%!       fputs(fd, "                            component, const, 0.,\n");
%!       fputs(fd, "                                       const, 0.,\n");
%!       fputs(fd, "                                       const, 0.;\n");
%!       fprintf(fd, "include: \"%s\";\n", elem_file);
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!       fd = -1;
%!     end_unwind_protect
%!     if (~options.verbose)
%!       opt_mbdyn.logfile = [filename, ".stdout"];
%!     endif
%!     opt_mbdyn.output_file = [filename, "_mbdyn"];
%!     info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!     [mesh_sol(i), sol(i)] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!     modal(i).mbdyn = mbdyn_post_load_output_eig(opt_mbdyn.output_file);
%!     if (options.post_proc_modes)
%!     for j=1:numel(modal(i).mbdyn.f)
%!       opt_modal.mode_index = j;
%!       opt_modal.scale = 100;
%!       mode_file = [opt_mbdyn.output_file, sprintf("_%02d_%02d", i, j)];
%!       mbdyn_post_eig_to_mov_file(opt_mbdyn.output_file, [mode_file, ".mov"], opt_modal, modal(i).mbdyn);
%!       [err, msg] = symlink([opt_mbdyn.output_file, ".log"], [mode_file, ".log"]);
%!       if (err ~= 0)
%!         error("symlink failed with status %d: %s", err, msg);
%!       endif
%!       opt_post.f_run_mbdyn = false;
%!       opt_post.f_run_mbdyn2easyanim = true;
%!       opt_post.f_runEasyAnim = false;
%!       opt_post.every = 1;
%!       opt_post.showAll = 1;
%!       info = mbdyn_solver_run(mode_file, opt_post);
%!     endfor
%!     endif
%!     fref(:, i) = modal(i).mbdyn.f(1:rows(fref));
%!   endfor
%!   f = zeros(options.number_of_modes, numel(sol_eig));
%!   for i=1:columns(f)
%!     f(:, i) = sort(sol_eig(i).f(:));
%!   endfor
%!   tol = 2e-2;
%!   assert_simple(f(floor(end/2+1):end,:), fref, tol * max(max(abs(fref))));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST19
%! pkg load mboct-fem-pkg;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_kelvin = 1;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   nodes_file = [filename, "_mbd.nod"];
%!   elem_file = [filename, "_mbd.elm"];
%!   set_file = [filename, "_mbd.set"];
%!   csl_file = [filename, "_mbd.csl"];
%!   control_file = [filename, "_mbd.con"];
%!   initial_value_file = [filename, "_mbd.inv"];
%!   input_file = [filename, "_mbd_inp.mbdyn"];
%!   output_file = [filename, "_mbd_out"];
%!   options.verbose = false;
%!   opt_mbd.output_file = output_file;
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     L = 100e-3 / SI_unit_meter;
%!     b = 1e-3 / SI_unit_meter;
%!     H = 1e-3 / SI_unit_meter;
%!     h = 1e-3 / SI_unit_meter;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "L = %g;\n", L);
%!     fprintf(fd, "b = %g;\n", b);
%!     fprintf(fd, "H = %g;\n", H);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!     fputs(fd, "Point(2) = {L,0.0,0.0,h};\n");
%!     fputs(fd, "Point(3) = {L,b,0.0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!     fputs(fd, "Line(1) = {4,3};\n");
%!     fputs(fd, "Line(2) = {3,2};\n");
%!     fputs(fd, "Line(3) = {2,1};\n");
%!     fputs(fd, "Line(4) = {1,4};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "Transfinite Surface(6) = {};\n");
%!     fputs(fd, "tmp1[] = Extrude {0,0.0,H} { Surface{6}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "tmp2[] = Extrude {0,0.0,H} { Surface{tmp1[0]}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!     fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume1\",1) = {tmp1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",2) = {tmp2[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {tmp1[4],tmp2[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",4) = {tmp1[2],tmp2[2]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=0;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.iso27 = zeros(rows(mesh.elements.iso27), 1, "int32");
%!   grp_id_mat1 = find([mesh.groups.iso27.id] == 1);
%!   grp_id_mat2 = find([mesh.groups.iso27.id] == 2);
%!   for i=1:numel(grp_id_mat1)
%!     mesh.materials.iso27(mesh.groups.iso27(grp_id_mat1(i)).elements(:)) = 1;
%!   endfor
%!   for i=1:numel(grp_id_mat2)
%!     mesh.materials.iso27(mesh.groups.iso27(grp_id_mat2(i)).elements(:)) = 2;
%!   endfor
%!   Ec = 210000e6 / SI_unit_pascal;
%!   Ei = 125000e6 / SI_unit_pascal;
%!   CTEc = 12.5e-6 / SI_unit_kelvin^-1;
%!   CTEi = 16.7e-6 / SI_unit_kelvin^-1;
%!   dT = 100 / SI_unit_kelvin;
%!   K1 = 14 + (Ec/Ei) + (Ei/Ec);
%!   Uz_ref = 3*(CTEc - CTEi)*dT*2*H*L^2/(H^2*K1);
%!   epsilon1 = [CTEc * dT * ones(3, 1); zeros(3, 1)];
%!   epsilon2 = [CTEi * dT * ones(3, 1); zeros(3, 1)];
%!   mesh.material_data(1).E = Ec;
%!   mesh.material_data(1).nu = 0.3;
%!   mesh.material_data(1).rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(1).type = "hookean linear elastic isotropic";
%!   mesh.material_data(1).extra_data = ", prestrain, reference, tpl_drive_id_epsilon1";
%!   mesh.material_data(2).E = Ei;
%!   mesh.material_data(2).nu = 0.35;
%!   mesh.material_data(2).rho = 8900 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(2).type = "hookean linear elastic isotropic";
%!   mesh.material_data(2).extra_data = ", prestrain, reference, tpl_drive_id_epsilon2";
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case_dof.locked_dof(mesh.groups.quad9(find([[mesh.groups.quad9].id] == 3)).nodes, 1:3) = true;
%!   load_case = struct();
%!   if (~options.verbose)
%!     opt_mbd.logfile = [opt_mbd.output_file, ".stdout"];
%!   endif
%!   opt_mbd.mbdyn_command = "mbdyn -C";
%!   opt_mbd.f_run_mbdyn = true;
%!   fd = -1;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(set_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", set_file, msg);
%!     endif
%!     fprintf(fd, "set: integer number_of_nodes = %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "set: integer number_of_nodes_hydraulic = %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!     fprintf(fd, "set: integer number_of_solids = %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "set: integer number_of_genels = %d;\n", opt_mbd_mesh.genels.number);
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon1 = 2001;\n");
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon2 = 2002;\n");
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon1_%d = %.16e;\n", i, epsilon1(i));
%!     endfor
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon2_%d = %.16e;\n", i, epsilon2(i));
%!     endfor
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(control_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", control_file, msg);
%!     endif
%!     fprintf(fd, "model: static;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(initial_value_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", initial_value_file, msg);
%!     endif
%!     fprintf(fd, "method: implicit euler;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen(input_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", input_file, msg);
%!     endif
%!     fprintf(fd, "include: \"%s\";\n", set_file);
%!     fprintf(fd, "begin: data;\n");
%!     fprintf(fd, "        problem: initial value; # the default\n");
%!     fprintf(fd, "end: data;\n");
%!     fprintf(fd, "begin: initial value;\n");
%!     fprintf(fd, "        initial time: 0;\n");
%!     fprintf(fd, "        final time: 1;\n");
%!     fprintf(fd, "        time step: 0.2;\n");
%!     fprintf(fd, "        max iterations: 100;\n");
%!     fprintf(fd, "        tolerance: 1e-5, test, norm, 1e-5, test, norm;\n");
%!     fprintf(fd, "        output: messages;\n");
%!     fprintf(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fprintf(fd, "        nonlinear solver: nox,\n");
%!     fprintf(fd, "                          modified, 30,\n");
%!     fprintf(fd, "                          keep jacobian matrix,\n");
%!     fprintf(fd, "                          use preconditioner as solver, no,\n");
%!     fprintf(fd, "                          linesearch method, backtrack,\n");
%!     fprintf(fd, "                          direction, newton,\n");
%!     fprintf(fd, "                          jacobian operator, newton krylov,\n");
%!     fprintf(fd, "                          forcing term, type2,\n");
%!     fprintf(fd, "                          linear solver tolerance, 1e-12,\n");
%!     fprintf(fd, "                          inner iterations before assembly, 15,\n");
%!     fprintf(fd, "                          linear solver max iterations, 300,\n");
%!     fprintf(fd, "                          krylov subspace size, 300,\n");
%!     fprintf(fd, "                          minimum step, 1e-12,\n");
%!     fprintf(fd, "                          recovery step type, constant,\n");
%!     fprintf(fd, "                          recovery step, 1e-12,\n");
%!     fprintf(fd, "                          verbose, 3,\n");
%!     fprintf(fd, "                          print convergence info, no;\n");
%!     fprintf(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fprintf(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!     fprintf(fd, "        derivatives tolerance: 1e-5, 1e-5;\n");
%!     fprintf(fd, "        derivatives max iterations: 10;\n");
%!     fprintf(fd, "        threads: assembly, 1;\n");
%!     fprintf(fd, "        threads: solver, 1;\n");
%!     fprintf(fd, "        output: cpu time;\n");
%!     fprintf(fd, "        include: \"%s\";\n", initial_value_file);
%!     fprintf(fd, "end: initial value;\n");
%!     fprintf(fd, "begin: control data;\n");
%!     fprintf(fd, "       output meter: closest next, 1, forever, 0.05;\n");
%!     fprintf(fd, "       skip initial joint assembly;\n");
%!     fprintf(fd, "       output precision: 16;\n");
%!     fprintf(fd, "       include: \"%s\";\n", control_file);
%!     fprintf(fd, "       default output: all, solids, accelerations;\n");
%!     fprintf(fd, "       structural nodes: number_of_nodes;\n");
%!     fprintf(fd, "       hydraulic nodes: number_of_nodes_hydraulic;\n");
%!     fprintf(fd, "       solids: number_of_solids;\n");
%!     fprintf(fd, "       genels: number_of_genels;\n");
%!     fprintf(fd, "       use automatic differentiation;\n");
%!     fprintf(fd, "end: control data;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon1, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon1_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon2, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon2_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fprintf(fd, "begin: nodes;\n");
%!     fprintf(fd, "       include: \"%s\";\n", nodes_file);
%!     fprintf(fd, "end: nodes;\n");
%!     fprintf(fd, "begin: elements;\n");
%!     fprintf(fd, "       include: \"%s\";\n", elem_file);
%!     fprintf(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.verbose)
%!     shell(sprintf("cat \"%s\" | nl", set_file));
%!     shell(sprintf("cat \"%s\" | nl", input_file));
%!     shell(sprintf("cat \"%s\" | nl", nodes_file));
%!     shell(sprintf("cat \"%s\" | nl", csl_file));
%!     shell(sprintf("cat \"%s\" | nl", elem_file));
%!   endif
%!   info = mbdyn_solver_run(input_file, opt_mbd);
%!   [mesh_sol, sol_stat] = mbdyn_post_load_output_sol(output_file);
%!   grp_id_load = find([mesh.groups.quad9.id] == 4);
%!   Uz = mean(sol_stat.def(mesh.groups.quad9(grp_id_load).nodes, 3));
%!   tol = 1e-2;
%!   assert_simple(Uz, Uz_ref, tol * abs(Uz_ref));
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST20
%! pkg load mboct-fem-pkg;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_kelvin = 1;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   nodes_file = [filename, "_mbd.nod"];
%!   elem_file = [filename, "_mbd.elm"];
%!   set_file = [filename, "_mbd.set"];
%!   csl_file = [filename, "_mbd.csl"];
%!   control_file = [filename, "_mbd.con"];
%!   initial_value_file = [filename, "_mbd.inv"];
%!   input_file = [filename, "_mbd_inp.mbdyn"];
%!   output_file = [filename, "_mbd_out"];
%!   options.verbose = false;
%!   opt_mbd.output_file = output_file;
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     L = 100e-3 / SI_unit_meter;
%!     b = 1e-3 / 2 / SI_unit_meter;
%!     H = 1e-3 / SI_unit_meter;
%!     h = 1e-3 / 2 / SI_unit_meter;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "L = %g;\n", L);
%!     fprintf(fd, "b = %g;\n", b);
%!     fprintf(fd, "H = %g;\n", H);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!     fputs(fd, "Point(2) = {L,0.0,0.0,h};\n");
%!     fputs(fd, "Point(3) = {L,b,0.0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!     fputs(fd, "Line(1) = {4,3};\n");
%!     fputs(fd, "Line(2) = {3,2};\n");
%!     fputs(fd, "Line(3) = {2,1};\n");
%!     fputs(fd, "Line(4) = {1,4};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "Transfinite Surface(6) = {};\n");
%!     fputs(fd, "tmp1[] = Extrude {0,0.0,H} { Surface{6}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "tmp2[] = Extrude {0,0.0,H} { Surface{tmp1[0]}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!     fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume1\",1) = {tmp1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",2) = {tmp2[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {tmp1[4],tmp2[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",4) = {tmp1[2],tmp2[2]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   grp_id_mat1 = find([mesh.groups.iso20.id] == 1);
%!   grp_id_mat2 = find([mesh.groups.iso20.id] == 2);
%!   for i=1:numel(grp_id_mat1)
%!     mesh.materials.iso20(mesh.groups.iso20(grp_id_mat1(i)).elements(:)) = 1;
%!   endfor
%!   for i=1:numel(grp_id_mat2)
%!     mesh.materials.iso20(mesh.groups.iso20(grp_id_mat2(i)).elements(:)) = 2;
%!   endfor
%!   Ec = 210000e6 / SI_unit_pascal;
%!   Ei = 125000e6 / SI_unit_pascal;
%!   CTEc = 12.5e-6 / SI_unit_kelvin^-1;
%!   CTEi = 16.7e-6 / SI_unit_kelvin^-1;
%!   dT = 100 / SI_unit_kelvin;
%!   K1 = 14 + (Ec/Ei) + (Ei/Ec);
%!   Uz_ref = 3*(CTEc - CTEi)*dT*2*H*L^2/(H^2*K1);
%!   epsilon1 = [CTEc * dT * ones(3, 1); zeros(3, 1)];
%!   epsilon2 = [CTEi * dT * ones(3, 1); zeros(3, 1)];
%!   mesh.material_data(1).E = Ec;
%!   mesh.material_data(1).nu = 0.3;
%!   mesh.material_data(1).rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(1).type = "hookean linear elastic isotropic";
%!   mesh.material_data(1).extra_data = ", prestrain, reference, tpl_drive_id_epsilon1";
%!   mesh.material_data(2).E = Ei;
%!   mesh.material_data(2).nu = 0.35;
%!   mesh.material_data(2).rho = 8900 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(2).type = "hookean linear elastic isotropic";
%!   mesh.material_data(2).extra_data = ", prestrain, reference, tpl_drive_id_epsilon2";
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case_dof.locked_dof(mesh.groups.quad8(find([[mesh.groups.quad8].id] == 3)).nodes, 1:3) = true;
%!   load_case = struct();
%!   if (~options.verbose)
%!     opt_mbd.logfile = [opt_mbd.output_file, ".stdout"];
%!   endif
%!   opt_mbd.mbdyn_command = "mbdyn -C";
%!   opt_mbd.f_run_mbdyn = true;
%!   fd = -1;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(set_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", set_file, msg);
%!     endif
%!     fprintf(fd, "set: integer number_of_nodes = %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "set: integer number_of_nodes_hydraulic = %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!     fprintf(fd, "set: integer number_of_solids = %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "set: integer number_of_genels = %d;\n", opt_mbd_mesh.genels.number);
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon1 = 2001;\n");
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon2 = 2002;\n");
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon1_%d = %.16e;\n", i, epsilon1(i));
%!     endfor
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon2_%d = %.16e;\n", i, epsilon2(i));
%!     endfor
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(control_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", control_file, msg);
%!     endif
%!     fprintf(fd, "model: static;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(initial_value_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", initial_value_file, msg);
%!     endif
%!     fprintf(fd, "method: implicit euler;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen(input_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", input_file, msg);
%!     endif
%!     fprintf(fd, "include: \"%s\";\n", set_file);
%!     fprintf(fd, "begin: data;\n");
%!     fprintf(fd, "        problem: initial value; # the default\n");
%!     fprintf(fd, "end: data;\n");
%!     fprintf(fd, "begin: initial value;\n");
%!     fprintf(fd, "        initial time: 0;\n");
%!     fprintf(fd, "        final time: 1;\n");
%!     fprintf(fd, "        time step: 0.2;\n");
%!     fprintf(fd, "        max iterations: 100;\n");
%!     fprintf(fd, "        tolerance: 1e-5, test, norm, 1e-5, test, norm;\n");
%!     fprintf(fd, "        output: messages;\n");
%!     fprintf(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fprintf(fd, "        nonlinear solver: nox,\n");
%!     fprintf(fd, "                          modified, 30,\n");
%!     fprintf(fd, "                          keep jacobian matrix,\n");
%!     fprintf(fd, "                          use preconditioner as solver, no,\n");
%!     fprintf(fd, "                          linesearch method, backtrack,\n");
%!     fprintf(fd, "                          direction, newton,\n");
%!     fprintf(fd, "                          jacobian operator, newton krylov,\n");
%!     fprintf(fd, "                          forcing term, type2,\n");
%!     fprintf(fd, "                          linear solver tolerance, 1e-12,\n");
%!     fprintf(fd, "                          inner iterations before assembly, 15,\n");
%!     fprintf(fd, "                          linear solver max iterations, 300,\n");
%!     fprintf(fd, "                          krylov subspace size, 300,\n");
%!     fprintf(fd, "                          minimum step, 1e-12,\n");
%!     fprintf(fd, "                          recovery step type, constant,\n");
%!     fprintf(fd, "                          recovery step, 1e-12,\n");
%!     fprintf(fd, "                          verbose, 3,\n");
%!     fprintf(fd, "                          print convergence info, no;\n");
%!     fprintf(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fprintf(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!     fprintf(fd, "        derivatives tolerance: 1e-5, 1e-5;\n");
%!     fprintf(fd, "        derivatives max iterations: 10;\n");
%!     fprintf(fd, "        threads: assembly, 1;\n");
%!     fprintf(fd, "        threads: solver, 1;\n");
%!     fprintf(fd, "        output: cpu time;\n");
%!     fprintf(fd, "        include: \"%s\";\n", initial_value_file);
%!     fprintf(fd, "end: initial value;\n");
%!     fprintf(fd, "begin: control data;\n");
%!     fprintf(fd, "       output meter: closest next, 1, forever, 0.05;\n");
%!     fprintf(fd, "       skip initial joint assembly;\n");
%!     fprintf(fd, "       output precision: 16;\n");
%!     fprintf(fd, "       include: \"%s\";\n", control_file);
%!     fprintf(fd, "       default output: all, solids, accelerations;\n");
%!     fprintf(fd, "       structural nodes: number_of_nodes;\n");
%!     fprintf(fd, "       hydraulic nodes: number_of_nodes_hydraulic;\n");
%!     fprintf(fd, "       solids: number_of_solids;\n");
%!     fprintf(fd, "       genels: number_of_genels;\n");
%!     fprintf(fd, "       use automatic differentiation;\n");
%!     fprintf(fd, "end: control data;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon1, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon1_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon2, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon2_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fprintf(fd, "begin: nodes;\n");
%!     fprintf(fd, "       include: \"%s\";\n", nodes_file);
%!     fprintf(fd, "end: nodes;\n");
%!     fprintf(fd, "begin: elements;\n");
%!     fprintf(fd, "       include: \"%s\";\n", elem_file);
%!     fprintf(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.verbose)
%!     shell(sprintf("cat \"%s\" | nl", set_file));
%!     shell(sprintf("cat \"%s\" | nl", input_file));
%!     shell(sprintf("cat \"%s\" | nl", nodes_file));
%!     shell(sprintf("cat \"%s\" | nl", csl_file));
%!     shell(sprintf("cat \"%s\" | nl", elem_file));
%!   endif
%!   info = mbdyn_solver_run(input_file, opt_mbd);
%!   [mesh_sol, sol_stat] = mbdyn_post_load_output_sol(output_file);
%!   grp_id_load = find([mesh.groups.quad8.id] == 4);
%!   Uz = mean(sol_stat.def(mesh.groups.quad8(grp_id_load).nodes, 3));
%!   tol = 1e-2;
%!   assert_simple(Uz, Uz_ref, tol * abs(Uz_ref));
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST21
%! pkg load mboct-fem-pkg;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_kelvin = 1;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   nodes_file = [filename, "_mbd.nod"];
%!   elem_file = [filename, "_mbd.elm"];
%!   set_file = [filename, "_mbd.set"];
%!   csl_file = [filename, "_mbd.csl"];
%!   control_file = [filename, "_mbd.con"];
%!   initial_value_file = [filename, "_mbd.inv"];
%!   input_file = [filename, "_mbd_inp.mbdyn"];
%!   output_file = [filename, "_mbd_out"];
%!   options.verbose = false;
%!   opt_mbd.output_file = output_file;
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     L = 100e-3 / SI_unit_meter;
%!     b = 1e-3 / SI_unit_meter;
%!     H = 1e-3 / SI_unit_meter;
%!     h = 1e-3 / SI_unit_meter;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "L = %g;\n", L);
%!     fprintf(fd, "b = %g;\n", b);
%!     fprintf(fd, "H = %g;\n", H);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!     fputs(fd, "Point(2) = {L,0.0,0.0,h};\n");
%!     fputs(fd, "Point(3) = {L,b,0.0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!     fputs(fd, "Line(1) = {4,3};\n");
%!     fputs(fd, "Line(2) = {3,2};\n");
%!     fputs(fd, "Line(3) = {2,1};\n");
%!     fputs(fd, "Line(4) = {1,4};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "Transfinite Surface(6) = {};\n");
%!     fputs(fd, "tmp1[] = Extrude {0,0.0,H} { Surface{6}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "tmp2[] = Extrude {0,0.0,H} { Surface{tmp1[0]}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!     fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume1\",1) = {tmp1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",2) = {tmp2[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {tmp1[4],tmp2[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",4) = {tmp1[2],tmp2[2]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_msh.elem_type = {"quad8r", "iso20r"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.iso20r = zeros(rows(mesh.elements.iso20r), 1, "int32");
%!   grp_id_mat1 = find([mesh.groups.iso20r.id] == 1);
%!   grp_id_mat2 = find([mesh.groups.iso20r.id] == 2);
%!   for i=1:numel(grp_id_mat1)
%!     mesh.materials.iso20r(mesh.groups.iso20r(grp_id_mat1(i)).elements(:)) = 1;
%!   endfor
%!   for i=1:numel(grp_id_mat2)
%!     mesh.materials.iso20r(mesh.groups.iso20r(grp_id_mat2(i)).elements(:)) = 2;
%!   endfor
%!   Ec = 210000e6 / SI_unit_pascal;
%!   Ei = 125000e6 / SI_unit_pascal;
%!   CTEc = 12.5e-6 / SI_unit_kelvin^-1;
%!   CTEi = 16.7e-6 / SI_unit_kelvin^-1;
%!   dT = 100 / SI_unit_kelvin;
%!   K1 = 14 + (Ec/Ei) + (Ei/Ec);
%!   Uz_ref = 3*(CTEc - CTEi)*dT*2*H*L^2/(H^2*K1);
%!   epsilon1 = [CTEc * dT * ones(3, 1); zeros(3, 1)];
%!   epsilon2 = [CTEi * dT * ones(3, 1); zeros(3, 1)];
%!   mesh.material_data(1).E = Ec;
%!   mesh.material_data(1).nu = 0.3;
%!   mesh.material_data(1).rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(1).type = "hookean linear elastic isotropic";
%!   mesh.material_data(1).extra_data = ", prestrain, reference, tpl_drive_id_epsilon1";
%!   mesh.material_data(2).E = Ei;
%!   mesh.material_data(2).nu = 0.35;
%!   mesh.material_data(2).rho = 8900 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(2).type = "hookean linear elastic isotropic";
%!   mesh.material_data(2).extra_data = ", prestrain, reference, tpl_drive_id_epsilon2";
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case_dof.locked_dof(mesh.groups.quad8r(find([[mesh.groups.quad8r].id] == 3)).nodes, 1:3) = true;
%!   load_case = struct();
%!   if (~options.verbose)
%!     opt_mbd.logfile = [opt_mbd.output_file, ".stdout"];
%!   endif
%!   opt_mbd.mbdyn_command = "mbdyn -C";
%!   opt_mbd.f_run_mbdyn = true;
%!   fd = -1;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(set_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", set_file, msg);
%!     endif
%!     fprintf(fd, "set: integer number_of_nodes = %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "set: integer number_of_nodes_hydraulic = %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!     fprintf(fd, "set: integer number_of_solids = %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "set: integer number_of_genels = %d;\n", opt_mbd_mesh.genels.number);
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon1 = 2001;\n");
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon2 = 2002;\n");
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon1_%d = %.16e;\n", i, epsilon1(i));
%!     endfor
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon2_%d = %.16e;\n", i, epsilon2(i));
%!     endfor
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(control_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", control_file, msg);
%!     endif
%!     fprintf(fd, "model: static;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(initial_value_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", initial_value_file, msg);
%!     endif
%!     fprintf(fd, "method: implicit euler;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen(input_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", input_file, msg);
%!     endif
%!     fprintf(fd, "include: \"%s\";\n", set_file);
%!     fprintf(fd, "begin: data;\n");
%!     fprintf(fd, "        problem: initial value; # the default\n");
%!     fprintf(fd, "end: data;\n");
%!     fprintf(fd, "begin: initial value;\n");
%!     fprintf(fd, "        initial time: 0;\n");
%!     fprintf(fd, "        final time: 1;\n");
%!     fprintf(fd, "        time step: 0.2;\n");
%!     fprintf(fd, "        max iterations: 100;\n");
%!     fprintf(fd, "        tolerance: 1e-5, test, norm, 1e-5, test, norm;\n");
%!     fprintf(fd, "        output: messages;\n");
%!     fprintf(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fprintf(fd, "        nonlinear solver: nox,\n");
%!     fprintf(fd, "                          modified, 30,\n");
%!     fprintf(fd, "                          keep jacobian matrix,\n");
%!     fprintf(fd, "                          use preconditioner as solver, no,\n");
%!     fprintf(fd, "                          linesearch method, backtrack,\n");
%!     fprintf(fd, "                          direction, newton,\n");
%!     fprintf(fd, "                          jacobian operator, newton krylov,\n");
%!     fprintf(fd, "                          forcing term, type2,\n");
%!     fprintf(fd, "                          linear solver tolerance, 1e-12,\n");
%!     fprintf(fd, "                          inner iterations before assembly, 15,\n");
%!     fprintf(fd, "                          linear solver max iterations, 300,\n");
%!     fprintf(fd, "                          krylov subspace size, 300,\n");
%!     fprintf(fd, "                          minimum step, 1e-12,\n");
%!     fprintf(fd, "                          recovery step type, constant,\n");
%!     fprintf(fd, "                          recovery step, 1e-12,\n");
%!     fprintf(fd, "                          verbose, 3,\n");
%!     fprintf(fd, "                          print convergence info, no;\n");
%!     fprintf(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fprintf(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!     fprintf(fd, "        derivatives tolerance: 1e-5, 1e-5;\n");
%!     fprintf(fd, "        derivatives max iterations: 10;\n");
%!     fprintf(fd, "        threads: assembly, 1;\n");
%!     fprintf(fd, "        threads: solver, 1;\n");
%!     fprintf(fd, "        output: cpu time;\n");
%!     fprintf(fd, "        include: \"%s\";\n", initial_value_file);
%!     fprintf(fd, "end: initial value;\n");
%!     fprintf(fd, "begin: control data;\n");
%!     fprintf(fd, "       output meter: closest next, 1, forever, 0.05;\n");
%!     fprintf(fd, "       skip initial joint assembly;\n");
%!     fprintf(fd, "       output precision: 16;\n");
%!     fprintf(fd, "       include: \"%s\";\n", control_file);
%!     fprintf(fd, "       default output: all, solids, accelerations;\n");
%!     fprintf(fd, "       structural nodes: number_of_nodes;\n");
%!     fprintf(fd, "       hydraulic nodes: number_of_nodes_hydraulic;\n");
%!     fprintf(fd, "       solids: number_of_solids;\n");
%!     fprintf(fd, "       genels: number_of_genels;\n");
%!     fprintf(fd, "       use automatic differentiation;\n");
%!     fprintf(fd, "end: control data;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon1, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon1_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon2, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon2_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fprintf(fd, "begin: nodes;\n");
%!     fprintf(fd, "       include: \"%s\";\n", nodes_file);
%!     fprintf(fd, "end: nodes;\n");
%!     fprintf(fd, "begin: elements;\n");
%!     fprintf(fd, "       include: \"%s\";\n", elem_file);
%!     fprintf(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.verbose)
%!     shell(sprintf("cat \"%s\" | nl", set_file));
%!     shell(sprintf("cat \"%s\" | nl", input_file));
%!     shell(sprintf("cat \"%s\" | nl", nodes_file));
%!     shell(sprintf("cat \"%s\" | nl", csl_file));
%!     shell(sprintf("cat \"%s\" | nl", elem_file));
%!   endif
%!   info = mbdyn_solver_run(input_file, opt_mbd);
%!   [mesh_sol, sol_stat] = mbdyn_post_load_output_sol(output_file);
%!   grp_id_load = find([mesh.groups.quad8r.id] == 4);
%!   Uz = mean(sol_stat.def(mesh.groups.quad8r(grp_id_load).nodes, 3));
%!   tol = 1e-2;
%!   assert_simple(Uz, Uz_ref, tol * abs(Uz_ref));
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST22
%! pkg load mboct-fem-pkg;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_kelvin = 1;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   nodes_file = [filename, "_mbd.nod"];
%!   elem_file = [filename, "_mbd.elm"];
%!   set_file = [filename, "_mbd.set"];
%!   csl_file = [filename, "_mbd.csl"];
%!   control_file = [filename, "_mbd.con"];
%!   initial_value_file = [filename, "_mbd.inv"];
%!   input_file = [filename, "_mbd_inp.mbdyn"];
%!   output_file = [filename, "_mbd_out"];
%!   options.verbose = false;
%!   opt_mbd.output_file = output_file;
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     L = 100e-3 / SI_unit_meter;
%!     b = 1e-3 / 4 / SI_unit_meter;
%!     H = 1e-3 / SI_unit_meter;
%!     h = 1e-3 / 4 / SI_unit_meter;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "L = %g;\n", L);
%!     fprintf(fd, "b = %g;\n", b);
%!     fprintf(fd, "H = %g;\n", H);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!     fputs(fd, "Point(2) = {L,0.0,0.0,h};\n");
%!     fputs(fd, "Point(3) = {L,b,0.0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!     fputs(fd, "Line(1) = {4,3};\n");
%!     fputs(fd, "Line(2) = {3,2};\n");
%!     fputs(fd, "Line(3) = {2,1};\n");
%!     fputs(fd, "Line(4) = {1,4};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "Transfinite Surface(6) = {};\n");
%!     fputs(fd, "tmp1[] = Extrude {0,0.0,H} { Surface{6}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "tmp2[] = Extrude {0,0.0,H} { Surface{tmp1[0]}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!     fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume1\",1) = {tmp1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",2) = {tmp2[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {tmp1[4],tmp2[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",4) = {tmp1[2],tmp2[2]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=0;\n");
%!     fputs(fd, "Mesh.ElementOrder = 1;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_msh.elem_type = {"iso4", "iso8"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!   grp_id_mat1 = find([mesh.groups.iso8.id] == 1);
%!   grp_id_mat2 = find([mesh.groups.iso8.id] == 2);
%!   for i=1:numel(grp_id_mat1)
%!     mesh.materials.iso8(mesh.groups.iso8(grp_id_mat1(i)).elements(:)) = 1;
%!   endfor
%!   for i=1:numel(grp_id_mat2)
%!     mesh.materials.iso8(mesh.groups.iso8(grp_id_mat2(i)).elements(:)) = 2;
%!   endfor
%!   Ec = 210000e6 / SI_unit_pascal;
%!   Ei = 125000e6 / SI_unit_pascal;
%!   CTEc = 12.5e-6 / SI_unit_kelvin^-1;
%!   CTEi = 16.7e-6 / SI_unit_kelvin^-1;
%!   dT = 100 / SI_unit_kelvin;
%!   K1 = 14 + (Ec/Ei) + (Ei/Ec);
%!   Uz_ref = 3*(CTEc - CTEi)*dT*2*H*L^2/(H^2*K1);
%!   epsilon1 = [CTEc * dT * ones(3, 1); zeros(3, 1)];
%!   epsilon2 = [CTEi * dT * ones(3, 1); zeros(3, 1)];
%!   mesh.material_data(1).E = Ec;
%!   mesh.material_data(1).nu = 0.3;
%!   mesh.material_data(1).rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(1).type = "hookean linear elastic isotropic";
%!   mesh.material_data(1).extra_data = ", prestrain, reference, tpl_drive_id_epsilon1";
%!   mesh.material_data(2).E = Ei;
%!   mesh.material_data(2).nu = 0.35;
%!   mesh.material_data(2).rho = 8900 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(2).type = "hookean linear elastic isotropic";
%!   mesh.material_data(2).extra_data = ", prestrain, reference, tpl_drive_id_epsilon2";
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case_dof.locked_dof(mesh.groups.iso4(find([[mesh.groups.iso4].id] == 3)).nodes, 1:3) = true;
%!   load_case = struct();
%!   if (~options.verbose)
%!     opt_mbd.logfile = [opt_mbd.output_file, ".stdout"];
%!   endif
%!   opt_mbd.mbdyn_command = "mbdyn -C";
%!   opt_mbd.f_run_mbdyn = true;
%!   fd = -1;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(set_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", set_file, msg);
%!     endif
%!     fprintf(fd, "set: integer number_of_nodes = %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "set: integer number_of_nodes_hydraulic = %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!     fprintf(fd, "set: integer number_of_solids = %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "set: integer number_of_genels = %d;\n", opt_mbd_mesh.genels.number);
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon1 = 2001;\n");
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon2 = 2002;\n");
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon1_%d = %.16e;\n", i, epsilon1(i));
%!     endfor
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon2_%d = %.16e;\n", i, epsilon2(i));
%!     endfor
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(control_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", control_file, msg);
%!     endif
%!     fprintf(fd, "model: static;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(initial_value_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", initial_value_file, msg);
%!     endif
%!     fprintf(fd, "method: implicit euler;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen(input_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", input_file, msg);
%!     endif
%!     fprintf(fd, "include: \"%s\";\n", set_file);
%!     fprintf(fd, "begin: data;\n");
%!     fprintf(fd, "        problem: initial value; # the default\n");
%!     fprintf(fd, "end: data;\n");
%!     fprintf(fd, "begin: initial value;\n");
%!     fprintf(fd, "        initial time: 0;\n");
%!     fprintf(fd, "        final time: 1;\n");
%!     fprintf(fd, "        time step: 0.2;\n");
%!     fprintf(fd, "        max iterations: 100;\n");
%!     fprintf(fd, "        tolerance: 1e-5, test, norm, 1e-5, test, norm;\n");
%!     fprintf(fd, "        output: messages;\n");
%!     fprintf(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fprintf(fd, "        nonlinear solver: nox,\n");
%!     fprintf(fd, "                          modified, 30,\n");
%!     fprintf(fd, "                          keep jacobian matrix,\n");
%!     fprintf(fd, "                          use preconditioner as solver, no,\n");
%!     fprintf(fd, "                          linesearch method, backtrack,\n");
%!     fprintf(fd, "                          direction, newton,\n");
%!     fprintf(fd, "                          jacobian operator, newton krylov,\n");
%!     fprintf(fd, "                          forcing term, type2,\n");
%!     fprintf(fd, "                          linear solver tolerance, 1e-12,\n");
%!     fprintf(fd, "                          inner iterations before assembly, 15,\n");
%!     fprintf(fd, "                          linear solver max iterations, 300,\n");
%!     fprintf(fd, "                          krylov subspace size, 300,\n");
%!     fprintf(fd, "                          minimum step, 1e-12,\n");
%!     fprintf(fd, "                          recovery step type, constant,\n");
%!     fprintf(fd, "                          recovery step, 1e-12,\n");
%!     fprintf(fd, "                          verbose, 3,\n");
%!     fprintf(fd, "                          print convergence info, no;\n");
%!     fprintf(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fprintf(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!     fprintf(fd, "        derivatives tolerance: 1e-5, 1e-5;\n");
%!     fprintf(fd, "        derivatives max iterations: 10;\n");
%!     fprintf(fd, "        threads: assembly, 1;\n");
%!     fprintf(fd, "        threads: solver, 1;\n");
%!     fprintf(fd, "        output: cpu time;\n");
%!     fprintf(fd, "        include: \"%s\";\n", initial_value_file);
%!     fprintf(fd, "end: initial value;\n");
%!     fprintf(fd, "begin: control data;\n");
%!     fprintf(fd, "       output meter: closest next, 1, forever, 0.05;\n");
%!     fprintf(fd, "       skip initial joint assembly;\n");
%!     fprintf(fd, "       output precision: 16;\n");
%!     fprintf(fd, "       include: \"%s\";\n", control_file);
%!     fprintf(fd, "       default output: all, solids, accelerations;\n");
%!     fprintf(fd, "       structural nodes: number_of_nodes;\n");
%!     fprintf(fd, "       hydraulic nodes: number_of_nodes_hydraulic;\n");
%!     fprintf(fd, "       solids: number_of_solids;\n");
%!     fprintf(fd, "       genels: number_of_genels;\n");
%!     fprintf(fd, "       use automatic differentiation;\n");
%!     fprintf(fd, "end: control data;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon1, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon1_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon2, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon2_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fprintf(fd, "begin: nodes;\n");
%!     fprintf(fd, "       include: \"%s\";\n", nodes_file);
%!     fprintf(fd, "end: nodes;\n");
%!     fprintf(fd, "begin: elements;\n");
%!     fprintf(fd, "       include: \"%s\";\n", elem_file);
%!     fprintf(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.verbose)
%!     shell(sprintf("cat \"%s\" | nl", set_file));
%!     shell(sprintf("cat \"%s\" | nl", input_file));
%!     shell(sprintf("cat \"%s\" | nl", nodes_file));
%!     shell(sprintf("cat \"%s\" | nl", csl_file));
%!     shell(sprintf("cat \"%s\" | nl", elem_file));
%!   endif
%!   info = mbdyn_solver_run(input_file, opt_mbd);
%!   [mesh_sol, sol_stat] = mbdyn_post_load_output_sol(output_file);
%!   grp_id_load = find([mesh.groups.iso4.id] == 4);
%!   Uz = mean(sol_stat.def(mesh.groups.iso4(grp_id_load).nodes, 3));
%!   tol = 1e-2;
%!   assert_simple(Uz, Uz_ref, tol * abs(Uz_ref));
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST23
%! pkg load mboct-fem-pkg;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_kelvin = 1;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   nodes_file = [filename, "_mbd.nod"];
%!   elem_file = [filename, "_mbd.elm"];
%!   set_file = [filename, "_mbd.set"];
%!   csl_file = [filename, "_mbd.csl"];
%!   control_file = [filename, "_mbd.con"];
%!   initial_value_file = [filename, "_mbd.inv"];
%!   input_file = [filename, "_mbd_inp.mbdyn"];
%!   output_file = [filename, "_mbd_out"];
%!   options.verbose = false;
%!   opt_mbd.output_file = output_file;
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     L = 100e-3 / SI_unit_meter;
%!     b = 1e-3 / 2 / SI_unit_meter;
%!     H = 1e-3 / SI_unit_meter;
%!     h = 1e-3 / 2 / SI_unit_meter;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "L = %g;\n", L);
%!     fprintf(fd, "b = %g;\n", b);
%!     fprintf(fd, "H = %g;\n", H);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!     fputs(fd, "Point(2) = {L,0.0,0.0,h};\n");
%!     fputs(fd, "Point(3) = {L,b,0.0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!     fputs(fd, "Line(1) = {4,3};\n");
%!     fputs(fd, "Line(2) = {3,2};\n");
%!     fputs(fd, "Line(3) = {2,1};\n");
%!     fputs(fd, "Line(4) = {1,4};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "Transfinite Surface(6) = {};\n");
%!     fputs(fd, "tmp1[] = Extrude {0,0.0,H} { Surface{6}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "tmp2[] = Extrude {0,0.0,H} { Surface{tmp1[0]}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(\"volume1\",1) = {tmp1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",2) = {tmp2[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {tmp1[4],tmp2[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",4) = {tmp1[2],tmp2[2]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!   grp_id_mat1 = find([mesh.groups.penta15.id] == 1);
%!   grp_id_mat2 = find([mesh.groups.penta15.id] == 2);
%!   for i=1:numel(grp_id_mat1)
%!     mesh.materials.penta15(mesh.groups.penta15(grp_id_mat1(i)).elements(:)) = 1;
%!   endfor
%!   for i=1:numel(grp_id_mat2)
%!     mesh.materials.penta15(mesh.groups.penta15(grp_id_mat2(i)).elements(:)) = 2;
%!   endfor
%!   Ec = 210000e6 / SI_unit_pascal;
%!   Ei = 125000e6 / SI_unit_pascal;
%!   CTEc = 12.5e-6 / SI_unit_kelvin^-1;
%!   CTEi = 16.7e-6 / SI_unit_kelvin^-1;
%!   dT = 100 / SI_unit_kelvin;
%!   K1 = 14 + (Ec/Ei) + (Ei/Ec);
%!   Uz_ref = 3*(CTEc - CTEi)*dT*2*H*L^2/(H^2*K1);
%!   epsilon1 = [CTEc * dT * ones(3, 1); zeros(3, 1)];
%!   epsilon2 = [CTEi * dT * ones(3, 1); zeros(3, 1)];
%!   mesh.material_data(1).E = Ec;
%!   mesh.material_data(1).nu = 0.3;
%!   mesh.material_data(1).rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(1).type = "hookean linear elastic isotropic";
%!   mesh.material_data(1).extra_data = ", prestrain, reference, tpl_drive_id_epsilon1";
%!   mesh.material_data(2).E = Ei;
%!   mesh.material_data(2).nu = 0.35;
%!   mesh.material_data(2).rho = 8900 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(2).type = "hookean linear elastic isotropic";
%!   mesh.material_data(2).extra_data = ", prestrain, reference, tpl_drive_id_epsilon2";
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case_dof.locked_dof(mesh.groups.quad8(find([[mesh.groups.quad8].id] == 3)).nodes, 1:3) = true;
%!   load_case = struct();
%!   if (~options.verbose)
%!     opt_mbd.logfile = [opt_mbd.output_file, ".stdout"];
%!   endif
%!   opt_mbd.mbdyn_command = "mbdyn -C";
%!   opt_mbd.f_run_mbdyn = true;
%!   fd = -1;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(set_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", set_file, msg);
%!     endif
%!     fprintf(fd, "set: integer number_of_nodes = %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "set: integer number_of_nodes_hydraulic = %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!     fprintf(fd, "set: integer number_of_solids = %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "set: integer number_of_genels = %d;\n", opt_mbd_mesh.genels.number);
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon1 = 2001;\n");
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon2 = 2002;\n");
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon1_%d = %.16e;\n", i, epsilon1(i));
%!     endfor
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon2_%d = %.16e;\n", i, epsilon2(i));
%!     endfor
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(control_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", control_file, msg);
%!     endif
%!     fprintf(fd, "model: static;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(initial_value_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", initial_value_file, msg);
%!     endif
%!     fprintf(fd, "method: implicit euler;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen(input_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", input_file, msg);
%!     endif
%!     fprintf(fd, "include: \"%s\";\n", set_file);
%!     fprintf(fd, "begin: data;\n");
%!     fprintf(fd, "        problem: initial value; # the default\n");
%!     fprintf(fd, "end: data;\n");
%!     fprintf(fd, "begin: initial value;\n");
%!     fprintf(fd, "        initial time: 0;\n");
%!     fprintf(fd, "        final time: 1;\n");
%!     fprintf(fd, "        time step: 0.2;\n");
%!     fprintf(fd, "        max iterations: 100;\n");
%!     fprintf(fd, "        tolerance: 1e-5, test, norm, 1e-5, test, norm;\n");
%!     fprintf(fd, "        output: messages;\n");
%!     fprintf(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fprintf(fd, "        nonlinear solver: nox,\n");
%!     fprintf(fd, "                          modified, 30,\n");
%!     fprintf(fd, "                          keep jacobian matrix,\n");
%!     fprintf(fd, "                          use preconditioner as solver, no,\n");
%!     fprintf(fd, "                          linesearch method, backtrack,\n");
%!     fprintf(fd, "                          direction, newton,\n");
%!     fprintf(fd, "                          jacobian operator, newton krylov,\n");
%!     fprintf(fd, "                          forcing term, type2,\n");
%!     fprintf(fd, "                          linear solver tolerance, 1e-12,\n");
%!     fprintf(fd, "                          inner iterations before assembly, 15,\n");
%!     fprintf(fd, "                          linear solver max iterations, 300,\n");
%!     fprintf(fd, "                          krylov subspace size, 300,\n");
%!     fprintf(fd, "                          minimum step, 1e-12,\n");
%!     fprintf(fd, "                          recovery step type, constant,\n");
%!     fprintf(fd, "                          recovery step, 1e-12,\n");
%!     fprintf(fd, "                          verbose, 3,\n");
%!     fprintf(fd, "                          print convergence info, no;\n");
%!     fprintf(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fprintf(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!     fprintf(fd, "        derivatives tolerance: 1e-5, 1e-5;\n");
%!     fprintf(fd, "        derivatives max iterations: 10;\n");
%!     fprintf(fd, "        threads: assembly, 1;\n");
%!     fprintf(fd, "        threads: solver, 1;\n");
%!     fprintf(fd, "        output: cpu time;\n");
%!     fprintf(fd, "        include: \"%s\";\n", initial_value_file);
%!     fprintf(fd, "end: initial value;\n");
%!     fprintf(fd, "begin: control data;\n");
%!     fprintf(fd, "       output meter: closest next, 1, forever, 0.05;\n");
%!     fprintf(fd, "       skip initial joint assembly;\n");
%!     fprintf(fd, "       output precision: 16;\n");
%!     fprintf(fd, "       include: \"%s\";\n", control_file);
%!     fprintf(fd, "       default output: all, solids, accelerations;\n");
%!     fprintf(fd, "       structural nodes: number_of_nodes;\n");
%!     fprintf(fd, "       hydraulic nodes: number_of_nodes_hydraulic;\n");
%!     fprintf(fd, "       solids: number_of_solids;\n");
%!     fprintf(fd, "       genels: number_of_genels;\n");
%!     fprintf(fd, "       use automatic differentiation;\n");
%!     fprintf(fd, "end: control data;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon1, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon1_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon2, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon2_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fprintf(fd, "begin: nodes;\n");
%!     fprintf(fd, "       include: \"%s\";\n", nodes_file);
%!     fprintf(fd, "end: nodes;\n");
%!     fprintf(fd, "begin: elements;\n");
%!     fprintf(fd, "       include: \"%s\";\n", elem_file);
%!     fprintf(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.verbose)
%!     shell(sprintf("cat \"%s\" | nl", set_file));
%!     shell(sprintf("cat \"%s\" | nl", input_file));
%!     shell(sprintf("cat \"%s\" | nl", nodes_file));
%!     shell(sprintf("cat \"%s\" | nl", csl_file));
%!     shell(sprintf("cat \"%s\" | nl", elem_file));
%!   endif
%!   info = mbdyn_solver_run(input_file, opt_mbd);
%!   [mesh_sol, sol_stat] = mbdyn_post_load_output_sol(output_file);
%!   grp_id_load = find([mesh.groups.quad8.id] == 4);
%!   Uz = mean(sol_stat.def(mesh.groups.quad8(grp_id_load).nodes, 3));
%!   tol = 1e-2;
%!   assert_simple(Uz, Uz_ref, tol * abs(Uz_ref));
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST24
%! pkg load mboct-fem-pkg;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_kelvin = 1;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   nodes_file = [filename, "_mbd.nod"];
%!   elem_file = [filename, "_mbd.elm"];
%!   set_file = [filename, "_mbd.set"];
%!   csl_file = [filename, "_mbd.csl"];
%!   control_file = [filename, "_mbd.con"];
%!   initial_value_file = [filename, "_mbd.inv"];
%!   input_file = [filename, "_mbd_inp.mbdyn"];
%!   output_file = [filename, "_mbd_out"];
%!   options.verbose = false;
%!   opt_mbd.output_file = output_file;
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     L = 100e-3 / SI_unit_meter;
%!     b = 1e-3 / SI_unit_meter;
%!     H = 1e-3 / SI_unit_meter;
%!     h = 1e-3 / SI_unit_meter;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "L = %g;\n", L);
%!     fprintf(fd, "b = %g;\n", b);
%!     fprintf(fd, "H = %g;\n", H);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!     fputs(fd, "Point(2) = {L,0.0,0.0,h};\n");
%!     fputs(fd, "Point(3) = {L,b,0.0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!     fputs(fd, "Line(1) = {4,3};\n");
%!     fputs(fd, "Line(2) = {3,2};\n");
%!     fputs(fd, "Line(3) = {2,1};\n");
%!     fputs(fd, "Line(4) = {1,4};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "Transfinite Surface(6) = {};\n");
%!     fputs(fd, "tmp1[] = Extrude {0,0.0,H} { Surface{6}; };\n");
%!     fputs(fd, "tmp2[] = Extrude {0,0.0,H} { Surface{tmp1[0]}; };\n");
%!     fputs(fd, "Physical Volume(\"volume1\",1) = {tmp1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",2) = {tmp2[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {tmp1[4],tmp2[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",4) = {tmp1[2],tmp2[2]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=0;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_msh.elem_type = {"tria6h", "tet10h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.tet10h = zeros(rows(mesh.elements.tet10h), 1, "int32");
%!   grp_id_mat1 = find([mesh.groups.tet10h.id] == 1);
%!   grp_id_mat2 = find([mesh.groups.tet10h.id] == 2);
%!   for i=1:numel(grp_id_mat1)
%!     mesh.materials.tet10h(mesh.groups.tet10h(grp_id_mat1(i)).elements(:)) = 1;
%!   endfor
%!   for i=1:numel(grp_id_mat2)
%!     mesh.materials.tet10h(mesh.groups.tet10h(grp_id_mat2(i)).elements(:)) = 2;
%!   endfor
%!   Ec = 210000e6 / SI_unit_pascal;
%!   Ei = 125000e6 / SI_unit_pascal;
%!   CTEc = 12.5e-6 / SI_unit_kelvin^-1;
%!   CTEi = 16.7e-6 / SI_unit_kelvin^-1;
%!   dT = 100 / SI_unit_kelvin;
%!   K1 = 14 + (Ec/Ei) + (Ei/Ec);
%!   Uz_ref = 3*(CTEc - CTEi)*dT*2*H*L^2/(H^2*K1);
%!   epsilon1 = [CTEc * dT * ones(3, 1); zeros(3, 1)];
%!   epsilon2 = [CTEi * dT * ones(3, 1); zeros(3, 1)];
%!   mesh.material_data(1).E = Ec;
%!   mesh.material_data(1).nu = 0.3;
%!   mesh.material_data(1).rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(1).type = "hookean linear elastic isotropic";
%!   mesh.material_data(1).extra_data = ", prestrain, reference, tpl_drive_id_epsilon1";
%!   mesh.material_data(2).E = Ei;
%!   mesh.material_data(2).nu = 0.35;
%!   mesh.material_data(2).rho = 8900 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(2).type = "hookean linear elastic isotropic";
%!   mesh.material_data(2).extra_data = ", prestrain, reference, tpl_drive_id_epsilon2";
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case_dof.locked_dof(mesh.groups.tria6h(find([[mesh.groups.tria6h].id] == 3)).nodes, 1:3) = true;
%!   load_case = struct();
%!   if (~options.verbose)
%!     opt_mbd.logfile = [opt_mbd.output_file, ".stdout"];
%!   endif
%!   opt_mbd.mbdyn_command = "mbdyn -C";
%!   opt_mbd.f_run_mbdyn = true;
%!   fd = -1;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(set_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", set_file, msg);
%!     endif
%!     fprintf(fd, "set: integer number_of_nodes = %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "set: integer number_of_nodes_hydraulic = %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!     fprintf(fd, "set: integer number_of_solids = %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "set: integer number_of_genels = %d;\n", opt_mbd_mesh.genels.number);
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon1 = 2001;\n");
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon2 = 2002;\n");
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon1_%d = %.16e;\n", i, epsilon1(i));
%!     endfor
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon2_%d = %.16e;\n", i, epsilon2(i));
%!     endfor
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(control_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", control_file, msg);
%!     endif
%!     fprintf(fd, "model: static;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(initial_value_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", initial_value_file, msg);
%!     endif
%!     fprintf(fd, "method: implicit euler;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen(input_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", input_file, msg);
%!     endif
%!     fprintf(fd, "include: \"%s\";\n", set_file);
%!     fprintf(fd, "begin: data;\n");
%!     fprintf(fd, "        problem: initial value; # the default\n");
%!     fprintf(fd, "end: data;\n");
%!     fprintf(fd, "begin: initial value;\n");
%!     fprintf(fd, "        initial time: 0;\n");
%!     fprintf(fd, "        final time: 1;\n");
%!     fprintf(fd, "        time step: 0.2;\n");
%!     fprintf(fd, "        max iterations: 100;\n");
%!     fprintf(fd, "        tolerance: 1e-5, test, norm, 1e-5, test, norm;\n");
%!     fprintf(fd, "        output: messages;\n");
%!     fprintf(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fprintf(fd, "        nonlinear solver: nox,\n");
%!     fprintf(fd, "                          modified, 30,\n");
%!     fprintf(fd, "                          keep jacobian matrix,\n");
%!     fprintf(fd, "                          use preconditioner as solver, no,\n");
%!     fprintf(fd, "                          linesearch method, backtrack,\n");
%!     fprintf(fd, "                          direction, newton,\n");
%!     fprintf(fd, "                          jacobian operator, newton krylov,\n");
%!     fprintf(fd, "                          forcing term, type2,\n");
%!     fprintf(fd, "                          linear solver tolerance, 1e-12,\n");
%!     fprintf(fd, "                          inner iterations before assembly, 15,\n");
%!     fprintf(fd, "                          linear solver max iterations, 300,\n");
%!     fprintf(fd, "                          krylov subspace size, 300,\n");
%!     fprintf(fd, "                          minimum step, 1e-12,\n");
%!     fprintf(fd, "                          recovery step type, constant,\n");
%!     fprintf(fd, "                          recovery step, 1e-12,\n");
%!     fprintf(fd, "                          verbose, 3,\n");
%!     fprintf(fd, "                          print convergence info, no;\n");
%!     fprintf(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fprintf(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!     fprintf(fd, "        derivatives tolerance: 1e-5, 1e-5;\n");
%!     fprintf(fd, "        derivatives max iterations: 10;\n");
%!     fprintf(fd, "        threads: assembly, 1;\n");
%!     fprintf(fd, "        threads: solver, 1;\n");
%!     fprintf(fd, "        output: cpu time;\n");
%!     fprintf(fd, "        include: \"%s\";\n", initial_value_file);
%!     fprintf(fd, "end: initial value;\n");
%!     fprintf(fd, "begin: control data;\n");
%!     fprintf(fd, "       output meter: closest next, 1, forever, 0.05;\n");
%!     fprintf(fd, "       skip initial joint assembly;\n");
%!     fprintf(fd, "       output precision: 16;\n");
%!     fprintf(fd, "       include: \"%s\";\n", control_file);
%!     fprintf(fd, "       default output: all, solids, accelerations;\n");
%!     fprintf(fd, "       structural nodes: number_of_nodes;\n");
%!     fprintf(fd, "       hydraulic nodes: number_of_nodes_hydraulic;\n");
%!     fprintf(fd, "       solids: number_of_solids;\n");
%!     fprintf(fd, "       genels: number_of_genels;\n");
%!     fprintf(fd, "       use automatic differentiation;\n");
%!     fprintf(fd, "end: control data;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon1, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon1_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon2, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon2_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fprintf(fd, "begin: nodes;\n");
%!     fprintf(fd, "       include: \"%s\";\n", nodes_file);
%!     fprintf(fd, "end: nodes;\n");
%!     fprintf(fd, "begin: elements;\n");
%!     fprintf(fd, "       include: \"%s\";\n", elem_file);
%!     fprintf(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.verbose)
%!     shell(sprintf("cat \"%s\" | nl", set_file));
%!     shell(sprintf("cat \"%s\" | nl", input_file));
%!     shell(sprintf("cat \"%s\" | nl", nodes_file));
%!     shell(sprintf("cat \"%s\" | nl", csl_file));
%!     shell(sprintf("cat \"%s\" | nl", elem_file));
%!   endif
%!   info = mbdyn_solver_run(input_file, opt_mbd);
%!   [mesh_sol, sol_stat] = mbdyn_post_load_output_sol(output_file);
%!   grp_id_load = find([mesh.groups.tria6h.id] == 4);
%!   Uz = mean(sol_stat.def(mesh.groups.tria6h(grp_id_load).nodes, 3));
%!   tol = 1e-2;
%!   assert_simple(Uz, Uz_ref, tol * abs(Uz_ref));
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST25
%! pkg load mboct-fem-pkg;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_kelvin = 1;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   nodes_file = [filename, "_mbd.nod"];
%!   elem_file = [filename, "_mbd.elm"];
%!   set_file = [filename, "_mbd.set"];
%!   csl_file = [filename, "_mbd.csl"];
%!   control_file = [filename, "_mbd.con"];
%!   initial_value_file = [filename, "_mbd.inv"];
%!   input_file = [filename, "_mbd_inp.mbdyn"];
%!   output_file = [filename, "_mbd_out"];
%!   options.verbose = false;
%!   opt_mbd.output_file = output_file;
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     L = 100e-3 / SI_unit_meter;
%!     b = 1e-3/2 / SI_unit_meter;
%!     H = 1e-3 / SI_unit_meter;
%!     h = 1e-3/2 / SI_unit_meter;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "L = %g;\n", L);
%!     fprintf(fd, "b = %g;\n", b);
%!     fprintf(fd, "H = %g;\n", H);
%!     fprintf(fd, "h = %g;\n", h);
%!     fputs(fd, "Point(1) = {0.0,0.0,0.0,h};\n");
%!     fputs(fd, "Point(2) = {L,0.0,0.0,h};\n");
%!     fputs(fd, "Point(3) = {L,b,0.0,h};\n");
%!     fputs(fd, "Point(4) = {0,b,0.0,h};\n");
%!     fputs(fd, "Line(1) = {4,3};\n");
%!     fputs(fd, "Line(2) = {3,2};\n");
%!     fputs(fd, "Line(3) = {2,1};\n");
%!     fputs(fd, "Line(4) = {1,4};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(L / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(b / h)) + 1;\n");
%!     fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!     fputs(fd, "Plane Surface(6) = {5};\n");
%!     fputs(fd, "Transfinite Surface(6) = {};\n");
%!     fputs(fd, "tmp1[] = Extrude {0,0.0,H} { Surface{6}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "tmp2[] = Extrude {0,0.0,H} { Surface{tmp1[0]}; Layers{Max(1, Round(H/h))}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{6, tmp1[0]};\n");
%!     fputs(fd, "Recombine Surface{tmp1[0], tmp2[0]};\n");
%!     fputs(fd, "Physical Volume(\"volume1\",1) = {tmp1[1]};\n");
%!     fputs(fd, "Physical Volume(\"volume2\",2) = {tmp2[1]};\n");
%!     fputs(fd, "Physical Surface(\"clamp\",3) = {tmp1[4],tmp2[4]};\n");
%!     fputs(fd, "Physical Surface(\"load\",4) = {tmp1[2],tmp2[2]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_msh.elem_type = {"quad8r", "iso20upcr"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   mesh.materials.iso20upcr = zeros(rows(mesh.elements.iso20upcr), 1, "int32");
%!   grp_id_mat1 = find([mesh.groups.iso20upcr.id] == 1);
%!   grp_id_mat2 = find([mesh.groups.iso20upcr.id] == 2);
%!   for i=1:numel(grp_id_mat1)
%!     mesh.materials.iso20upcr(mesh.groups.iso20upcr(grp_id_mat1(i)).elements(:)) = 1;
%!   endfor
%!   for i=1:numel(grp_id_mat2)
%!     mesh.materials.iso20upcr(mesh.groups.iso20upcr(grp_id_mat2(i)).elements(:)) = 2;
%!   endfor
%!   Ec = 210000e6 / SI_unit_pascal;
%!   Ei = 125000e6 / SI_unit_pascal;
%!   CTEc = 12.5e-6 / SI_unit_kelvin^-1;
%!   CTEi = 16.7e-6 / SI_unit_kelvin^-1;
%!   dT = 100 / SI_unit_kelvin;
%!   K1 = 14 + (Ec/Ei) + (Ei/Ec);
%!   Uz_ref = 3*(CTEc - CTEi)*dT*2*H*L^2/(H^2*K1);
%!   epsilon1 = [CTEc * dT * ones(3, 1); zeros(3, 1)];
%!   epsilon2 = [CTEi * dT * ones(3, 1); zeros(3, 1)];
%!   mesh.material_data(1).E = Ec;
%!   mesh.material_data(1).nu = 0.3;
%!   mesh.material_data(1).rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(1).type = "hookean linear elastic isotropic";
%!   mesh.material_data(1).extra_data = ", prestrain, reference, tpl_drive_id_epsilon1";
%!   mesh.material_data(2).E = Ei;
%!   mesh.material_data(2).nu = 0.35;
%!   mesh.material_data(2).rho = 8900 / (SI_unit_kilogram / SI_unit_meter^3);
%!   mesh.material_data(2).type = "hookean linear elastic isotropic";
%!   mesh.material_data(2).extra_data = ", prestrain, reference, tpl_drive_id_epsilon2";
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case_dof.locked_dof(mesh.groups.quad8r(find([[mesh.groups.quad8r].id] == 3)).nodes, 1:3) = true;
%!   load_case = struct();
%!   if (~options.verbose)
%!     opt_mbd.logfile = [opt_mbd.output_file, ".stdout"];
%!   endif
%!   opt_mbd.mbdyn_command = "mbdyn -C";
%!   opt_mbd.f_run_mbdyn = true;
%!   fd = -1;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(set_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", set_file, msg);
%!     endif
%!     fprintf(fd, "set: integer number_of_nodes = %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "set: integer number_of_nodes_hydraulic = %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!     fprintf(fd, "set: integer number_of_solids = %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "set: integer number_of_genels = %d;\n", opt_mbd_mesh.genels.number);
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon1 = 2001;\n");
%!     fprintf(fd, "set: integer tpl_drive_id_epsilon2 = 2002;\n");
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon1_%d = %.16e;\n", i, epsilon1(i));
%!     endfor
%!     for i=1:6
%!       fprintf(fd, "set: real epsilon2_%d = %.16e;\n", i, epsilon2(i));
%!     endfor
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(control_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", control_file, msg);
%!     endif
%!     fprintf(fd, "model: static;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   unwind_protect
%!     [fd, msg] = fopen(initial_value_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", initial_value_file, msg);
%!     endif
%!     fprintf(fd, "method: implicit euler;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen(input_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", input_file, msg);
%!     endif
%!     fprintf(fd, "include: \"%s\";\n", set_file);
%!     fprintf(fd, "begin: data;\n");
%!     fprintf(fd, "        problem: initial value; # the default\n");
%!     fprintf(fd, "end: data;\n");
%!     fprintf(fd, "begin: initial value;\n");
%!     fprintf(fd, "        initial time: 0;\n");
%!     fprintf(fd, "        final time: 1;\n");
%!     fprintf(fd, "        time step: 0.2;\n");
%!     fprintf(fd, "        max iterations: 100;\n");
%!     fprintf(fd, "        tolerance: 1e-5, test, norm, 1e-5, test, norm;\n");
%!     fprintf(fd, "        output: messages;\n");
%!     fprintf(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fprintf(fd, "        nonlinear solver: nox,\n");
%!     fprintf(fd, "                          modified, 30,\n");
%!     fprintf(fd, "                          keep jacobian matrix,\n");
%!     fprintf(fd, "                          use preconditioner as solver, no,\n");
%!     fprintf(fd, "                          linesearch method, backtrack,\n");
%!     fprintf(fd, "                          direction, newton,\n");
%!     fprintf(fd, "                          jacobian operator, newton krylov,\n");
%!     fprintf(fd, "                          forcing term, type2,\n");
%!     fprintf(fd, "                          linear solver tolerance, 1e-12,\n");
%!     fprintf(fd, "                          inner iterations before assembly, 15,\n");
%!     fprintf(fd, "                          linear solver max iterations, 300,\n");
%!     fprintf(fd, "                          krylov subspace size, 300,\n");
%!     fprintf(fd, "                          minimum step, 1e-12,\n");
%!     fprintf(fd, "                          recovery step type, constant,\n");
%!     fprintf(fd, "                          recovery step, 1e-12,\n");
%!     fprintf(fd, "                          verbose, 3,\n");
%!     fprintf(fd, "                          print convergence info, no;\n");
%!     fprintf(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!     fprintf(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!     fprintf(fd, "        derivatives tolerance: 1e-5, 1e-5;\n");
%!     fprintf(fd, "        derivatives max iterations: 10;\n");
%!     fprintf(fd, "        threads: assembly, 1;\n");
%!     fprintf(fd, "        threads: solver, 1;\n");
%!     fprintf(fd, "        output: cpu time;\n");
%!     fprintf(fd, "        include: \"%s\";\n", initial_value_file);
%!     fprintf(fd, "end: initial value;\n");
%!     fprintf(fd, "begin: control data;\n");
%!     fprintf(fd, "       output meter: closest next, 1, forever, 0.05;\n");
%!     fprintf(fd, "       skip initial joint assembly;\n");
%!     fprintf(fd, "       output precision: 16;\n");
%!     fprintf(fd, "       include: \"%s\";\n", control_file);
%!     fprintf(fd, "       default output: all, solids, accelerations;\n");
%!     fprintf(fd, "       structural nodes: number_of_nodes;\n");
%!     fprintf(fd, "       hydraulic nodes: number_of_nodes_hydraulic;\n");
%!     fprintf(fd, "       solids: number_of_solids;\n");
%!     fprintf(fd, "       genels: number_of_genels;\n");
%!     fprintf(fd, "       use automatic differentiation;\n");
%!     fprintf(fd, "end: control data;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon1, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon1_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "template drive caller: tpl_drive_id_epsilon2, 6,\n");
%!     fprintf(fd, " green lagrange strain, single,\n");
%!     for i=1:6
%!       fprintf(fd, "  epsilon2_%d,\n", i);
%!     endfor
%!     fprintf(fd, "  time;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fprintf(fd, "begin: nodes;\n");
%!     fprintf(fd, "       include: \"%s\";\n", nodes_file);
%!     fprintf(fd, "end: nodes;\n");
%!     fprintf(fd, "begin: elements;\n");
%!     fprintf(fd, "       include: \"%s\";\n", elem_file);
%!     fprintf(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (options.verbose)
%!     shell(sprintf("cat \"%s\" | nl", set_file));
%!     shell(sprintf("cat \"%s\" | nl", input_file));
%!     shell(sprintf("cat \"%s\" | nl", nodes_file));
%!     shell(sprintf("cat \"%s\" | nl", csl_file));
%!     shell(sprintf("cat \"%s\" | nl", elem_file));
%!   endif
%!   info = mbdyn_solver_run(input_file, opt_mbd);
%!   [mesh_sol, sol_stat] = mbdyn_post_load_output_sol(output_file);
%!   grp_id_load = find([mesh.groups.quad8r.id] == 4);
%!   Uz = mean(sol_stat.def(mesh.groups.quad8r(grp_id_load).nodes, 3));
%!   tol = 1e-2;
%!   assert_simple(Uz, Uz_ref, tol * abs(Uz_ref));
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
