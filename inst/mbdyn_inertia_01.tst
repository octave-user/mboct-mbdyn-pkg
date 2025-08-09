%!test
%! try
%!   ## TEST1
%!   status = pkg("list", "netcdf");
%!   if (~status{1}.loaded)
%!     pkg load netcdf;
%!   endif
%!   param.Vx = 0;
%!   param.Vy = 0;
%!   param.Vz = 0;
%!   param.Wx = 2. * pi;
%!   param.Wy = pi;
%!   param.Wz = pi/2.;
%!   param.a = 30e-3;
%!   param.b = 10e-3;
%!   param.c = 5e-3;
%!   param.x0 = 10e-3;
%!   param.y0 = 12e-3;
%!   param.z0 = 15e-3;
%!   param.ox = 0.5 * param.a;
%!   param.oy = 0.;
%!   param.oz = 0.;
%!   param.dx = 50e-3;
%!   param.dy = 60e-3;
%!   param.dz = 70e-3;
%!   param.phix = 30. * pi / 180.;
%!   param.phiy = 45. * pi / 180.;
%!   param.phiz = 70. * pi / 180.;
%!   param.rho = 7850.;
%!   param.m = param.rho * param.a * param.b * param.c;
%!   param.Jxx = param.m * (param.c^2+param.b^2)/12.;
%!   param.Jyy = param.m * (param.c^2+param.a^2)/12.;
%!   param.Jzz = param.m * (param.b^2+param.a^2)/12.;
%!   node_types = {"modal", "dynamic"};
%!   imethods = {"impliciteuler", "cranknicolson", "ms2,0.6", "ms3,0.6", "ms4,0.6", "ss2,0.6", "ss3,0.6", "ss4,0.6", "hope,0.6", "Bathe,0.6", "msstc3,0.6", "msstc4,0.6", "msstc5,0.6", "mssth3,0.6", "mssth4,0.6", "mssth5,0.6", "DIRK33", "DIRK43", "DIRK54", "hybrid,ms,0.6"};
%!   stages  = [              1,               1,         1,         1,         1,         1,         1,         1,          1,           2,            3,            4,            5,            3,            4,            5,        3,        4,        5,               1];
%!   info = repmat(struct(), numel(node_types), numel(imethods));
%!   err_gamma = err_beta = zeros(numel(node_types), numel(imethods));
%!   for idx_node_type=1:numel(node_types)
%!     for idx_method=1:numel(imethods)
%!       fd = -1;
%! %unwind_protect
%!         unwind_protect
%!           [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_inertia_XXXXXX"));
%!           if (fd == -1)
%!             error("failed to open temporary file");
%!           endif
%!           mbdyn_pre_write_param_file(fd, param);
%!           fputs(fd, " begin: data;\n");
%!           fputs(fd, "         problem: initial value; # the default\n");
%!           fputs(fd, " end: data;\n");
%!           fputs(fd, " begin: initial value;\n");
%!           fputs(fd, "     threads: assembly, 1;\n");
%!           fputs(fd, "         initial time: 0;\n");
%!           fputs(fd, "         final time: 1;\n");
%!           fprintf(fd, "         time step: 1e-5 * %d;\n", stages(idx_method));
%!           fputs(fd, "     max iterations: 10;\n");
%!           fputs(fd, "     linear solver: naive, colamd;\n");
%!           fputs(fd, "     nonlinear solver: nox,\n");
%!           fputs(fd, "               modified, 100,\n");
%!           fputs(fd, "               keep jacobian matrix,\n");
%!           fputs(fd, "               use preconditioner as solver, no,\n");
%!           fputs(fd, "               linesearch method, backtrack,\n");
%!           fputs(fd, "               direction, newton,\n");
%!           fputs(fd, "               jacobian operator, newton krylov,\n");
%!           fputs(fd, "               forcing term, type 2,\n");
%!           fputs(fd, "               forcing term max tolerance, 1e-8,\n");
%!           fputs(fd, "               forcing term min tolerance, 1e-10,\n");
%!           fputs(fd, "               linear solver tolerance, 1e-8,\n");
%!           fputs(fd, "               inner iterations before assembly, 25,\n");
%!           fputs(fd, "               linear solver max iterations, 50,\n");
%!           fputs(fd, "               linear solver, gmres,\n");
%!           fputs(fd, "               krylov subspace size, 50,\n");
%!           fputs(fd, "               minimum step, 1e-6,\n");
%!           fputs(fd, "               recovery step type, constant,\n");
%!           fputs(fd, "               recovery step, 1e-6,\n");
%!           fputs(fd, "               verbose, 0,\n");
%!           fputs(fd, "               print convergence info, no;\n");
%!           fputs(fd, "         tolerance: 1e-8, test, norm, 1e-8, test, norm;\n");
%!           fputs(fd, "         derivatives tolerance: 1e-6, 1e-6;\n");
%!           fputs(fd, "         derivatives max iterations: 20;\n");
%!           fputs(fd, "         derivatives coefficient: 1e-6;\n");
%!           fprintf(fd, "     method: %s;\n", imethods{idx_method});
%!           fputs(fd, " end: initial value;\n");
%!           fputs(fd, " begin: control data;\n");
%!           fputs(fd, "     output meter: closest next, 0., forever, 0.05;\n");
%!           fputs(fd, "     print: dof description;\n");
%!           fputs(fd, "     print: equation description;\n");
%!           fputs(fd, "     output precision: 16;\n");
%!           fputs(fd, "     default orientation: orientation vector;\n");
%!           fputs(fd, "     default output: all;\n");
%!           fputs(fd, "     output results: netcdf, text;\n");
%!           fputs(fd, "     use automatic differentiation;\n");
%!           fputs(fd, "        tolerance: 1e-6;\n");
%!           fputs(fd, "        max iterations: 10;\n");
%!           fputs(fd, "         structural nodes: 2;\n");
%!           fputs(fd, "     rigid bodies: 1;\n");
%!           fputs(fd, "     inertia: 1;\n");
%!           fputs(fd, " end: control data;\n");
%!           fputs(fd, " reference: 0,\n");,
%!           fputs(fd, "     position, reference, global, x0, y0, z0,\n");,
%!           fputs(fd, "     orientation, reference, global, eye,\n");,
%!           fputs(fd, "     velocity, reference, global, null,\n");,
%!           fputs(fd, "     angular velocity, reference, global, null;\n");
%!           fputs(fd, " reference: 1,\n");,
%!           fputs(fd, "     position, reference, 0, null,\n");,
%!           fputs(fd, "     orientation, reference, 0, eye,\n");,
%!           fputs(fd, "     velocity, reference, 0, Vx, Vy, Vz,\n");,
%!           fputs(fd, "     angular velocity, reference, 0, Wx, Wy, Wz;\n");
%!           fputs(fd, " reference: 2,\n");,
%!           fputs(fd, "     position, reference, 1, ox, oy, oz,\n");,
%!           fputs(fd, "     orientation, reference, 1, eye,\n");,
%!           fputs(fd, "     velocity, reference, 1, null,\n");,
%!           fputs(fd, "     angular velocity, reference, 1, null;\n");
%!           fputs(fd, " reference: 3,\n");,
%!           fputs(fd, "     position, reference, 1, dx, dy, dz,\n");,
%!           fputs(fd, "     orientation, reference, 1, euler123, phix, phiy, phiz,\n");,
%!           fputs(fd, "     velocity, reference, 1, null,\n");,
%!           fputs(fd, "     angular velocity, reference, 1, null;\n");
%!           fputs(fd, " begin: nodes;\n");
%!           fprintf(fd, "         structural: 1, %s,\n", node_types{idx_node_type});,
%!           fputs(fd, "                 position, reference, 3, null,\n");,
%!           fputs(fd, "                 orientation, reference, 3, eye,\n");,
%!           fputs(fd, "                 velocity, reference, 3, null,\n");,
%!           fputs(fd, "                 angular velocity, reference, 3, null;\n");
%!           fputs(fd, "         structural: 2, dummy, 1, offset,\n");,
%!           fputs(fd, "                 reference, 2, null, reference, 2, eye;\n");
%!           fputs(fd, " end: nodes;\n");
%!           fputs(fd, " begin: elements;\n");
%!           fputs(fd, "        body: 1, 1, m, reference, 2, null, diag, Jxx, Jyy, Jzz, orientation, reference, 2, eye;\n");
%!           fputs(fd, "        inertia: 1, name, \"bodies1\", position, reference, 2, null, orientation, reference, 2, eye, body, 1, output, always, output, both;\n");
%!           fputs(fd, "        drive caller: 1, element, 1, body, string, \"E\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 100, element, 1, inertia, string, \"X[1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 200, element, 1, inertia, string, \"X[2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 300, element, 1, inertia, string, \"X[3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 400, element, 1, inertia, string, \"V[1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 500, element, 1, inertia, string, \"V[2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 600, element, 1, inertia, string, \"V[3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 700, element, 1, inertia, string, \"Omega[1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 800, element, 1, inertia, string, \"Omega[2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 900, element, 1, inertia, string, \"Omega[3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 1000, element, 1, inertia, string, \"JP[1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 1100, element, 1, inertia, string, \"JP[2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 1200, element, 1, inertia, string, \"JP[3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 1300, element, 1, inertia, string, \"J[1,1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 1400, element, 1, inertia, string, \"J[2,1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 1500, element, 1, inertia, string, \"J[3,1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 1600, element, 1, inertia, string, \"J[1,2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 1700, element, 1, inertia, string, \"J[2,2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 1800, element, 1, inertia, string, \"J[3,2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 1900, element, 1, inertia, string, \"J[1,3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 2000, element, 1, inertia, string, \"J[2,3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 2100, element, 1, inertia, string, \"J[3,3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 2200, element, 1, inertia, string, \"m\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 2300, element, 1, inertia, string, \"beta[1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 2400, element, 1, inertia, string, \"beta[2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 2500, element, 1, inertia, string, \"beta[3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 2600, element, 1, inertia, string, \"gamma[1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 2700, element, 1, inertia, string, \"gamma[2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 2800, element, 1, inertia, string, \"gamma[3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 2900, element, 1, inertia, string, \"Jcg[1,1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 3000, element, 1, inertia, string, \"Jcg[2,1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 3100, element, 1, inertia, string, \"Jcg[3,1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 3200, element, 1, inertia, string, \"Jcg[1,2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 3300, element, 1, inertia, string, \"Jcg[2,2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 3400, element, 1, inertia, string, \"Jcg[3,2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 3500, element, 1, inertia, string, \"Jcg[1,3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 3600, element, 1, inertia, string, \"Jcg[2,3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 3700, element, 1, inertia, string, \"Jcg[3,3]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 3800, element, 1, inertia, string, \"Phi[1]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 3900, element, 1, inertia, string, \"Phi[2]\", direct, output, yes;\n");
%!           fputs(fd, "        drive caller: 4000, element, 1, inertia, string, \"Phi[3]\", direct, output, yes;\n");
%!           fputs(fd, " end: elements;\n");
%!         unwind_protect_cleanup
%!           if (fd ~= -1)
%!             fclose(fd);
%!           endif
%!         end_unwind_protect
%!         options.output_file = fname;
%!         options.verbose = false;
%!         if (~options.verbose)
%!           options.logfile = [fname, ".stdout"];
%!         endif
%!         info(idx_node_type, idx_method) = mbdyn_solver_run(fname, options);
%!         ncfile = [options.output_file, ".nc"];
%!         beta = ncread(ncfile, "elem.inertia.1.B");
%!         gamma = ncread(ncfile, "elem.inertia.1.G_cm");
%!         t0 = [param.ox; param.oy; param.oz];
%!         V = ncread(ncfile, "node.struct.2.XP");
%!         W = ncread(ncfile, "node.struct.2.Omega");
%!         W0 = [param.Wx;
%!               param.Wy;
%!               param.Wz];
%!         V0 = [param.Vx;
%!               param.Vy;
%!               param.Vz];
%!         J0 = diag([param.Jxx;
%!                    param.Jyy;
%!                    param.Jzz]);
%!         gamma0 = J0 * W0;
%!         beta0 = (cross(W0, t0) + V0) * param.m;
%!         l = max([param.a, param.b, param.c]);
%!         err_gamma(idx_node_type, idx_method) = max(norm(gamma - gamma0, "cols")) / norm(gamma0);
%!         err_beta(idx_node_type, idx_method) = max(norm(beta - beta0, "cols")) / norm(beta0);
%! %unwind_protect_cleanup
%!         if (fd ~= -1)
%!           unlink(fname);
%!           files = dir([fname, "*"]);
%!           for i=1:numel(files)
%!             unlink(fullfile(files(i).folder, files(i).name));
%!           endfor
%!         endif
%! %end_unwind_protect
%!     endfor
%!   endfor
%!   status = false(numel(idx_node_type), numel(imethods));
%!   for idx_node_type=1:numel(node_types)
%!     for idx_method=1:numel(imethods)
%!       switch (imethods{idx_method})
%!         case {"msstc3,0.6", "msstc4,0.6", "msstc5,0.6"}
%!           tol = 1e-7;
%!         case {"mssth3,0.6", "mssth4,0.6", "mssth5,0.6"}
%!           tol = 1e-7;
%!         case {"Bathe,0.6"}
%!           tol = 1e-7;
%!         case {"DIRK33", "DIRK43", "DIRK54"}
%!           tol = 1e-8;
%!         case "impliciteuler"
%!           tol = 1e-2;
%!         otherwise
%!           tol = 1e-7;
%!       endswitch
%!       status_msg = "failed";
%!       if (err_gamma(idx_node_type, idx_method) < tol && err_beta(idx_node_type, idx_method) < tol)
%!         status(idx_method) = true;
%!         status_msg = "passed";
%!       endif
%!       fprintf(stderr, "%s: %-20s: (%6d/%4d/%4d/%3.2e/%3.3f) %8.3e / %8.3e < %5.3e - %s\n", node_types{idx_node_type}, imethods{idx_method}, info(idx_node_type, idx_method).total_steps, info(idx_node_type, idx_method).total_iter, info(idx_node_type, idx_method).total_jac, info(idx_node_type, idx_method).total_err, info(idx_node_type, idx_method).total_cpu, err_beta(idx_node_type, idx_method), err_gamma(idx_node_type, idx_method), tol, status_msg);
%!     endfor
%!   endfor
%!   assert_simple(all(status));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
