%! mbdyn_pre_solid_write_elements.tst:01
%!test
%! try
%!   ## TEST 1
%!   close all;
%!   pkg load mboct-fem-pkg;
%!   ## Define the unit system
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1;
%!   SI_unit_kilogram = 1e3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   h0 = [0.4e-3, 1.4e-3, 1.8e-3] / SI_unit_meter;
%!   options.verbose = true;
%!   options.plot = false;
%!   for j=1:numel(h0)
%!   param.E = 210000e6 / SI_unit_pascal;
%!   param.nu = 0.3;
%!   param.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.t = 1e-3 / SI_unit_meter;
%!   param.t_ = param.t;
%!   param.Di = 20e-3 / SI_unit_meter;
%!   param.De = 30e-3 / SI_unit_meter;
%!   param.h = param.t;
%!   param.h0 = h0(j);
%!   param.s = param.h0;
%!   param.alpha = atan(2 * param.h0 / (param.De - param.Di))
%!   param.l0 = param.h0 + param.t * cos(param.alpha);
%!   #unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     fd = -1;
%!     mbdyn_file = [filename, ".mbdyn"];
%!     elem_file = [filename, ".elm"];
%!     nodes_file = [filename, ".nod"];
%!     csl_file = [filename, ".csl"];
%!     geometry_file = [filename, ".geo"];
%!     unwind_protect
%!       [fd, msg] = fopen(geometry_file, "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s\"", geometry_file);
%!       endif
%!       fprintf(fd, "SetFactory(\"Built-in\");\n");
%!       fprintf(fd, "t = %g;\n", param.t);
%!       fprintf(fd, "h0 = %g;\n", param.h0);
%!       fprintf(fd, "De = %g;\n", param.De);
%!       fprintf(fd, "Di = %g;\n", param.Di);
%!       fprintf(fd, "h = %g;\n", param.h);
%!       fprintf(fd, "alpha = %.16e;\n", param.alpha);
%!       fprintf(fd, "Point(1) = {h0, 0.5 * Di, 0};\n");
%!       fprintf(fd, "Point(2) = {h0 + t * Cos(alpha), 0.5 * Di + t * Sin(alpha), 0};\n");
%!       fprintf(fd, "Point(3) = {h0 + t * Cos(alpha) - (De - Di)/2 * Sin(alpha), 0.5 * De, 0};\n");
%!       fprintf(fd, "Point(4) = {0, 1/2 * De - t * Sin(alpha), 0};\n");
%!       fprintf(fd, "Line(1) = {1,2};\n");
%!       fprintf(fd, "Line(2) = {2,3};\n");
%!       fprintf(fd, "Line(3) = {3,4};\n");
%!       fprintf(fd, "Line(4) = {4,1};\n");
%!       fprintf(fd, "Transfinite Curve(1) = Max(1, Round(t / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(2) = Max(1, Round((De - Di) / (2 * h))) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(3) = Max(1, Round(t / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(4) = Max(1, Round((De - Di) / (2 * h))) + 1;\n");
%!       fprintf(fd, "Line Loop(1) = {1,2,3,4};\n");
%!       fprintf(fd, "Plane Surface(1) = {1};\n");
%!       fprintf(fd, "Transfinite Surface(1) = {PointsOf{Surface{1};}};\n");
%!       fprintf(fd, "v1[] = Extrude {{1,0,0},{0,0,0},Pi/2}{ Surface{1}; Layers{Round(Pi / 2 * (Di + De) / 4 / h)}; Recombine;};\n");
%!       fprintf(fd, "Recombine Surface{1, v1[0]};\n");
%!       fprintf(fd, "Physical Volume(\"volume\", 1) = {v1[1]};\n");
%!       fprintf(fd, "Physical Surface(\"xz\", 2) = {26};\n");
%!       fprintf(fd, "Physical Surface(\"xy\", 3) = {1};\n");
%!       fprintf(fd, "Physical Curve(\"support\", 4) = {20};\n");
%!       fprintf(fd, "Physical Curve(\"load\", 5) = {12};\n");
%!       fprintf(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!       fprintf(fd, "Mesh.ElementOrder = 2;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!       fd = -1;
%!     end_unwind_protect
%!     [~] = unlink([filename, ".msh"]);
%!     ## spawn_wait(spawn("gmsh", {[filename, ".geo"]})); return;
%!     pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       warning("gmsh failed with status %d", status);
%!     endif
%!     [~] = unlink([filename, ".geo"]);
%!     param.material = "hookean linear elastic isotropic";
%!     opt_msh.elem_type = {"iso20r", "quad8r", "line3", "point1"};
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!     [~] = unlink([filename, ".msh"]);
%!     grp_idx_volume = find([mesh.groups.iso20r.id] == 1);
%!     grp_idx_xz = find([mesh.groups.quad8r.id] == 2);
%!     grp_idx_xy = find([mesh.groups.quad8r.id] == 3);
%!     grp_idx_support = find([mesh.groups.line3.id] == 4);
%!     grp_idx_load = find([mesh.groups.line3.id] == 5);
%!     load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!     load_case_dof.locked_dof(mesh.groups.quad8r(grp_idx_xz).nodes, 2) = true;
%!     load_case_dof.locked_dof(mesh.groups.quad8r(grp_idx_xy).nodes, 3) = true;
%!     load_case_dof.locked_dof(mesh.groups.line3(grp_idx_support).nodes, 1) = true;
%!     mesh.materials.iso20r = zeros(rows(mesh.elements.iso20r), 1, "int32");
%!     mesh.materials.iso20r(mesh.groups.iso20r(grp_idx_volume).elements) = 1;
%!     mesh.material_data(1).E = param.E;
%!     mesh.material_data(1).nu = param.nu;
%!     mesh.material_data(1).rho = param.rho;
%!     load_case = struct();
%!     opt_mbd_mesh = struct();
%!     opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!     unwind_protect
%!       [fd, msg] = fopen(mbdyn_file, "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s\"", mbdyn_file);
%!       endif
%!       fprintf(fd, "set: real t1 = %g;\n", 1 / SI_unit_second);
%!       fprintf(fd, "set: real h0 = %g;\n", param.h0);
%!       fprintf(fd, "set: real t = %g;\n", param.t);
%!       fprintf(fd, "set: real alpha = %.16e;\n", param.alpha);
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "    problem: initial value; # the default\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "    initial time: 0;\n");
%!       fputs(fd, "    final time: t1;\n");
%!       fputs(fd, "    time step: t1 / 50;\n");
%!       fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!       fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!       fputs(fd, "    max iterations: 10000;\n");
%!       fputs(fd, "    tolerance: 1.e-6, test, norm, 1e-6, test, norm;\n");
%!       fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!       fputs(fd, "    method: implicit euler;\n");
%!       fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!       fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!       fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 20, linear solver max iterations, 40, minimum step, 1e-12, recovery step, 1e-12, verbose, 3;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "    model: static;\n");
%!       fputs(fd, "    use automatic differentiation;\n");
%!       fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!       fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!       fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!       fprintf(fd, "    genels: %d;\n", opt_mbd_mesh.genels.number + numel(mesh.groups.line3(grp_idx_load).nodes));
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", nodes_file);
%!       fputs(fd, " end: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", csl_file);
%!       fputs(fd, " begin: elements;\n");
%!       fprintf(fd, "include: \"%s\";\n", elem_file);
%!       for i=1:numel(mesh.groups.line3(grp_idx_load).nodes)
%!         genel_id = opt_mbd_mesh.genels.number + i;
%!         fprintf(fd, "genel: %d, clamp, %d, structural, 1, algebraic,\n", ...
%!                 genel_id, mesh.groups.line3(grp_idx_load).nodes(i));
%!         fprintf(fd, "  piecewise linear, 2,\n");
%!         fprintf(fd, "  0., h0 + t * cos(alpha),\n");
%!         fprintf(fd, " t1, t;\n");
%!       endfor
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
%!     [mesh_sol, sol_stat] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!     [genel_id, genel_data] = mbdyn_post_load_output([opt_mbdyn.output_file, ".gen"],1);
%!     s = F = zeros(numel(sol_stat.t), 1);
%!     for i=1:numel(mesh.groups.line3(grp_idx_load).nodes)
%!       genel_idx = find(opt_mbd_mesh.genels.number + i == genel_id);
%!       F += genel_data{genel_idx};
%!       s -= sol_stat.def(mesh.groups.line3(grp_idx_load).nodes(i), 1, :)(:);
%!     endfor
%!     s /= numel(mesh.groups.line3(grp_idx_load).nodes);
%!     ## Roloff/Matek, Maschinenelemente, 2005, page 293, formula (10.24)-(10.25)
%!     delta = param.De / param.Di;
%!     c1 = (param.t_ / param.t)^2 / ((0.25 * param.l0 / param.t - param.t_/param.t + 0.75) * (0.625 * param.l0 / param.t - param.t_ / param.t + 0.375));
%!     c2 = (0.156 * (param.l0 / param.t - 1)^2 + 1) * c1 / (param.t_/param.t)^3;
%!     K1 = 1 / pi * ((delta - 1) / delta)^2 / ((delta + 1) / (delta - 1) - 2 / log(delta));
%!     K4 = sqrt(-0.5 * c1 + sqrt((0.5 * c1)^2 + c2));
%!     F_ref = 4 * param.E / (1 - param.nu^2) * param.t^4 / (K1 * param.De^2) * K4^2 * s / param.t .* (K4^2 * (param.h0 / param.t - s / param.t) .* (param.h0 / param.t - s / (2 * param.t)) + 1);
%!     if (options.plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(s * SI_unit_meter, 4 * F * SI_unit_newton, "-;F;r");
%!     plot(s * SI_unit_meter, F_ref * SI_unit_newton, "-;Fref;k");
%!     xlabel("s [m]");
%!     ylabel("F [N]");
%!     grid minor on;
%!     title(sprintf("reaction force versus displacement h0/t=%.2f", param.h0 / param.t));
%!     endif
%!     tol = 7e-2;
%!     assert_simple(4 * F, F_ref, tol * max(abs(F_ref)));
%!   #unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%!   #end_unwind_protect
%!   endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
