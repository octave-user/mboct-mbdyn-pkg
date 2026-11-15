%!test
%! try
%!   pkg load mbdyn_util_oct;
%!   pkg load mboct-fem-pkg;
%!   filename = "";
%!                                 %unwind_protect
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
%!   a = 30e-3 / SI_unit_meter;
%!   b = 20e-3 / SI_unit_meter;
%!   c = 10e-3 / SI_unit_meter;
%!   d = 0e-3 / SI_unit_meter;
%!   cgx = 5e-3 / SI_unit_meter;
%!   cgy = 2e-3 / SI_unit_meter;
%!   cgz = 3e-3 / SI_unit_meter;
%!   h = c;
%!   rho = 2700 / (SI_unit_kilogram / SI_unit_meter^3);
%!   m = rho * a * b * c;
%!   Jcg = diag(m * [b^2 + c^2, a^2 + c^2, a^2 + b^2] / 12);
%!   X0 = [3e-3; 7e-3; 5e-3] / SI_unit_meter;
%!   V0 = [1; 2; 3];
%!   W0 = [pi/2; pi/4; pi];
%!   cg = [cgx; cgy; cgz];
%!   alpha = [pi / 4; pi / 2; pi / 3];
%!   R0 = euler123_to_rotation_matrix(alpha);
%!   V0cg = cross(W0, R0 * cg) + V0;
%!   S = m * cg;
%!   J = Jcg - skew(cg) * skew(cg) * m;
%!   beta0 = m * V0 + cross(W0, R0 * S);
%!   gamma0 = R0 * Jcg * R0.' * W0;
%!   tol = 1e-3;
%!   tolm = 1e-12;
%!   options.interactive = false;
%!   options.plot = false;
%!   options.verbose = false;
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
%!     fprintf(fd, "cgx = %.16e;\n", cgx);
%!     fprintf(fd, "cgy = %.16e;\n", cgy);
%!     fprintf(fd, "cgz = %.16e;\n", cgz);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {cgx - 0.5 * a, cgy - 0.5 * b, cgz - 0.5 * c};\n");
%!     fputs(fd, "Point(2) = {cgx - 0.5 * a, cgy + 0.5 * b, cgz - 0.5 * c};\n");
%!     fputs(fd, "Point(3) = {cgx + 0.5 * a, cgy + 0.5 * b, cgz - 0.5 * c};\n");
%!     fputs(fd, "Point(4) = {cgx + 0.5 * a, cgy - 0.5 * b, cgz - 0.5 * c};\n");
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
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
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
%!                        filename_geo, ...
%!                        "-o", [filename, ".msh"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!   mesh.material_data.E = 10000 * 70000e6 / SI_unit_pascal;
%!   mesh.material_data.nu = 0.3;
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.alpha = 0 / (1 / SI_unit_second);
%!   mesh.material_data.beta = 0 / (SI_unit_second);
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   grp_idx_beam = find([[mesh.groups.iso20].id] == 1);
%!   grp_idx_clamp = find([[mesh.groups.quad8].id] == 1);
%!   mesh.materials.iso20(mesh.groups.iso20(grp_idx_beam).elements) = 1;
%!   cms_opt.number_of_threads = options.number_of_threads;
%!   cms_opt.algorithm = "shift-invert";
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
%!   cms_opt.modes.number = 2;
%!   cms_opt.tolerance_tau = -1;
%!   cms_opt.element.name = "elem_id_modal";
%!   cms_opt.create_binary = true;
%!   cms_opt.use_binary = true;
%!   cms_opt.update_binary = true;
%!   cms_opt.invariants = true;
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_dof.locked_dof(cms_opt.nodes.modal.number, :) = true;
%!   load_case_dof.locked_dof(mesh.groups.quad8(grp_idx_clamp).nodes, :) = true;
%!   [cms_data.mesh, cms_data.mat_ass, cms_data.dof_map, cms_data.sol_eig, cms_data.cms_opt, cms_data.sol_tau] = fem_cms_create2(mesh, load_case_dof, cms_opt);
%!   assert(cms_data.mat_ass.dm, m, m * tolm);
%!   assert(cms_data.mat_ass.J, J, norm(J) * tolm);
%!   assert(cms_data.mat_ass.S, S, norm(S) * tolm);
%!   fem_cms_export(filename, cms_data.mesh, cms_data.dof_map, cms_data.mat_ass, cms_data.cms_opt);
%!   for j=1:2
%!   unwind_protect
%!     filename_mbdyn = sprintf("%s_%d.mbdyn", filename, j);
%!     [fd, msg] = fopen(filename_mbdyn, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!     endif
%!     fprintf(fd, "set: real a = %.16e;\n", a);
%!     fprintf(fd, "set: real b = %.16e;\n", b);
%!     fprintf(fd, "set: real c = %.16e;\n", c);
%!     fprintf(fd, "set: real d = %.16e;\n", d);
%!     for i=1:3
%!       fprintf(fd, "set: real alpha%d = %.16e;\n", i, alpha(i));
%!     endfor
%!     for i=1:3
%!       fprintf(fd, "set: real X0%d = %.16e;\n", i, X0(i));
%!     endfor
%!     fprintf(fd, "set: real m = %.16e;\n", m);
%!     for i=1:3
%!       fprintf(fd, "set: real cg%s = %.16e;\n", {"x","y","z"}{i}, cg(i));
%!     endfor
%!     for i=1:3
%!       fprintf(fd, "set: real J%s = %.16e;\n", {"xx","yy","zz"}{i}, Jcg(i, i));
%!     endfor
%!     for i=1:3
%!       fprintf(fd, "set: real V0%s = %.16e;\n", {"x", "y", "z"}{i}, V0(i));
%!     endfor
%!     for i=1:3
%!       fprintf(fd, "set: real W0%s = %.16e;\n", {"x", "y", "z"}{i}, W0(i));
%!     endfor
%!     fputs(fd, "set: integer ref_id_ground = 1001;\n");
%!     fputs(fd, "set: integer ref_id_modal = 1002;\n");
%!     fputs(fd, "set: integer ref_id_interface1 = 1003;\n");
%!     fputs(fd, "set: integer node_id_modal = 2001;\n");
%!     fputs(fd, "set: integer node_id_interface1 = 2002;\n");
%!     fputs(fd, "set: integer node_id_cg = 2003;\n");
%!     fputs(fd, "set: integer elem_id_modal = 3001;\n");
%!     fputs(fd, "set: integer elem_id_inertia = 3002;\n");
%!     fputs(fd, "set: integer drive_id_beta = 4000;\n");
%!     fputs(fd, "set: integer drive_id_gamma = 5000;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: 1;\n");
%!     fputs(fd, "        time step: 1e-4;\n");
%!     fputs(fd, "        method: ms, 0.6;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-6, test, norm;\n");
%!     fputs(fd, "        max iterations: 1000;\n");
%!     fputs(fd, "        derivatives max iterations: 50;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 100;\n");
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
%!     fputs(fd, "             linear solver max iterations, 300,\n");
%!     fputs(fd, "             minimum step, 1e-12,\n");
%!     fputs(fd, "             krylov subspace size, 300;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: disable;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "        output meter: closest next, 0., forever, 1e-2;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       default orientation: euler123;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "        structural nodes: 3;\n");
%!     switch (j)
%!     case 1
%!     fputs(fd, "        joints: 1;\n");
%!     case 2
%!     fputs(fd, "        rigid bodies: 1;\n");
%!     endswitch
%!     fputs(fd, "        inertia: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_ground,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_modal,\n");
%!     fputs(fd, "        reference, ref_id_ground, X01, X02, X03,\n");
%!     fputs(fd, "        reference, ref_id_ground, euler123, alpha1, alpha2, alpha3,\n");
%!     fputs(fd, "        reference, ref_id_ground, V0x, V0y, V0z,\n");
%!     fputs(fd, "        reference, ref_id_ground, W0x, W0y, W0z;\n");
%!     fputs(fd, "reference: ref_id_interface1,\n");
%!     fputs(fd, "        reference, ref_id_modal, a + d,  0., 0.,\n");
%!     fputs(fd, "        reference, ref_id_modal, eye,\n");
%!     fputs(fd, "        reference, ref_id_modal, null,\n");
%!     fputs(fd, "        reference, ref_id_modal, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     switch (j)
%!     case 1
%!     fputs(fd, "        structural: node_id_modal, modal,\n");
%!     case 2
%!     fputs(fd, "        structural: node_id_modal, dynamic,\n");
%!     endswitch
%!     fputs(fd, "                reference, ref_id_modal, null,\n");
%!     fputs(fd, "                reference, ref_id_modal, eye,\n");
%!     fputs(fd, "                reference, ref_id_modal, null,\n");
%!     fputs(fd, "                reference, ref_id_modal, null, accelerations, yes;\n");
%!     fputs(fd, "        structural: node_id_cg, dummy, node_id_modal, offset,\n");
%!     fputs(fd, "                reference, ref_id_modal, cgx, cgy, cgz,\n");
%!     fputs(fd, "                reference, ref_id_modal, eye;\n");
%!     switch (j)
%!     case 1
%!     fputs(fd, "        structural: node_id_interface1, static,\n");
%!     fputs(fd, "                reference, ref_id_interface1, null,\n");
%!     fputs(fd, "                reference, ref_id_interface1, eye,\n");
%!     fputs(fd, "                reference, ref_id_interface1, null,\n");
%!     fputs(fd, "                reference, ref_id_interface1, null;\n");
%!     case 2
%!     fputs(fd, "        structural: node_id_interface1, dummy, node_id_modal, offset,\n");
%!     fputs(fd, "                reference, ref_id_interface1, null,\n");
%!     fputs(fd, "                reference, ref_id_interface1, eye;\n");
%!     endswitch
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     switch (j)
%!     case 1
%!     fprintf(fd, "        include: \"%s.elm\";\n", filename);
%!     case 2
%!     fputs(fd, "         body: elem_id_modal, node_id_modal, m, reference, ref_id_modal, cgx, cgy, cgz, diag, Jxx, Jyy, Jzz, orientation, reference, ref_id_modal, eye;\n");
%!     endswitch
%!     fputs(fd, "        inertia: elem_id_inertia,\n");
%!     fputs(fd, "                 position, reference, ref_id_ground, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_ground, eye,\n");
%!     switch (j)
%!     case 1
%!     fputs(fd, "                 joint, all,\n");
%!     case 2
%!     fputs(fd, "                 body, all,\n");
%!     endswitch
%!     fputs(fd, "                 output, always;\n");
%!     fputs(fd, "        drive caller: drive_id_beta + 1, element, elem_id_inertia, inertia, string, \"beta[1]\", direct, output, yes;\n");
%!     fputs(fd, "        drive caller: drive_id_beta + 2, element, elem_id_inertia, inertia, string, \"beta[2]\", direct, output, yes;\n");
%!     fputs(fd, "        drive caller: drive_id_beta + 3, element, elem_id_inertia, inertia, string, \"beta[3]\", direct, output, yes;\n");
%!     fputs(fd, "        drive caller: drive_id_gamma + 1, element, elem_id_inertia, inertia, string, \"gamma[1]\", direct, output, yes;\n");
%!     fputs(fd, "        drive caller: drive_id_gamma + 2, element, elem_id_inertia, inertia, string, \"gamma[2]\", direct, output, yes;\n");
%!     fputs(fd, "        drive caller: drive_id_gamma + 3, element, elem_id_inertia, inertia, string, \"gamma[3]\", direct, output, yes;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   options_mbd.output_file = sprintf("%s_%d_mbd", filename, j);
%!   if (options.verbose)
%!     shell(sprintf("cat %s | nl", filename_mbdyn));
%!   endif
%!   res(j).info = mbdyn_solver_run(filename_mbdyn, options_mbd);
%!   [res(j).t, ...
%!    res(j).trajectory, ...
%!    res(j).deformation, ...
%!    res(j).velocity, ...
%!    res(j).acceleration, ...
%!    res(j).node_id] = mbdyn_post_load_output_struct(options_mbd.output_file);
%!   res(j).log_dat = mbdyn_post_load_log(options_mbd.output_file);
%!   res(j).bodies = mbdyn_post_load_log_body(options_mbd.output_file);
%!   switch (j)
%!   case 1
%!     res(j).sol_dyn = fem_post_cms_sol_import(options_mbd.output_file, cms_data);
%!     options_cms.scale_type = "least square";
%!     res(j).sol_tot = fem_post_cms_expand(res(j).sol_dyn, cms_data, 1:numel(res(j).t), options_cms);
%!   endswitch
%!   [res(j).drive_id, res(j).drive_data] = mbdyn_post_load_output_drv(options_mbd.output_file);
%!   res(j).beta = zeros(3, numel(res(j).t));
%!   for i=1:3
%!     res(j).beta(i, :) = res(j).drive_data{res(j).drive_id == res(j).log_dat.vars.drive_id_beta + i};
%!   endfor
%!   res(j).gamma = zeros(3, numel(res(j).t));
%!   for i=1:3
%!     res(j).gamma(i, :) = res(j).drive_data{res(j).drive_id == res(j).log_dat.vars.drive_id_gamma + i};
%!   endfor
%!   assert(res(j).bodies.dm, m, tolm * m);
%!   assert(res(j).bodies.J, J, tolm * norm(J));
%!   assert(res(j).bodies.Xgc, cg, tolm * norm(cg));
%!   assert(max(norm(res(j).beta - beta0, "cols")) / norm(beta0) < tol);
%!   assert(max(norm(res(j).gamma - gamma0, "cols")) / norm(gamma0) < tol);
%!   endfor
%!   for j=1:3
%!     assert(max(norm(res(1).trajectory{j} - res(2).trajectory{j}, "rows")) < tol * max(norm(res(2).trajectory{j}, "rows")));
%!     assert(max(norm(res(1).velocity{j} - res(2).velocity{j}, "rows")) < tol * max(norm(res(2).velocity{j}, "rows")));
%!     switch (j)
%!     case 1
%!     assert(max(norm(res(1).acceleration{j}(10:end,:) - res(2).acceleration{j}(10:end,:), "rows")) < tol * max(norm(res(2).acceleration{j}, "rows")));
%!     endswitch
%!   endfor
%!   for j=1:2
%!     assert(max(norm(res(j).deformation{3}(:, 1:3) - res(j).t * V0cg.', "rows")) < tol * norm(cg));
%!     assert(max(norm(res(j).velocity{3}(:, 1:3) - V0cg.', "rows")) < tol * norm(V0));
%!   endfor
%!                                 %unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%!                                 %end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
