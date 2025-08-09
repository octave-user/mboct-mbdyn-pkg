## mbdyn_pre_write_fem_data.tst:25
%!test
%! try
%! ## TEST25
%! pkg load mboct-fem-pkg;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_kelvin = 1;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! filename = "";
%! %unwind_protect
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
%! # opt_mbd.logfile = [opt_mbd.output_file, ".stdout"];
%!   endif
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
%! %unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! %end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
