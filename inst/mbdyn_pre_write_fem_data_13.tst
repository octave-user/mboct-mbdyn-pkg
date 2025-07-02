## mbdyn_pre_write_fem_data.tst:13
%!test
%! try
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
%! # options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!   endif
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
