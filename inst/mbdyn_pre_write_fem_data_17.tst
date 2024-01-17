## mbdyn_pre_write_fem_data.tst:17
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
