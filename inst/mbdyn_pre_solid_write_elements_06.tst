## mbdyn_pre_solid_write_elements.tst:06
%!test
%! try
%! ## TEST 6
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! SI_unit_joule = SI_unit_newton * SI_unit_meter;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 10e6 / SI_unit_pascal;
%! param.nu = 0.3;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.beta = 1e-3 / SI_unit_second;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.h = 0.625e-3/2 / SI_unit_meter;
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.d2 = param.d1 - 2 * param.h;
%! param.l1 = param.h;
%! param.z0 = 5e-3 / SI_unit_meter;
%! param.vz0 = 0 / (SI_unit_meter / SI_unit_second);
%! param.g = 9.81 / (SI_unit_meter / SI_unit_second^2);
%! param.Zeta = -80 * pi / 180;
%! param.t1 = 0.1 / SI_unit_second;
%! param.N = 100;
%! param.mus = 0.5;
%! param.muc = 0.5;
%! param.vs = 1 / (SI_unit_meter / SI_unit_second);
%! param.kv = 0 / (SI_unit_pascal / (SI_unit_meter / SI_unit_second));
%! param.i = 1;
%! param.sigma0x = 1e4 / SI_unit_meter^-1;
%! param.sigma0y = 1e4 / SI_unit_meter^-1;
%! param.sigma1x = 0 / (SI_unit_second * SI_unit_meter^-1);
%! param.sigma1y = 0 / (SI_unit_second * SI_unit_meter^-1);
%! param.m = (param.d1^2 - param.d2^2) * pi / 4 * param.l1 * param.rho;
%! param.gamma = 1e-16;
%! param.sn = param.m * param.g / (param.gamma * param.d1);
%! param.z0 -= param.m * param.g * sin(abs(param.Zeta)) / (4 * param.sn);
%! options.verbose = false;
%! options.do_plot = false;
%! options.var_step_size = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
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
%!     fputs(fd, "Point(1) = {0,  0, z0 + 0.5 * d2};\n");
%!     fputs(fd, "Point(2) = {0, l1, z0 + 0.5 * d2};\n");
%!     fputs(fd, "Point(3) = {0, l1, z0 + 0.5 * d1};\n");
%!     fputs(fd, "Point(4) = {0,  0, z0 + 0.5 * d1};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Curve Loop(1) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(l1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round((d1 - d2) / 2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(l1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round((d1 - d2) / 2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Surface(1) = {1,2,3,4};\n");
%!     fputs(fd, "v1[] = Extrude {{0, 1, 0}, {0, 0, z0}, Pi}{ Surface{1}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{0, 1, 0}, {0, 0, z0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Recombine Surface{v1[0],v2[0]};\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {4,9};\n");
%!     fputs(fd, "Mesh.ElementOrder = 1;\n");
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
%!   opt_msh.elem_type = {"iso8", "iso4", "tet4", "penta6", "tria3"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_id_ground = rows(mesh.nodes) + 1;
%!   mesh.nodes(node_id_ground, 1:3) = zeros(1, 3);
%!   mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!   mesh.materials.iso8([mesh.groups.iso8(find([[mesh.groups.iso8.id] == 1])).elements]) = 1;
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean viscoelastic";
%!   mesh.material_data(1).beta = param.beta;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_id_ground) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_sphere";
%!   opt_mbd_mesh.joints.number = 1;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   node_id_cont = mesh.groups.iso4(1).nodes;
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fputs(fd, "set: integer ref_id_sphere = 1000;\n");
%!     fputs(fd, "set: integer drive_id_E = 2000;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fx = 2001;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fy = 2002;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fz = 2003;\n");
%!     fprintf(fd, " set: integer node_id_ground = %d;\n", node_id_ground);
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     if (options.var_step_size)
%!       fputs(fd, "    min time step: t1 / N / 100;\n");
%!       fputs(fd, "    max time step: t1 / N;\n");
%!       fputs(fd, " strategy: factor, 0.5, 1, 1.25, 3, 3, 5;\n");
%!     endif
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-6, test, norm, 1e-6, test, norm;\n");
%!     fputs(fd, "    linear solver: pardiso, grad, max iterations, 100;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: mcp newton min fb;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    #output meter: closest next, 0., forever, t1 / 100;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "    loadable elements: %d;\n", numel(node_id_cont));
%!     fputs(fd, "      gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, "reference: ref_id_sphere,\n");
%!     fputs(fd, "  position, reference, global, null,\n");
%!     fputs(fd, "  orientation, reference, global, eye,\n");
%!     fputs(fd, "  velocity, reference, global, 0., 0., vz0,\n");
%!     fputs(fd, "  angular velocity, reference, global, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     for i=1:numel(node_id_cont)
%!       fprintf(fd, "user defined: %d, unilateral disp in plane,\n", i);
%!       fprintf(fd, "node1, %d,\n", node_id_cont(i));
%!       fputs(fd, "    offset, 1,\n");
%!       fputs(fd, "    reference, node, null,\n");
%!       fputs(fd, "    stiffness, sn,\n");
%!       fputs(fd, "    enable mcp, yes,\n");
%!       fputs(fd, " node2, node_id_ground,\n");
%!       fputs(fd, " offset, reference, node, null,\n");
%!       fputs(fd, " hinge, 3, 0, 0, 1,\n");
%!       fputs(fd, "        1, 1, 0, 0,\n");
%!       fputs(fd, " coulomb friction coefficient, muc,\n");
%!       fputs(fd, " static friction coefficient, mus,\n");
%!       fputs(fd, " sliding velocity coefficient, vs,\n");
%!       fputs(fd, " sliding velocity exponent, i,\n");
%!       fputs(fd, " micro slip stiffness x, sigma0x,\n");
%!       fputs(fd, " micro slip stiffness y, sigma0y,\n");
%!       fputs(fd, " micro slip damping x, sigma1x,\n");
%!       fputs(fd, " micro slip damping y, sigma1y,\n");
%!       fputs(fd, " viscous friction coefficient, kv;\n");
%!     endfor
%!     fputs(fd, "  joint: 1, clamp, node_id_ground, node, node;\n");
%!     fputs(fd, "  gravity: uniform, component, g * cos(Zeta), 0., g * sin(Zeta);\n");
%!     fprintf(fd, "drive caller: drive_id_E, array, %d", rows(mesh.elements.iso8));
%!     for i=1:rows(mesh.elements.iso8)
%!       fprintf(fd, ",\n    element, %d, solid, string, \"E\", direct", i);
%!     endfor
%!     fputs(fd, ",\n  output, yes;\n\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fx, element, 1, joint, string, \"Fx\", direct, output, yes;\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fy, element, 1, joint, string, \"Fy\", direct, output, yes;\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fz, element, 1, joint, string, \"Fz\", direct, output, yes;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [drive_id, drive_data] = mbdyn_post_load_output_drv(opt_mbdyn.output_file);
%!   Wkin_res = drive_data{1};
%!   Fx_res = drive_data{2};
%!   Fy_res = drive_data{3};
%!   Fz_res = drive_data{4};
%!   ra = 0.5 * param.d1;
%!   ri = 0.5 * param.d2;
%!   h = param.l1;
%!   m = param.rho * pi * (ra^2 - ri^2) * h;
%!   J = m * (ra^2 + ri^2) / 2;
%!   mred = m + J / ra^2;
%!   Phi = param.Zeta + pi / 2;
%!   qddot = m * param.g * sin(Phi) / mred;
%!   qdot = qddot * sol.t;
%!   q = 0.5 * qddot * sol.t.^2;
%!   Wkin_ref = 0.5 * mred * qdot.^2;
%!   Fz_ref = repmat(-m * param.g * cos(Phi), size(sol.t));
%!   Fx_ref = repmat(m * param.g * sin(Phi) * (1 - m / mred), size(sol.t));
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Wkin_ref * SI_unit_joule, "-x;Wkin_r_e_f;k");
%!     plot(sol.t * SI_unit_second, Wkin_res * SI_unit_joule, "-o;Wkin_r_e_s;r");
%!     xlabel("t [s]");
%!     ylabel("W [J]");
%!     grid on;
%!     grid minor on;
%!     title("kinetic energy versus time");
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Fz_ref * SI_unit_newton, "-x;Fz_r_e_f;k");
%!     plot(sol.t * SI_unit_second, Fz_res * SI_unit_newton, "-o;Fz_r_e_s;r");
%!     xlabel("t [s]");
%!     ylabel("F [N]");
%!     grid on;
%!     grid minor on;
%!     title("normal force");
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Fx_ref * SI_unit_newton, "-x;Fx_r_e_f;k");
%!     plot(sol.t * SI_unit_second, Fx_res * SI_unit_newton, "-o;Fx_r_e_s;r");
%!     xlabel("t [s]");
%!     ylabel("F [N]");
%!     grid on;
%!     grid minor on;
%!     title("tangent force");
%!   endif
%!   tol = 2e-2;
%!   assert_simple(Wkin_res, Wkin_ref, tol * max(abs(Wkin_res)));
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
