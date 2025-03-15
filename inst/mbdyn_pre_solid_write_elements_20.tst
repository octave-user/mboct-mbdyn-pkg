## mbdyn_pre_solid_write_elements.tst:06
%!demo
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
%! param.m1 = 20 / SI_unit_kilogram;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 10e6 / SI_unit_pascal;
%! param.nu = 0.3;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.beta = 1e-3 / SI_unit_second;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.h = 5e-3 / SI_unit_meter;
%! param.ri = 15e-3 / SI_unit_meter;
%! param.ra = 20e-3 / SI_unit_meter;
%! param.R = 30e-3 / SI_unit_meter;
%! param.Phi = 45 * pi / 180;
%! param.z0 = param.R + param.ra + 20e-3 / SI_unit_meter;
%! param.vz0 = 0 / (SI_unit_meter / SI_unit_second);
%! param.g = 9.81 / (SI_unit_meter / SI_unit_second^2);
%! param.Zeta = -90 * pi / 180;
%! param.t1 = 0.5 / SI_unit_second;
%! param.tp = 0.1 * param.t1;
%! param.pref = 5e5 / SI_unit_pascal;
%! param.N = 500;
%! param.mus = 0.5;
%! param.muc = 0.5;
%! param.vs = 1 / (SI_unit_meter / SI_unit_second);
%! param.kv = 0 / (SI_unit_pascal / (SI_unit_meter / SI_unit_second));
%! param.i = 1;
%! param.sigma0x = 1e4 / SI_unit_meter^-1;
%! param.sigma0y = 1e4 / SI_unit_meter^-1;
%! param.sigma1x = 0 / (SI_unit_second * SI_unit_meter^-1);
%! param.sigma1y = 0 / (SI_unit_second * SI_unit_meter^-1);
%! options.verbose = true;
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
%!     fputs(fd, "Point(1) = {0,  0, R};\n");
%!     fputs(fd, "Point(2) = {0, ri * Sin(Phi/2), R - ri * Cos(Phi/2)};\n");
%!     fputs(fd, "Point(3) = {0, 0, R + ri};\n");
%!     fputs(fd, "Point(4) = {0, -ri * Sin(Phi/2), R - ri * Cos(Phi/2)};\n");
%!     fputs(fd, "Point(5) = {0, ra * Sin(Phi/2), R - ra * Cos(Phi/2)};\n");
%!     fputs(fd, "Point(6) = {0, 0, R + ra};\n");
%!     fputs(fd, "Point(7) = {0, -ra * Sin(Phi/2), R - ra * Cos(Phi/2)};\n");
%!     fputs(fd, "Circle(1) = {2, 1, 3};\n");
%!     fputs(fd, "Circle(2) = {3, 1, 4};\n");
%!     fputs(fd, "Line(3) = {4, 7};\n");
%!     fputs(fd, "Circle(4) = {7,1,6};\n");
%!     fputs(fd, "Circle(5) = {6,1,5};\n");
%!     fputs(fd, "Line(6) = {5,2};\n");
%!     fputs(fd, "Curve Loop(1) = {1,2,3,4,5,6};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(0.5*(ri+ra) * (Pi - Phi/2) / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round(0.5*(ri+ra) * (Pi - Phi/2) / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round(0.5*(ri+ra) * (Pi - Phi/2) / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(5) = Max(1, Round(0.5*(ri+ra) * (Pi - Phi/2) / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round((ra - ri) / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(6) = Max(1, Round((ra - ri) / h)) + 1;\n");
%!     fputs(fd, "v1[] = Extrude {{0, 1, 0}, {0, 0, 0}, Pi}{ Surface{1}; Layers{Ceil(Pi * R / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{0, 1, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * R / h)}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Recombine Surface{v1[0],v2[0]};\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(\"pressure\", 1) = {9, 10, 3, 2};\n");
%!     fputs(fd, "Physical Surface(\"run\", 2) = {12, 5, 13, 6};\n");
%!     fputs(fd, "Physical Surface(\"hub\", 3) = {7, 14, 11, 4};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "ReorientMesh Volume{v1[1],v2[1]};\n");
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
%!   opt_msh.elem_type = {"iso20r", "quad8r"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_id_ground = rows(mesh.nodes) + 1;
%!   mesh.nodes(node_id_ground, 1:3) = [0, 0, -param.z0];
%!   mesh.materials.iso20r = zeros(rows(mesh.elements.iso20r), 1, "int32");
%!   mesh.materials.iso20r([mesh.groups.iso20r(find([[mesh.groups.iso20r.id] == 1])).elements]) = 1;
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean viscoelastic";
%!   mesh.material_data(1).beta = param.beta;
%!   opt_mbd_mesh.surface_loads.time_function = {"piecewise linear, 2, 0., 0., tp, pref"};
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_id_ground) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_tire";
%!   opt_mbd_mesh.joints.number = 1;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   grp_idx_run = find([mesh.groups.quad8r.id] == 2);
%!   node_id_cont = mesh.groups.quad8r(grp_idx_run).nodes;
%!   grp_idx_pressure = find([mesh.groups.quad8r.id] == 1);
%!   load_case.pressure.quad8r.elements = mesh.elements.quad8r(mesh.groups.quad8r(grp_idx_pressure).elements, :);
%!   load_case.pressure.quad8r.p = ones(size(load_case.pressure.quad8r.elements));
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!   grp_idx_hub = find([mesh.groups.quad8r.id] == 3);
%!   node_idx_hub = mesh.groups.quad8r(grp_idx_hub).nodes;
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fputs(fd, "set: integer ref_id_tire = 1000;\n");
%!     fputs(fd, "set: integer ref_id_hub = 1100;\n");
%!     fputs(fd, "set: integer drive_id_E = 2000;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fx = 2001;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fy = 2002;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fz = 2003;\n");
%!     fputs(fd, "set: integer node_id_hub = 1000000;\n");
%!     fputs(fd, "set: integer node_id_mass = 1100000;\n");
%!     fputs(fd, "set: integer joint_id_hub = 2000000;\n");
%!     fputs(fd, "set: integer joint_id_axle = 3000000;\n");
%!     fputs(fd, "set: integer joint_id_mass = 3100000;\n");
%!     fputs(fd, "set: integer body_id_mass = 4000000;\n");
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
%!     fprintf(fd, "    structural nodes: %d + 2;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number + numel(node_idx_hub) + 2);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "    surface loads: %d;\n", opt_mbd_mesh.surface_loads.number);
%!     fprintf(fd, "    loadable elements: %d;\n", numel(node_id_cont));
%!     fprintf(fd, "    rigid bodies: %d;\n", 1);
%!     fputs(fd, "      gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, "reference: ref_id_tire,\n");
%!     fputs(fd, "  position, reference, global, 0., 0., z0,\n");
%!     fputs(fd, "  orientation, reference, global, eye,\n");
%!     fputs(fd, "  velocity, reference, global, 0., 0., vz0,\n");
%!     fputs(fd, "  angular velocity, reference, global, null;\n");
%!     fputs(fd, "  reference: ref_id_hub,\n");
%!     fputs(fd, "  position, reference, ref_id_tire, null,\n");
%!     fputs(fd, "  orientation, reference, ref_id_tire, eye,\n");
%!     fputs(fd, "  velocity, reference, ref_id_tire, null,\n");
%!     fputs(fd, "  angular velocity, reference, ref_id_tire, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "  structural: node_id_hub, static,\n");
%!     fputs(fd, "  position, reference, ref_id_hub, null,\n");
%!     fputs(fd, "  orientation, reference, ref_id_hub, eye,\n");
%!     fputs(fd, "  velocity, reference, ref_id_hub, null,\n");
%!     fputs(fd, "  angular velocity, reference, ref_id_hub, null;\n");

%!     fputs(fd, "  structural: node_id_mass, dynamic,\n");
%!     fputs(fd, "  position, reference, ref_id_hub, null,\n");
%!     fputs(fd, "  orientation, reference, ref_id_hub, eye,\n");
%!     fputs(fd, "  velocity, reference, ref_id_hub, null,\n");
%!     fputs(fd, "  angular velocity, reference, ref_id_hub, null;\n");

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
%!       #fputs(fd, "    stiffness, sn,\n");
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
%!     fputs(fd, "  gravity: uniform, cos(Zeta), 0., sin(Zeta), string, \"(Time >= tp) * g\";\n");
%!     fprintf(fd, "drive caller: drive_id_E, array, %d", rows(mesh.elements.iso20r));
%!     for i=1:rows(mesh.elements.iso20r)
%!       fprintf(fd, ",\n    element, %d, solid, string, \"E\", direct", i);
%!     endfor
%!     fputs(fd, ",\n  output, yes;\n\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fx, element, 1, joint, string, \"Fx\", direct, output, yes;\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fy, element, 1, joint, string, \"Fy\", direct, output, yes;\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fz, element, 1, joint, string, \"Fz\", direct, output, yes;\n");
%!     for i=1:numel(node_idx_hub)
%!       fprintf(fd, "joint: joint_id_hub + %d - 1, offset displacement joint, node_id_hub, %d;\n", i, node_idx_hub(i));
%!     endfor
%!     fputs(fd, "body: body_id_mass, node_id_mass, m1, null, diag, 0., 0., 0.;\n");
%!     fputs(fd, "joint: joint_id_axle, total joint,\n");
%!     fputs(fd, "node_id_hub,\n");
%!     fputs(fd, "position, reference, ref_id_hub, null,\n");
%!     fputs(fd, "position orientation, reference, ref_id_hub, eye,\n");
%!     fputs(fd, "rotation orientation, reference, ref_id_hub, eye,\n");
%!     fputs(fd, "node_id_mass,\n");
%!     fputs(fd, "position, reference, ref_id_hub, null,\n");
%!     fputs(fd, "position orientation, reference, ref_id_hub, eye,\n");
%!     fputs(fd, "rotation orientation, reference, ref_id_hub, eye,\n");
%!     fputs(fd, "position constraint,\n");
%!     fputs(fd, "active, active, active, null,\n");
%!     fputs(fd, "orientation constraint,\n");
%!     fputs(fd, "active, inactive, active, null;\n");
%!
%!     fputs(fd, "joint: joint_id_mass, total pin joint,\n");
%!     fputs(fd, "node_id_mass,\n");
%!     fputs(fd, "position, reference, ref_id_hub, null,\n");
%!     fputs(fd, "position orientation, reference, ref_id_hub, eye,\n");
%!     fputs(fd, "rotation orientation, reference, ref_id_hub, eye,\n");
%!     fputs(fd, "position, reference, ref_id_hub, null,\n");
%!     fputs(fd, "position orientation, reference, ref_id_hub, eye,\n");
%!     fputs(fd, "rotation orientation, reference, ref_id_hub, eye,\n");
%!     fputs(fd, "position constraint,\n");
%!     fputs(fd, "inactive, active, inactive, null,\n");
%!     fputs(fd, "orientation constraint,\n");
%!     fputs(fd, "active, active, active, null;\n");
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
%!   if (options.verbose)
%!     shell(["nl ", mbdyn_file]);
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [drive_id, drive_data] = mbdyn_post_load_output_drv(opt_mbdyn.output_file);
%!   Wkin_res = drive_data{1};
%!   Fx_res = drive_data{2};
%!   Fy_res = drive_data{3};
%!   Fz_res = drive_data{4};
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
