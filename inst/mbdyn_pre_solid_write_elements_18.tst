## mbdyn_pre_solid_write_elements.tst:10
%!test
%! try
%! ## TEST 18
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.114)
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
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! M_alpha = 0.25;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.kappa = param.E / (3 * (1 - 2 * param.nu));
%! param.C1 = param.G / (2 * (1 + M_alpha));
%! param.C2 = M_alpha * param.C1;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 2.5e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! %unwind_protect
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
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
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
%!   opt_msh.elem_type = {"quad8", "iso20fr", "penta15f", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"quad8", "tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"quad8", "tria6h"});
%!   if (isfield(mesh.elements, "iso20fr"))
%!     mesh.materials.iso20fr = zeros(rows(mesh.elements.iso20fr), 1, "int32");
%!     mesh.materials.iso20fr([mesh.groups.iso20fr(find([[mesh.groups.iso20fr.id] == 1])).elements]) = 1;
%!   endif
%!   if (isfield(mesh.elements, "penta15f"))
%!     mesh.materials.penta15f = zeros(rows(mesh.elements.penta15f), 1, "int32");
%!     mesh.materials.penta15f([mesh.groups.penta15f(find([[mesh.groups.penta15f.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).C1 = param.C1;
%!   mesh.material_data(1).C2 = param.C2;
%!   mesh.material_data(1).kappa = param.kappa;
%!   mesh.material_data(1).type = "mooney rivlin elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%! # opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * (param.C1 + param.C2) * D; ## (4.114)
%!   Nref = -1/2 * pi * Ra^4 * (param.C1 + 2 * param.C2) * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-3;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! %unwind_protect_cleanup
%!   if (numel(filename))
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
