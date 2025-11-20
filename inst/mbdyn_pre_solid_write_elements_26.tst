## fem_pre_mesh_import.m:449
%!test
%! try
%!   pkg load mboct-fem-pkg
%!   ## MBDyn rigid body dynamics
%!   n = 5;
%!   fref = [1.27458549458623;	4.73428636065456];
%!   param.h = 100e-3;
%!   param.d = 10e-3;
%!   param.D = 280e-3;
%!   param.t = 30e-3;
%!   param.dx = 10e-3;
%!   param.omega = 2 * pi * n;
%!   param.g = 9.81;
%!   param.num_modes = 10;
%!   param.T = 1;
%!   options.verbose = false;
%!   filename = "";
%!   unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     fd = -1;
%!     mbdyn_file = [filename, ".mbdyn"];
%!     elem_file = [filename, ".elm"];
%!     nodes_file = [filename, ".nod"];
%!     csl_file = [filename, ".csl"];
%!     geo_file = [filename, ".geo"];
%!     unwind_protect
%!       [fd, msg] = fopen(geo_file, "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s\"", geo_file);
%!       endif
%!       fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!       fprintf(fd, "h=%g;\n", param.h);
%!       fprintf(fd, "d=%g;\n", param.d);
%!       fprintf(fd, "D=%g;\n", param.D);
%!       fprintf(fd, "t=%g;\n", param.t);
%!       fprintf(fd, "dx=%g;\n", param.dx);
%!       fputs(fd, "Point(1) = {0, 0, 0};\n");
%!       fputs(fd, "Point(2) = {0.5 * d, 0, 0};\n");
%!       fputs(fd, "Point(3) = {0.5 * d, 0, 0.5 * (h - t)};\n");
%!       fputs(fd, "Point(4) = {0.5 * D, 0, 0.5 * (h - t)};\n");
%!       fputs(fd, "Point(5) = {0.5 * D, 0, 0.5 * (h + t)};\n");
%!       fputs(fd, "Point(6) = {0.5 * d, 0, 0.5 * (h + t)};\n");
%!       fputs(fd, "Point(7) = {0.5 * d, 0, h};\n");
%!       fputs(fd, "Point(8) = {0, 0, h};\n");
%!       fputs(fd, "Point(9) = {0, 0, (h+t)/2};\n");
%!       fputs(fd, "Point(10) = {0, 0, (h-t)/2};\n");
%!       fputs(fd, "Line(1) = {1, 2};\n");
%!       fputs(fd, "Line(2) = {2, 3};\n");
%!       fputs(fd, "Line(3) = {3, 10};\n");
%!       fputs(fd, "Line(4) = {10, 1};\n");
%!       fputs(fd, "Line(5) = {3, 6};\n");
%!       fputs(fd, "Line(6) = {6, 9};\n");
%!       fputs(fd, "Line(7) = {9, 10};\n");
%!       fputs(fd, "Line(8) = {3, 4};\n");
%!       fputs(fd, "Line(9) = {4,5};\n");
%!       fputs(fd, "Line(10) = {5,6};\n");
%!       fputs(fd, "Line(11) = {6,7};\n");
%!       fputs(fd, "Line(12) = {7,8};\n");
%!       fputs(fd, "Line(13) = {8,9};\n");
%!       fputs(fd, "Transfinite Curve(1) = Round(d/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(2) = Round((h-t)/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(3) = Round(d/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(4) = Round((h-t)/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(5) = Round(t/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(6) = Round(d/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(7) = Round(t/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(8) = Round((D-d)/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(9) = Round(t/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(10) = Round((D-d)/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(11) = Round((h-t)/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(12) = Round(d/2/dx)+1;\n");
%!       fputs(fd, "Transfinite Curve(13) = Round((h-t)/2/dx)+1;\n");
%!       fputs(fd, "Line Loop(14) = {1,2,3,4};\n");
%!       fputs(fd, "Line Loop(15) = {-3,5,6,7};\n");
%!       fputs(fd, "Line Loop(16) = {8,9,10,-5};\n");
%!       fputs(fd, "Line Loop(17) = {11,12,13,6};\n");
%!       fputs(fd, "Plane Surface(18) = {14};\n");
%!       fputs(fd, "Plane Surface(19) = {15};\n");
%!       fputs(fd, "Plane Surface(20) = {16};\n");
%!       fputs(fd, "Plane Surface(21) = {17};\n");
%!       fputs(fd, "Transfinite Surface(18)={1,2,3,10};\n");
%!       fputs(fd, "Transfinite Surface(19)={10,3,6,9};\n");
%!       fputs(fd, "Transfinite Surface(20)={3,4,5,6};\n");
%!       fputs(fd, "Transfinite Surface(21)={6,7,8,9};\n");
%!       fputs(fd, "Recombine Surface{18};\n");
%!       fputs(fd, "Recombine Surface{19};\n");
%!       fputs(fd, "Recombine Surface{20};\n");
%!       fputs(fd, "Recombine Surface{21};\n");
%!       fputs(fd, "v1[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{18};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v2[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{19};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v3[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{20};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v4[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{21};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "Recombine Surface{18,v1[0]};\n");
%!       fputs(fd, "Recombine Surface{19,v2[0]};\n");
%!       fputs(fd, "Recombine Surface{20,v3[0]};\n");
%!       fputs(fd, "Recombine Surface{21,v4[0]};\n");
%!       fputs(fd, "v5[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v1[0]};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v6[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v2[0]};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v7[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v3[0]};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "v8[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v4[0]};Layers{Round(d/2*Pi/dx)};Recombine;};\n");
%!       fputs(fd, "Recombine Surface{v1[0],v5[0]};\n");
%!       fputs(fd, "Recombine Surface{v2[0],v6[0]};\n");
%!       fputs(fd, "Recombine Surface{v3[0],v7[0]};\n");
%!       fputs(fd, "Recombine Surface{v4[0],v8[0]};\n");
%!       fputs(fd, "Physical Volume(\"volume\",1) = {v1[1],v2[1],v3[1],v4[1],v5[1],v6[1],v7[1],v8[1]};\n");
%!       fputs(fd, "Physical Surface(\"support\", 2) = {22, 39};\n");
%!       fputs(fd, "Mesh.ElementOrder = 2;\n");
%!       fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!       fputs(fd, "Mesh.Algorithm = 3;\n");
%!       fputs(fd, "Mesh 3;\n");
%!       fputs(fd, "Coherence Mesh;\n");
%!       fputs(fd, "Mesh.Format = 1;\n");
%!       fprintf(fd, "Save \"%s.msh\";\n", filename);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!     end_unwind_protect
%!     pid = spawn("gmsh", {"-format", "msh2", "-0", "-order", "2", geo_file});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       warning("gmsh failed with status %d", status);
%!     endif
%!     [~] = unlink(geo_file);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     mesh.material_data.E = 210000e6;
%!     mesh.material_data.nu = 0.3;
%!     mesh.material_data.rho = 7850;
%!     mesh.materials.penta18 = ones(rows(mesh.elements.penta18),1, "int32");
%!     mesh.materials.iso27 = ones(rows(mesh.elements.iso27), 1, "int32");
%!     node_idx_tip = rows(mesh.nodes) + 1;
%!     mesh.nodes(node_idx_tip, :) = zeros(1, 6);
%!     mesh.elements.rbe2 = fem_pre_mesh_rbe2_from_surf(mesh, 2, node_idx_tip);
%!     load_case_dof.locked_dof = false(size(mesh.nodes));
%!     load_case = struct();
%!     opt_mbd_mesh = struct();
%!     opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!     opt_mbd_mesh.struct_nodes.type(node_idx_tip) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!     opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!     unwind_protect
%!       [fd, msg] = fopen(mbdyn_file, "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s\"", mbdyn_file);
%!       endif
%!       fprintf(fd, " set: real T = %g.;\n", param.T);
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "    problem: initial value; # the default\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "    initial time: 0;\n");
%!       fputs(fd, "    final time: T;\n");
%!       fputs(fd, "    time step: T / 10.;\n");
%!       fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!       fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!       fputs(fd, "    max iterations: 10000;\n");
%!       fputs(fd, "    tolerance: 1.e-4, test, norm, 1e-6, test, norm;\n");
%!       fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 10;\n");
%!       fputs(fd, "    method: implicit euler;\n");
%!       fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!       fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!       fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 20, linear solver max iterations, 40, minimum step, 1e-12, recovery step, 1e-12, verbose, 3;\n");
%!       fputs(fd, "    eigenanalysis: T,\n");
%!       fputs(fd, "          mode, largest imaginary part,\n");
%!       fputs(fd, "            suffix format, \"%02d\",\n");
%!       fputs(fd, "            output eigenvectors, output geometry,\n");
%!       fputs(fd, "            lower frequency limit, 1.,\n");
%!       fputs(fd, "            use arpack, 10, 200, 0;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "      use automatic differentiation;\n");
%!       fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!       fprintf(fd, "    hydraulic nodes: %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!       fprintf(fd, "    surface loads: %d;\n", opt_mbd_mesh.surface_loads.number);
%!       fprintf(fd, "    joints: %d + 1;\n", opt_mbd_mesh.joints.number);
%!       fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!       fprintf(fd, "    genels: %d;\n", opt_mbd_mesh.genels.number);
%!       fprintf(fd, "    gravity;\n");
%!       fprintf(fd, "    rigid body kinematics: drive, angular velocity, component, const, 0., const, 0.,\n");
%!       fprintf(fd, "piecewise linear, 2, 0., 0.,\n");
%!       fprintf(fd, "                  T, %g;\n", param.omega);
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", nodes_file);
%!       fputs(fd, " end: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", csl_file);
%!       fputs(fd, " begin: elements;\n");
%!       fprintf(fd, "include: \"%s\";\n", elem_file);
%!       fprintf(fd, "joint: 1000, total pin joint,\n");
%!       fprintf(fd, "%d,\n", node_idx_tip);
%!       fputs(fd, "position, reference, global, null,\n");
%!       fputs(fd, "position orientation, reference, global, eye,\n");
%!       fputs(fd, "rotation orientation, reference, global, eye,\n");
%!       fputs(fd, "position, reference, global, null,\n");
%!       fputs(fd, "position orientation, reference, global, eye,\n");
%!       fputs(fd, "rotation orientation, reference, global, eye,     \n");
%!       fputs(fd, "position constraint, \n");
%!       fputs(fd, "		active, \n");
%!       fputs(fd, "		active, \n");
%!       fputs(fd, "		active,\n");
%!       fputs(fd, "	component, \n");
%!       fputs(fd, "		const, 0.,\n");
%!       fputs(fd, "           const, 0.,\n");
%!       fputs(fd, "		const, 0.,\n");
%!       fputs(fd, "orientation constraint,\n");
%!       fputs(fd, "		inactive,\n");
%!       fputs(fd, "		inactive,\n");
%!       fputs(fd, "		inactive,\n");
%!       fputs(fd, "	component,\n");
%!       fputs(fd, "           const, 0.,\n");
%!       fputs(fd, "           const, 0.,\n");
%!       fputs(fd, "		const, 0.;\n");
%!       fprintf(fd, "gravity: uniform, 0., 0., -1.,\n");
%!       fprintf(fd, "piecewise linear, 2, 0., 0.,\n");
%!       fprintf(fd, "                  T, %g;\n", param.g);
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!       fd = -1;
%!     end_unwind_protect
%!     if (~options.verbose)
%!       ## opt_mbdyn.logfile = [filename, ".stdout"];
%!     endif
%!     opt_mbdyn.output_file = [filename, "_mbdyn"];
%!     info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!     [mesh_sol, sol_stat] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!     [mesh_eig, sol_eig] = mbdyn_post_load_output_eig_sol(opt_mbdyn.output_file);
%!     tol = 2e-2;
%!     assert_simple(sol_eig.f(1:2), fref, tol * max(fref));
%!   unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!       endfor
%!     endif
%!   end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
