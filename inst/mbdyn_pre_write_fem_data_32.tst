%! fem_pre_mesh_import.m:449
%!test
%! try
%!   pkg load mboct-fem-pkg;
%!   pkg load mbdyn_util_oct;
%!   ## MBDyn rigid body dynamics
%!   param.n = 10;
%!   param.h = 100e-3;
%!   param.d = 50e-3;
%!   param.D = 280e-3;
%!   param.t = 30e-3;
%!   param.dx = 50e-3;
%!   param.g = 9.81;
%!   options.post_proc = false;
%!   filename = "";
%!   ## unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     fd = -1;
%!     unwind_protect
%!       [fd, msg] = fopen([filename, ".geo"], "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s.geo\"", filename);
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
%!       fputs(fd, "Transfinite Curve(1) = Max(2,Round(d/2/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(2) = Max(2,Round((h-t)/2/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(3) = Max(2,Round(d/2/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(4) = Max(2,Round((h-t)/2/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(5) = Max(2,Round(t/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(6) = Max(2,Round(d/2/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(7) = Max(2,Round(t/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(8) = Max(2,Round((D-d)/2/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(9) = Max(2,Round(t/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(10) = Max(2,Round((D-d)/2/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(11) = Max(2,Round((h-t)/2/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(12) = Max(2,Round(d/2/dx)+1);\n");
%!       fputs(fd, "Transfinite Curve(13) = Max(2,Round((h-t)/2/dx)+1);\n");
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
%!       fputs(fd, "v1[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{18};Layers{Max(2,Round(d/2*Pi/dx))};Recombine;};\n");
%!       fputs(fd, "v2[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{19};Layers{Max(2,Round(d/2*Pi/dx))};Recombine;};\n");
%!       fputs(fd, "v3[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{20};Layers{Max(2,Round(d/2*Pi/dx))};Recombine;};\n");
%!       fputs(fd, "v4[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{21};Layers{Max(2,Round(d/2*Pi/dx))};Recombine;};\n");
%!       fputs(fd, "Recombine Surface{18,v1[0]};\n");
%!       fputs(fd, "Recombine Surface{19,v2[0]};\n");
%!       fputs(fd, "Recombine Surface{20,v3[0]};\n");
%!       fputs(fd, "Recombine Surface{21,v4[0]};\n");
%!       fputs(fd, "v5[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v1[0]};Layers{Max(2,Round(d/2*Pi/dx))};Recombine;};\n");
%!       fputs(fd, "v6[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v2[0]};Layers{Max(2,Round(d/2*Pi/dx))};Recombine;};\n");
%!       fputs(fd, "v7[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v3[0]};Layers{Max(2,Round(d/2*Pi/dx))};Recombine;};\n");
%!       fputs(fd, "v8[] = Extrude {{0,0,1},{0,0,0},Pi}{Surface{v4[0]};Layers{Max(2,Round(d/2*Pi/dx))};Recombine;};\n");
%!       fputs(fd, "Recombine Surface{v1[0],v5[0]};\n");
%!       fputs(fd, "Recombine Surface{v2[0],v6[0]};\n");
%!       fputs(fd, "Recombine Surface{v3[0],v7[0]};\n");
%!       fputs(fd, "Recombine Surface{v4[0],v8[0]};\n");
%!       fputs(fd, "Physical Volume(\"volume\",1) = {v1[1],v2[1],v3[1],v4[1],v5[1],v6[1],v7[1],v8[1]};\n");
%!       fputs(fd, "Physical Surface(\"support\", 2) = {22, 39};\n");
%!       fputs(fd, "Physical Surface(\"center\", 3) = {48, 31};\n");
%!       fputs(fd, "Mesh.ElementOrder = 2;\n");
%!       fputs(fd, "Mesh.SecondOrderIncomplete=0;\n");
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
%!     ## spawn_wait(spawn("gmsh", {[filename, ".geo"]})); return;
%!     pid = spawn("gmsh", {"-format", "msh2", "-0", [filename, ".geo"]});
%!     status = spawn_wait(pid);
%!     if (status ~= 0)
%!       warning("gmsh failed with status %d", status);
%!     endif
%!     [~] = unlink([filename, ".geo"]);
%!     mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh"));
%!     mesh.material_data.E = 210000e6;
%!     mesh.material_data.nu = 0.3;
%!     mesh.material_data.rho = 7850;
%!     mesh.materials.penta18 = ones(rows(mesh.elements.penta18),1, "int32");
%!     mesh.materials.iso27 = ones(rows(mesh.elements.iso27), 1, "int32");
%!     node_idx_tip = rows(mesh.nodes) + 1;
%!     node_idx_center = rows(mesh.nodes) + 2;
%!     mesh.nodes(node_idx_tip, :) = zeros(1, 6);
%!     mesh.nodes(node_idx_center, 1:3) = [0, 0, 0.5 * param.h];
%!     mesh.elements.rbe2 = fem_pre_mesh_rbe2_from_surf(mesh, 2, node_idx_tip);
%!     mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, 3, node_idx_center, "quad9");
%!     cms_data.load_case.locked_dof = false(size(mesh.nodes));
%!     cms_data.cms_opt.nodes.modal.number = node_idx_tip;
%!     cms_data.cms_opt.nodes.modal.name = "node_id_modal";
%!     cms_data.cms_opt.nodes.interfaces.number = node_idx_center;
%!     cms_data.cms_opt.nodes.interfaces.name = "node_id_interface1";
%!     cms_data.cms_opt.modes.number = 20;
%!     cms_data.cms_opt.element.name = "elem_id_pegtop";
%!     cms_data.cms_opt.solver = "umfpack";
%!     cms_data.cms_opt.pre_scaling = true;
%!     cms_data.cms_opt.refine_max_iter = int32(10);
%!     ## cms_data.cms_opt.tolerance_tau = 1e-6;
%!     [cms_data.mesh, cms_data.mat_ass, cms_data.dof_map, cms_data.sol_eig, cms_data.cms_opt, cms_data.sol_tau] = fem_cms_create2(mesh, cms_data.load_case, cms_data.cms_opt);
%!     ##cms_data.mat_ass = rmfield(cms_data.mat_ass, "KTAU0red");
%!     modal_file = [filename, "_modal"];
%!     nodes_file = [filename, "_solid.nod"];
%!     csl_file = [filename, "_solid.csl"];
%!     elem_file = [filename, "_solid.elm"];
%!     fem_cms_export(modal_file, cms_data.mesh, cms_data.dof_map, cms_data.mat_ass, cms_data.cms_opt);
%!     opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_pegtop";
%!     opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!     opt_mbd_mesh.struct_nodes.type(node_idx_tip) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!     opt_mbd_mesh.struct_nodes.type(node_idx_center) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!     opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, cms_data.load_case, cms_data.load_case, elem_file, opt_mbd_mesh);
%!     for i=1:numel(param.n)
%!       for j=1:3
%!         fd = -1;
%!         unwind_protect
%!           filename_mbdyn = sprintf("%s_%02d_%d.mbdyn", filename, i, j);
%!           [fd, msg] = fopen(filename_mbdyn, "w");
%!           if (fd == -1)
%!             error("failed to open file \"%s\": %s", filename_mbdyn, msg);
%!           endif
%!           fputs(fd, "set: integer ref_id_ground = 100;\n");
%!           fputs(fd, "set: integer ref_id_pegtop = 101;\n");
%!           switch (j)
%!           case {1, 2}
%!             fputs(fd, "set: integer ref_id_interface1 = 102;\n");
%!             fputs(fd, "set: integer node_id_modal = 1001;\n");
%!             fputs(fd, "set: integer node_id_interface1 = 1002;\n");
%!           endswitch
%!           switch (j)
%!           case {1, 2}
%!             fputs(fd, "set: integer elem_id_pegtop = 2001;\n");
%!           endswitch
%!           fputs(fd, "set: integer joint_id_ground = 2002;\n");
%!           fputs(fd, "set: real k = 2;\n");
%!           fputs(fd, "set: real omega0x = 2. * pi * 0.1;\n");
%!           fputs(fd, "set: real omega0y = 0.;\n");
%!           fprintf(fd, "set: real omega0z = 2. * pi * %g;\n", param.n(i));
%!           fputs(fd, "set: real g = 9.81;\n");
%!           fputs(fd, "set: real t1 = 2 * pi * k / abs(omega0z);\n");
%!           fprintf(fd, "set: real h = %.16g;\n", param.h);
%!           fputs(fd, "begin: data;\n");
%!           fputs(fd, "        problem: initial value;\n");
%!           fputs(fd, "end: data;\n");
%!           fputs(fd, "begin: initial value;\n");
%!           fputs(fd, "        initial time: 0;\n");
%!           fputs(fd, "        final time: t1;\n");
%!           fputs(fd, "        time step: t1 / k / 360.;\n");
%!           fputs(fd, "        method: DIRK33;\n");
%!           fputs(fd, "        tolerance: 1e-4, test, minmax, 1e-6, test,minmax;\n");
%!           fputs(fd, "        max iterations: 1000;\n");
%!           fputs(fd, "        derivatives max iterations: 50;\n");
%!           fputs(fd, "        derivatives coefficient: 1e-8, auto;\n");
%!           fputs(fd, "        derivatives tolerance: 1e-6, 1e-6;\n");
%!           fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!           fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 100;\n");
%!           fputs(fd, "        nonlinear solver: nox, modified, 10,\n");
%!           fputs(fd, "             inner iterations before assembly, 6,\n");
%!           fputs(fd, "             jacobian operator, newton krylov,\n");
%!           fputs(fd, "             solver, line search based,\n");
%!           fputs(fd, "             line search method, backtrack,\n");
%!           fputs(fd, "             recovery step type, constant,\n");
%!           fputs(fd, "             recovery step, 1e-6,\n");
%!           fputs(fd, "             verbose, yes,\n");
%!           fputs(fd, "             forcing term, type 2,\n");
%!           fputs(fd, "             direction, newton,\n");
%!           fputs(fd, "             weighted rms absolute tolerance, 0*1e-3,\n");
%!           fputs(fd, "             weighted rms relative tolerance, 0*1e-3,\n");
%!           fputs(fd, "             linear solver, gmres,\n");
%!           fputs(fd, "             linear solver max iterations, 300,\n");
%!           fputs(fd, "             minimum step, 1e-12,\n");
%!           fputs(fd, "             krylov subspace size, 300;\n");
%!           fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!           fputs(fd, "        threads: assembly, 1;\n");
%!           fputs(fd, "end: initial value;\n");
%!           fputs(fd, "begin: control data;\n");
%!           fputs(fd, "       skip initial joint assembly;\n");
%!           fputs(fd, "       use automatic differentiation;\n");
%!           fputs(fd, "       default output: structural nodes, drive callers;\n");
%!           fputs(fd, "       default orientation: euler123;\n");
%!           fputs(fd, "       output precision: 16;\n");
%!           fputs(fd, "       max iterations: 0;\n");
%!           switch (j)
%!             case 1
%!               fputs(fd, "    structural nodes: 2;\n");
%!               fputs(fd, "    joints: 2;\n");
%!             case 2
%!               fputs(fd, "    structural nodes: 2;\n");
%!               fputs(fd, "    joints: 1;\n");
%!               fputs(fd, "    rigid bodies: 1;\n");
%!             case 3
%!               fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!               fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number + 1);
%!               fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!           endswitch
%!           fputs(fd, "        gravity;\n");
%!           fputs(fd, "end: control data;\n");
%!           fputs(fd, "reference: ref_id_ground,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, eye,\n");
%!           fputs(fd, "        reference, global, null,\n");
%!           fputs(fd, "        reference, global, null;\n");
%!           fputs(fd, "reference: ref_id_pegtop,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, eye,\n");
%!           fputs(fd, "        reference, ref_id_ground, null,\n");
%!           fputs(fd, "        reference, ref_id_ground, omega0x, omega0y, omega0z;\n");
%!           switch (j)
%!           case {1, 2}
%!             fputs(fd, "reference: ref_id_interface1,\n");
%!             fputs(fd, "        reference, ref_id_pegtop, 0., 0., h/2.,\n");
%!             fputs(fd, "        reference, ref_id_pegtop, eye,\n");
%!             fputs(fd, "        reference, ref_id_pegtop, null,\n");
%!             fputs(fd, "        reference, ref_id_pegtop, null;\n");
%!           endswitch
%!           fprintf(fd, "      include: \"%s\";\n", csl_file);
%!           fputs(fd, "begin: nodes;\n")
%!           switch (j)
%!           case {1, 2}
%!             fputs(fd, "        structural: node_id_modal, modal,\n");
%!             fputs(fd, "                reference, ref_id_pegtop, null,\n");
%!             fputs(fd, "                reference, ref_id_pegtop, eye,\n");
%!             fputs(fd, "                reference, ref_id_pegtop, null,\n");
%!             fputs(fd, "                reference, ref_id_pegtop, null, accelerations, yes;\n");
%!           case 3
%!             fprintf(fd, "      include: \"%s\";\n", nodes_file);
%!           endswitch
%!           switch (j)
%!             case 1
%!               fputs(fd, "        structural: node_id_interface1, static,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, eye,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null;\n");
%!             case 2
%!               fputs(fd, "        structural: node_id_interface1, dummy, node_id_modal, offset,\n");
%!               fputs(fd, "                reference, ref_id_interface1, null,\n");
%!               fputs(fd, "                reference, ref_id_interface1, eye;\n");
%!           endswitch
%!           fputs(fd, "end: nodes;\n");
%!           fputs(fd, "begin: elements;\n");
%!           fputs(fd, "        joint: joint_id_ground, total pin joint,\n");
%!           switch (j)
%!           case {1, 2}
%!             fputs(fd, "                node_id_modal,\n");
%!           case 3
%!             fprintf(fd, "              %d,\n", node_idx_tip);
%!           endswitch
%!           fputs(fd, "                        position, reference, ref_id_ground, null,\n");
%!           fputs(fd, "                        position orientation, reference, ref_id_ground, eye,\n");
%!           fputs(fd, "                        rotation orientation, reference, ref_id_ground, eye,\n");
%!           fputs(fd, "                        position, reference, ref_id_ground, null,\n");
%!           fputs(fd, "                        position orientation, reference, ref_id_ground, eye,\n");
%!           fputs(fd, "                        rotation orientation, reference, ref_id_ground, eye,\n");
%!           fputs(fd, "               position constraint,\n");
%!           fputs(fd, "                          active, active, active,\n");
%!           fputs(fd, "                          null,\n");
%!           fputs(fd, "               orientation constraint,\n");
%!           fputs(fd, "                        inactive, inactive, inactive,\n");
%!           fputs(fd, "                          null;\n");
%!           switch (j)
%!             case 1
%!               fprintf(fd, "include: \"%s.elm\";\n", modal_file);
%!             case 2
%!               dm = cms_data.mat_ass.dm;
%!               Xcg = cms_data.mat_ass.S / dm;
%!               Jcg = cms_data.mat_ass.J + skew(Xcg) * skew(Xcg) * dm;
%!               fputs(fd, "body: elem_id_pegtop, node_id_modal,\n");
%!               fprintf(fd, "%.16e,\n", cms_data.mat_ass.dm);
%!               fputs(fd, "reference, node,");
%!               fprintf(fd, "%.16e, ", Xcg);
%!               fputs(fd, "\nmatr,\n")
%!               for k=1:3
%!                 fprintf(fd, "%.16e, ", Jcg(k, :));
%!                 fputs(fd, "\n");
%!               endfor
%!               fputs(fd, "orientation, reference, node, eye;\n");
%!             case 3
%!               fprintf(fd, "include: \"%s\";\n", elem_file);
%!           endswitch
%!           fputs(fd, "        gravity: uniform, 0., 0., -1., const, g;\n");
%!           fputs(fd, "end: elements;\n");
%!         unwind_protect_cleanup
%!           if (fd ~= -1)
%!             fclose(fd);
%!           endif
%!           fd = -1;
%!         end_unwind_protect
%!         opt_mbdyn.output_file = sprintf("%s_%02d_%d", filename, i, j);
%!         mbdyn_solver_run(filename_mbdyn, opt_mbdyn);
%!         res(j).log_dat = mbdyn_post_load_log(opt_mbdyn.output_file);
%!         [res(j).t, ...
%!          res(j).trajectory, ...
%!          res(j).deformation, ...
%!          res(j).velocity, ...
%!          res(j).acceleration, ...
%!          res(j).node_id] = mbdyn_post_load_output_struct([opt_mbdyn.output_file]);
%!         [res(j).joint_id, ...
%!          res(j).local_reaction, ...
%!          res(j).global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!         switch (j)
%!         case 1
%!         if (options.post_proc)
%!           opt_scale.scale_type = "least square";
%!           opt_scale.scale = 300;
%!           opt_post.print_and_exit = false;
%!           opt_scale.output_stress = FEM_SCA_STRESS_VMIS;
%!           opt_post.elem_types = {"iso27", "penta18"};
%!           fem_post_cms_sol_export(cms_data, opt_mbdyn.output_file, opt_mbdyn.output_file, opt_scale, opt_post);
%!           spawn_wait(spawn("gmsh", {[opt_mbdyn.output_file, "_struct.geo"]}));
%!         endif
%!         case 3
%!           [res(j).mesh, res(j).sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!         endswitch
%!       endfor
%!     endfor
%!     tol = 1e-2;
%!     for i=1:2
%!       assert_simple(res(1).trajectory{i}, res(2).trajectory{i}, tol * max(max(abs(res(2).trajectory{i}))));
%!       assert_simple(res(1).velocity{i}, res(2).velocity{i}, tol * max(max(abs(res(2).velocity{i}))));
%!     endfor
%!     tol = 2e-3;
%!     assert_simple(res(1).trajectory{1}, res(3).trajectory{node_idx_tip}, tol * max(max(abs(res(3).trajectory{node_idx_tip}))));
%!     assert_simple(res(1).trajectory{2}, res(3).trajectory{node_idx_center}, tol * max(max(abs(res(3).trajectory{node_idx_center}))));
%!     assert_simple(res(1).velocity{1}, res(3).velocity{node_idx_tip}, tol * max(max(abs(res(3).velocity{node_idx_tip}))));
%!     assert_simple(res(1).velocity{2}, res(3).velocity{node_idx_center}, tol * max(max(abs(res(3).velocity{node_idx_center}))));
%!     tol = 3e-2;
%!     a3_tip = (res(3).velocity{node_idx_tip}(2:end, :) - res(3).velocity{node_idx_tip}(1:end-1, :)) / (res(3).t(2) - res(3).t(1));
%!     a3_tip = interp1(0.5 * (res(3).t(2:end) + res(3).t(1:end-1)), a3_tip, res(3).t, "extrap");
%!     assert_simple(res(1).acceleration{1}(floor(0.1 * end):end, :), a3_tip(floor(0.1 * end):end, :), tol * max(max(abs(a3_tip))));
%!     tol = 1e-2;
%!     assert_simple(res(1).global_reaction{1}(ceil(0.1 * end):end, :), res(3).global_reaction{end}(ceil(0.1 * end):end,:), tol * max(max(abs(res(3).global_reaction{end}))));
%!  ## unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!       endfor
%!     endif
%!  ## end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
