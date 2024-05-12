## mbdyn_post_load_output_eig.tst:06
%!test
%! try
%! ## TEST 330
%! ## VIBRATIONS OF COMPLETE SPHERICAL SHELLS WITH IMPERFECTIONS
%! ## Thomas A. Duffey
%! ## Jason E. Pepin
%! ## Amy N. Robertson
%! ## Michael L. Steinzig
%! ## Internatial Modal Analysis Conference (IMAC-XXIII)
%! ## Orlando, Florida
%! ## January 31-February 3, 2005
%! ## Los Alamos
%! ## NATIONAL LABORATORY
%! pkg load mboct-fem-pkg;
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     options.verbose = true;
%!     SI_unit_meter = 1;
%!     SI_unit_second = 1;
%!     SI_unit_kilogram = 1;
%!     SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!     SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!     SI_unit_rad = 1;
%!     E = 28e6 * 6895 / SI_unit_pascal;
%!     nu = 0.28;
%!     rho = 0.000751 * 4.4482 / (25.4e-3^4) / (SI_unit_kilogram / SI_unit_meter^3);
%!     R = 4.4688 * 25.4e-3 / SI_unit_meter;
%!     t = 0.0625 * 25.4e-3 / SI_unit_meter;
%!     W0 = [0; 0; 0] / (SI_unit_rad / SI_unit_second);
%!     fmin = 1. / SI_unit_second^-1;
%!     fmax = 10000 / SI_unit_second^-1;
%!     number_of_modes = 39*2;
%!     t1 = 1 / SI_unit_second;
%!     N = 1;
%!     mesh_size = 2 * R * pi / 36 / SI_unit_meter;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "R = %g;\n", R);
%!     fprintf(fd, "t = %g;\n", t);
%!     fprintf(fd, "h = %g;\n", mesh_size);
%!     fputs(fd, "Point(1) = {0,0,-R - 0.5 * t,h};\n");
%!     fputs(fd, "Point(2) = {0,0,-R + 0.5 * t,h};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "A[] = Extrude{{0,1,0},{0,0,0},Pi}{Line{1}; Layers{Ceil(R * Pi / h)};Recombine;};\n");
%!     fputs(fd, "V1[] = Extrude{{0,0,1},{0,0,0},Pi}{Surface{A[1]};Layers{Ceil(R * Pi/h)};Recombine;};\n");
%!     fputs(fd, "V2[] = Extrude{{0,0,1},{0,0,0},Pi}{Surface{V1[0]};Layers{Ceil(R * Pi/h)};Recombine;};\n");
%!     fputs(fd, "Physical Volume(\"volume\", 1) = {V1[1],V2[1]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   [~] = unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-0", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_msh.elem_type = {"iso20r", "penta15"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   mesh.materials.iso20r = ones(rows(mesh.elements.iso20r), 1, "int32");
%!   mesh.materials.penta15 = ones(rows(mesh.elements.penta15), 1, "int32");
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   elem_file = [filename, ".elem"];
%!   mbdyn_file = [filename, ".mbdyn"];
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_center";
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fputs(fd, "set: integer ref_id_center = 1001;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fprintf(fd, "    final time: %g;\n", t1);
%!     fprintf(fd, "    time step: %g;\n", t1 / N);
%!     fputs(fd, "    max iterations: 100;\n");
%!     fprintf(fd, "    tolerance: %g, test, norm, %g, test, norm;\n", 1e-4 / SI_unit_newton, 1e-6 / SI_unit_meter);
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 100;\n");
%!     fputs(fd, "    method: ms4, 0.6;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30, linear solver max iterations, 60, verbose, 3;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    eigenanalysis: list, 1, 0.,\n");
%!     fputs(fd, "    output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fputs(fd, "        results output precision, 16,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", fmin, fmax);
%!     fprintf(fd, "    use arpack,%d,%d,0.,suffix format,\"%%02d\";\n", number_of_modes, 2 * number_of_modes + 1);
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, "     output results: netcdf, text;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " reference: ref_id_center,\n");
%!     fputs(fd, "   position, reference, global, null,\n");
%!     fputs(fd, "   orientation, reference, global, eye,\n");
%!     fputs(fd, "   velocity, reference, global, null,\n");
%!     fprintf(fd, "   angular velocity, reference, global%s;\n", sprintf(", %g", W0));
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!     options_mbd.output_file = sprintf("%s_mbd", filename);
%!     options_eig.positive_frequencies = false;
%!     if (~options.verbose)
%!       options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!     endif
%!     info = mbdyn_solver_run(mbdyn_file, options_mbd);
%!     [mesh_sol, sol] = mbdyn_post_load_output_sol(options_mbd.output_file);
%!     options_eig.use_netcdf = false;
%!     modal = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig);
%!     options_eig.use_netcdf = true;
%!     modalnc = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig);
%!     fref = [5078; 6005; 6378; 6729] / SI_unit_second^-1;
%!     assert_simple(modal.f([7, 12, 19, 39] - 6), fref, 5e-4 * max(fref));
%!     assert_simple(modalnc.f([7, 12, 19, 39] - 6), fref, 5e-4 * max(fref));
%!     fn = fieldnames(modal);
%!     for idxfn=1:numel(fn)
%!       assert_simple(getfield(modal, fn{idxfn}), getfield(modalnc, fn{idxfn}));
%!     endfor
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
