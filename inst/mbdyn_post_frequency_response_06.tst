## mbdyn_post_frequency_response.tst:06
%!test
%! try
%! ## TEST 4
%! pkg load mboct-fem-pkg;
%! close all;
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.number_of_modes = 50;
%! param.Udyn = eye(3) / SI_unit_meter;
%! param.fmin = 0 / (SI_unit_second^-1);
%! param.fmax = 1500 / (SI_unit_second^-1);
%! param.num_freq = 250;
%! param.t1 = 1000 / SI_unit_second;
%! param.num_steps = 10;
%! param.damping = "beta";
%! geometry(1).user_data.helspr.L = 25.8e-3 / SI_unit_meter;
%! geometry(1).user_data.helspr.Di = 12.12e-3 / SI_unit_meter;
%! geometry(1).user_data.helspr.d = 1.3e-3 / SI_unit_meter;
%! geometry(1).user_data.helspr.n = 1.3;
%! geometry(1).user_data.helspr.ni = 3;
%! geometry(1).user_data.helspr.ng = 0.75;
%! geometry(1).user_data.color = "r";
%! geometry(2).user_data.helspr.L = 27.7e-3 / SI_unit_meter;
%! geometry(2).user_data.helspr.Di = 12.12e-3 / SI_unit_meter;
%! geometry(2).user_data.helspr.d = 1.3e-3 / SI_unit_meter;
%! geometry(2).user_data.helspr.n = 1.7;
%! geometry(2).user_data.helspr.ni = 2.7;
%! geometry(2).user_data.helspr.ng = 0.75;
%! geometry(2).user_data.color = "g";
%! geometry(3).user_data.helspr.L = 28.63e-3 / SI_unit_meter;
%! geometry(3).user_data.helspr.Di = 12.12e-3 / SI_unit_meter;
%! geometry(3).user_data.helspr.d = 1.25e-3 / SI_unit_meter;
%! geometry(3).user_data.helspr.n = 2;
%! geometry(3).user_data.helspr.ni = 2.7;
%! geometry(3).user_data.helspr.ng = 0.75;
%! geometry(3).user_data.color = "b";
%! geometry = geometry(1);
%! material.E = 206000e6 / SI_unit_pascal;
%! material.G = 81500e6 / SI_unit_pascal;
%! material.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%! damp.D = [1e-2; 1e-2];
%! damp.f = [20, 1500] / (SI_unit_second^-1);
%! param.Fz = -7.2*0 / SI_unit_newton;
%! switch (param.damping)
%! case "none"
%! case "Rayleigh"
%!   [material.alpha, material.beta] = fem_pre_mat_rayleigh_damping(damp.D, damp.f); ## alpha damping not yet supported by MBDyn
%! case "beta"
%!   material.beta = 2 * damp.D(end) / (2 * pi * damp.f(end));
%! endswitch
%! param.omega = 2 * pi * linspace(param.fmin, param.fmax, param.num_freq);
%! param.f_plot_response = false;
%! function [x, y, z, R, Phi] = helspr_geo(geo, r, s, t)
%!   Phi = 2 * pi * (r * geo.n + geo.ni);
%!   r1 = 0.5 * geo.d * s;
%!   Theta = 2 * pi * t;
%!   x1 = [0.5 * geo.D * cos(Phi);
%!         0.5 * geo.D * sin(Phi);
%!         (geo.L - geo.d * (2 * (geo.ni - geo.ng) + 1)) * r + geo.d * (geo.ni - geo.ng + 0.5)];
%!   x2 = [0;
%!         0;
%!         x1(3)];
%!   e2 = x2 - x1;
%!   e1 = [-0.5 * geo.D * sin(Phi) * 2 * pi * geo.n;
%!          0.5 * geo.D * cos(Phi) * 2 * pi * geo.n;
%!          (geo.L - geo.d * (2 * (geo.ni - geo.ng) + 1))];
%!   e3 = cross(e1, e2);
%!   e2 = cross(e3, e1);
%!   R1 = [e1, e2, e3];
%!   R1 *= diag(1 ./ norm(R1, "cols"));
%!   assert_simple(R1.' * R1, eye(3), eps^0.9);
%!   assert_simple(R1 * R1.', eye(3), eps^0.9);
%!   x3 = R1 * [0; r1 * cos(Theta); r1 * sin(Theta)] + x1;
%!   x = x3(1);
%!   y = x3(2);
%!   z = x3(3);
%! endfunction
%!
%! function [F, locked] = boundary_cond(r, s, t, geo, loads)
%!   locked = [];
%!   F = [];
%! endfunction
%!
%! function p = pressure_boundary_cond(r, s, t, geometry, load, perm_idx)
%!   p = [];
%! endfunction
%!
%! function [modal, load_case_dof] = create_modal_data_mbdyn(mesh, dof_map, load_case_dof, load_case_stat, options, param)
%!   for i=1:numel(mesh.elements.joints)
%!     [ridx, cidx, constr] = find(mesh.elements.joints(i).C);
%!     if (numel(mesh.elements.joints(i).nodes) > 1)
%!       error("joints with multiple nodes not yet implemented");
%!     endif
%!     for j=1:numel(cidx)
%!       load_case_dof.locked_dof(mesh.elements.joints(i).nodes, cidx(j)) = true;
%!     endfor
%!   endfor
%!   mesh.elements = rmfield(mesh.elements, "joints");
%!   filename = "";
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%! %unwind_protect
%!     nodes_file = [filename, ".nod"];
%!     csl_file = [filename, ".csl"];
%!     elem_file = [filename, ".elem"];
%!     mbdyn_file = [filename, ".mbdyn"];
%!     opt_mbd_mesh = struct();
%!     opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_center";
%!     opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!     opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_stat, elem_file, opt_mbd_mesh);
%!     unwind_protect
%!       [fd, msg] = fopen(mbdyn_file, "w");
%!       if (fd == -1)
%!         error("failed to open file \"%s\"", mbdyn_file);
%!       endif
%!       fputs(fd, "set: integer ref_id_center = 1001;\n");
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "    problem: initial value; # the default\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "    initial time: 0;\n");
%!       fprintf(fd, "    final time: %g;\n", param.t1);
%!       fprintf(fd, "    time step: %g;\n", param.t1 / param.num_steps);
%!       fputs(fd, "    max iterations: 100;\n");
%!       fprintf(fd, "    tolerance: %g, test, norm, %g, test, norm;\n", 1e-4, 1e-6);
%!       fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 100;\n");
%!       fputs(fd, "    method: implicit euler;\n");
%!       fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!       fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!       fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30, linear solver max iterations, 60, verbose, 3;\n");
%!       fputs(fd, "    threads: assembly, 1;\n");
%!       fputs(fd, "    threads: solver, 1;\n");
%!       fputs(fd, "    eigenanalysis: 0.,\n");
%!       fputs(fd, "    output eigenvectors,\n");
%!       fputs(fd, "    output sparse matrices,\n");
%!       fputs(fd, "        output geometry,\n");
%!       fputs(fd, "        results output precision, 16,\n");
%!       fputs(fd, "        mode, largest imaginary part,\n");
%!       fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", options.fmin, options.fmax);
%!       fprintf(fd, "    use arpack,%d,%d,0.,suffix format,\"%%02d\";\n", param.number_of_modes, 2 * param.number_of_modes + 1);
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "        use automatic differentiation;\n");
%!       fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!       fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!       fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!       fprintf(fd, "    genels: %d;\n", opt_mbd_mesh.genels.number);
%!       fputs(fd, "     default output: none;\n");
%!       fputs(fd, "     output results: netcdf, no text;\n");
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " reference: ref_id_center,\n");
%!       fputs(fd, "   position, reference, global, null,\n");
%!       fputs(fd, "   orientation, reference, global, eye,\n");
%!       fputs(fd, "   velocity, reference, global, null,\n");
%!       fprintf(fd, "   angular velocity, reference, global, null;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", nodes_file);
%!       fputs(fd, " end: nodes;\n");
%!       fprintf(fd, "include: \"%s\";\n", csl_file);
%!       fputs(fd, " begin: elements;\n");
%!       fprintf(fd, "include: \"%s\";\n", elem_file);
%!       fputs(fd, " end: elements;\n");
%!       options_mbd.output_file = sprintf("%s_mbd", filename);
%!       options_eig.positive_frequencies = false;
%!       if (~options.verbose)
%!         options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!       endif
%!       info = mbdyn_solver_run(mbdyn_file, options_mbd);
%!       options_eig.use_netcdf = true;
%!       modal = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig);
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!       fd = -1;
%!     end_unwind_protect
%! %unwind_protect_cleanup
%!     if (~isempty(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!       endfor
%!     endif
%! %end_unwind_protect
%! endfunction
%!
%! function [Freact, Freact_mbd, X1_mbd, mesh, mesh_def, sol_stat, sol_eig, sol_eig_def, modal_mbd, load_case_mbd] = transmissibility(geometry, material, param)
%!   material.nu = material.E / (2 * material.G) - 1;
%!   geometry.user_data.helspr.D = geometry.user_data.helspr.Di + geometry.user_data.helspr.d;
%!   kz = material.G * geometry.user_data.helspr.d^4 / (8 * geometry.user_data.helspr.n * geometry.user_data.helspr.D^3);
%!   Ustat = [0; 0; param.Fz / kz];
%!   h = geometry.user_data.helspr.d * pi;
%!   geometry.user_data.helspr.nPhi = max([2, round(sqrt((geometry.user_data.helspr.D * pi * geometry.user_data.helspr.n)^2 + geometry.user_data.helspr.L^2) / h)]) + 1;
%!   geometry.user_data.helspr.nr = max([1, round(0.5 * geometry.user_data.helspr.d / h)]) + 1;
%!   geometry.user_data.helspr.nTheta = max([3, round(geometry.user_data.helspr.d * pi / h)]) + 1;
%!   geometry.mesh_size.r = linspace(0, 1, geometry.user_data.helspr.nPhi);
%!   geometry.mesh_size.s = linspace(0, 1, geometry.user_data.helspr.nr);
%!   geometry.mesh_size.t = linspace(0, 1, geometry.user_data.helspr.nTheta);
%!   geometry.sewing.tolerance = sqrt(eps) * geometry.user_data.helspr.D;
%!   geometry.spatial_coordinates = @(r, s, t, geo) feval("helspr_geo", geometry.user_data.helspr, r, s, t);
%!   geometry.material_selector = @(r, s, t, geo) 1;
%!   geometry.boundary_condition = @(r, s, t, geo, load) feval("boundary_cond", r, s, t, geo, load);
%!   geometry.pressure_boundary_condition = @(r, s, t, geo, load, perm_idx) feval("pressure_boundary_cond", r, s, t, geo, load, perm_idx);
%!   options.elem_type = "iso20";
%!   loads = struct();
%!   [mesh, load_case_dof] = fem_pre_mesh_struct_create(geometry, loads, material, options);
%!   idx_node_bottom = unique(mesh.structured.inode_idx(1, :, :)(:));
%!   idx_node_top = unique(mesh.structured.inode_idx(end, :, :)(:));
%!   idx_node_bottom = idx_node_bottom(idx_node_bottom > 0);
%!   idx_node_top = idx_node_top(idx_node_top > 0);
%!   idx_node_joint = [idx_node_bottom; idx_node_top];
%!   empty_cell = cell(1, numel(idx_node_joint));
%!   mesh.elements.joints = struct("nodes", mat2cell(idx_node_joint, ones(numel(idx_node_joint), 1, "int32"), 1), "C", repmat({[eye(3), zeros(3, 3)]}, numel(idx_node_joint), 1));
%!   [dof_map] = fem_ass_dof_map(mesh, load_case_dof);
%!   load_case_stat = fem_pre_load_case_create_empty(1);
%!   for i=1:numel(load_case_stat)
%!     load_case_stat(i).joints = struct("U", repmat({zeros(3, 1)}, numel(idx_node_joint), 1));
%!     for j=1:numel(idx_node_top)
%!       load_case_stat(i).joints(numel(idx_node_bottom) + j).U = Ustat(:, i);
%!     endfor
%!   endfor
%!   load_case_dyn = fem_pre_load_case_create_empty(columns(param.Udyn));
%!   for i=1:numel(load_case_dyn)
%!     load_case_dyn(i).joints = struct("U", repmat({zeros(3, 1)}, numel(idx_node_joint), 1));
%!     for j=1:numel(idx_node_top)
%!       load_case_dyn(i).joints(numel(idx_node_bottom) + j).U = param.Udyn(:, i);
%!     endfor
%!   endfor
%!   [mat_ass.M, ...
%!    mat_ass.K, ...
%!    mat_ass.R, ...
%!    mat_ass.mat_info, ...
%!    mat_ass.mesh_info] = fem_ass_matrix(mesh, ...
%!                                        dof_map, ...
%!                                        [FEM_MAT_MASS, ...
%!                                         FEM_MAT_STIFFNESS, ...
%!                                         FEM_VEC_LOAD_CONSISTENT], ...
%!                                        load_case_stat);
%!   [sol_stat] = fem_sol_static(mesh, dof_map, mat_ass);
%!   [sol_eig] = fem_sol_modal(mesh, dof_map, mat_ass, param.number_of_modes);
%!   options_mbd.fmin = 0.8 * min(sol_eig.f);
%!   options_mbd.fmax = 1.2 * max(sol_eig.f);
%!   options_mbd.verbose = true;
%!   [modal_mbd, load_case_mbd] = create_modal_data_mbdyn(mesh, dof_map, load_case_dof, load_case_stat, options_mbd, param);
%!   sol_stat.stress = fem_ass_matrix(mesh, ...
%!                                    dof_map, ...
%!                                    [FEM_VEC_STRESS_CAUCH], ...
%!                                    load_case_stat, ...
%!                                    sol_stat);
%!   for i=1:numel(load_case_dyn)
%!     load_case_dyn(i).tau0 = sol_stat.stress.tau;
%!   endfor
%!   mesh_def = mesh;
%!   mesh_def.nodes += sol_stat.def;
%!   [mat_ass_def.M, ...
%!    mat_ass_def.D, ...
%!    mat_ass_def.K, ...
%!    mat_ass_def.KTAU0, ...
%!    mat_ass_def.R, ...
%!    mat_ass_def.mat_info, ...
%!    mat_ass_def.mesh_info] = fem_ass_matrix(mesh_def, ...
%!                                            dof_map, ...
%!                                            [FEM_MAT_MASS, ...
%!                                             FEM_MAT_DAMPING, ...
%!                                             FEM_MAT_STIFFNESS, ...
%!                                             FEM_MAT_STIFFNESS_TAU0, ...
%!                                             FEM_VEC_LOAD_CONSISTENT], ...
%!                                            load_case_dyn);
%!   mat_ass_def.K += mat_ass_def.KTAU0;
%!   [sol_eig_def] = fem_sol_modal(mesh_def, dof_map, mat_ass_def, param.number_of_modes, struct("solver", "pardiso", "refine_max_iter", int32(10)));
%!   Freact = complex(zeros(3, columns(mat_ass_def.R), numel(param.omega)));
%!   opt_factor.number_of_threads = mbdyn_solver_num_threads_default();
%!   opt_factor.verbose = int32(0);
%!   opt_factor.solver = "pardiso";
%!   opt_factor.refine_max_iter = int32(20);
%!   for i=1:numel(param.omega)
%!     fprintf(stderr, "%3d:%.2fHz\n", i, param.omega(i) / (2 * pi));
%!     A = -param.omega(i)^2 * mat_ass_def.M + 1j * param.omega(i) * mat_ass_def.D + mat_ass_def.K;
%!     Uij = fem_sol_factor(A, opt_factor) \ mat_ass_def.R;
%!     for j=1:size(Freact, 2)
%!       Freact(:, j, i) = sum(Uij(:, j)(dof_map.edof.joints(1:numel(idx_node_bottom), :)), 1)(:) * mat_ass_def.mat_info.beta(3);
%!     endfor
%!   endfor
%!   X1_mbd = Freact_mbd = complex(zeros(3, 3, numel(param.omega)));
%!   df_dy_dot = (modal_mbd.Aplus + modal_mbd.Aminus) / 2;
%!   df_dy = (modal_mbd.Aplus - modal_mbd.Aminus) / (2 * modal_mbd.dCoef);
%!   f = zeros(rows(df_dy), 3);
%!   [nidx, dofidx] = find(load_case_mbd.locked_dof);
%!   idx_node_bottom = unique(mesh.structured.inode_idx(1, :, :)(:));
%!   idx_node_top = unique(mesh.structured.inode_idx(end, :, :)(:));
%!   idx_node_bottom = idx_node_bottom(idx_node_bottom > 0);
%!   idx_node_top = idx_node_top(idx_node_top > 0);
%!   for i=1:numel(idx_node_top)
%!     for j=1:3
%!       idx_joint = find(nidx == idx_node_top(i) & dofidx == j);
%!       if (isempty(idx_joint))
%!         continue;
%!       endif
%!       idx_genel = find(idx_joint == modal_mbd.genel_labels);
%!       if (isempty(idx_genel))
%!         continue;
%!       endif
%!       f(modal_mbd.genel_idx(idx_genel) + 1, j) += param.Udyn(j, j);
%!     endfor
%!   endfor
%!   idx_X1_mbd = modal_mbd.idx(find(idx_node_top(1) == modal_mbd.labels)) + (1:3);
%!   opt_factor.symmetric = false;
%!   for i=2:numel(param.omega) ## FIXME: A will be become singular if omega == 0
%!     fprintf(stderr, "%3d:%.2fHz\n", i, param.omega(i) / (2 * pi));
%!     A = df_dy + 1j * param.omega(i) * df_dy_dot;
%!     y = fem_sol_factor(A, opt_factor) \ (1j * param.omega(i) * f);
%!     for j=1:columns(y)
%!       X1_mbd(:, j, i) = y(idx_X1_mbd, j);
%!       for k=1:numel(idx_node_bottom)
%!         idx_joint = find(nidx == idx_node_bottom(k));
%!         if (isempty(idx_joint))
%!           continue;
%!         endif
%!         for l=1:numel(idx_joint)
%!           idx_genel = find(idx_joint(l) == modal_mbd.genel_labels);
%!           if (isempty(idx_genel))
%!             continue;
%!           endif
%!           Freact_mbd(dofidx(idx_joint(l)), j, i) += 1j * param.omega(i) * y(modal_mbd.genel_idx(idx_genel) + 1, j);
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%! endfunction
%! Freact = Freact_mbd = X1_mbd = cell(1, numel(geometry));
%! for i=1:numel(geometry)
%!   [Freact{i}, Freact_mbd{i}, X1_mbd{i}, mesh, mesh_def, sol_stat, sol_eig, sol_eig_def, modal_mbd, load_case_mbd] = transmissibility(geometry(i), material, param);
%! endfor
%! if (param.f_plot_response)
%! figure("visible", "off");
%! hold on;
%! for i=1:numel(geometry)
%!   plot(param.omega / (2 * pi) * SI_unit_second, 20 * log10(SI_unit_newton / SI_unit_meter * sqrt(Freact{i}(1, 1, :).^2 + Freact{i}(2, 2, :).^2 + Freact{i}(3, 3, :).^2)), sprintf("-;fem-%d;%s", i, geometry(i).user_data.color));
%!   plot(param.omega / (2 * pi) * SI_unit_second, 20 * log10(SI_unit_newton / SI_unit_meter * sqrt(Freact_mbd{i}(1, 1, :).^2 + Freact_mbd{i}(2, 2, :).^2 + Freact_mbd{i}(3, 3, :).^2)), sprintf("-.;mbd-%d;%s", i, geometry(i).user_data.color));
%! endfor
%! xlabel("f [Hz]");
%! ylabel("kdyn [dB/(1N/m)]");
%! title("overall transmissibility");
%! grid minor on;
%! for j=1:columns(param.Udyn)
%!   figure("visible", "off");
%!   hold on;
%!   for i=1:numel(geometry)
%!     plot(param.omega / (2 * pi) * SI_unit_second, 20 * log10(SI_unit_newton / SI_unit_meter * abs(Freact{i}(i, i, :))), sprintf("-;fem-%d;%s", i, geometry(i).user_data.color));
%!     plot(param.omega / (2 * pi) * SI_unit_second, 20 * log10(SI_unit_newton / SI_unit_meter * abs(Freact{i}(i, i, :))), sprintf("--;mbd-%d;%s", i, geometry(i).user_data.color));
%!   endfor
%!   xlabel("f [Hz]");
%!   ylabel(sprintf("kdyn%s [dB/(1N/m)]", {"x", "y", "z"}{j}));
%!   title(sprintf("transmissibility %s-direction", {"x", "y", "z"}{j}));
%!   grid minor on;
%! endfor
%! endif
%! for i=1:numel(Freact)
%!   assert_simple(max(max(max(abs(Freact_mbd{i}(:, :, 2:end) ./ Freact{i}(:, :, 2:end) - 1)))) < 1e-4);
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
