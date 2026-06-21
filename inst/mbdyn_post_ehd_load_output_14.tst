## mbdyn_post_ehd_load_output.tst:01
%!test
%! try
%!   pkg load mbdyn_util_oct;
%!   ## References: Elastohydrodynamic Lubrication in a Cylindrical Journal Bearing
%!   ##             COMSOL Multiphysics 6.4
%!   SI_unit_meter = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   omega = [100, 316.23, 1000, 3162.3, 10000] * pi / 30 / SI_unit_second^-1;
%!   omega = omega(1);
%!   res = repmat(struct(), 1, numel(omega));
%!   for j=1:numel(omega)
%!     param.F1 = 20e3 / SI_unit_newton;
%!     param.omega = omega(j);
%!     param.d = 200e-3 / SI_unit_meter;
%!     param.h0 = 200e-6 / SI_unit_meter;
%!     param.D = param.d + 2 * param.h0;
%!     param.B = 67e-3 / SI_unit_meter;
%!     param.wg = 10e-3 / SI_unit_meter;
%!     param.eta0 = 75e-3 / (SI_unit_pascal * SI_unit_second);
%!     param.alpha = 25e-9 / SI_unit_pascal^-1;
%!     param.rho0 = 850 / (SI_unit_kilogram / SI_unit_meter^3);
%!     param.Arho1 = 0*0.6e-9 / SI_unit_pascal^-1;
%!     param.Arho2 = 0*1.7e-9 / SI_unit_pascal^-1;
%!     param.pref = 100 / SI_unit_pascal;
%!     param.p_side = param.pref;
%!     param.p_in = 2e5 / SI_unit_pascal;
%!     param.pmax = 100e5 / SI_unit_pascal;
%!     param.M = int32(20);
%!     param.N = int32(85);
%!     param.output_bearing_data = true;
%!     param.cavitation_model = "mass conserving";
%!     param.hydraulic_nodes = true;
%!     param.jacobian_check = false;
%!     param.nonlinear_solver = "nox";
%!     fd = -1;
%!     output_file = "";
%!     unwind_protect
%!       unwind_protect
%!         output_dir = tempdir();
%!         [fd, output_file] = mkstemp(fullfile(output_dir, "oct-mbdyn_post_ehd_load_output_XXXXXX"));
%!         mbdyn_pre_write_param_file(fd, param);
%!         fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!         fputs(fd, "set: integer node_id_stator = 1002;\n");
%!         if (param.hydraulic_nodes)
%!           fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!           fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!           fputs(fd, "set: integer hyd_node_id_inlet = 2003;\n");
%!         endif
%!         fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!         fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!         fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!         fputs(fd, "set: integer force_id_load = 4002;\n");
%!         fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!         fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!         fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!         if (param.hydraulic_nodes)
%!           fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!           fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!           fputs(fd, "set: integer genel_id_inlet = 4005;\n");
%!         endif
%!         fputs(fd, "set: real n = 3;\n");
%!         fputs(fd, "set: real t1 = abs(2 * pi * n / omega);\n");
%!         fputs(fd, "set: real dt = t1 / 100;\n");
%!         fputs(fd, "begin: data;\n");
%!         fputs(fd, "        problem: initial value; # the default\n");
%!         fputs(fd, "end: data;\n");
%!         fputs(fd, "begin: initial value;\n");
%!         fputs(fd, "        initial time: 0;\n");
%!         fputs(fd, "        final time: t1;\n");
%!         fputs(fd, "        time step: dt;\n");
%!         fputs(fd, "        max iterations: 50;\n");
%!         fputs(fd, "        tolerance: 1e-6, test, norm, 1e-4, test, norm;\n");
%!         fputs(fd, "        linear solver: umfpack, grad, colamd, scale, row max column max, always;\n");
%!         fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e5;\n");
%!         fputs(fd, "        threads: assembly, 1;\n");
%!         fputs(fd, "        method: hybrid, ms, 0.;\n");
%!         fputs(fd, "        output: messages;\n");
%!         fputs(fd, "        derivatives tolerance: 3e-3, 1e-3;\n");
%!         fputs(fd, "        derivatives max iterations: 20;\n");
%!         fputs(fd, "        derivatives coefficient: 1e-3, auto;\n");
%!         fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!         switch (param.nonlinear_solver)
%!         case "nox"
%!         fputs(fd, "        nonlinear solver: nox,\n");
%!         fputs(fd, "             jacobian operator, newton,\n");
%!         fputs(fd, "             solver, line search based,\n");
%!         fputs(fd, "             forcing term, type 2,\n");
%!         fputs(fd, "             direction, newton,\n");
%!         fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!         fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!         fputs(fd, "             linear solver, gmres,\n");
%!         fputs(fd, "             linear solver max iterations, 100,\n");
%!         fputs(fd, "             krylov subspace size, 100;\n");
%!         case "mcp"
%!         fputs(fd, "        nonlinear solver: mcp newton min fb;\n");
%!         otherwise
%!         fputs(fd, "        nonlinear solver: linesearch, default solver options, heavy nonlinear, divergence check, no, lambda min, 0., tolerance x, 1e-4;\n");
%!         endswitch
%!         fputs(fd, "end: initial value;\n");
%!         fputs(fd, "begin: control data;\n");
%!         if (param.jacobian_check)
%!           fputs(fd, "    finite difference jacobian meter: const, 1., iterations, string, \"Var == 1\", coefficient, 1e-8, order, 1, output, none, statistics iteration, yes;\n");
%!         endif
%!         fputs(fd, "    use automatic differentiation;\n");
%!         fputs(fd, "    skip initial joint assembly;\n");
%!         fputs(fd, "    #output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!         fputs(fd, "    use: loadable elements, in assembly;\n");
%!         fputs(fd, "    default orientation: euler123;\n");
%!         fputs(fd, "        structural nodes: 2;\n");
%!         fputs(fd, "        joints: 2;\n");
%!         fputs(fd, "        loadable elements: 1;\n");
%!         fputs(fd, "        forces: 1;\n");
%!         if (param.hydraulic_nodes)
%!           fputs(fd, "        hydraulic nodes: 3;\n");
%!           fputs(fd, "        genels: 3;\n");
%!         endif
%!         fputs(fd, "        print: dof stats, to file;\n");
%!         fputs(fd, "        print: equation description, to file;\n");
%!         fputs(fd, "        print: dof description, to file;\n");
%!         fputs(fd, "        print: element connection, to file;\n");
%!         fputs(fd, "        print: node connection, to file;\n");
%!         fputs(fd, "    output precision: 8;\n");
%!         fputs(fd, "end: control data;\n");
%!         fputs(fd, "reference: ref_id_bearing,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, eye,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, null;");
%!         fputs(fd, "reference: ref_id_rotor,\n");
%!         fputs(fd, "        reference, ref_id_bearing, null,\n");
%!         fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!         fputs(fd, "        reference, ref_id_bearing, null,\n");
%!         fputs(fd, "        reference, ref_id_bearing, 0., 0., omega;\n");
%!         fputs(fd, "reference: ref_id_stator,\n");
%!         fputs(fd, "        reference, ref_id_bearing, null,\n");
%!         fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!         fputs(fd, "        reference, ref_id_bearing, null,\n");
%!         fputs(fd, "        reference, ref_id_bearing, null;\n");
%!         fputs(fd, "begin: nodes;\n");
%!         fputs(fd, "        structural: node_id_rotor, static,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null;\n");
%!         fputs(fd, "        structural: node_id_stator, static,\n");
%!         fputs(fd, "                reference, ref_id_stator, null,\n");
%!         fputs(fd, "                reference, ref_id_stator, eye,\n");
%!         fputs(fd, "                reference, ref_id_stator, null,\n");
%!         fputs(fd, "                reference, ref_id_stator, null;\n");
%!         if (param.hydraulic_nodes)
%!           fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_side;\n");
%!           fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_side;\n");
%!           fputs(fd, "    hydraulic: hyd_node_id_inlet, p_in;\n");
%!         endif
%!         fputs(fd, "end: nodes;\n");
%!         fputs(fd, "begin: elements;\n");
%!         fputs(fd, "        force: force_id_load, absolute, node_id_rotor,\n");
%!         fputs(fd, "               position, reference, ref_id_rotor, null,\n");
%!         fputs(fd, "               0., -1., 0., string, \"F1 * (sin(pi / 2. * Time / (t1 / 2))^2 * (Time <= (t1 / 2)) + (Time > (t1 / 2)))\";\n");
%!         fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!         fputs(fd, "                node_id_rotor,\n");
%!         fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!         fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!         fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!         fputs(fd, "        position constraint,\n");
%!         fputs(fd, "            inactive,\n");
%!         fputs(fd, "            inactive,\n");
%!         fputs(fd, "            active,\n");
%!         fputs(fd, "        component,\n");
%!         fputs(fd, "            const, 0.,\n");
%!         fputs(fd, "            const, 0.,\n");
%!         fputs(fd, "            const, 0.,\n");
%!         fputs(fd, "        orientation constraint,\n");
%!         fputs(fd, "            active,\n");
%!         fputs(fd, "            active,\n");
%!         fputs(fd, "            angular velocity,\n");
%!         fputs(fd, "        component,\n");
%!         fputs(fd, "            const, 0.,\n");
%!         fputs(fd, "            const, 0.,\n");
%!         fputs(fd, "            omega;\n");
%!         fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!         fputs(fd, "                node_id_stator,\n");
%!         fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!         fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!         fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!         fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!         fputs(fd, "        position constraint,\n");
%!         fputs(fd, "            active,\n");
%!         fputs(fd, "            active,\n");
%!         fputs(fd, "            active,\n");
%!         fputs(fd, "        component,\n");
%!         fputs(fd, "            const, 0.,\n");
%!         fputs(fd, "            const, 0.,\n");
%!         fputs(fd, "            const, 0.,\n");
%!         fputs(fd, "        orientation constraint,\n");
%!         fputs(fd, "            active,\n");
%!         fputs(fd, "            active,\n");
%!         fputs(fd, "            active,\n");
%!         fputs(fd, "        component,\n");
%!         fputs(fd, "            const, 0.,\n");
%!         fputs(fd, "            const, 0.,\n");
%!         fputs(fd, "            const, 0.;\n");
%!         if (param.hydraulic_nodes)
%!           fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!           fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_side;\n");
%!           fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!           fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_side;\n");
%!           fputs(fd, "    genel: genel_id_inlet, clamp,\n");
%!           fputs(fd, "        hyd_node_id_inlet, hydraulic, p_in;\n");
%!         endif
%!         fputs(fd, "    user defined: elem_id_bearing,\n");
%!         fputs(fd, "        hydrodynamic plain bearing2,\n");
%!         switch(param.cavitation_model)
%!           case "non mass conserving"
%!             fputs(fd, "            hydraulic fluid, incompressible,\n");
%!             fputs(fd, "                density, rho0,\n");
%!             fputs(fd, "                viscosity, eta0,\n");
%!             fputs(fd, "                pressure, pref,\n");
%!             fputs(fd, "                temperature, 0,\n");
%!             fputs(fd, "                alpha, alpha,\n");
%!             fputs(fd, "                Arho1, Arho1,\n");
%!             fputs(fd, "                Arho2, Arho2,\n");
%!           case "mass conserving"
%!             fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!             fputs(fd, "                density, rho0,\n");
%!             fputs(fd, "                1.,\n");
%!             fputs(fd, "                0.,\n");
%!             fputs(fd, "                viscosity, eta0,\n");
%!             fputs(fd, "                temperature, 0,\n");
%!             fputs(fd, "                alpha, alpha,\n");
%!             fputs(fd, "                Arho1, Arho1,\n");
%!             fputs(fd, "                Arho2, Arho2,\n");
%!             fputs(fd, "            viscosity vapor, factor, 1e-3,\n");
%!         endswitch
%!         fputs(fd, "            mesh, linear finite difference,\n");
%!         switch (param.nonlinear_solver)
%!         case "mcp"
%!           fputs(fd, "          enable mcp, yes,\n");
%!         endswitch
%!         fputs(fd, "            geometry, cylindrical,\n");
%!         fputs(fd, "                mesh position, at bearing,\n");
%!         fputs(fd, "                bearing width, B,\n");
%!         fputs(fd, "                shaft diameter, d,\n");
%!         fputs(fd, "                bearing diameter, D,\n");
%!         fputs(fd, "                shaft node,\n");
%!         fputs(fd, "                    node_id_rotor,\n");
%!         fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                bearing node,\n");
%!         fputs(fd, "                    node_id_stator,\n");
%!         fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!         fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!         fputs(fd, "            number of nodes z, M,\n");
%!         fputs(fd, "            number of nodes Phi, N,\n");
%!         if (param.hydraulic_nodes)
%!           fputs(fd, "            pressure coupling conditions axial,\n");
%!           fputs(fd, "                hyd_node_id_outlet1,\n");
%!           fputs(fd, "                hyd_node_id_outlet2,\n");
%!           fputs(fd, "            pressure coupling conditions radial, 2,\n");
%!           fputs(fd, "                position, 0., 0.,\n");
%!           fputs(fd, "                rectangle, width, wg, height, B,\n");
%!           fputs(fd, "                pressure node, hyd_node_id_inlet,\n");
%!           fputs(fd, "                position, D * pi / 2., 0.,\n");
%!           fputs(fd, "                rectangle, width, wg, height, B,\n");
%!           fputs(fd, "                pressure node, hyd_node_id_inlet,\n");
%!         else
%!           fputs(fd, "         boundary conditions, const, p_side, const, p_side,\n");
%!           fputs(fd, "         lubrication grooves, 2,\n");
%!           fputs(fd, "                at bearing,\n");
%!           fputs(fd, "                pressure, const, p_in,\n");
%!           fputs(fd, "                position, 0., 0.,\n");
%!           fputs(fd, "                rectangle, width, wg, height, B,\n");
%!           fputs(fd, "                at bearing,\n");
%!           fputs(fd, "                pressure, const, p_in,\n");
%!           fputs(fd, "                position, D * pi / 2., 0.,\n");
%!           fputs(fd, "                rectangle, width, wg, height, B,\n");
%!         endif
%!         fputs(fd, "            pressure dof scale, pmax,\n");
%!         fputs(fd, "            output pressure, yes,\n");
%!         fputs(fd, "            output density, yes,\n");
%!         fputs(fd, "            output clearance, yes,\n");
%!         fputs(fd, "            output clearance derivative, no,\n");
%!         fputs(fd, "            output velocity, no,\n");
%!         fputs(fd, "            output stress, no,\n");
%!         fputs(fd, "            output reaction force, yes,\n");
%!         fputs(fd, "            output friction loss, yes,\n");
%!         fputs(fd, "            output mesh, yes,\n");
%!         fputs(fd, "            output, output_bearing_data;\n");
%!         fputs(fd, "end: elements;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!           fd = -1;
%!         endif
%!       end_unwind_protect
%!       opt_sol.output_file = output_file;
%!       shell(sprintf("nl %s", output_file));
%!       info = mbdyn_solver_run(output_file, opt_sol);
%!       res(j).log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!       [res(j).t, ...
%!        res(j).trajectory, ...
%!        res(j).deformation, ...
%!        res(j).velocity, ...
%!        res(j).acceleration, ...
%!        res(j).node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!       if (param.hydraulic_nodes)
%!         [res(j).genel_id, res(j).genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!       endif
%!       res(j).log_dat.vars = mbdyn_post_id_to_index(res(j), res(j).log_dat.vars);
%!       opt_load.verbose = false;
%!       opt_load.num_steps = numel(res.t);
%!       opt_load.output_index = 1:numel(res.t);
%!       opt_load.loaded_fields = {};
%!       opt_load.interpolate_mesh = true;
%!       res(j).bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res(j).log_dat, opt_load);
%!     unwind_protect_cleanup
%!       if (numel(output_file))
%!         fn = dir([output_file, "*"]);
%!         for i=1:numel(fn)
%!           fn_i = fullfile(fn(i).folder, fn(i).name);
%!           status = unlink(fn_i);
%!           if (status ~= 0)
%!             warning("failed to remove file \"%s\"", fn_i);
%!           endif
%!         endfor
%!       endif
%!     end_unwind_protect
%!   endfor
%!   figure("visible", "off");
%!   hold on;
%!   for j=1:numel(res)
%!     plot(1e3 * res(j).bearings.xi(1,:) * SI_unit_meter, 1e-5 * res(j).bearings.columns.p(floor(end/2),:, end) * SI_unit_pascal, sprintf("-;%.2f;", 30 / pi * omega(j) * SI_unit_second^-1));
%!   endfor
%!   xlabel("x [mm]");
%!   ylabel("p [bar]");
%!   title("midplane pressure");
%!   grid minor on;
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
