%!test
%! ## TEST 1
%! ## References:
%! ## Hans Juergen Butenschoen
%! ## Das hydrodynamische, zylindrische Gleitlager endlicher Breite unter instationaerer Belastung, Karlsruhe 1976
%! ## Test case for rotation
%! B_d_r = [1, 1/2, 1/2.5, 1/3, 1/4, 1/5, 1/6, 1/8];
%!
%! epsilon_r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99];
%!
%! So_r =  [0.1196 0.0368  0.0243  0.0171  0.0098  0.0063  0.0044  0.0025
%!          0.2518 0.0783  0.0518  0.0366  0.0210  0.0136  0.0095  0.0054
%!          0.4091 0.1304  0.0867  0.0615  0.0354  0.0229  0.0160  0.0091
%!          0.6108 0.2026  0.1357  0.0968  0.0560  0.0364  0.0254  0.0144
%!          0.8903 0.3124  0.2117  0.1522  0.0888  0.0579  0.0406  0.0231
%!          1.3146 0.4982  0.3435  0.2496  0.1476  0.0969  0.0682  0.0390
%!          2.0432 0.8595  0.6075  0.4492  0.2708  0.1797  0.1274  0.0732
%!          3.5663 1.7339  1.2756  0.9687  0.6043  0.4085  0.2930  0.1706
%!          8.4352 5.0881  4.0187  3.2201  2.1595  1.5261  1.1263  0.6776
%!          18.7895 13.3083 11.2225 9.4993  6.9248  5.1833  3.9831  2.5202
%!          24.1172 17.7934 15.2823 13.1525 9.8517  7.5257  5.8712  3.7907
%!          33.1297 25.5920 22.4503 19.7034 15.2697 11.9804 9.5425  6.3397
%!          51.4774 41.9230 37.7412 33.9523 27.5040 22.3556 18.3874 12.7695
%!          107.7868 93.7881 87.2906 81.1597 70.0359 60.3874 52.1425 39.2568];
%!
%! beta_r = [79.410 81.767  82.119  82.299  82.481  82.561  82.608  82.653
%!           73.858 75.142  75.283  75.344  75.287  75.406  75.414  75.422
%!           68.264 68.493  68.423  68.368  68.204  68.261  68.238  68.211
%!           62.571 61.778  61.544  61.382  61.208  61.115  61.062  61.007
%!           56.705 54.993  54.603  54.348  54.069  53.922  53.863  53.784
%!           50.536 48.049  47.521  47.191  46.825  46.647  46.554  46.449
%!           43.859 40.803  40.156  39.756  39.326  39.108  38.983  38.865
%!           36.235 32.938  32.216  31.761  31.249  30.988  30.840  30.692
%!           26.482 23.566  22.849  22.368  21.790  21.476  21.289  21.089
%!           19.450 17.265  16.648  16.207  15.632  15.291  15.075  14.832
%!           17.609 15.660  15.089  14.669  14.109  13.768  13.547  13.262
%!           15.484 13.818  13.313  12.929  12.299  12.062  11.838  11.570
%!           12.903 11.598  11.178  10.850  10.375  10.057  9.835   9.558
%!            9.416  8.587   8.301   8.066   7.703   7.440   7.242   6.975] * pi / 180;
%!
%! Q_r = [0.1603  0.0939  0.0768  0.0648  0.0492  0.0396  0.0331  0.0249
%!        0.3196  0.1878  0.1536  0.1296  0.0984  0.0792  0.0662  0.0498
%!        0.4765  0.2816  0.2304  0.1944  0.1476  0.1188  0.0993  0.0747
%!        0.6318  0.3755  0.3072  0.2592  0.1968  0.1583  0.1324  0.0996
%!        0.7852  0.4694  0.384   0.324   0.246   0.198   0.1655  0.1245
%!        0.9374  0.5634  0.461   0.3889  0.2953  0.2376  0.1986  0.1494
%!        1.0888  0.6578  0.5383  0.454   0.3446  0.2772  0.2317  0.1743
%!        1.2392  0.7529  0.6159  0.5193  0.394   0.3169  0.2648  0.1992
%!        1.3898  0.8453  0.6944  0.5852  0.4436  0.3567  0.298   0.2241
%!        1.4655  0.8984  0.7343  0.6186  0.4687  0.3767  0.3147  0.2366
%!        1.4808  0.9084  0.7425  0.6254  0.4737  0.3807  0.3181  0.2391
%!        1.4968  0.9185  0.7507  0.6322  0.4788  0.3848  0.3214  0.2417
%!        1.5124  0.9287  0.7585  0.6391  0.484   0.3889  0.3248  0.2442
%!        1.5277  0.9351  0.7674  0.6461  0.4892  0.393   0.3282  0.2467];
%!
%! So = beta = mu = Q = dQ = mu_r = nan(numel(epsilon_r), numel(B_d_r));
%! test_freq = 4;
%! verbose = false;
%! for j=1:test_freq:numel(epsilon_r)
%! for k=1:test_freq:numel(B_d_r)
%! param.d1 = 10e-3;
%! param.Psi = 1e-3;
%! param.M = int32(20);
%! param.N = int32(85);
%! param.output_bearing_data = true;
%! param.B = param.d1 * B_d_r(k);
%! param.epsilon = epsilon_r(j);
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "oct-mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer hyd_node_id_inlet = 2003;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: integer genel_id_inlet = 4005;\n");
%!     fputs(fd, "set: real omega1z = 2 * pi * 50;\n");
%!     fputs(fd, "set: real omega2z = 2 * pi * 15;\n");
%!     fputs(fd, "set: real v1z = 0;\n");
%!     fputs(fd, "set: real v2z = 0;\n");
%!     fputs(fd, "set: real n = 1e-3;\n");
%!     fputs(fd, "set: integer K = 36000;\n");
%!     fputs(fd, "set: real t1 = abs(2 * pi * n / omega1z);\n");
%!     fputs(fd, "set: real dt = t1 / (n * K);\n");
%!     fputs(fd, "set: real D2 = d1 / (1. - Psi);\n");
%!     fputs(fd, "set: real eta = 1e-3;\n");
%!     fputs(fd, "set: real rho = 3600. * 1000.; # [l/h]\n\n");
%!     fputs(fd, "set: real epsilon_t0 = epsilon;\n");
%!     fputs(fd, "set: real epsilon_t1 = epsilon;\n");
%!     fputs(fd, "set: real delta_t0 = pi;\n");
%!     fputs(fd, "set: real delta_t1 = pi + omega2z * t1;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real p_pB2 = 0;\n");
%!     fputs(fd, "set: real p_mB2 = 0;\n");
%!     fputs(fd, "set: real p_in = 0;\n");
%!     fputs(fd, "set: real pmax = 500000.;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 50;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-4, test, norm;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e5;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: hybrid, ms, 0.;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        nonlinear solver: nox,\n");
%!     fputs(fd, "             jacobian operator, newton,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    use automatic differentiation;\n");
%!     fputs(fd, "    skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "        structural nodes: 2;\n");
%!     fputs(fd, "        joints: 2;\n");
%!     fputs(fd, "        loadable elements: 1;\n");
%!     fputs(fd, "        hydraulic nodes: 3;\n");
%!     fputs(fd, "        genels: 3;\n");
%!     fputs(fd, "        print: dof stats, to file;\n");
%!     fputs(fd, "        print: equation description, to file;\n");
%!     fputs(fd, "        print: dof description, to file;\n");
%!     fputs(fd, "        print: element connection, to file;\n");
%!     fputs(fd, "        print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 8;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        v1z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing, null,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        v2z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_inlet, p_in;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            mult, const, v1z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, v2z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!     fputs(fd, "    genel: genel_id_inlet, clamp,\n");
%!     fputs(fd, "        hyd_node_id_inlet, hydraulic, p_in;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!     fputs(fd, "                density, rho,\n");
%!     fputs(fd, "                1.,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                viscosity, eta,\n");
%!     fputs(fd, "                temperature, 0,\n");
%!     fputs(fd, "            viscosity vapor, eta,\n");
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at bearing,\n");
%!     fputs(fd, "                bearing width, B,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure coupling conditions radial, 1,\n");
%!     fputs(fd, "                position, 0.5 * D2 * (delta_t0 + pi), 0.,\n");
%!     fputs(fd, "                rectangle, width, D2 * pi / (N - 1), height, B,\n");
%!     fputs(fd, "                pressure node, hyd_node_id_inlet,\n");
%!     fputs(fd, "            pressure dof scale, pmax,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, no,\n");
%!     fputs(fd, "            output velocity, no,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   if (~opt_sol.verbose)
%!     opt_sol.logfile = [output_file, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%!   d = res.log_dat.bearings.cylindrical.d;
%!   D = res.log_dat.bearings.cylindrical.D;
%!   B = res.log_dat.bearings.cylindrical.B;
%!   s = (D - d);
%!   Psi = (D - d) / D;
%!   eta = res.log_dat.bearings.eta;
%!   node_idx_2 = find(res.node_id == res.log_dat.bearings.cylindrical.nodes(2).label);
%!   Phi2 = res.trajectory{node_idx_2}(:, 4:6).';
%!   R2 = euler123_to_rotation_matrix(Phi2);
%!   Rb1 = res.log_dat.bearings.cylindrical.nodes(1).Rb;
%!   Rb2 = res.log_dat.bearings.cylindrical.nodes(2).Rb;
%!   F1 = -res.bearings.columns.F1.';
%!   M1 = -res.bearings.columns.M1.';
%!   omega1z = res.bearings.cylindrical.omega1z;
%!   omega2z = res.bearings.cylindrical.omega2z;
%!   omega_res = res.bearings.cylindrical.omega_res;
%!   epsilon = res.bearings.cylindrical.epsilon;
%!   delta = res.bearings.cylindrical.delta;
%!   Q_out1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end) / res.log_dat.vars.rho;
%!   Q_out2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end) / res.log_dat.vars.rho;
%!   Q_in = -res.genel_data{res.log_dat.vars.genel_idx_inlet}(end) / res.log_dat.vars.rho;
%!   P = norm(Rb2(:, 1:2).' * R2.' * F1);
%!   alpha = atan2(Rb2(:, 2).' * R2.' * F1, Rb2(:, 1).' * R2.' * F1);
%!   beta(j, k) = sign(omega_res) * (delta - alpha);
%!   beta(j, k) = mod(beta(j, k) + pi, 2 * pi) - pi;
%!   mu(j, k) = 2 * Rb2(:, 3).' * R2.' * M1 / (P * d) * sign(omega1z - omega2z);
%!   So(j, k) = P * Psi^2 / (B * D * eta * abs(omega_res));
%!   Q(j, k) = (Q_out1 + Q_out2) / ((0.5 * D)^3 * Psi * abs(omega_res));
%!   dQ(j, k) = (Q_out1 + Q_out2 + Q_in) / Q_in;
%!   mu_r(j, k) = Psi * (abs((omega1z - omega2z) / omega_res) * pi / (sqrt(1 - epsilon^2) * So_r(j, k)) + sin(beta_r(j, k)) * abs(epsilon) / 2);
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
%! endfor
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   semilogy(epsilon_r(1:test_freq:end), So_r(1:test_freq:end, i), "-;So_r(epsilon);0");
%!   semilogy(epsilon_r(1:test_freq:end), So(1:test_freq:end, i), "-;So(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("So [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Sommerfeld number for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta_r(1:test_freq:end, i), "-;beta_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta(1:test_freq:end, i), "-;beta(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("beta [deg]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Angular offset of reaction force for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), mu_r(1:test_freq:end, i), "-;mu_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), mu(1:test_freq:end, i), "-;mu(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("mu [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Coefficient of friction for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), Q_r(1:test_freq:end, i), "-;Q_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), Q(1:test_freq:end, i), "-;Q(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("Q [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Non-dimensional oil flow pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! assert_simple(mean(mean(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.03);
%! assert_simple(mean(mean(abs(beta(1:test_freq:end, 1:test_freq:end) ./ beta_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.02);
%! assert_simple(mean(mean(abs(mu(1:test_freq:end, 1:test_freq:end) ./ mu_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.02);
%! assert_simple(mean(mean(abs(Q(1:test_freq:end, 1:test_freq:end) ./ Q_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.07);
%! assert_simple(max(max(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.05);
%! assert_simple(max(max(abs(beta(1:test_freq:end, 1:test_freq:end) ./ beta_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.04);
%! assert_simple(max(max(abs(mu(1:test_freq:end, 1:test_freq:end) ./ mu_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.05);
%! assert_simple(max(max(abs(Q(1:test_freq:end, 1:test_freq:end) ./ Q_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.08);

%!test
%! ## TEST 2
%! ## References:
%! ## Hans Juergen Butenschoen
%! ## Das hydrodynamische, zylindrische Gleitlager endlicher Breite unter instationaerer Belastung, Karlsruhe 1976
%! ## Test case for squeeze flow
%!
%! B_d_r = [1, 1/2, 1/2.5, 1/3, 1/4, 1/5, 1/6, 1/8];
%! epsilon_r = [-0.99, -0.98, -0.97, -0.96, -0.95, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, ...
%!                0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99];

%! So_r = [0.4925  0.1303  0.0842  0.0588  0.0333  0.0214  0.0149  0.0084
%!         0.4986  0.1319  0.0853  0.0596  0.0338  0.0217  0.0151  0.0085
%!         0.5049  0.1337  0.0864  0.0604  0.0342  0.022  0.0153  0.0086
%!         0.5113  0.1354  0.0875  0.0612  0.0347  0.0223  0.0155  0.0087
%!         0.5177  0.1372  0.0887  0.062  0.0351  0.0226  0.0157  0.0085
%!         0.5518  0.1466  0.0948  0.0663  0.0376  0.0241  0.0168  0.0095
%!         0.6293  0.1685  0.1085  0.0761  0.0432  0.0277  0.0193  0.0109
%!         0.7242  0.1944  0.126  0.0881  0.05  0.0321  0.0223  0.0126
%!         0.8342  0.2265  0.1469  0.1028  0.0583  0.0375  0.0261  0.0147
%!         0.971  0.2664  0.173  0.1211  0.0688  0.0442  0.0308  0.0174
%!         1.1391  0.3164  0.2058  0.1443  0.082  0.0527  0.0367  0.0207
%!         1.3494  0.3803  0.2479  0.1739  0.0989  0.0637  0.0444  0.025
%!         1.6134  0.4632  0.3027  0.2127  0.1212  0.078  0.0544  0.0307
%!         1.9579  0.5732  0.3757  0.2645  0.1509  0.0973  0.0678  0.0383
%!         2.4076  0.7224  0.4754  0.3355  0.1919  0.1238  0.0863  0.0488
%!         3.0122  0.9814  0.6159  0.4357  0.2499  0.1614  0.1127  0.0637
%!         3.8485  1.2236  0.8205  0.5826  0.3354  0.2171  0.1517  0.0859
%!         5.0479  1.6904  1.1329  0.8081  0.4675  0.3033  0.2122  0.1203
%!         6.8503  2.4202  1.6392  1.1756  0.6845  0.4455  0.3123  0.1774
%!         9.7319  3.675  2.5206  1.8236  1.0715  0.7005  0.4923  0.2803
%!         14.769  6.0649  4.2362  3.1002  1.8455  1.2147  0.857  0.4899
%!         24.842  11.373  8.156  6.0733  3.6891  2.4543  1.7424  1.0027
%!         50.349  26.676  15.94  15.282  9.6161  6.5227  4.6863  2.7312
%!         160.25  104.75  84.39  68.484  46.559  33.102  24.5  14.77
%!         490.18  371.98  320.59  275.76  205.15  155.23  119.99  76.315
%!         699.59  549.66  492.49  422.2  323.24  250.15  196.59  127.78
%!         1098.6  900.95  806.73  720.35  571.88  455.42  366.1  245.48
%!         2073.5  1774.3  1629.4  1490.8  1239.4  1027.2  853.42  600.66
%!         6057.9  5486.9  5171.8  4884.9  4335.3  3813.7  2331  2560.1];
%!
%! Q_r = [3.7975	1.9601	1.5761	1.3172	0.9917	0.7954	0.6638	0.4988
%!        3.7949	1.9601	1.5761	1.3172	0.9917	0.7954	0.6638	0.4988
%!        3.7927	1.9600	1.5760	1.3172	0.9917	0.7954	0.6638	0.4988
%!        3.7901	1.9600	1.5760	1.3172	0.9916	0.7954	0.6638	0.4988
%!        3.7878	1.9599	1.5759	1.3171	0.9916	0.7953	0.6638	0.4988
%!        3.7749	1.9594	1.5755	1.3171	0.9915	0.7953	0.6638	0.4988
%!        3.7498	1.9577	1.5743	1.3170	0.9914	0.7952	0.6637	0.4988
%!        3.7216	1.9549	1.5723	1.3164	0.9912	0.7951	0.6637	0.4988
%!        3.6937	1.9512	1.5706	1.3156	0.9909	0.7950	0.6636	0.4988
%!        3.6648	1.9470	1.5688	1.3144	0.9906	0.7949	0.6635	0.4988
%!        3.6338	1.9422	1.5664	1.3133	0.9903	0.7948	0.6634	0.4987
%!        3.5994	1.9369	1.5639	1.3119	0.9899	0.7947	0.6634	0.4987
%!        3.5585	1.9308	1.5611	1.3104	0.9895	0.7945	0.6633	0.4987
%!        3.5124	1.9231	1.5581	1.3090	0.9890	0.7943	0.6632	0.4987
%!        3.4640	1.9152	1.5549	1.3072	0.9885	0.7941	0.6632	0.4986
%!        3.4002	1.9063	1.5512	1.3051	0.9879	0.7939	0.6631	0.4986
%!        3.3313	1.8961	1.5476	1.3031	0.9873	0.7937	0.6630	0.4986
%!        3.2523	1.8859	1.5431	1.3010	0.9866	0.7935	0.6629	0.4985
%!        3.1603	1.8738	1.5381	1.2987	0.9860	0.7932	0.6627	0.4985
%!        3.0839	1.8603	1.5328	1.2959	0.9852	0.7929	0.6626	0.4985
%!        2.9292	1.8455	1.5266	1.2927	0.9843	0.7925	0.6625	0.4984
%!        2.8114	1.8257	1.5188	1.2889	0.9832	0.7921	0.6623	0.4984
%!        2.6815	1.8013	1.5086	1.2841	0.9819	0.7917	0.6622	0.4984
%!        2.5363	1.7723	1.4934	1.2770	0.9803	0.7911	0.6620	0.4983
%!        2.4560	1.7543	1.4843	1.2735	0.9793	0.7908	0.6619	0.4983
%!        2.4407	1.7503	1.4822	1.2725	0.9791	0.7907	0.6618	0.4983
%!        2.4244	1.7463	1.4799	1.2715	0.9788	0.7906	0.6618	0.4982
%!        2.4089	1.7416	1.4773	1.2702	0.9785	0.7905	0.6618	0.4982
%!        2.3923	1.7362	1.4746	1.2690	0.9782	0.7904	0.6617	0.4982];
%! So = beta = mu = Q = dQ = nan(numel(epsilon_r), numel(B_d_r));
%! cavitation = "mass conserving";
%! test_freq = 4;
%! verbose = false;
%! for j=1:test_freq:numel(epsilon_r)
%! for k=1:test_freq:numel(B_d_r)
%! param.d1 = 10e-3;
%! param.D2 = 10.01e-3;
%! param.M = int32(20);
%! param.N = int32(85);
%! param.output_bearing_data = true;
%! param.B = param.d1 * B_d_r(k);
%! param.epsilon = epsilon_r(j);
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "oct-mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer hyd_node_id_inlet = 2003;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: integer genel_id_inlet = 4005;\n");
%!     fputs(fd, "set: real omega1z = 2 * pi * 0;\n");
%!     fputs(fd, "set: real omega2z = 2 * pi * 0;\n");
%!     fputs(fd, "set: real v1z = 0;\n");
%!     fputs(fd, "set: real v2z = 0;\n");
%!     fputs(fd, "set: real n = 1e-3;\n");
%!     fputs(fd, "set: integer K = 36;\n");
%!     fputs(fd, "set: real t1 = abs(2 * pi * n / max(1.,abs(omega1z)));\n");
%!     fputs(fd, "set: real dt = t1 / K;\n");
%!     fputs(fd, "set: real eta = 1e-3;\n");
%!     fputs(fd, "set: real rho = 3600. * 1000.; # [l/h]\n\n");
%!     fputs(fd, "set: real epsilon_t0 = (abs(epsilon) * (1 - sign(epsilon) * n)) * (abs(epsilon) > 0) + n * (abs(epsilon) == 0);\n");
%!     fputs(fd, "set: real epsilon_t1 = abs(epsilon);\n");
%!     fputs(fd, "set: real delta_t0 = pi;\n");
%!     fputs(fd, "set: real delta_t1 = pi + (omega1z + omega2z) / 2 * t1;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real p_pB2 = 0;\n");
%!     fputs(fd, "set: real p_mB2 = 0;\n");
%!     fputs(fd, "set: real p_in = 0;\n");
%!     fputs(fd, "set: real pmax = 500000.;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 50;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-4, test, norm;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e7;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: hybrid, ms, 0.;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        nonlinear solver: nox,\n");
%!     fputs(fd, "             jacobian operator, newton,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    use automatic differentiation;\n");
%!     fputs(fd, "    skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "    structural nodes: 2;\n");
%!     fputs(fd, "    joints: 2;\n");
%!     fputs(fd, "    loadable elements: 1;\n");
%!     fputs(fd, "    hydraulic nodes: 3;\n");
%!     fputs(fd, "    genels: 3;\n");
%!     fputs(fd, "    print: dof stats, to file;\n");
%!     fputs(fd, "    print: equation description, to file;\n");
%!     fputs(fd, "    print: dof description, to file;\n");
%!     fputs(fd, "    print: element connection, to file;\n");
%!     fputs(fd, "    print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 16;\n");
%!     fputs(fd, "    default output: reference frames;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.5 * (D2 - d1) * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        0.5 * (D2 - d1) * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        v1z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing, null,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        v2z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_inlet, p_in;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            mult, const, v1z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, v2z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!     fputs(fd, "    genel: genel_id_inlet, clamp,\n");
%!     fputs(fd, "        hyd_node_id_inlet, hydraulic, p_in;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     switch (cavitation)
%!     case "non mass conserving"
%!       fputs(fd, "            hydraulic fluid, incompressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                pressure, 0.,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!     case "mass conserving"
%!       fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                1.,\n");
%!       fputs(fd, "                0.,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!       fputs(fd, "            viscosity vapor, eta,\n");
%!     otherwise
%!       error("unknown cavitation model: %s", cavitation);
%!     endswitch
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at bearing,\n");
%!     fputs(fd, "                bearing width, B,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure coupling conditions radial, 1,\n");
%!     fputs(fd, "                position, 0.5 * D2 * (delta_t0 + pi * (epsilon_dot > 0)), 0.,\n");
%!     fputs(fd, "                rectangle, width, 1.5 * D2 * pi / (N - 1), height, B,\n");
%!     fputs(fd, "                pressure node, hyd_node_id_inlet,\n");
%!     fputs(fd, "            pressure dof scale, pmax,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, yes,\n");
%!     fputs(fd, "            output velocity, yes,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   if (~opt_sol.verbose)
%!     opt_sol.logfile = [output_file, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%!   d = res.log_dat.bearings.cylindrical.d;
%!   D = res.log_dat.bearings.cylindrical.D;
%!   B = res.log_dat.bearings.cylindrical.B;
%!   s = (D - d);
%!   Psi = (D - d) / D;
%!   eta = res.log_dat.bearings.eta;
%!   node_idx_2 = find(res.node_id == res.log_dat.bearings.cylindrical.nodes(2).label);
%!   Phi2 = res.trajectory{node_idx_2}(end, 4:6).';
%!   R2 = euler123_to_rotation_matrix(Phi2);
%!   Rb1 = res.log_dat.bearings.cylindrical.nodes(1).Rb;
%!   Rb2 = res.log_dat.bearings.cylindrical.nodes(2).Rb;
%!   F1 = -res.bearings.columns.F1(end, :).';
%!   M1 = -res.bearings.columns.M1(end, :).';
%!   omega1z = res.bearings.cylindrical.omega1z(end, :);
%!   omega2z = res.bearings.cylindrical.omega2z(end, :);
%!   omega_res = res.bearings.cylindrical.omega_res(end);
%!   epsilon = res.bearings.cylindrical.epsilon(end);
%!   delta = res.bearings.cylindrical.delta(end);
%!   e1_dot = res.bearings.cylindrical.e_dot_R2(end, :);
%!   epsilon_dot = res.bearings.cylindrical.epsilon_dot(end);
%!   Q_out1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end) / res.log_dat.vars.rho;
%!   Q_out2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end) / res.log_dat.vars.rho;
%!   Q_in = -res.genel_data{res.log_dat.vars.genel_idx_inlet}(end) / res.log_dat.vars.rho;
%!   P = norm(Rb2(:, 1:2).' * R2.' * F1);
%!   alpha = atan2(Rb2(:, 2).' * R2.' * F1, Rb2(:, 1).' * R2.' * F1);
%!   beta(j, k) = sign(omega_res) * (delta - alpha);
%!   beta(j, k) = mod(beta(j, k) + pi, 2 * pi) - pi;
%!   mu(j, k) = 2 * Rb2(:, 3).' * R2.' * M1 / (P * d) * sign(omega1z - omega2z);
%!   So(j, k) = P * Psi^2 / (B * D * eta * abs(epsilon_dot));
%!   Q(j, k) = (Q_out1 + Q_out2) / ((0.5 * D)^3 * Psi * abs(epsilon_dot));
%!   dQ(j, k) = (Q_out1 + Q_out2 + Q_in) / Q_in;
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
%! endfor
%! endfor
%! for i=1:columns(beta)
%!   beta(beta(:, i) < -0.9 * pi, i) += 2 * pi;
%! endfor
%! mu_r = zeros(numel(epsilon_r), numel(B_d_r));
%! beta_r = zeros(numel(epsilon_r), numel(B_d_r));
%! beta_r(find(epsilon_r < 0), :) = pi;
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   semilogy(epsilon_r(1:test_freq:end), So_r(1:test_freq:end, i), "-;So_r(epsilon);0");
%!   semilogy(epsilon_r(1:test_freq:end), So(1:test_freq:end, i), "-;So(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("So [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Sommerfeld number for pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta_r(1:test_freq:end, i), "-;beta_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta(1:test_freq:end, i), "-;beta(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("beta [deg]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Angular offset of reaction force for pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), mu_r(1:test_freq:end, i), "-;mu_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), mu(1:test_freq:end, i), "-;mu(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("mu [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Coefficient of friction for pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), Q_r(1:test_freq:end, i), "-;Q_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), Q(1:test_freq:end, i), "-;Q(epsilon);1");
%!   ylim([0, ylim()(2)]);
%!   xlabel("epsilon [1]");
%!   ylabel("Q [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Non-dimensional oil flow pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! assert_simple(mean(mean(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.02);
%! assert_simple(mean(mean(abs(beta(1:test_freq:end,1:test_freq:end) - beta_r(1:test_freq:end, 1:test_freq:end)))) < 1e-6 * pi / 180);
%! assert_simple(mean(mean(abs(mu(1:test_freq:end,1:test_freq:end)))) < 1e-8);
%! assert_simple(mean(mean(abs(Q(1:test_freq:end,1:test_freq:end) ./ Q_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.06);
%! assert_simple(max(max(abs(So(1:test_freq:end,1:test_freq:end) ./ So_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.11);
%! assert_simple(max(max(abs(beta(1:test_freq:end,1:test_freq:end) - beta_r(1:test_freq:end,1:test_freq:end)))) < 1e-6 * pi / 180);
%! assert_simple(max(max(abs(mu(1:test_freq:end,1:test_freq:end)))) < 1e-8);
%! assert_simple(max(max(abs(Q(1:test_freq:end,1:test_freq:end) ./ Q_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.07);

%!test
%! ## TEST 3
%! ## References:
%! ## Hans Juergen Butenschoen
%! ## Das hydrodynamische, zylindrische Gleitlager endlicher Breite unter instationaerer Belastung, Karlsruhe 1976
%! ## Test case for rotation
%! B_d_r = [1, 1/2, 1/2.5, 1/3, 1/4, 1/5, 1/6, 1/8];
%!
%! epsilon_r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99];
%!
%! So_r =  [0.1196 0.0368  0.0243  0.0171  0.0098  0.0063  0.0044  0.0025
%!          0.2518 0.0783  0.0518  0.0366  0.0210  0.0136  0.0095  0.0054
%!          0.4091 0.1304  0.0867  0.0615  0.0354  0.0229  0.0160  0.0091
%!          0.6108 0.2026  0.1357  0.0968  0.0560  0.0364  0.0254  0.0144
%!          0.8903 0.3124  0.2117  0.1522  0.0888  0.0579  0.0406  0.0231
%!          1.3146 0.4982  0.3435  0.2496  0.1476  0.0969  0.0682  0.0390
%!          2.0432 0.8595  0.6075  0.4492  0.2708  0.1797  0.1274  0.0732
%!          3.5663 1.7339  1.2756  0.9687  0.6043  0.4085  0.2930  0.1706
%!          8.4352 5.0881  4.0187  3.2201  2.1595  1.5261  1.1263  0.6776
%!          18.7895 13.3083 11.2225 9.4993  6.9248  5.1833  3.9831  2.5202
%!          24.1172 17.7934 15.2823 13.1525 9.8517  7.5257  5.8712  3.7907
%!          33.1297 25.5920 22.4503 19.7034 15.2697 11.9804 9.5425  6.3397
%!          51.4774 41.9230 37.7412 33.9523 27.5040 22.3556 18.3874 12.7695
%!          107.7868 93.7881 87.2906 81.1597 70.0359 60.3874 52.1425 39.2568];
%!
%! beta_r = [79.410 81.767  82.119  82.299  82.481  82.561  82.608  82.653
%!           73.858 75.142  75.283  75.344  75.287  75.406  75.414  75.422
%!           68.264 68.493  68.423  68.368  68.204  68.261  68.238  68.211
%!           62.571 61.778  61.544  61.382  61.208  61.115  61.062  61.007
%!           56.705 54.993  54.603  54.348  54.069  53.922  53.863  53.784
%!           50.536 48.049  47.521  47.191  46.825  46.647  46.554  46.449
%!           43.859 40.803  40.156  39.756  39.326  39.108  38.983  38.865
%!           36.235 32.938  32.216  31.761  31.249  30.988  30.840  30.692
%!           26.482 23.566  22.849  22.368  21.790  21.476  21.289  21.089
%!           19.450 17.265  16.648  16.207  15.632  15.291  15.075  14.832
%!           17.609 15.660  15.089  14.669  14.109  13.768  13.547  13.262
%!           15.484 13.818  13.313  12.929  12.299  12.062  11.838  11.570
%!           12.903 11.598  11.178  10.850  10.375  10.057  9.835   9.558
%!            9.416  8.587   8.301   8.066   7.703   7.440   7.242   6.975] * pi / 180;
%!
%! Q_r = [0.1603  0.0939  0.0768  0.0648  0.0492  0.0396  0.0331  0.0249
%!        0.3196  0.1878  0.1536  0.1296  0.0984  0.0792  0.0662  0.0498
%!        0.4765  0.2816  0.2304  0.1944  0.1476  0.1188  0.0993  0.0747
%!        0.6318  0.3755  0.3072  0.2592  0.1968  0.1583  0.1324  0.0996
%!        0.7852  0.4694  0.384   0.324   0.246   0.198   0.1655  0.1245
%!        0.9374  0.5634  0.461   0.3889  0.2953  0.2376  0.1986  0.1494
%!        1.0888  0.6578  0.5383  0.454   0.3446  0.2772  0.2317  0.1743
%!        1.2392  0.7529  0.6159  0.5193  0.394   0.3169  0.2648  0.1992
%!        1.3898  0.8453  0.6944  0.5852  0.4436  0.3567  0.298   0.2241
%!        1.4655  0.8984  0.7343  0.6186  0.4687  0.3767  0.3147  0.2366
%!        1.4808  0.9084  0.7425  0.6254  0.4737  0.3807  0.3181  0.2391
%!        1.4968  0.9185  0.7507  0.6322  0.4788  0.3848  0.3214  0.2417
%!        1.5124  0.9287  0.7585  0.6391  0.484   0.3889  0.3248  0.2442
%!        1.5277  0.9351  0.7674  0.6461  0.4892  0.393   0.3282  0.2467];
%!
%! So = beta = mu = Q = dQ = mu_r = nan(numel(epsilon_r), numel(B_d_r));
%! test_freq = 4;
%! verbose = false;
%! for j=1:test_freq:numel(epsilon_r)
%! for k=1:test_freq:numel(B_d_r)
%! param.d1 = 10e-3;
%! param.Psi = 1e-3;
%! param.M = int32(20);
%! param.N = int32(85);
%! param.output_bearing_data = true;
%! param.B = param.d1 * B_d_r(k);
%! param.epsilon = epsilon_r(j);
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "oct-mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer hyd_node_id_inlet = 2003;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: integer genel_id_inlet = 4005;\n");
%!     fputs(fd, "set: real omega1z = 2 * pi * 50;\n");
%!     fputs(fd, "set: real omega2z = 2 * pi * 15;\n");
%!     fputs(fd, "set: real v1z = 0;\n");
%!     fputs(fd, "set: real v2z = 0;\n");
%!     fputs(fd, "set: real n = 1e-3;\n");
%!     fputs(fd, "set: integer K = 36000;\n");
%!     fputs(fd, "set: real t1 = abs(2 * pi * n / omega1z);\n");
%!     fputs(fd, "set: real dt = t1 / (n * K);\n");
%!     fputs(fd, "set: real D2 = d1 / (1. - Psi);\n");
%!     fputs(fd, "set: real eta = 1e-3;\n");
%!     fputs(fd, "set: real rho = 3600. * 1000.; # [l/h]\n\n");
%!     fputs(fd, "set: real epsilon_t0 = epsilon;\n");
%!     fputs(fd, "set: real epsilon_t1 = epsilon;\n");
%!     fputs(fd, "set: real delta_t0 = pi;\n");
%!     fputs(fd, "set: real delta_t1 = pi + omega2z * t1;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real p_pB2 = 0;\n");
%!     fputs(fd, "set: real p_mB2 = 0;\n");
%!     fputs(fd, "set: real p_in = 0;\n");
%!     fputs(fd, "set: real pmax = 500000.;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 50;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-4, test, norm;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e5;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: hybrid, ms, 0.;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4, 1e-3;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    use automatic differentiation;\n");
%!     fputs(fd, "    skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "        structural nodes: 2;\n");
%!     fputs(fd, "        joints: 2;\n");
%!     fputs(fd, "        loadable elements: 1;\n");
%!     fputs(fd, "        hydraulic nodes: 3;\n");
%!     fputs(fd, "        genels: 3;\n");
%!     fputs(fd, "        print: dof stats, to file;\n");
%!     fputs(fd, "        print: equation description, to file;\n");
%!     fputs(fd, "        print: dof description, to file;\n");
%!     fputs(fd, "        print: element connection, to file;\n");
%!     fputs(fd, "        print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 8;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        v1z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing, null,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        v2z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_inlet, p_in;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            mult, const, v1z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, v2z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!     fputs(fd, "    genel: genel_id_inlet, clamp,\n");
%!     fputs(fd, "        hyd_node_id_inlet, hydraulic, p_in;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!     fputs(fd, "                density, rho,\n");
%!     fputs(fd, "                1.,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                viscosity, eta,\n");
%!     fputs(fd, "                temperature, 0,\n");
%!     fputs(fd, "            viscosity vapor, eta,\n");
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at bearing,\n");
%!     fputs(fd, "                bearing width, B,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure coupling conditions radial, 1,\n");
%!     fputs(fd, "                position, 0.5 * D2 * (delta_t0 + pi), 0.,\n");
%!     fputs(fd, "                rectangle, width, D2 * pi / (N - 1), height, B,\n");
%!     fputs(fd, "                pressure node, hyd_node_id_inlet,\n");
%!     fputs(fd, "            pressure dof scale, pmax,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, no,\n");
%!     fputs(fd, "            output velocity, no,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   if (~opt_sol.verbose)
%!     opt_sol.logfile = [output_file, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%!   d = res.log_dat.bearings.cylindrical.d;
%!   D = res.log_dat.bearings.cylindrical.D;
%!   B = res.log_dat.bearings.cylindrical.B;
%!   s = (D - d);
%!   Psi = (D - d) / D;
%!   eta = res.log_dat.bearings.eta;
%!   node_idx_2 = find(res.node_id == res.log_dat.bearings.cylindrical.nodes(2).label);
%!   Phi2 = res.trajectory{node_idx_2}(:, 4:6).';
%!   R2 = euler123_to_rotation_matrix(Phi2);
%!   Rb1 = res.log_dat.bearings.cylindrical.nodes(1).Rb;
%!   Rb2 = res.log_dat.bearings.cylindrical.nodes(2).Rb;
%!   F1 = -res.bearings.columns.F1.';
%!   M1 = -res.bearings.columns.M1.';
%!   omega1z = res.bearings.cylindrical.omega1z;
%!   omega2z = res.bearings.cylindrical.omega2z;
%!   omega_res = res.bearings.cylindrical.omega_res;
%!   epsilon = res.bearings.cylindrical.epsilon;
%!   delta = res.bearings.cylindrical.delta;
%!   Q_out1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end) / res.log_dat.vars.rho;
%!   Q_out2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end) / res.log_dat.vars.rho;
%!   Q_in = -res.genel_data{res.log_dat.vars.genel_idx_inlet}(end) / res.log_dat.vars.rho;
%!   P = norm(Rb2(:, 1:2).' * R2.' * F1);
%!   alpha = atan2(Rb2(:, 2).' * R2.' * F1, Rb2(:, 1).' * R2.' * F1);
%!   beta(j, k) = sign(omega_res) * (delta - alpha);
%!   beta(j, k) = mod(beta(j, k) + pi, 2 * pi) - pi;
%!   mu(j, k) = 2 * Rb2(:, 3).' * R2.' * M1 / (P * d) * sign(omega1z - omega2z);
%!   So(j, k) = P * Psi^2 / (B * D * eta * abs(omega_res));
%!   Q(j, k) = (Q_out1 + Q_out2) / ((0.5 * D)^3 * Psi * abs(omega_res));
%!   dQ(j, k) = (Q_out1 + Q_out2 + Q_in) / Q_in;
%!   mu_r(j, k) = Psi * (abs((omega1z - omega2z) / omega_res) * pi / (sqrt(1 - epsilon^2) * So_r(j, k)) + sin(beta_r(j, k)) * abs(epsilon) / 2);
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
%! endfor
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   semilogy(epsilon_r(1:test_freq:end), So_r(1:test_freq:end, i), "-;So_r(epsilon);0");
%!   semilogy(epsilon_r(1:test_freq:end), So(1:test_freq:end, i), "-;So(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("So [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Sommerfeld number for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta_r(1:test_freq:end, i), "-;beta_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta(1:test_freq:end, i), "-;beta(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("beta [deg]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Angular offset of reaction force for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), mu_r(1:test_freq:end, i), "-;mu_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), mu(1:test_freq:end, i), "-;mu(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("mu [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Coefficient of friction for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), Q_r(1:test_freq:end, i), "-;Q_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), Q(1:test_freq:end, i), "-;Q(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("Q [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Non-dimensional oil flow pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! assert_simple(mean(mean(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.03);
%! assert_simple(mean(mean(abs(beta(1:test_freq:end, 1:test_freq:end) ./ beta_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.02);
%! assert_simple(mean(mean(abs(mu(1:test_freq:end, 1:test_freq:end) ./ mu_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.02);
%! assert_simple(mean(mean(abs(Q(1:test_freq:end, 1:test_freq:end) ./ Q_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.07);
%! assert_simple(max(max(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.05);
%! assert_simple(max(max(abs(beta(1:test_freq:end, 1:test_freq:end) ./ beta_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.04);
%! assert_simple(max(max(abs(mu(1:test_freq:end, 1:test_freq:end) ./ mu_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.05);
%! assert_simple(max(max(abs(Q(1:test_freq:end, 1:test_freq:end) ./ Q_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.08);

%!test
%! ## TEST 4
%! ## References:
%! ## Hans Juergen Butenschoen
%! ## Das hydrodynamische, zylindrische Gleitlager endlicher Breite unter instationaerer Belastung, Karlsruhe 1976
%! ## Test case for squeeze flow
%!
%! B_d_r = [1, 1/2, 1/2.5, 1/3, 1/4, 1/5, 1/6, 1/8];
%! epsilon_r = [-0.99, -0.98, -0.97, -0.96, -0.95, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, ...
%!                0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99];

%! So_r = [0.4925  0.1303  0.0842  0.0588  0.0333  0.0214  0.0149  0.0084
%!         0.4986  0.1319  0.0853  0.0596  0.0338  0.0217  0.0151  0.0085
%!         0.5049  0.1337  0.0864  0.0604  0.0342  0.022  0.0153  0.0086
%!         0.5113  0.1354  0.0875  0.0612  0.0347  0.0223  0.0155  0.0087
%!         0.5177  0.1372  0.0887  0.062  0.0351  0.0226  0.0157  0.0085
%!         0.5518  0.1466  0.0948  0.0663  0.0376  0.0241  0.0168  0.0095
%!         0.6293  0.1685  0.1085  0.0761  0.0432  0.0277  0.0193  0.0109
%!         0.7242  0.1944  0.126  0.0881  0.05  0.0321  0.0223  0.0126
%!         0.8342  0.2265  0.1469  0.1028  0.0583  0.0375  0.0261  0.0147
%!         0.971  0.2664  0.173  0.1211  0.0688  0.0442  0.0308  0.0174
%!         1.1391  0.3164  0.2058  0.1443  0.082  0.0527  0.0367  0.0207
%!         1.3494  0.3803  0.2479  0.1739  0.0989  0.0637  0.0444  0.025
%!         1.6134  0.4632  0.3027  0.2127  0.1212  0.078  0.0544  0.0307
%!         1.9579  0.5732  0.3757  0.2645  0.1509  0.0973  0.0678  0.0383
%!         2.4076  0.7224  0.4754  0.3355  0.1919  0.1238  0.0863  0.0488
%!         3.0122  0.9814  0.6159  0.4357  0.2499  0.1614  0.1127  0.0637
%!         3.8485  1.2236  0.8205  0.5826  0.3354  0.2171  0.1517  0.0859
%!         5.0479  1.6904  1.1329  0.8081  0.4675  0.3033  0.2122  0.1203
%!         6.8503  2.4202  1.6392  1.1756  0.6845  0.4455  0.3123  0.1774
%!         9.7319  3.675  2.5206  1.8236  1.0715  0.7005  0.4923  0.2803
%!         14.769  6.0649  4.2362  3.1002  1.8455  1.2147  0.857  0.4899
%!         24.842  11.373  8.156  6.0733  3.6891  2.4543  1.7424  1.0027
%!         50.349  26.676  15.94  15.282  9.6161  6.5227  4.6863  2.7312
%!         160.25  104.75  84.39  68.484  46.559  33.102  24.5  14.77
%!         490.18  371.98  320.59  275.76  205.15  155.23  119.99  76.315
%!         699.59  549.66  492.49  422.2  323.24  250.15  196.59  127.78
%!         1098.6  900.95  806.73  720.35  571.88  455.42  366.1  245.48
%!         2073.5  1774.3  1629.4  1490.8  1239.4  1027.2  853.42  600.66
%!         6057.9  5486.9  5171.8  4884.9  4335.3  3813.7  2331  2560.1];
%!
%! Q_r = [3.7975	1.9601	1.5761	1.3172	0.9917	0.7954	0.6638	0.4988
%!        3.7949	1.9601	1.5761	1.3172	0.9917	0.7954	0.6638	0.4988
%!        3.7927	1.9600	1.5760	1.3172	0.9917	0.7954	0.6638	0.4988
%!        3.7901	1.9600	1.5760	1.3172	0.9916	0.7954	0.6638	0.4988
%!        3.7878	1.9599	1.5759	1.3171	0.9916	0.7953	0.6638	0.4988
%!        3.7749	1.9594	1.5755	1.3171	0.9915	0.7953	0.6638	0.4988
%!        3.7498	1.9577	1.5743	1.3170	0.9914	0.7952	0.6637	0.4988
%!        3.7216	1.9549	1.5723	1.3164	0.9912	0.7951	0.6637	0.4988
%!        3.6937	1.9512	1.5706	1.3156	0.9909	0.7950	0.6636	0.4988
%!        3.6648	1.9470	1.5688	1.3144	0.9906	0.7949	0.6635	0.4988
%!        3.6338	1.9422	1.5664	1.3133	0.9903	0.7948	0.6634	0.4987
%!        3.5994	1.9369	1.5639	1.3119	0.9899	0.7947	0.6634	0.4987
%!        3.5585	1.9308	1.5611	1.3104	0.9895	0.7945	0.6633	0.4987
%!        3.5124	1.9231	1.5581	1.3090	0.9890	0.7943	0.6632	0.4987
%!        3.4640	1.9152	1.5549	1.3072	0.9885	0.7941	0.6632	0.4986
%!        3.4002	1.9063	1.5512	1.3051	0.9879	0.7939	0.6631	0.4986
%!        3.3313	1.8961	1.5476	1.3031	0.9873	0.7937	0.6630	0.4986
%!        3.2523	1.8859	1.5431	1.3010	0.9866	0.7935	0.6629	0.4985
%!        3.1603	1.8738	1.5381	1.2987	0.9860	0.7932	0.6627	0.4985
%!        3.0839	1.8603	1.5328	1.2959	0.9852	0.7929	0.6626	0.4985
%!        2.9292	1.8455	1.5266	1.2927	0.9843	0.7925	0.6625	0.4984
%!        2.8114	1.8257	1.5188	1.2889	0.9832	0.7921	0.6623	0.4984
%!        2.6815	1.8013	1.5086	1.2841	0.9819	0.7917	0.6622	0.4984
%!        2.5363	1.7723	1.4934	1.2770	0.9803	0.7911	0.6620	0.4983
%!        2.4560	1.7543	1.4843	1.2735	0.9793	0.7908	0.6619	0.4983
%!        2.4407	1.7503	1.4822	1.2725	0.9791	0.7907	0.6618	0.4983
%!        2.4244	1.7463	1.4799	1.2715	0.9788	0.7906	0.6618	0.4982
%!        2.4089	1.7416	1.4773	1.2702	0.9785	0.7905	0.6618	0.4982
%!        2.3923	1.7362	1.4746	1.2690	0.9782	0.7904	0.6617	0.4982];
%! So = beta = mu = Q = dQ = nan(numel(epsilon_r), numel(B_d_r));
%! cavitation = "mass conserving";
%! test_freq = 4;
%! verbose = false;
%! for j=1:test_freq:numel(epsilon_r)
%! for k=1:test_freq:numel(B_d_r)
%! param.d1 = 10e-3;
%! param.D2 = 10.01e-3;
%! param.M = int32(20);
%! param.N = int32(85);
%! param.output_bearing_data = true;
%! param.B = param.d1 * B_d_r(k);
%! param.epsilon = epsilon_r(j);
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "oct-mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer hyd_node_id_inlet = 2003;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: integer genel_id_inlet = 4005;\n");
%!     fputs(fd, "set: real omega1z = 2 * pi * 0;\n");
%!     fputs(fd, "set: real omega2z = 2 * pi * 0;\n");
%!     fputs(fd, "set: real v1z = 0;\n");
%!     fputs(fd, "set: real v2z = 0;\n");
%!     fputs(fd, "set: real n = 1e-3;\n");
%!     fputs(fd, "set: integer K = 36;\n");
%!     fputs(fd, "set: real t1 = abs(2 * pi * n / max(1.,abs(omega1z)));\n");
%!     fputs(fd, "set: real dt = t1 / K;\n");
%!     fputs(fd, "set: real eta = 1e-3;\n");
%!     fputs(fd, "set: real rho = 3600. * 1000.; # [l/h]\n\n");
%!     fputs(fd, "set: real epsilon_t0 = (abs(epsilon) * (1 - sign(epsilon) * n)) * (abs(epsilon) > 0) + n * (abs(epsilon) == 0);\n");
%!     fputs(fd, "set: real epsilon_t1 = abs(epsilon);\n");
%!     fputs(fd, "set: real delta_t0 = pi;\n");
%!     fputs(fd, "set: real delta_t1 = pi + (omega1z + omega2z) / 2 * t1;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real p_pB2 = 0;\n");
%!     fputs(fd, "set: real p_mB2 = 0;\n");
%!     fputs(fd, "set: real p_in = 0;\n");
%!     fputs(fd, "set: real pmax = 500000.;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 50;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-4, test, norm;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e7;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: hybrid, ms, 0.;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    use automatic differentiation;\n");
%!     fputs(fd, "    skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "    structural nodes: 2;\n");
%!     fputs(fd, "    joints: 2;\n");
%!     fputs(fd, "    loadable elements: 1;\n");
%!     fputs(fd, "    hydraulic nodes: 3;\n");
%!     fputs(fd, "    genels: 3;\n");
%!     fputs(fd, "    print: dof stats, to file;\n");
%!     fputs(fd, "    print: equation description, to file;\n");
%!     fputs(fd, "    print: dof description, to file;\n");
%!     fputs(fd, "    print: element connection, to file;\n");
%!     fputs(fd, "    print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 16;\n");
%!     fputs(fd, "    default output: reference frames;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.5 * (D2 - d1) * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        0.5 * (D2 - d1) * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        v1z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing, null,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        v2z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_inlet, p_in;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            mult, const, v1z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, v2z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!     fputs(fd, "    genel: genel_id_inlet, clamp,\n");
%!     fputs(fd, "        hyd_node_id_inlet, hydraulic, p_in;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     switch (cavitation)
%!     case "non mass conserving"
%!       fputs(fd, "            hydraulic fluid, incompressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                pressure, 0.,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!     case "mass conserving"
%!       fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                1.,\n");
%!       fputs(fd, "                0.,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!       fputs(fd, "            viscosity vapor, eta,\n");
%!     otherwise
%!       error("unknown cavitation model: %s", cavitation);
%!     endswitch
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at bearing,\n");
%!     fputs(fd, "                bearing width, B,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure coupling conditions radial, 1,\n");
%!     fputs(fd, "                position, 0.5 * D2 * (delta_t0 + pi * (epsilon_dot > 0)), 0.,\n");
%!     fputs(fd, "                rectangle, width, 1.5 * D2 * pi / (N - 1), height, B,\n");
%!     fputs(fd, "                pressure node, hyd_node_id_inlet,\n");
%!     fputs(fd, "            pressure dof scale, pmax,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, yes,\n");
%!     fputs(fd, "            output velocity, yes,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   if (~opt_sol.verbose)
%!     opt_sol.logfile = [output_file, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%!   d = res.log_dat.bearings.cylindrical.d;
%!   D = res.log_dat.bearings.cylindrical.D;
%!   B = res.log_dat.bearings.cylindrical.B;
%!   s = (D - d);
%!   Psi = (D - d) / D;
%!   eta = res.log_dat.bearings.eta;
%!   node_idx_2 = find(res.node_id == res.log_dat.bearings.cylindrical.nodes(2).label);
%!   Phi2 = res.trajectory{node_idx_2}(end, 4:6).';
%!   R2 = euler123_to_rotation_matrix(Phi2);
%!   Rb1 = res.log_dat.bearings.cylindrical.nodes(1).Rb;
%!   Rb2 = res.log_dat.bearings.cylindrical.nodes(2).Rb;
%!   F1 = -res.bearings.columns.F1(end, :).';
%!   M1 = -res.bearings.columns.M1(end, :).';
%!   omega1z = res.bearings.cylindrical.omega1z(end, :);
%!   omega2z = res.bearings.cylindrical.omega2z(end, :);
%!   omega_res = res.bearings.cylindrical.omega_res(end);
%!   epsilon = res.bearings.cylindrical.epsilon(end);
%!   delta = res.bearings.cylindrical.delta(end);
%!   e1_dot = res.bearings.cylindrical.e_dot_R2(end, :);
%!   epsilon_dot = res.bearings.cylindrical.epsilon_dot(end);
%!   Q_out1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end) / res.log_dat.vars.rho;
%!   Q_out2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end) / res.log_dat.vars.rho;
%!   Q_in = -res.genel_data{res.log_dat.vars.genel_idx_inlet}(end) / res.log_dat.vars.rho;
%!   P = norm(Rb2(:, 1:2).' * R2.' * F1);
%!   alpha = atan2(Rb2(:, 2).' * R2.' * F1, Rb2(:, 1).' * R2.' * F1);
%!   beta(j, k) = sign(omega_res) * (delta - alpha);
%!   beta(j, k) = mod(beta(j, k) + pi, 2 * pi) - pi;
%!   mu(j, k) = 2 * Rb2(:, 3).' * R2.' * M1 / (P * d) * sign(omega1z - omega2z);
%!   So(j, k) = P * Psi^2 / (B * D * eta * abs(epsilon_dot));
%!   Q(j, k) = (Q_out1 + Q_out2) / ((0.5 * D)^3 * Psi * abs(epsilon_dot));
%!   dQ(j, k) = (Q_out1 + Q_out2 + Q_in) / Q_in;
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
%! endfor
%! endfor
%! for i=1:columns(beta)
%!   beta(beta(:, i) < -0.9 * pi, i) += 2 * pi;
%! endfor
%! mu_r = zeros(numel(epsilon_r), numel(B_d_r));
%! beta_r = zeros(numel(epsilon_r), numel(B_d_r));
%! beta_r(find(epsilon_r < 0), :) = pi;
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   semilogy(epsilon_r(1:test_freq:end), So_r(1:test_freq:end, i), "-;So_r(epsilon);0");
%!   semilogy(epsilon_r(1:test_freq:end), So(1:test_freq:end, i), "-;So(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("So [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Sommerfeld number for pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta_r(1:test_freq:end, i), "-;beta_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta(1:test_freq:end, i), "-;beta(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("beta [deg]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Angular offset of reaction force for pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), mu_r(1:test_freq:end, i), "-;mu_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), mu(1:test_freq:end, i), "-;mu(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("mu [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Coefficient of friction for pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), Q_r(1:test_freq:end, i), "-;Q_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), Q(1:test_freq:end, i), "-;Q(epsilon);1");
%!   ylim([0, ylim()(2)]);
%!   xlabel("epsilon [1]");
%!   ylabel("Q [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Non-dimensional oil flow pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! assert_simple(mean(mean(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.02);
%! assert_simple(mean(mean(abs(beta(1:test_freq:end,1:test_freq:end) - beta_r(1:test_freq:end, 1:test_freq:end)))) < 1e-6 * pi / 180);
%! assert_simple(mean(mean(abs(mu(1:test_freq:end,1:test_freq:end)))) < 1e-8);
%! assert_simple(mean(mean(abs(Q(1:test_freq:end,1:test_freq:end) ./ Q_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.06);
%! assert_simple(max(max(abs(So(1:test_freq:end,1:test_freq:end) ./ So_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.11);
%! assert_simple(max(max(abs(beta(1:test_freq:end,1:test_freq:end) - beta_r(1:test_freq:end,1:test_freq:end)))) < 1e-6 * pi / 180);
%! assert_simple(max(max(abs(mu(1:test_freq:end,1:test_freq:end)))) < 1e-8);
%! assert_simple(max(max(abs(Q(1:test_freq:end,1:test_freq:end) ./ Q_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.07);

%!test
%! ## TEST 5
%! function [p, mdotz] = reference_solution(h0, B, dh_dz, p0, p1, eta, U1z, U2z, rho, dx, dh_dt, z)
%!    p = "(dh_dz*h0^2*p1*z^2+2*B*dh_dz^2*h0*p1*z^2+B^2*dh_dz^3*p1*z^2-dh_dz*h0^2*p0*z^2+6*B*U2z*dh_dz*eta*z^2+6*B*U1z*dh_dz*eta*z^2+12*B*dh_dt*eta*z^2+2*h0^3*p1*z+4*B*dh_dz*h0^2*p1*z+2*B^2*dh_dz^2*h0*p1*z-2*h0^3*p0*z-6*B^2*U2z*dh_dz*eta*z-6*B^2*U1z*dh_dz*eta*z-12*B^2*dh_dt*eta*z+2*B*h0^3*p0+B^2*dh_dz*h0^2*p0)/(B*(2*h0+B*dh_dz)*(dh_dz*z+h0)^2)";
%!    mdotz = "-((12*B*dh_dt*eta*h0+6*B^2*dh_dt*dh_dz*eta)*rho*z+((h0^4+2*B*dh_dz*h0^3+B^2*dh_dz^2*h0^2)*p1+(-h0^4-2*B*dh_dz*h0^3-B^2*dh_dz^2*h0^2)*p0+(-6*B*U2z-6*B*U1z)*eta*h0^2+((-6*B^2*U2z-6*B^2*U1z)*dh_dz-6*B^2*dh_dt)*eta*h0)*rho)/(12*B*eta*h0+6*B^2*dh_dz*eta)";
%!    p = eval(vectorize(p));
%!    mdotz = eval(mdotz) * dx;
%! endfunction
%! param.d1 = 10e-3;
%! param.h = 1e-5;
%! param.M = int32(750);
%! param.N = int32(4);
%! param.output_bearing_data = true;
%! param.B = 25e-3;
%! cavitation = "mass conserving";
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "oct-mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: real omega1z = 0;\n");
%!     fputs(fd, "set: real omega2z = 0;\n");
%!     fputs(fd, "set: real v1z = 15;\n");
%!     fputs(fd, "set: real v2z = -7;\n");
%!     fputs(fd, "set: integer K = 36;\n");
%!     fputs(fd, "set: real t1 = 1e-2;\n");
%!     fputs(fd, "set: real dt = t1 / K;\n");
%!     fputs(fd, "set: real D2 = d1 + 2 * h;\n");
%!     fputs(fd, "set: real dy = 1.5 * h;\n");
%!     fputs(fd, "set: real eta = 1e-3;\n");
%!     fputs(fd, "set: real rho = 3600. * 1000.; # [l/h]\n\n");
%!     fputs(fd, "set: real epsilon_t0 = 0;\n");
%!     fputs(fd, "set: real epsilon_t1 = 0;\n");
%!     fputs(fd, "set: real delta_t0 = 0;\n");
%!     fputs(fd, "set: real delta_t1 = 0;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real p_pB2 = 10e5;\n");
%!     fputs(fd, "set: real p_mB2 = 0.8e5;\n");
%!     fputs(fd, "set: real pmax = 500000.;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 50;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-4, test, norm;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e5;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: hybrid, ms, 0.;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        nonlinear solver: nox,\n");
%!     fputs(fd, "             jacobian operator, newton,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             forcing term min tolerance, 1e-10,\n");
%!     fputs(fd, "             forcing term max tolerance, 1e-3,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 200,\n");
%!     fputs(fd, "             krylov subspace size, 200;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    use automatic differentiation;\n");
%!     fputs(fd, "    skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "        structural nodes: 2;\n");
%!     fputs(fd, "        joints: 2;\n");
%!     fputs(fd, "        loadable elements: 1;\n");
%!     fputs(fd, "        hydraulic nodes: 2;\n");
%!     fputs(fd, "        genels: 2;\n");
%!     fputs(fd, "        print: dof stats, to file;\n");
%!     fputs(fd, "        print: equation description, to file;\n");
%!     fputs(fd, "        print: dof description, to file;\n");
%!     fputs(fd, "        print: element connection, to file;\n");
%!     fputs(fd, "        print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 8;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        -v1z * t1,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        v1z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        -v2z * t1,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        v2z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            array, 2, -v1z * t1,\n");
%!     fputs(fd, "                      mult, const, v1z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            array, 2, -v2z * t1, mult, const, v2z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     switch (cavitation)
%!     case "non mass conserving"
%!       fputs(fd, "            hydraulic fluid, incompressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                pressure, 0.,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!     case "mass conserving"
%!       fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                1.,\n");
%!       fputs(fd, "                0.,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!       fputs(fd, "            viscosity vapor, eta,\n");
%!     otherwise
%!       error("unknown cavitation model: %s", cavitation);
%!     endswitch
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at bearing,\n");
%!     fputs(fd, "                bearing width, B,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                pockets,\n");
%!     fputs(fd, "                  shaft, 1,\n");
%!     fputs(fd, "                    position, d1 * pi * 0.5, 0.,\n");
%!     fputs(fd, "                    complete surface,\n");
%!     fputs(fd, "                    pocket height,\n");
%!     fputs(fd, "                    linear,\n");
%!     fputs(fd, "                      x, 0., d1 * pi,\n");
%!     fputs(fd, "                      z, -0.5 * B, 0.5 * B,\n");
%!     fputs(fd, "                      delta y, 0., -dy,\n");
%!     fputs(fd, "                               0., -dy,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure dof scale, pmax,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, no,\n");
%!     fputs(fd, "            output velocity, no,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   if (~opt_sol.verbose)
%!     opt_sol.logfile = [output_file, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%!   mdot1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end);
%!   mdot2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end);
%!   d = res.log_dat.bearings.cylindrical.d;
%!   D = res.log_dat.bearings.cylindrical.D;
%!   B = res.log_dat.bearings.cylindrical.B;
%!   h0 = 0.5 * (D - d);
%!   dh_dz = res.log_dat.vars.dy / res.log_dat.vars.B;
%!   p0 = res.log_dat.vars.p_mB2;
%!   p1 = res.log_dat.vars.p_pB2;
%!   eta = res.log_dat.vars.eta;
%!   U1z = res.log_dat.vars.v1z - res.log_dat.vars.v2z;
%!   U2z = 0;
%!   rho = res.log_dat.vars.rho;
%!   dx = res.log_dat.bearings.cylindrical.dm * pi;
%!   dh_dt = -dh_dz * (U1z - U2z);
%!   z = res.bearings.zi(:, 1) + 0.5 * B;
%!   [p_ref, mdotz_ref] = reference_solution(h0, B, dh_dz, p0, p1, eta, U1z, U2z, rho, dx, dh_dt, z);
%!   figure("visible", "off");
%!   hold on;
%!   plot(res.bearings.zi(:, 1) + 0.5 * B, res.bearings.columns.p(:, 1), "-;p(z);1");
%!   plot(z, p_ref, "-;p_r;0");
%!   xlabel("z [m]");
%!   ylabel("p [Pa]");
%!   grid on;
%!   grid minor on;
%!   title("axial pressure distribution");
%!   figure("visible", "off");
%!   hold on;
%!   plot(res.bearings.zi(:, 1) + 0.5 * B, res.bearings.columns.h(:, 1), "-;h(z);1");
%!   plot(z, h0 + dh_dz * z, "-;h_r(z);0");
%!   xlabel("z [m]");
%!   ylabel("h [m]");
%!   grid on;
%!   grid minor on;
%!   title("radial clearance versus time");
%!   assert_simple(res.bearings.columns.p(:, 1), p_ref, 1e-4 * max(abs(p_ref)));
%!   assert_simple(-mdot1, mdotz_ref(1), 0.5e-2 * abs(mdotz_ref(1)));
%!   assert_simple(mdot2, mdotz_ref(end), 0.5e-2 * abs(mdotz_ref(end)));
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 6
%! function [p, mdotz] = reference_solution(h0, B, dh_dz, p0, p1, eta, U1z, U2z, rho, dx, dh_dt, z)
%!    p = "(dh_dz*h0^2*p1*z^2+2*B*dh_dz^2*h0*p1*z^2+B^2*dh_dz^3*p1*z^2-dh_dz*h0^2*p0*z^2+6*B*U2z*dh_dz*eta*z^2+6*B*U1z*dh_dz*eta*z^2+12*B*dh_dt*eta*z^2+2*h0^3*p1*z+4*B*dh_dz*h0^2*p1*z+2*B^2*dh_dz^2*h0*p1*z-2*h0^3*p0*z-6*B^2*U2z*dh_dz*eta*z-6*B^2*U1z*dh_dz*eta*z-12*B^2*dh_dt*eta*z+2*B*h0^3*p0+B^2*dh_dz*h0^2*p0)/(B*(2*h0+B*dh_dz)*(dh_dz*z+h0)^2)";
%!    mdotz = "-((12*B*dh_dt*eta*h0+6*B^2*dh_dt*dh_dz*eta)*rho*z+((h0^4+2*B*dh_dz*h0^3+B^2*dh_dz^2*h0^2)*p1+(-h0^4-2*B*dh_dz*h0^3-B^2*dh_dz^2*h0^2)*p0+(-6*B*U2z-6*B*U1z)*eta*h0^2+((-6*B^2*U2z-6*B^2*U1z)*dh_dz-6*B^2*dh_dt)*eta*h0)*rho)/(12*B*eta*h0+6*B^2*dh_dz*eta)";
%!    p = eval(vectorize(p));
%!    mdotz = eval(mdotz) * dx;
%! endfunction
%! param.d1 = 10e-3;
%! param.h = 1e-5;
%! param.M = int32(750);
%! param.N = int32(4);
%! param.output_bearing_data = true;
%! param.B = 25e-3;
%! cavitation = "mass conserving";
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "oct-mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: real omega1z = 0;\n");
%!     fputs(fd, "set: real omega2z = 0;\n");
%!     fputs(fd, "set: real v1z = 15;\n");
%!     fputs(fd, "set: real v2z = -7;\n");
%!     fputs(fd, "set: integer K = 36;\n");
%!     fputs(fd, "set: real t1 = 1e-2;\n");
%!     fputs(fd, "set: real dt = t1 / K;\n");
%!     fputs(fd, "set: real D2 = d1 + 2 * h;\n");
%!     fputs(fd, "set: real dy = 1.5 * h;\n");
%!     fputs(fd, "set: real eta = 1e-3;\n");
%!     fputs(fd, "set: real rho = 3600. * 1000.; # [l/h]\n\n");
%!     fputs(fd, "set: real epsilon_t0 = 0;\n");
%!     fputs(fd, "set: real epsilon_t1 = 0;\n");
%!     fputs(fd, "set: real delta_t0 = 0;\n");
%!     fputs(fd, "set: real delta_t1 = 0;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real p_pB2 = 10e5;\n");
%!     fputs(fd, "set: real p_mB2 = 0.8e5;\n");
%!     fputs(fd, "set: real pmax = 500000.;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 50;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-4, test, norm;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e5;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: hybrid, ms, 0.;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    use automatic differentiation;\n");
%!     fputs(fd, "    skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "        structural nodes: 2;\n");
%!     fputs(fd, "        joints: 2;\n");
%!     fputs(fd, "        loadable elements: 1;\n");
%!     fputs(fd, "        hydraulic nodes: 2;\n");
%!     fputs(fd, "        genels: 2;\n");
%!     fputs(fd, "        print: dof stats, to file;\n");
%!     fputs(fd, "        print: equation description, to file;\n");
%!     fputs(fd, "        print: dof description, to file;\n");
%!     fputs(fd, "        print: element connection, to file;\n");
%!     fputs(fd, "        print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 8;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        -v1z * t1,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        v1z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        -v2z * t1,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        v2z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            array, 2, -v1z * t1,\n");
%!     fputs(fd, "                      mult, const, v1z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            array, 2, -v2z * t1, mult, const, v2z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     switch (cavitation)
%!     case "non mass conserving"
%!       fputs(fd, "            hydraulic fluid, incompressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                pressure, 0.,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!     case "mass conserving"
%!       fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                1.,\n");
%!       fputs(fd, "                0.,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!       fputs(fd, "            viscosity vapor, eta,\n");
%!     otherwise
%!       error("unknown cavitation model: %s", cavitation);
%!     endswitch
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at bearing,\n");
%!     fputs(fd, "                bearing width, B,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                pockets,\n");
%!     fputs(fd, "                  shaft, 1,\n");
%!     fputs(fd, "                    position, d1 * pi * 0.5, 0.,\n");
%!     fputs(fd, "                    complete surface,\n");
%!     fputs(fd, "                    pocket height,\n");
%!     fputs(fd, "                    linear,\n");
%!     fputs(fd, "                      x, 0., d1 * pi,\n");
%!     fputs(fd, "                      z, -0.5 * B, 0.5 * B,\n");
%!     fputs(fd, "                      delta y, 0., -dy,\n");
%!     fputs(fd, "                               0., -dy,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure dof scale, pmax,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, no,\n");
%!     fputs(fd, "            output velocity, no,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   if (~opt_sol.verbose)
%!     opt_sol.logfile = [output_file, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%!   mdot1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end);
%!   mdot2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end);
%!   d = res.log_dat.bearings.cylindrical.d;
%!   D = res.log_dat.bearings.cylindrical.D;
%!   B = res.log_dat.bearings.cylindrical.B;
%!   h0 = 0.5 * (D - d);
%!   dh_dz = res.log_dat.vars.dy / res.log_dat.vars.B;
%!   p0 = res.log_dat.vars.p_mB2;
%!   p1 = res.log_dat.vars.p_pB2;
%!   eta = res.log_dat.vars.eta;
%!   U1z = res.log_dat.vars.v1z - res.log_dat.vars.v2z;
%!   U2z = 0;
%!   rho = res.log_dat.vars.rho;
%!   dx = res.log_dat.bearings.cylindrical.dm * pi;
%!   dh_dt = -dh_dz * (U1z - U2z);
%!   z = res.bearings.zi(:, 1) + 0.5 * B;
%!   [p_ref, mdotz_ref] = reference_solution(h0, B, dh_dz, p0, p1, eta, U1z, U2z, rho, dx, dh_dt, z);
%!   figure("visible", "off");
%!   hold on;
%!   plot(res.bearings.zi(:, 1) + 0.5 * B, res.bearings.columns.p(:, 1), "-;p(z);1");
%!   plot(z, p_ref, "-;p_r;0");
%!   xlabel("z [m]");
%!   ylabel("p [Pa]");
%!   grid on;
%!   grid minor on;
%!   title("axial pressure distribution");
%!   figure("visible", "off");
%!   hold on;
%!   plot(res.bearings.zi(:, 1) + 0.5 * B, res.bearings.columns.h(:, 1), "-;h(z);1");
%!   plot(z, h0 + dh_dz * z, "-;h_r(z);0");
%!   xlabel("z [m]");
%!   ylabel("h [m]");
%!   grid on;
%!   grid minor on;
%!   title("radial clearance versus time");
%!   assert_simple(res.bearings.columns.p(:, 1), p_ref, 1e-4 * max(abs(p_ref)));
%!   assert_simple(-mdot1, mdotz_ref(1), 0.5e-2 * abs(mdotz_ref(1)));
%!   assert_simple(mdot2, mdotz_ref(end), 0.5e-2 * abs(mdotz_ref(end)));
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 7
%! ## References
%! ## Erich Hansen et. al. 2022
%! ## An EHL extension of the unsteady FBNS Algorithm
%! param.output_bearing_data = true;
%! close all;
%! cavitation = "mass conserving";
%! param.h_min = 1e-6;
%! param.M = int32(201);
%! param.N = int32(4);
%! param.p_cav = 0;
%! param.beta = 2.4e9;
%! param.output_bearing_data = true;
%! param.L_x1 = 1e-2;
%! param.u_up = 0;
%! param.u_low = 1;
%! param.h_max = 1.05e-6;
%! param.a = 2e-3;
%! param.b = 3e-3;
%! param.h_p = 1e-6;
%! param.d1 = 10e-3;
%! param.D2 = param.d1 + 2 * param.h_min;
%! param.mu_0 = 1e-2;
%! param.rho_0 = 850;
%! param.p_amb = 1e5;
%! param.pref = 1e6;
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "oct-mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: integer drive_id_time_step_cfl = 5001;\n");
%!     fputs(fd, "set: integer drive_id_time_step = 5002;\n");
%!     fputs(fd, "set: real omega1z = 0;\n");
%!     fputs(fd, "set: real omega2z = 0;\n");
%!     fputs(fd, "set: real CFL = 20;\n");
%!     fputs(fd, "set: real t1 = 3;\n");
%!     fputs(fd, "set: real dt = CFL * L_x1 / M / (0.5 * abs(u_low + u_up));\n");
%!     fputs(fd, "set: real epsilon_t0 = 0;\n");
%!     fputs(fd, "set: real epsilon_t1 = 0;\n");
%!     fputs(fd, "set: real delta_t0 = 0;\n");
%!     fputs(fd, "set: real delta_t1 = 0;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / L_x1);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / L_x1);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        dummy steps tolerance: 1e-4;\n");
%!     fputs(fd, "        dummy steps max iterations: 100;\n");
%!     fputs(fd, "        dummy steps number: 0;\n");
%!     fputs(fd, "        dummy steps ratio: 1e-6;\n");
%!     fputs(fd, "        dummy steps method: implicit euler;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 500;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, norm, 1e-4, test, norm;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4, 1e5*1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 200;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-8;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        nonlinear solver: mcp newton min fb, lambda min, 0, tolerance x, 1e-8, mcp sigma, 0.9, mcp rho, 0, print convergence info, no;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    max iterations: 10;\n");
%!     fputs(fd, "    use automatic differentiation;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "        structural nodes: 2;\n");
%!     fputs(fd, "        joints: 2;\n");
%!     fputs(fd, "        loadable elements: 1;\n");
%!     fputs(fd, "        hydraulic nodes: 2;\n");
%!     fputs(fd, "        genels: 2;\n");
%!     fputs(fd, "        print: dof stats, to file;\n");
%!     fputs(fd, "        print: equation description, to file;\n");
%!     fputs(fd, "        print: dof description, to file;\n");
%!     fputs(fd, "        print: element connection, to file;\n");
%!     fputs(fd, "        print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 8;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        -u_up * t1,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        u_up,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        -u_low * t1,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        u_low,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_amb;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_amb;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            array, 2, -u_up * t1,\n");
%!     fputs(fd, "                      mult, const, u_up,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            array, 2, -u_low * t1, mult, const, u_low,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_amb;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_amb;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     switch (cavitation)
%!     case "non mass conserving"
%!       fputs(fd, "            hydraulic fluid, incompressible,\n");
%!       fputs(fd, "                density, rho_0,\n");
%!       fputs(fd, "                viscosity, mu_0,\n");
%!       fputs(fd, "                pressure, 0.,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!     case "mass conserving"
%!       fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!       fputs(fd, "                density, rho_0,\n");
%!       fputs(fd, "                beta,\n");
%!       fputs(fd, "                p_cav,\n");
%!       fputs(fd, "                viscosity, mu_0,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!       fputs(fd, "            viscosity vapor, mu_0,\n");
%!     otherwise
%!       error("unknown cavitation model: %s", cavitation);
%!     endswitch
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            enable mcp, yes,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at shaft,\n");
%!     fputs(fd, "                bearing width, L_x1,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                pockets,\n");
%!     fputs(fd, "                  shaft, 2,\n");
%!     fputs(fd, "                    position, D2 * pi * 0.5, a + 0.5 * b - 0.5 * L_x1,\n");
%!     fputs(fd, "                    rectangle, width, D2 * pi, height, b,\n");
%!     fputs(fd, "                    pocket height, const, -h_p,\n");
%!     fputs(fd, "                    position, D2 * pi * 0.5, 0.,\n");
%!     fputs(fd, "                    complete surface,\n");
%!     fputs(fd, "                    pocket height,\n");
%!     fputs(fd, "                    linear,\n");
%!     fputs(fd, "                      x, 0., d1 * pi,\n");
%!     fputs(fd, "                      z, -0.5 * L_x1, 0.5 * L_x1,\n");
%!     fputs(fd, "                      delta y, -(h_max - h_min), 0.,\n");
%!     fputs(fd, "                               -(h_max - h_min), 0.,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure dof scale, pref,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, no,\n");
%!     fputs(fd, "            output velocity, yes,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   opt_sol.mbdyn_command = "mbdyn";
%!   if (~opt_sol.verbose)
%!     opt_sol.logfile = [output_file, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%! q1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end) / (param.rho_0 * param.D2 * pi);
%! q2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end) / (param.rho_0 * param.D2 * pi);
%! %! Comparison to analytical solution of Fowell et al., 2007
%! %! Erik Hansen, 04.03.2022
%! %! Analytical solution:
%! analy.a         = param.a;
%! analy.b         = param.b;
%! analy.L_x1      = param.L_x1;
%! analy.h_max     = param.h_max;
%! analy.h_min     = param.h_min;
%! analy.h_p       = param.h_p;
%! analy.p_amb     = param.p_amb;
%! analy.mu_l      = param.mu_0;
%! analy.u_up_nc   = param.u_low;
%! analy.u_up_wc   = param.u_low;
%! analy.p_cav     = param.p_cav;
%! analy.K         = (analy.h_max - analy.h_min)/analy.h_min;
%! analy.p_cf_nc   = 6*analy.u_up_nc*analy.mu_l*analy.L_x1/analy.K/analy.h_min;
%! analy.p_cf_wc   = 6*analy.u_up_wc*analy.mu_l*analy.L_x1/analy.K/analy.h_min;
%! analy.Nx1       = rows(res.bearings.zi);
%! analy.x1        = linspace(0,analy.L_x1,analy.Nx1);
%! analy.h         = linspace(analy.h_max,analy.h_min,analy.Nx1);
%! analy.cond.p    = analy.x1 > analy.a     &    analy.x1 < analy.a + analy.b;
%! analy.h(analy.cond.p) = analy.h(analy.cond.p) + analy.h_p;
%! analy.h_2       = analy.h(analy.x1 == analy.a);
%! analy.h_2p      = analy.h_2 + analy.h_p;
%! analy.h_3       = analy.h(analy.x1 == analy.a + analy.b);
%! analy.h_3p      = analy.h_3 + analy.h_p;
%! analy.q_u1_nc   = (1/analy.h_2^2 - 1/analy.h_max^2) + (1/analy.h_3p^2 - 1/analy.h_2p^2) - (1/analy.h_3^2 - 1/analy.h_min^2);
%! analy.q_u1_nc   =((1/analy.h_2   - 1/analy.h_max)   + (1/analy.h_3p   - 1/analy.h_2p)   - (1/analy.h_3   - 1/analy.h_min)) / analy.q_u1_nc;
%! analy.A_1       = (analy.p_amb - analy.p_cav)/analy.p_cf_wc;
%! analy.q_u1_wc   = analy.h_max*analy.h_2/(analy.h_max + analy.h_2)*(analy.A_1*analy.h_max*analy.h_2/(analy.h_max - analy.h_2) + 1);
%! analy.A_3       = analy.A_1 + (1/analy.h_3 - 1/analy.h_min - 1/analy.h_3p) - analy.q_u1_wc*(1/analy.h_3^2 - 1/analy.h_min^2 - 1/analy.h_3p^2);
%! analy.h_bp      = (-1 - sqrt(1 + 4*analy.q_u1_wc*analy.A_3))/(2*analy.A_3);
%! analy.p_nc      = zeros(1,analy.Nx1);
%! analy.p_wc      = zeros(1,analy.Nx1);
%! analy.cond.entr = analy.x1 <= analy.a;
%! analy.p_nc(analy.cond.entr) = analy.p_amb + analy.p_cf_nc*((1./analy.h(analy.cond.entr) - 1/analy.h_max) - analy.q_u1_nc*(1./analy.h(analy.cond.entr).^2 - 1/analy.h_max^2));
%! analy.p_wc(analy.cond.entr) = analy.p_amb + analy.p_cf_wc*((1./analy.h(analy.cond.entr) - 1/analy.h_max) - analy.q_u1_wc*(1./analy.h(analy.cond.entr).^2 - 1/analy.h_max^2));
%! analy.p_nc_2    = analy.p_nc(analy.x1 == analy.a);
%! analy.p_nc(analy.cond.p) = analy.p_nc_2 + analy.p_cf_nc*((1./analy.h(analy.cond.p) - 1/analy.h_2p) - analy.q_u1_nc*(1./analy.h(analy.cond.p).^2 - 1/analy.h_2p^2));
%! analy.cond.cav = analy.x1 > analy.a & analy.h > analy.h_bp;
%! analy.p_wc(analy.cond.cav) = analy.p_cav;
%! analy.cond.no_cav = analy.x1 > analy.a & analy.x1 < analy.a + analy.b & analy.h <= analy.h_bp;
%! analy.p_wc(analy.cond.no_cav) = analy.p_cav + analy.p_cf_wc*((1./analy.h(analy.cond.no_cav) - 1/analy.h_bp) - analy.q_u1_wc*(1./analy.h(analy.cond.no_cav).^2 - 1/analy.h_bp^2));
%! analy.cond.exit = analy.x1 >= analy.a + analy.b;
%! analy.p_nc(analy.cond.exit) = analy.p_amb + analy.p_cf_nc*((1./analy.h(analy.cond.exit) - 1/analy.h_min) - analy.q_u1_nc*(1./analy.h(analy.cond.exit).^2 - 1/analy.h_min^2));
%! analy.p_wc(analy.cond.exit) = analy.p_amb + analy.p_cf_wc*((1./analy.h(analy.cond.exit) - 1/analy.h_min) - analy.q_u1_wc*(1./analy.h(analy.cond.exit).^2 - 1/analy.h_min^2));
%! z = res.bearings.zi(:, 1) + 0.5 * param.L_x1;
%! figure("visible", "off");
%! hold on;
%! plot(res.bearings.zi(:, 1) + 0.5 * param.L_x1, res.bearings.columns.p(:, 1), "-x;p(z);1");
%! plot(analy.x1, analy.p_wc, "-;Fowell;0");
%! xlabel("z [m]");
%! ylabel("p [Pa]");
%! grid on;
%! grid minor on;
%! title("axial pressure distribution");
%! figure("visible", "off");
%! hold on;
%! plot(res.bearings.zi(:, 1) + 0.5 * param.L_x1, res.bearings.columns.h(:, 1), "-+;h(z);1");
%! plot([0, param.L_x1], [param.h_max, param.h_min], "-;h_r(z);0");
%! xlabel("z [m]");
%! ylabel("h [m]");
%! grid on;
%! grid minor on;
%! title("radial clearance versus time");
%! tol_q = 1e-3;
%! tol_p = 2e-2;
%! assert_simple(-q1, analy.q_u1_wc, tol_q * analy.q_u1_wc);
%! assert_simple(q2, analy.q_u1_wc, tol_q * analy.q_u1_wc);
%! assert_simple(res.bearings.columns.p(:, 1), analy.p_wc(:), tol_p * max(abs(analy.p_wc)));
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%!demo
%! close all;
%! ## References:
%! ## Hans Juergen Butenschoen
%! ## Das hydrodynamische, zylindrische Gleitlager endlicher Breite unter instationaerer Belastung, Karlsruhe 1976
%!
%! B_d_r = [1, 1/2, 1/2.5, 1/3, 1/4, 1/5, 1/6, 1/8];
%!
%! epsilon_r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999];
%!
%! So_r =  [0.1196 0.0368  0.0243  0.0171  0.0098  0.0063  0.0044  0.0025
%!          0.2518 0.0783  0.0518  0.0366  0.0210  0.0136  0.0095  0.0054
%!          0.4091 0.1304  0.0867  0.0615  0.0354  0.0229  0.0160  0.0091
%!          0.6108 0.2026  0.1357  0.0968  0.0560  0.0364  0.0254  0.0144
%!          0.8903 0.3124  0.2117  0.1522  0.0888  0.0579  0.0406  0.0231
%!          1.3146 0.4982  0.3435  0.2496  0.1476  0.0969  0.0682  0.0390
%!          2.0432 0.8595  0.6075  0.4492  0.2708  0.1797  0.1274  0.0732
%!          3.5663 1.7339  1.2756  0.9687  0.6043  0.4085  0.2930  0.1706
%!          8.4352 5.0881  4.0187  3.2201  2.1595  1.5261  1.1263  0.6776
%!          18.7895 13.3083 11.2225 9.4993  6.9248  5.1833  3.9831  2.5202
%!          24.1172 17.7934 15.2823 13.1525 9.8517  7.5257  5.8712  3.7907
%!          33.1297 25.5920 22.4503 19.7034 15.2697 11.9804 9.5425  6.3397
%!          51.4774 41.9230 37.7412 33.9523 27.5040 22.3556 18.3874 12.7695
%!          107.7868 93.7881 87.2906 81.1597 70.0359 60.3874 52.1425 39.2568
%!          223.8850 203.2450 193.3490 183.8040 166.1540 149.8690 134.8910 109.6090
%!          1174.5400 1124.6200 1102.0700 1078.3100 1032.8700 989.1500 945.6700 864.7400];
%!
%! beta_r = [79.410 81.767  82.119  82.299  82.481  82.561  82.608  82.653
%!           73.858 75.142  75.283  75.344  75.287  75.406  75.414  75.422
%!           68.264 68.493  68.423  68.368  68.204  68.261  68.238  68.211
%!           62.571 61.778  61.544  61.382  61.208  61.115  61.062  61.007
%!           56.705 54.993  54.603  54.348  54.069  53.922  53.863  53.784
%!           50.536 48.049  47.521  47.191  46.825  46.647  46.554  46.449
%!           43.859 40.803  40.156  39.756  39.326  39.108  38.983  38.865
%!           36.235 32.938  32.216  31.761  31.249  30.988  30.840  30.692
%!           26.482 23.566  22.849  22.368  21.790  21.476  21.289  21.089
%!           19.450 17.265  16.648  16.207  15.632  15.291  15.075  14.832
%!           17.609 15.660  15.089  14.669  14.109  13.768  13.547  13.262
%!           15.484 13.818  13.313  12.929  12.299  12.062  11.838  11.570
%!           12.903 11.598  11.178  10.850  10.375  10.057  9.835   9.558
%!            9.416  8.587   8.301   8.066   7.703   7.440   7.242   6.975
%!            6.829  6.325   6.143   5.987   5.733   5.526   5.380   5.151
%!            3.196  3.048   2.989   2.940   2.848   2.769   2.708   2.599] * pi / 180;
%!
%! Q_r = [0.1603  0.0939  0.0768  0.0648  0.0492  0.0396  0.0331  0.0249
%!        0.3196  0.1878  0.1536  0.1296  0.0984  0.0792  0.0662  0.0498
%!        0.4765  0.2816  0.2304  0.1944  0.1476  0.1188  0.0993  0.0747
%!        0.6318  0.3755  0.3072  0.2592  0.1968  0.1583  0.1324  0.0996
%!        0.7852  0.4694  0.384   0.324   0.246   0.198   0.1655  0.1245
%!        0.9374  0.5634  0.461   0.3889  0.2953  0.2376  0.1986  0.1494
%!        1.0888  0.6578  0.5383  0.454   0.3446  0.2772  0.2317  0.1743
%!        1.2392  0.7529  0.6159  0.5193  0.394   0.3169  0.2648  0.1992
%!        1.3898  0.8453  0.6944  0.5852  0.4436  0.3567  0.298   0.2241
%!        1.4655  0.8984  0.7343  0.6186  0.4687  0.3767  0.3147  0.2366
%!        1.4808  0.9084  0.7425  0.6254  0.4737  0.3807  0.3181  0.2391
%!        1.4968  0.9185  0.7507  0.6322  0.4788  0.3848  0.3214  0.2417
%!        1.5124  0.9287  0.7585  0.6391  0.484   0.3889  0.3248  0.2442
%!        1.5277  0.9351  0.7674  0.6461  0.4892  0.393   0.3282  0.2467
%!        1.5354  0.9441  0.7716  0.6496  0.4918  0.3951  0.33    0.248
%!        1.5409  0.9482  0.7745  0.6525  0.494   0.3968  0.3314  0.2491];
%!
%! param.d1 = 10e-3;
%! param.Psi = 1e-3;
%! param.M = int32(20);
%! param.N = int32(75);
%! param.output_bearing_data = true;
%! param.B = 7e-3;
%! param.epsilon = 0.6;
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "oct-mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer hyd_node_id_inlet = 2003;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: integer genel_id_inlet = 4005;\n");
%!     fputs(fd, "set: real omega1z = 2 * pi * 50;\n");
%!     fputs(fd, "set: real omega2z = 2 * pi * 15;\n");
%!     fputs(fd, "set: real v1z = 0;\n");
%!     fputs(fd, "set: real v2z = 0;\n");
%!     fputs(fd, "set: real n = 3;\n");
%!     fputs(fd, "set: integer K = 36;\n");
%!     fputs(fd, "set: real t1 = abs(2 * pi * n / omega1z);\n");
%!     fputs(fd, "set: real dt = t1 / (n * K);\n");
%!     fputs(fd, "set: real D2 = d1 / (1. - Psi);\n");
%!     fputs(fd, "set: real eta = 1e-3;\n");
%!     fputs(fd, "set: real rho = 3600. * 1000.; # [l/h]\n\n");
%!     fputs(fd, "set: real epsilon_t0 = epsilon;\n");
%!     fputs(fd, "set: real epsilon_t1 = epsilon;\n");
%!     fputs(fd, "set: real delta_t0 = pi;\n");
%!     fputs(fd, "set: real delta_t1 = pi + omega2z * t1;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real p_pB2 = 0;\n");
%!     fputs(fd, "set: real p_mB2 = 0;\n");
%!     fputs(fd, "set: real p_in = 0;\n");
%!     fputs(fd, "set: real pmax = 500000.;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 50;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-4, test, norm;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e5;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: hybrid, ms, 0.;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        nonlinear solver: nox,\n");
%!     fputs(fd, "             jacobian operator, newton,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    use automatic differentiation;\n");
%!     fputs(fd, "    #skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "        structural nodes: 2;\n");
%!     fputs(fd, "        joints: 2;\n");
%!     fputs(fd, "        loadable elements: 1;\n");
%!     fputs(fd, "        hydraulic nodes: 3;\n");
%!     fputs(fd, "        genels: 3;\n");
%!     fputs(fd, "        print: dof stats, to file;\n");
%!     fputs(fd, "        print: equation description, to file;\n");
%!     fputs(fd, "        print: dof description, to file;\n");
%!     fputs(fd, "        print: element connection, to file;\n");
%!     fputs(fd, "        print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 8;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        v1z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing, null,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        v2z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_inlet, p_in;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            mult, const, v1z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, v2z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!     fputs(fd, "    genel: genel_id_inlet, clamp,\n");
%!     fputs(fd, "        hyd_node_id_inlet, hydraulic, p_in;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!     fputs(fd, "                density, rho,\n");
%!     fputs(fd, "                1.,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                viscosity, eta,\n");
%!     fputs(fd, "                temperature, 0,\n");
%!     fputs(fd, "            viscosity vapor, eta,\n");
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at bearing,\n");
%!     fputs(fd, "                bearing width, B,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure coupling conditions radial, 1,\n");
%!     fputs(fd, "                position, 0.5 * D2 * (delta_t0 + pi), 0.,\n");
%!     fputs(fd, "                rectangle, width, D2 * pi / (N - 1), height, B,\n");
%!     fputs(fd, "                pressure node, hyd_node_id_inlet,\n");
%!     fputs(fd, "            pressure dof scale, pmax,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, no,\n");
%!     fputs(fd, "            output velocity, no,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   if (~opt_sol.verbose)
%!     opt_sol.logfile = [output_file, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%!   d = res.log_dat.bearings.cylindrical.d;
%!   D = res.log_dat.bearings.cylindrical.D;
%!   B = res.log_dat.bearings.cylindrical.B;
%!   s = (D - d);
%!   Psi = (D - d) / D;
%!   eta = res.log_dat.bearings.eta;
%!   node_idx_2 = find(res.node_id == res.log_dat.bearings.cylindrical.nodes(2).label);
%!   Phi2 = res.trajectory{node_idx_2}(:, 4:6).';
%!   R2 = euler123_to_rotation_matrix(Phi2);
%!   Rb1 = res.log_dat.bearings.cylindrical.nodes(1).Rb;
%!   Rb2 = res.log_dat.bearings.cylindrical.nodes(2).Rb;
%!   F1 = -res.bearings.columns.F1.';
%!   M1 = -res.bearings.columns.M1.';
%!   omega1z = res.bearings.cylindrical.omega1z;
%!   omega2z = res.bearings.cylindrical.omega2z;
%!   omega_res = res.bearings.cylindrical.omega_res;
%!   epsilon = res.bearings.cylindrical.epsilon;
%!   delta = res.bearings.cylindrical.delta;
%!   Q_out1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end) / res.log_dat.vars.rho;
%!   Q_out2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end) / res.log_dat.vars.rho;
%!   Q_in = -res.genel_data{res.log_dat.vars.genel_idx_inlet}(end) / res.log_dat.vars.rho;
%!   P = norm(Rb2(:, 1:2).' * R2.' * F1);
%!   alpha = atan2(Rb2(:, 2).' * R2.' * F1, Rb2(:, 1).' * R2.' * F1);
%!   beta = sign(omega_res) * (delta - alpha);
%!   beta = mod(beta + pi, 2 * pi) - pi;
%!   mu = 2 * Rb2(:, 3).' * R2.' * M1 / (P * d) * sign(omega1z - omega2z);
%!   So = P * Psi^2 / (B * D * eta * abs(omega_res));
%!   Q = (Q_out1 + Q_out2) / ((0.5 * D)^3 * Psi * abs(omega_res));
%!   dQ = (Q_out1 + Q_out2 + Q_in) / Q_in;
%!   So_ref = interp2(B_d_r, epsilon_r, So_r, B / d, epsilon, "linear");
%!   beta_ref = interp2(B_d_r, epsilon_r, beta_r, B / d, epsilon, "linear");
%!   mu_ref = Psi * (abs((omega1z - omega2z) / omega_res) * pi / (sqrt(1 - epsilon^2) * So_ref) + sin(beta_ref) * abs(epsilon) / 2);
%!   Q_ref = interp2(B_d_r, epsilon_r, Q_r, B / d, epsilon, "linear");
%!   figure("visible", "off");
%!   set(plot(res.bearings.xi(floor(end/2), :) / (0.5 * D) * 180 / pi, res.bearings.columns.p(floor(end/2), :)), "linewidth", 3);
%!   grid on;
%!   grid minor on;
%!   xticks([0:30:360]);
%!   xlabel("Phi [deg]");
%!   ylabel("p [Pa]");
%!   title("midplane pressure distribution");
%!   figure("visible", "off");
%!   set(plot(res.bearings.xi(floor(end/2), :) / (0.5 * D) * 180 / pi, res.bearings.columns.h(floor(end/2), :) / (0.5 * (D - d))), "linewidth", 3);
%!   grid on;
%!   grid minor on;
%!   xlabel("Phi [deg]");
%!   ylabel("h/h0 [1]");
%!   title("relative radial clearance");
%!   figure("visible", "off");
%!   set(plot(res.bearings.xi(floor(end/2), :) / (0.5 * D) * 180 / pi, res.bearings.columns.rho(floor(end/2), :) / res.log_dat.vars.rho), "linewidth", 3);
%!   grid on;
%!   grid minor on;
%!   xlabel("Phi [deg]");
%!   ylabel("rho/rho0 [1]");
%!   title("volumetric filling ratio");
%!   figure_list();
%!   fprintf(stderr, "So / So_ref - 1 = %.2f\n", So / So_ref - 1);
%!   fprintf(stderr, "beta / beta_ref - 1 = %.2f\n", beta / beta_ref - 1);
%!   fprintf(stderr, "mu / mu_ref - 1 = %.2f\n", mu / mu_ref - 1);
%!   fprintf(stderr, "Q / Q_ref - 1 = %.2f\n", Q / Q_ref - 1);
%!   assert_simple(So, So_ref, 0.03 * So_ref);
%!   assert_simple(beta, beta_ref, 0.02 * beta_ref);
%!   assert_simple(mu, mu_ref, 0.03 * mu_ref);
%!   assert_simple(Q, Q_ref, 0.07 * Q_ref);
%!   assert_simple(dQ, 0, 1e-3);
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect


%!test
%! ############################################################################################################################
%! ## EHD TEST CASE according
%! ## Freund, Norman Owen, A thermo-elasto-hydrodynamic study of journal bearings, Doctor of Philosophy thesis, University of
%! ## Wollongong. Dept. of Mechanical Engineering, University of Wollongong, 1995. http://ro.uow.edu.au/theses/1568
%! ############################################################################################################################
%!
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Figure 46b, page 175
%! ref_data_w = [  0, 195;
%!                10, 122;
%!                20,  65;
%!                30,  28;
%!                40,  18;
%!                50,  20;
%!                60,  30;
%!                70,  40;
%!                80,  48;
%!                90,  50;
%!               100,  60;
%!               110,  45;
%!               120,  30;
%!               130,   0];
%! ref_p = [  0,     0;
%!           10,   100;
%!           20,   200;
%!           30,   300;
%!           40,   375;
%!           50,   425;
%!           60,   525;
%!           70,   650;
%!           80,   750;
%!           90,   875;
%!          100,   975;
%!          110,  1100;
%!          120,  1175;
%!          130,  1175;
%!          140,  1050;
%!          150,   925;
%!          160,   700;
%!          170,   450;
%!          180,   200;
%!          190,    50;
%!          200,     0];
%! output_file = "";
%! have_mesh_size_binary = false;
%! unwind_protect
%!   output_file = tempname();
%!   if (ispc())
%!     output_file(output_file == "\\") = "/";
%!   endif
%!   [status, output] = shell("which fem_pre_mesh_size", true);
%!   if (status == 0)
%!     have_mesh_size_binary = true;
%!   endif
%!   if (~have_mesh_size_binary)
%!     error("fem_pre_mesh_size was not installed\nrun ./configure && make install inside the src directory");
%!   endif
%!   SI_unit_meter = 1;
%!   SI_unit_kilogram = 1;
%!   SI_unit_second = 1;
%!   SI_unit_kelvin = 1;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   ## Table 6, page 154
%!   ## Table 7, page 156
%!   param.E = 200000e6 / SI_unit_pascal;
%!   param.nu = 0.3;
%!   param.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.Di = 200e-3 / SI_unit_meter;       %! inner diameter shell (named D according to Freund)
%!   param.Do = 299e-3 / SI_unit_meter;       %! outer diameter bearing
%!   param.Wo = 200e-3 / SI_unit_meter;       %! total widht of bearing (named L according to Freund)
%!   param.d1 = 1e-3 / SI_unit_meter;         %! diaphragm thickness (1mm-5mm)
%!   param.d2 = 64e-3 / SI_unit_meter;        %! diaphragm half free edge width
%!   param.d3 = 132 * pi / 180 * param.Di / 2;# diaphragm length
%!   param.d4 = 7.5e-3 / SI_unit_meter;       %! diaphragm half end width
%!   param.d5 = 0e-3 / SI_unit_meter;         %! diaphragm straight extent
%!   param.cr = 150e-6 / SI_unit_meter;       %! bearing radial clearance
%!   param.etal = 0.00929 / (SI_unit_pascal * SI_unit_second); %! dynamic viscosity lubricant
%!   param.fact_etav = 1;                     %! etav/etal
%!   param.rhol = 844 / (SI_unit_kilogram / SI_unit_meter^3); %! density lubricant
%!   param.betal = 2.4e9 / SI_unit_pascal;    %! bulk modulus
%!   param.pcl = 0 / SI_unit_pascal;          %! cavitation pressure
%!   param.pin = 0 / SI_unit_pascal;        %! oil feed pressure
%!   param.pside = 0 / SI_unit_pascal;      %! side pressure
%!   param.Tl = 40 / SI_unit_kelvin;          %! liquid temperature
%!   param.d7 = 1e-3 / SI_unit_meter;         %! circumferential length of oil supply slot
%!   param.U1 = 28 / (SI_unit_meter / SI_unit_second); %! fluid velocity at the journal surface (named U2 according to Freund)
%!   param.U2 = 0 / (SI_unit_meter / SI_unit_second); %! fluid velocity at the shell surface (named U1 according to Freund)
%!   param.delta = 180 * pi / 180;            %! position of eccentricity
%!   param.epsilon = 0.3;                     %! relative eccentricity
%!   param.diaph_sec_n = 72;                  %! number of cross sections for diaphragm
%!   param.h = 15e-3 / SI_unit_meter;         %! mesh size for hydraulic mesh
%!   param.h1 = 4e-3 / SI_unit_meter;         %! mesh size in the area of the diaphragm
%!   param.h2 = 20e-3 / SI_unit_meter;        %! mesh size outside the are of the diaphragm
%!   param.ht = 10e-3 / SI_unit_meter;        %! mesh transition region
%!   param.pref = 1e6 / SI_unit_pascal;       %! reference pressure
%!   param.damp_alpha = 0 / (1 / SI_unit_second);  %! mass damping factor
%!   param.damp_beta = 0 / SI_unit_second;         %! stiffness damping factor
%!   n = 2;
%!   k = 72;
%!   param.t1 = n * param.Di * pi / param.U1;
%!   param.dt = param.t1 / (n * k);
%!   if (~have_mesh_size_binary)
%!     opt_mesh.mesh.element_size = param.h;
%!   endif
%!   opt_mesh.mesh.jacobian_range = [0.5, 1.5];
%!   opt_mesh.verbose = false;
%!   opt_mesh.output_file = [output_file, "_msh"];
%!   options.geo_tol = sqrt(eps);
%!   options_mbdyn.mbdyn_command = "mbdyn";
%!   empty_cell = cell(1, 2);
%!   group_defs = struct("id", empty_cell, ...
%!                       "name", empty_cell, ...
%!                       "R", empty_cell, ...
%!                       "X0", empty_cell, ...
%!                       "Xi", empty_cell, ...
%!                       "type", empty_cell, ...
%!                       "geometry", empty_cell, ...
%!                       "compliance_matrix", empty_cell);
%!   group_defs(1).id = 1;
%!   group_defs(1).name = "node_id_shell_support";
%!   group_defs(1).R = eye(3);
%!   group_defs(1).X0 = zeros(3, 1);
%!   group_defs(1).Xi = zeros(3, 1);
%!   group_defs(1).type = "cylinder";
%!   group_defs(1).geometry.rmin = 0.5 * param.Do;
%!   group_defs(1).geometry.rmax = 0.5 * param.Do;
%!   group_defs(1).geometry.zmin = -0.5 * param.Wo;
%!   group_defs(1).geometry.zmax = 0.5 * param.Wo;
%!   group_defs(1).compliance_matrix.matrix_type = "none";

%!   group_defs(2).id = 2;
%!   group_defs(2).name = "node_id_shell_bearing";
%!   group_defs(2).R = eye(3);
%!   group_defs(2).X0 = zeros(3, 1);
%!   group_defs(2).Xi = zeros(3, 1);
%!   group_defs(2).type = "cylinder";
%!   group_defs(2).geometry.rmin = 0.5 * param.Di;
%!   group_defs(2).geometry.rmax = 0.5 * param.Di;
%!   group_defs(2).geometry.zmin = -0.5 * param.Wo;
%!   group_defs(2).geometry.zmax = 0.5 * param.Wo;
%!   group_defs(2).compliance_matrix.matrix_type = "nodal substruct";
%!   group_defs(2).compliance_matrix.bearing_type = "shell";
%!   group_defs(2).compliance_matrix.bearing_model = "EHD/FD";
%!   group_defs(2).compliance_matrix.reference_pressure = param.pref;
%!   group_defs(2).compliance_matrix.mesh_size = param.h;
%!   group_defs(2).bearing = "elem_id_bearing";
%!   fd = -1;
%!   unwind_protect
%!     fd = fopen([output_file, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", [output_file, ".geo"]);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Geometry.OCCUnionUnify = 0;\n");
%!     fputs(fd, "Wd_0 = 2 * d2;\n");
%!     fputs(fd, "Wd_1 = 2 * d4;\n");
%!     fputs(fd, "td = d1;\n");
%!     fputs(fd, "dxo = d7;\n");
%!     fputs(fd, "Phid_0 = 0 / (0.5 * Di);\n");
%!     fputs(fd, "Phid_1 = d5 / (0.5 * Di);\n");
%!     fputs(fd, "Phid_2 = d3 / (0.5 * Di);\n");
%!     fputs(fd, "hd = (Do - Di) / 2. - td;\n");
%!     fputs(fd, "Phio = dxo / (0.5 * Di);\n");
%!     fputs(fd, "diaph_sec_m = 4;\n");
%!     fputs(fd, "diaph_sec_p[] = {};\n");
%!     fputs(fd, "For i In {0:diaph_sec_n - 1}\n");
%!     fputs(fd, "    j = i * diaph_sec_m;\n");
%!     fputs(fd, "    alpha = i / (diaph_sec_n - 1);\n");
%!     fputs(fd, "    Phii = (Phid_0 + 0.5 * Phio + Phid_2 * alpha);   \n");
%!     fputs(fd, "    If (Phid_2 * alpha >= Phid_1)\n");
%!     fputs(fd, "        Wdi = Wd_0 + (Wd_1 - Wd_0) * (Phid_2 * alpha - Phid_1) / (Phid_2 - Phid_1);\n");
%!     fputs(fd, "    Else\n");
%!     fputs(fd, "        Wdi = Wd_0;\n");
%!     fputs(fd, "    EndIf\n");
%!     fputs(fd, "    diaph_sec_p[j + 0] = newp;\n");
%!     fputs(fd, "    Point(diaph_sec_p[j + 0]) = {(0.5 * Di + td) * Cos(Phii), (0.5 * Di + td) * Sin(Phii), -0.5 * Wdi};\n");
%!     fputs(fd, "    diaph_sec_p[j + 1] = newp;\n");
%!     fputs(fd, "    Point(diaph_sec_p[j + 1]) = {(0.5 * Do + hd) * Cos(Phii), (0.5 * Do + hd) * Sin(Phii), -0.5 * Wdi};\n");
%!     fputs(fd, "    diaph_sec_p[j + 2] = newp;\n");
%!     fputs(fd, "    Point(diaph_sec_p[j + 2]) = {(0.5 * Do + hd) * Cos(Phii), (0.5 * Do + hd) * Sin(Phii), 0.5 * Wdi};\n");
%!     fputs(fd, "    diaph_sec_p[j + 3] = newp;\n");
%!     fputs(fd, "    Point(diaph_sec_p[j + 3]) = {(0.5 * Di + td) * Cos(Phii), (0.5 * Di + td) * Sin(Phii), 0.5 * Wdi};\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "bearing_sec_p[] = {};\n");
%!     fputs(fd, "bearing_sec_p[0] = newp;\n");
%!     fputs(fd, "Point(bearing_sec_p[0]) = {0.5 * Di, 0., -0.5 * Wo};\n");
%!     fputs(fd, "bearing_sec_p[1] = newp;\n");
%!     fputs(fd, "Point(bearing_sec_p[1]) = {0.5 * Do, 0., -0.5 * Wo};\n");
%!     fputs(fd, "bearing_sec_p[2] = newp;\n");
%!     fputs(fd, "Point(bearing_sec_p[2]) = {0.5 * Do, 0., 0.5 * Wo};\n");
%!     fputs(fd, "bearing_sec_p[3] = newp;\n");
%!     fputs(fd, "Point(bearing_sec_p[3]) = {0.5 * Di, 0., 0.5 * Wo};\n");
%!     fputs(fd, "Phi = Phid_0 - 0.5 * Phio;\n");
%!     fputs(fd, "ossl_sec_p[] = {};\n");
%!     fputs(fd, "ossl_sec_p[0] = newp;\n");
%!     fputs(fd, "Point(ossl_sec_p[0]) = {(0.5 * Di - td) * Cos(Phi), (0.5 * Di - td) * Sin(Phi), -0.5 * Wd_0};\n");
%!     fputs(fd, "ossl_sec_p[1] = newp;\n");
%!     fputs(fd, "Point(ossl_sec_p[1]) = {(0.5 * Do + hd) * Cos(Phi), (0.5 * Do + hd) * Sin(Phi), -0.5 * Wd_0};\n");
%!     fputs(fd, "ossl_sec_p[2] = newp;\n");
%!     fputs(fd, "Point(ossl_sec_p[2]) = {(0.5 * Do + hd) * Cos(Phi), (0.5 * Do + hd) * Sin(Phi), 0.5 * Wd_0};\n");
%!     fputs(fd, "ossl_sec_p[3] = newp;\n");
%!     fputs(fd, "Point(ossl_sec_p[3]) = {(0.5 * Di - td) * Cos(Phi), (0.5 * Di - td) * Sin(Phi), 0.5 * Wd_0};\n");
%!     fputs(fd, "diaph_sec_l[] = {};\n");
%!     fputs(fd, "For i In {0:diaph_sec_n - 1}\n");
%!     fputs(fd, "    For j In {0:diaph_sec_m - 1}\n");
%!     fputs(fd, "        k0 = i * diaph_sec_m + j;\n");
%!     fputs(fd, "        If (j == diaph_sec_m - 1)\n");
%!     fputs(fd, "           k1 = i * diaph_sec_m;\n");
%!     fputs(fd, "        Else\n");
%!     fputs(fd, "           k1 = k0 + 1;\n");
%!     fputs(fd, "        EndIf\n");
%!     fputs(fd, "        diaph_sec_l[k0] = newreg;\n");
%!     fputs(fd, "        Line(diaph_sec_l[k0]) = {diaph_sec_p[k0], diaph_sec_p[k1]};\n");
%!     fputs(fd, "    EndFor\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "bearing_sec_l[] = {};\n");
%!     fputs(fd, "For i In {0:#bearing_sec_p[] - 1} \n");
%!     fputs(fd, "    bearing_sec_l[i] = newreg;\n");
%!     fputs(fd, "    If (i == #bearing_sec_p[] - 1)\n");
%!     fputs(fd, "       j = 0;\n");
%!     fputs(fd, "    Else\n");
%!     fputs(fd, "       j = i + 1;\n");
%!     fputs(fd, "    EndIf\n");
%!     fputs(fd, "    Line(bearing_sec_l[i]) = {bearing_sec_p[i], bearing_sec_p[j]};\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "ossl_sec_l[] = {};\n");
%!     fputs(fd, "For i In {0:#ossl_sec_p[] - 1} \n");
%!     fputs(fd, "    ossl_sec_l[i] = newreg;\n");
%!     fputs(fd, "    If (i == #ossl_sec_p[] - 1)\n");
%!     fputs(fd, "       j = 0;\n");
%!     fputs(fd, "    Else\n");
%!     fputs(fd, "       j = i + 1;\n");
%!     fputs(fd, "    EndIf\n");
%!     fputs(fd, "    Line(ossl_sec_l[i]) = {ossl_sec_p[i], ossl_sec_p[j]};\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "diaph_sec_ll[] = {};\n");
%!     fputs(fd, "For i In {0:diaph_sec_n - 1}\n");
%!     fputs(fd, "    diaph_sec_ll[i] = newreg;\n");
%!     fputs(fd, "    Line Loop(diaph_sec_ll[i]) = {diaph_sec_l[i * diaph_sec_m],\n");
%!     fputs(fd, "                                  diaph_sec_l[i * diaph_sec_m + 1],\n");
%!     fputs(fd, "                                  diaph_sec_l[i * diaph_sec_m + 2],\n");
%!     fputs(fd, "                                  diaph_sec_l[i * diaph_sec_m + 3]};\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "bearing_sec_ll = newreg;\n");
%!     fputs(fd, "Line Loop(bearing_sec_ll) = {bearing_sec_l[]};\n");
%!     fputs(fd, "ossl_sec_ll = newreg;\n");
%!     fputs(fd, "Line Loop(ossl_sec_ll) = {ossl_sec_l[]};\n");
%!     fputs(fd, "bearing_sec_s = newreg;\n");
%!     fputs(fd, "Plane Surface(bearing_sec_s) = {bearing_sec_ll};\n");
%!     fputs(fd, "ossl_sec_s = newreg;\n");
%!     fputs(fd, "Plane Surface(ossl_sec_s) = {ossl_sec_ll};\n");
%!     fputs(fd, "diaph_v = newv;\n");
%!     fputs(fd, "ThruSections(diaph_v) = {diaph_sec_ll[]};\n");
%!     fputs(fd, "For i In {0:#diaph_sec_ll[] - 1}\n");
%!     fputs(fd, "    Recursive Delete { Curve{diaph_sec_ll[i]}; }\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "For i In {0:#diaph_sec_l[] - 1}\n");
%!     fputs(fd, "    Recursive Delete { Curve{diaph_sec_l[i]}; }\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "For i In {0:#diaph_sec_p[] - 1}\n");
%!     fputs(fd, "    Recursive Delete { Curve{diaph_sec_p[i]}; }\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "bearing_rot_v1[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{bearing_sec_s}; };\n");
%!     fputs(fd, "bearing_rot_v2[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{bearing_rot_v1[0]}; };\n");
%!     fputs(fd, "bearing_rot_v3[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{bearing_rot_v2[0]}; };\n");
%!     fputs(fd, "bearing_rot_v4[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} { Surface{bearing_rot_v3[0]}; };\n");
%!     fputs(fd, "ossl_rot_v[] = Extrude {{0, 0, 1}, {0, 0, 0}, Phio} { Surface{ossl_sec_s}; };\n");
%!     fputs(fd, "bearing_v5 = newv;\n");
%!     fputs(fd, "BooleanUnion(bearing_v5) = {Volume{bearing_rot_v1[1]}; Delete; }{ Volume{bearing_rot_v2[1], bearing_rot_v3[1], bearing_rot_v4[1]}; Delete; };\n");
%!     fputs(fd, "bearing_v = newv;\n");
%!     fputs(fd, "BooleanDifference(bearing_v) = {Volume{bearing_v5}; Delete; }{ Volume{diaph_v, ossl_rot_v[1]}; Delete; };\n");
%!     fputs(fd, "bearing_bnd[] = Unique(Abs(Boundary{Volume{bearing_v};}));\n");
%!     fputs(fd, "Physical Volume(\"bearing\", 1) = {bearing_v};\n");
%!     fputs(fd, "For i In {0:#bearing_bnd[] - 1}\n");
%!     fputs(fd, "    Physical Surface(i) = {bearing_bnd[i]};\n");
%!     fputs(fd, "EndFor\n");
%!     if (have_mesh_size_binary)
%!       fputs(fd, "iCurrField = 0;\n");
%!       fputs(fd, "iCurrField++;\n");
%!       fputs(fd, "Field[iCurrField] = ExternalProcess;\n");
%!       fputs(fd, "Field[iCurrField].CommandLine = Sprintf(\"fem_pre_mesh_size diaphragm Wo=%g Di=%g Do=%g d1=%g d2=%g d3=%g d4=%g d7=%g h1=%g h2=%g ht=%g\", Wo, Di, Do, 1.5 * d1, d2, d3, d4, d7, h1, h2, ht);\n");
%!       fputs(fd, "Background Field = iCurrField;\n");
%!     endif
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect

%!   mesh = fem_pre_mesh_unstruct_create([output_file, ".geo"], param, opt_mesh);
%!   mesh.groups.tria6 = fem_pre_mesh_groups_create(mesh, group_defs, options.geo_tol).tria6;
%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.materials.tet10(mesh.groups.tet10(find([mesh.groups.tet10.id == 1])).elements) = 1;
%!   mesh.material_data.rho = param.rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(param.E, param.nu);

%!   cms_opt.invariants = true;
%!   cms_opt.refine_max_iter = int32(0);
%!   cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!   cms_opt.verbose = false;
%!   cms_opt.modes.number = 0;
%!   cms_opt.element.name = "elem_id_diaphragm_cms";

%!   node_set = int32(rows(mesh.nodes) + (1:numel(group_defs)));
%!   node_names = {group_defs.name};
%!   cms_opt.nodes.modal.number = node_set(1);
%!   cms_opt.nodes.modal.name = node_names{1};

%!   for j=1:numel(node_set) - 1
%!     cms_opt.nodes.interfaces(j).number = node_set(j + 1);
%!     cms_opt.nodes.interfaces(j).name = node_names{j + 1};
%!   endfor

%!   idx_grp_itf = find([group_defs.id] > 1);
%!   idx_grp_modal = find([group_defs.id] == 1);

%!   mesh.nodes([cms_opt.nodes.interfaces.number], 1:3) = [group_defs(idx_grp_itf).Xi].';
%!   mesh.nodes(cms_opt.nodes.modal.number, 1:3) = group_defs(idx_grp_modal).Xi.';

%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, [group_defs(idx_grp_itf).id], node_set(idx_grp_itf));

%!   idx_grp_outer = find([mesh.groups.tria6.id] == 1);

%!   load_case.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
%!   load_case.locked_dof(mesh.groups.tria6(idx_grp_outer).nodes, 1:3) = true; %! Norman Owen Freund 1995, page 8
%!   load_case.locked_dof(cms_opt.nodes.modal.number, 1:6) = true;

%!   bearing_surf = repmat(struct("group_idx", [], "X0", [], "R", [], "options", [], "name", [], "bearing", []), 1, numel(group_defs));
%!   num_comp_mat = int32(0);

%!   for j=1:numel(group_defs)
%!     switch group_defs(j).compliance_matrix.matrix_type
%!       case "none"
%!       otherwise
%!         ++num_comp_mat;
%!         bearing_surf(num_comp_mat).group_idx = find([mesh.groups.tria6.id] == group_defs(j).id);
%!         bearing_surf(num_comp_mat).name = group_defs(j).name;
%!         bearing_surf(num_comp_mat).bearing = group_defs(j).bearing;
%!         bearing_surf(num_comp_mat).X0 = group_defs(j).X0;
%!         bearing_surf(num_comp_mat).R = group_defs(j).R;
%!         bearing_surf(num_comp_mat).options = group_defs(j).compliance_matrix;
%!         bearing_surf(num_comp_mat).master_node_no = node_set(j);

%!         switch group_defs(j).type
%!           case "cylinder"
%!             bearing_surf(num_comp_mat).r = mean([group_defs(j).geometry.rmax, group_defs(j).geometry.rmin]);
%!             bearing_surf(num_comp_mat).w = group_defs(j).geometry.zmax - group_defs(j).geometry.zmin;
%!           otherwise
%!             error("bearing geometry type \"%s\" not implemented", group_defs(j).type);
%!         endswitch

%!         bearing_surf(num_comp_mat).nodes = mesh.groups.tria6(find([[mesh.groups.tria6].id] == group_defs(j).id)).nodes;
%!     endswitch
%!   endfor

%!   bearing_surf = bearing_surf(1:num_comp_mat);

%!   [load_case_pressure, bearing_surf] = fem_ehd_pre_comp_mat_load_case(mesh, bearing_surf);

%!   load_case = fem_pre_load_case_merge(load_case, load_case_pressure);

%!   [mesh, ...
%!    mat_ass, ...
%!    dof_map, ...
%!    sol_eig, ...
%!    cms_opt] = fem_cms_create(mesh, load_case, cms_opt);

%!   mat_ass.Dred = param.damp_alpha * mat_ass.Mred + param.damp_beta * mat_ass.Kred;

%!   fem_cms_export([output_file, "_cms"], mesh, dof_map, mat_ass, cms_opt);

%!   comp_mat = fem_ehd_pre_comp_mat_unstruct(mesh, ...
%!                                            mat_ass, ...
%!                                            dof_map, ...
%!                                            cms_opt, ...
%!                                            bearing_surf);

%!   for j=1:numel(comp_mat)
%!     comp_mat_file = [output_file, "_", bearing_surf(j).bearing, ".dat"];
%!     fem_ehd_pre_comp_mat_export(comp_mat(j), bearing_surf(j).options, comp_mat_file);
%!   endfor

%!   unwind_protect
%!     fd = -1;

%!     [fd, msg] = fopen([output_file, "_shell.set"], "w");

%!     if (fd == -1)
%!       error("failed to open file \"%s\": %s", [output_file, "_shell.set"], msg);
%!     endif

%!     for j=1:numel(comp_mat)
%!       fprintf(fd, "set: number_of_nodes_x = %d;\n", numel(comp_mat(j).bearing_surf.grid_x) + 1);
%!       fprintf(fd, "set: number_of_nodes_z = %d;\n", numel(comp_mat(j).bearing_surf.grid_z));
%!     endfor
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect

%!   mbdyn_pre_write_param_file([output_file, ".set"], param);

%!   options_mbdyn.output_file = [output_file, "_mbd"];
%!   options_mbdyn.f_run_mbdyn2easyanim = false;
%!   fd = -1;
%!   unwind_protect
%!     fd = fopen([output_file, "_mbd.mbdyn"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", [output_file, "_mbd.mbdyn"]);
%!     endif
%!     fputs(fd, "set: integer ref_id_assembly = 1001;\n");
%!     fputs(fd, "set: integer ref_id_shell_bearing = 1002;\n");
%!     fputs(fd, "set: integer ref_id_shell_support = 1003;\n");
%!     fputs(fd, "set: integer ref_id_journal_bearing = 1004;\n");
%!     fputs(fd, "set: integer ref_id_journal_support = 1005;\n");
%!     fputs(fd, "set: integer node_id_shell_support = 2001;\n");
%!     fputs(fd, "set: integer node_id_shell_bearing = 2002;\n");
%!     fputs(fd, "set: integer node_id_journal_bearing = 2004;\n");
%!     fputs(fd, "set: integer joint_id_shell_support = 3001;\n");
%!     fputs(fd, "set: integer joint_id_journal_support = 3002;\n");
%!     fputs(fd, "set: integer elem_id_diaphragm_cms = 3003;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 3005;\n");
%!     fputs(fd, "set: integer number_of_nodes_x;\n");
%!     fputs(fd, "set: integer number_of_nodes_z;\n");
%!     fprintf(fd, "include: \"%s.set\";\n", output_file);
%!     fprintf(fd, "include: \"%s_shell.set\";\n", output_file);
%!     fputs(fd, "set: real omega1 = U1 / (0.5 * Di - cr);\n");
%!     fputs(fd, "set: real omega2 = U2 / (0.5 * Di);\n");
%!     fputs(fd, "set: real omega = omega1 - omega2;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        linear solver: umfpack, grad, scale, row max column max, always;\n");
%!     fputs(fd, "        nonlinear solver: nox,\n");
%!     fputs(fd, "             jacobian operator, newton,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, minmax, 1e-4, test, norm;\n");
%!     fputs(fd, "        max iterations: 100;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-3;\n");
%!     fputs(fd, "        derivatives max iterations: 5;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!     fputs(fd, "        output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       print: dof stats, to file;\n");
%!     fputs(fd, "       print: dof description, to file;\n");
%!     fputs(fd, "       print: equation description, to file;\n");
%!     fputs(fd, "       structural nodes: 3;\n");
%!     fputs(fd, "       joints: 3;\n");
%!     fputs(fd, "       loadable elements: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_assembly,\n");
%!     fputs(fd, "        position, reference, global, null,\n");
%!     fputs(fd, "        orientation, reference, global, eye,\n");
%!     fputs(fd, "        velocity, reference, global, null,\n");
%!     fputs(fd, "        angular velocity, reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_shell_bearing,\n");
%!     fputs(fd, "        position, reference, ref_id_assembly, null,\n");
%!     fputs(fd, "        orientation, reference, ref_id_assembly, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_assembly, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_assembly,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                omega2;\n");
%!     fputs(fd, "reference: ref_id_shell_support,\n");
%!     fputs(fd, "        position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "        orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_shell_bearing, null;\n");
%!     fputs(fd, "reference: ref_id_journal_bearing,\n");
%!     fputs(fd, "        position, reference, ref_id_assembly,\n");
%!     fputs(fd, "                  epsilon * cr * cos(delta),\n");
%!     fputs(fd, "                  epsilon * cr * sin(delta),\n");
%!     fputs(fd, "                  0.,\n");
%!     fputs(fd, "        orientation, reference, ref_id_assembly, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_assembly, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_assembly,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                omega1;\n");
%!     fputs(fd, "reference: ref_id_journal_support,\n");
%!     fputs(fd, "        position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "        orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_journal_bearing, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_shell_support, modal,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_support, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shell_support, eye,\n");
%!     fputs(fd, "                velocity, reference, ref_id_shell_support, null,\n");
%!     fputs(fd, "                angular velocity, reference, ref_id_shell_support, null;\n");
%!     fputs(fd, "        structural: node_id_shell_bearing, static,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                velocity, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                angular velocity, reference, ref_id_shell_bearing, null;\n");
%!     fputs(fd, "        structural: node_id_journal_bearing, static,\n");
%!     fputs(fd, "                position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                velocity, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                angular velocity, reference, ref_id_journal_bearing, null;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_shell_support, total pin joint,\n");
%!     fputs(fd, "                node_id_shell_support,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                position constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        component,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                orientation constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        component,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                mult, time, omega2;\n");
%!     fputs(fd, "        joint: joint_id_journal_support, total pin joint,\n");
%!     fputs(fd, "                node_id_journal_bearing,\n");
%!     fputs(fd, "                position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                position constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        null,\n");
%!     fputs(fd, "                orientation constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        component,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                mult, time, omega1;\n");
%!     fprintf(fd, "        include: \"%s_cms.elm\";\n", output_file);
%!     fputs(fd, "        user defined: elem_id_bearing, hydrodynamic plain bearing2,\n");
%!     fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!     fputs(fd, "                density, rhol,\n");
%!     fputs(fd, "                betal,\n");
%!     fputs(fd, "                pressure, pcl,\n");
%!     fputs(fd, "                viscosity, etal,\n");
%!     fputs(fd, "                temperature, Tl,\n");
%!     fputs(fd, "            viscosity vapor, factor, fact_etav,\n");
%!     fputs(fd, "                mesh, linear finite difference,\n");
%!     fputs(fd, "                geometry, cylindrical,\n");
%!     fputs(fd, "                        mesh position, at bearing,\n");
%!     fputs(fd, "                        bearing width, Wo,\n");
%!     fputs(fd, "                        shaft diameter, Di - 2 * cr,\n");
%!     fputs(fd, "                        bearing diameter, Di,\n");
%!     fputs(fd, "                shaft node, node_id_journal_bearing,\n");
%!     fputs(fd, "                offset, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                bearing node, node_id_shell_bearing,\n");
%!     fputs(fd, "                offset, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                        number of nodes z, number_of_nodes_z,\n");
%!     fputs(fd, "                        number of nodes Phi, number_of_nodes_x,\n");
%!     fputs(fd, "                        boundary conditions,\n");
%!     fputs(fd, "                                        pressure, pside,\n");
%!     fputs(fd, "                                        pressure, pside,\n");
%!     fputs(fd, "                        lubrication grooves, 1,\n");
%!     fputs(fd, "                                at bearing,\n");
%!     fputs(fd, "                                        pressure, pin,\n");
%!     fputs(fd, "                                        position, 0., 0.,\n");
%!     fputs(fd, "                                        rectangle, width, d7, height, Wo,\n");
%!     fputs(fd, "                compliance model,\n");
%!     fprintf(fd, "                        matrix, 1, from file, \"%s_elem_id_bearing.dat\",\n", output_file);
%!     fputs(fd, "                        E1, E,\n");
%!     fputs(fd, "                        nu1, nu,\n");
%!     fputs(fd, "                        modal element, elem_id_diaphragm_cms,\n");
%!     fputs(fd, "                pressure dof scale, pref,\n");
%!     fputs(fd, "                reynolds equation scale, dt / (rhol * cr),\n");
%!     fputs(fd, "                elasticity equation scale, dt / cr,\n");
%!     fputs(fd, "                output pressure, yes,\n");
%!     fputs(fd, "                output stress, yes,\n");
%!     fputs(fd, "                output density, yes,\n");
%!     fputs(fd, "                output friction loss, yes,\n");
%!     fputs(fd, "                output clearance, yes,\n");
%!     fputs(fd, "                output reaction force, yes,\n");
%!     fputs(fd, "                output mesh, yes,\n");
%!     fputs(fd, "                output total deformation, yes,\n");
%!     fputs(fd, "                output, yes;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   mbdyn_solver_run([output_file, "_mbd.mbdyn"], options_mbdyn);

%!   res.log_dat = mbdyn_post_load_log(options_mbdyn.output_file);

%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbdyn.output_file);

%!   [res.elem_id, res.q, res.qdot, res.qddot] = mbdyn_post_load_output_mod(options_mbdyn.output_file, numel(res.t));

%!   res.bearings = mbdyn_post_ehd_load_output(options_mbdyn.output_file, res.log_dat);

%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);

%!   res.elem_idx_cms = find(res.elem_id == res.log_dat.vars.elem_id_diaphragm_cms);
%!   res.elem_idx_bearing = find(res.log_dat.vars.elem_id_bearing == [res.log_dat.bearings.label]);

%!   cms_data.mesh = mesh;
%!   cms_data.dof_map = dof_map;
%!   cms_data.cms_opt = cms_opt;
%!   res.sol_dyn = fem_post_cms_sol_import(options_mbdyn.output_file, cms_data);

%!   figure("visible", "off");
%!   plot(res.t * SI_unit_second, 1e-3 * max(res.bearings.columns.p_n, [], 2) * SI_unit_pascal, "-;max(p(t));1");
%!   xlabel("t [s]");
%!   ylabel("p [kPa]");
%!   grid on;
%!   grid minor on;
%!   title("convergence history of peak pressure");

%!   figure("visible", "off");
%!   plot(res.t * SI_unit_second, 1e6 * max(res.bearings.columns.wtot_n, [], 2) * SI_unit_pascal, "-;max(w(t));1");
%!   xlabel("t [s]");
%!   ylabel("w [um]");
%!   grid on;
%!   grid minor on;
%!   title("convergence history of peak deformation");

%!   figure("visible", "off");
%!   Phi = res.bearings.xi(1, :) / (0.5 * param.Di);
%!   p = interp2(res.bearings.xi, res.bearings.zi, res.bearings.columns.p(:, :, end), res.bearings.xi(1, :), 0);
%!   w = interp2(res.bearings.xi, res.bearings.zi, res.bearings.columns.wtot(:, :, end), res.bearings.xi(1, :), 0);
%!   ax = plotyy(180 / pi * Phi, 1e-3 * p * SI_unit_pascal, 180 / pi * Phi, 1e6 * w * SI_unit_meter);
%!   xlabel("Phi [deg]");
%!   ylabel(ax(1), "p [kPa]");
%!   ylabel(ax(2), "w [um]");
%!   grid on;
%!   grid minor on;
%!   for i=1:2
%!     xlim(ax(i), [0, 360]);
%!   endfor
%!   xticks(0:30:360);
%!   title("midplane pressure and deformation");

%!   figure("visible","off");
%!   hold on;
%!   set(plot(180 / pi * Phi, 1e6 * w * SI_unit_meter, "-;w(Phi) [um];1"), "linewidth", 5);
%!   set(plot(ref_data_w(:, 1), ref_data_w(:, 2), "-;reference w(Phi) [um];0"), "linewidth", 3);
%!   xlabel("Phi [deg]");
%!   ylabel("w [um]");
%!   grid on;
%!   grid minor on;
%!   xlim([0,360]);
%!   xticks(0:30:360);
%!   title("midplane deformation");

%!   figure("visible","off");
%!   hold on;
%!   set(plot(180 / pi * Phi, 1e-3 * p * SI_unit_pascal, "-;p(Phi) [kPa];1"), "linewidth", 5);
%!   set(plot(ref_p(:, 1), ref_p(:, 2), "-;reference p(Phi) [kPa];0"), "linewidth", 3);
%!   xlabel("Phi [deg]");
%!   ylabel("p* [kPa]");
%!   grid on;
%!   grid minor on;
%!   xlim([0,360]);
%!   xticks(0:30:360);
%!   title("midplane pressure");

%!   figure("visible","off");
%!   set(plot(180 / pi * Phi, p * param.cr^2 / (0.5 * param.Di * param.U1 * param.etal), "-;p*(Phi) [1];1"), "linewidth", 5);
%!   xlabel("Phi [deg]");
%!   ylabel("p* [1]");
%!   grid on;
%!   grid minor on;
%!   xlim([0,360]);
%!   xticks(0:30:360);
%!   title("midplane nondimensional pressure");

%!   figure("visible", "off");
%!   hold on;
%!   h = interp2(res.bearings.xi, res.bearings.zi, res.bearings.columns.h(:, :, end), res.bearings.xi(1, :), 0);
%!   set(plot(180 / pi * Phi, 1e6 * h * SI_unit_meter, "-;h(Phi) [um];1"), "linewidth", 5);
%!   set(plot(180 / pi * Phi, 1e6 * w * SI_unit_meter, "-;w(Phi) [um];3"), "linewidth", 5);
%!   xlabel("Phi [deg]");
%!   ylabel(ax(2), "h, w [um]");
%!   grid on;
%!   grid minor on;
%!   xlim([0, 360]);
%!   title("midplane clearance and deformation");

%!   figure("visible","off");
%!   contourf(180/pi * res.bearings.xi / (0.5 * param.Di), 1e3 * res.bearings.zi * SI_unit_meter, 1e-3 * res.bearings.columns.p(:, :, end) * SI_unit_pascal);
%!   colormap jet;
%!   colorbar;
%!   xlabel("Phi [deg]");
%!   ylabel("z [mm]");
%!   grid on;
%!   grid minor on;
%!   title("pressure distribution p [kPa]");

%!   figure("visible","off");
%!   contourf(180/pi * res.bearings.xi / (0.5 * param.Di), 1e3 * res.bearings.zi * SI_unit_meter, 1e6 * res.bearings.columns.wtot(:, :, end) * SI_unit_meter);
%!   colormap jet;
%!   colorbar;
%!   xlabel("Phi [deg]");
%!   ylabel("z [mm]");
%!   grid on;
%!   grid minor on;
%!   title("radial deformation [um]");

%!   figure("visible","off");
%!   contourf(180/pi * res.bearings.xi / (0.5 * param.Di), 1e3 * res.bearings.zi * SI_unit_meter, res.bearings.columns.rho(:, :, end) / param.rhol);
%!   colormap jet;
%!   colorbar;
%!   xlabel("Phi [deg]");
%!   ylabel("z [mm]");
%!   grid on;
%!   grid minor on;
%!   title("fractional film content [1]");

%!   figure_list();
%!   p_int = interp1(180 / pi * Phi, 1e-3 * p * SI_unit_pascal, ref_p(:, 1), "linear");
%!   w_int = interp1(180 / pi * Phi, 1e6 * w * SI_unit_meter, ref_data_w(:, 1), "linear");
%!   assert_simple(p_int, ref_p(:, 2), 0.15 * max(abs(ref_p(:, 2))));
%!   assert_simple(mean(abs(p_int - ref_p(:, 2))) < 0.05 * max(abs(ref_p(:, 2))));
%!   ## Don't check the first point because of issues related to the interpolation of the compliance matrix near to the slot
%!   assert_simple(w_int(2:end), ref_data_w(2:end, 2), 0.07 * max(abs(ref_data_w(:, 2))));
%!   assert_simple(mean(abs(w_int(2:end) - ref_data_w(2:end, 2))) < 0.04 * max(abs(ref_data_w(:, 2))));
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       status = unlink(fullfile(fn(i).folder, fn(i).name));
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
