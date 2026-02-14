## mbdyn_post_ehd_load_output.tst:01
%!test
%! try
%! pkg load mbdyn_util_oct;
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
%! f_plot = false;
%! for j=1:test_freq:numel(epsilon_r)
%! for k=1:test_freq:numel(B_d_r)
%! param.d1 = 10e-3;
%! param.Psi = 1e-3;
%! param.M = int32(20);
%! param.N = int32(85);
%! param.output_bearing_data = true;
%! param.B = param.d1 * B_d_r(k);
%! param.epsilon = epsilon_r(j);
%! param.cavitation_model = "mass conserving";
%! param.hydraulic_nodes = true;
%! param.jacobian_check = false;
%! fd = -1;
%! output_file = "";
%! %unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "oct-mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     if (param.hydraulic_nodes)
%!       fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!       fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!       fputs(fd, "set: integer hyd_node_id_inlet = 2003;\n");
%!     endif
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     if (param.hydraulic_nodes)
%!       fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!       fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!       fputs(fd, "set: integer genel_id_inlet = 4005;\n");
%!     endif
%!     fputs(fd, "set: real omega1z = 2 * pi * 50;\n");
%!     fputs(fd, "set: real omega2z = 2 * pi * 0;\n");
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
%!     fputs(fd, "        derivatives coefficient: 1e-3, auto;\n");
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
%!     if (param.jacobian_check)
%!       fputs(fd, "    finite difference jacobian meter: const, 1., iterations, string, \"Var == 1\", coefficient, 1e-8, order, 1, output, none, statistics iteration, yes;\n");
%!     endif
%!     fputs(fd, "    use automatic differentiation;\n");
%!     fputs(fd, "    skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "        structural nodes: 2;\n");
%!     fputs(fd, "        joints: 2;\n");
%!     fputs(fd, "        loadable elements: 1;\n");
%!     if (param.hydraulic_nodes)
%!       fputs(fd, "        hydraulic nodes: 3;\n");
%!       fputs(fd, "        genels: 3;\n");
%!     endif
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
%!     if (param.hydraulic_nodes)
%!       fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!       fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!       fputs(fd, "    hydraulic: hyd_node_id_inlet, p_in;\n");
%!     endif
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
%!     if (param.hydraulic_nodes)
%!       fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!       fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!       fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!       fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!       fputs(fd, "    genel: genel_id_inlet, clamp,\n");
%!       fputs(fd, "        hyd_node_id_inlet, hydraulic, p_in;\n");
%!     endif
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     switch(param.cavitation_model)
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
%!     if (param.hydraulic_nodes)
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure coupling conditions radial, 1,\n");
%!     fputs(fd, "                position, 0.5 * D2 * (delta_t0 + pi), 0.,\n");
%!     fputs(fd, "                rectangle, width, D2 * pi / (N - 1), height, B,\n");
%!     fputs(fd, "                pressure node, hyd_node_id_inlet,\n");
%!     else
%!       fputs(fd, "         boundary conditions, const, p_mB2, const, p_pB2,\n");
%!       fputs(fd, "         lubrication grooves, 1,\n");
%!       fputs(fd, "                at bearing,\n");
%!       fputs(fd, "                pressure, const, p_in,\n");
%!       fputs(fd, "                position, 0.5 * D2 * (delta_t0 + pi), 0.,\n");
%!       fputs(fd, "                rectangle, width, D2 * pi / (N - 1), height, B,\n");
%!     endif
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
%! # opt_sol.logfile = [output_file, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   if (param.hydraulic_nodes)
%!     [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   endif
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
%!   M2 = -res.bearings.columns.M2.';
%!   omega1z = res.bearings.cylindrical.omega1z;
%!   omega2z = res.bearings.cylindrical.omega2z;
%!   omega_res = res.bearings.cylindrical.omega_res;
%!   epsilon = res.bearings.cylindrical.epsilon;
%!   delta = res.bearings.cylindrical.delta;
%!   if (param.hydraulic_nodes)
%!     Q_out1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end) / res.log_dat.vars.rho;
%!     Q_out2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end) / res.log_dat.vars.rho;
%!     Q_in = -res.genel_data{res.log_dat.vars.genel_idx_inlet}(end) / res.log_dat.vars.rho;
%!   endif
%!   P = norm(Rb2(:, 1:2).' * R2.' * F1);
%!   alpha = atan2(Rb2(:, 2).' * R2.' * F1, Rb2(:, 1).' * R2.' * F1);
%!   beta(j, k) = sign(omega_res) * (delta - alpha);
%!   beta(j, k) = mod(beta(j, k) + pi, 2 * pi) - pi;
%!   mu(j, k) = 2 * Rb2(:, 3).' * R2.' * M1 / (P * d) * sign(omega1z - omega2z);
%!   Pf1(j, k) = -M1(3) * (omega1z - omega2z);
%!   Pf2(j, k) = M2(3) * (omega1z - omega2z);
%!   Pf3(j, k) = res.bearings.columns.Pff + res.bearings.columns.Pfc;
%!   So(j, k) = P * Psi^2 / (B * D * eta * abs(omega_res));
%!   if (param.hydraulic_nodes)
%!     Q(j, k) = (Q_out1 + Q_out2) / ((0.5 * D)^3 * Psi * abs(omega_res));
%!     dQ(j, k) = (Q_out1 + Q_out2 + Q_in) / Q_in;
%!   endif
%!   mu_r(j, k) = Psi * (abs((omega1z - omega2z) / omega_res) * pi / (sqrt(1 - epsilon^2) * So_r(j, k)) + sin(beta_r(j, k)) * abs(epsilon) / 2);
%! %unwind_protect_cleanup
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
%! %end_unwind_protect
%! endfor
%! endfor
%! if (f_plot)
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   semilogy(epsilon_r(1:test_freq:end), So_r(1:test_freq:end, i), "-;So_r(epsilon);k");
%!   semilogy(epsilon_r(1:test_freq:end), So(1:test_freq:end, i), "-;So(epsilon);r");
%!   xlabel("epsilon [1]");
%!   ylabel("So [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Sommerfeld number for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta_r(1:test_freq:end, i), "-;beta_r(epsilon);k");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta(1:test_freq:end, i), "-;beta(epsilon);r");
%!   xlabel("epsilon [1]");
%!   ylabel("beta [deg]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Angular offset of reaction force for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), mu_r(1:test_freq:end, i), "-;mu_r(epsilon);k");
%!   plot(epsilon_r(1:test_freq:end), mu(1:test_freq:end, i), "-;mu(epsilon);r");
%!   xlabel("epsilon [1]");
%!   ylabel("mu [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Coefficient of friction for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! if (param.hydraulic_nodes)
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), Q_r(1:test_freq:end, i), "-;Q_r(epsilon);k");
%!   plot(epsilon_r(1:test_freq:end), Q(1:test_freq:end, i), "-;Q(epsilon);r");
%!   xlabel("epsilon [1]");
%!   ylabel("Q [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Non-dimensional oil flow pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! endif
%! endif
%! assert_simple(mean(mean(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.03);
%! assert_simple(mean(mean(abs(beta(1:test_freq:end, 1:test_freq:end) ./ beta_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.02);
%! assert_simple(mean(mean(abs(mu(1:test_freq:end, 1:test_freq:end) ./ mu_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.02);
%! if (param.hydraulic_nodes)
%!   assert_simple(mean(mean(abs(Q(1:test_freq:end, 1:test_freq:end) ./ Q_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.07);
%! endif
%! assert_simple(max(max(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.05);
%! assert_simple(max(max(abs(beta(1:test_freq:end, 1:test_freq:end) ./ beta_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.04);
%! assert_simple(max(max(abs(mu(1:test_freq:end, 1:test_freq:end) ./ mu_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.05);
%! if (param.hydraulic_nodes)
%!   assert_simple(max(max(abs(Q(1:test_freq:end, 1:test_freq:end) ./ Q_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.08);
%! endif
%! assert_simple(max(max(abs(Pf2 - Pf3))) / max(max(abs(Pf2))) < 1e-2);
%! #assert_simple(max(max(abs(Pf1 - Pf3))) / max(max(abs(Pf1))) < 1e-2);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
