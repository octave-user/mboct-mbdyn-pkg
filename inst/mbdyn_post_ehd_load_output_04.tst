## mbdyn_post_ehd_load_output.tst:04
%!test
%! try
%! pkg load mbdyn_util_oct;
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
%! f_plot = false;
%! fd = -1;
%! output_file = "";
%! %unwind_protect
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
%! for i=1:columns(beta)
%!   beta(beta(:, i) < -0.9 * pi, i) += 2 * pi;
%! endfor
%! mu_r = zeros(numel(epsilon_r), numel(B_d_r));
%! beta_r = zeros(numel(epsilon_r), numel(B_d_r));
%! beta_r(find(epsilon_r < 0), :) = pi;
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
%!   title(sprintf("Sommerfeld number for pure displacement B/d=%.2f", B_d_r(i)));
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
%!   title(sprintf("Angular offset of reaction force for pure displacement B/d=%.2f", B_d_r(i)));
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
%!   title(sprintf("Coefficient of friction for pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), Q_r(1:test_freq:end, i), "-;Q_r(epsilon);k");
%!   plot(epsilon_r(1:test_freq:end), Q(1:test_freq:end, i), "-;Q(epsilon);r");
%!   ylim([0, ylim()(2)]);
%!   xlabel("epsilon [1]");
%!   ylabel("Q [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Non-dimensional oil flow pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! endif
%! assert_simple(mean(mean(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.02);
%! assert_simple(mean(mean(abs(beta(1:test_freq:end,1:test_freq:end) - beta_r(1:test_freq:end, 1:test_freq:end)))) < 1e-6 * pi / 180);
%! assert_simple(mean(mean(abs(mu(1:test_freq:end,1:test_freq:end)))) < 1e-8);
%! assert_simple(mean(mean(abs(Q(1:test_freq:end,1:test_freq:end) ./ Q_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.06);
%! assert_simple(max(max(abs(So(1:test_freq:end,1:test_freq:end) ./ So_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.11);
%! assert_simple(max(max(abs(beta(1:test_freq:end,1:test_freq:end) - beta_r(1:test_freq:end,1:test_freq:end)))) < 1e-6 * pi / 180);
%! assert_simple(max(max(abs(mu(1:test_freq:end,1:test_freq:end)))) < 1e-8);
%! assert_simple(max(max(abs(Q(1:test_freq:end,1:test_freq:end) ./ Q_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.07);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
