## mbdyn_post_ehd_load_output.tst:07
%!test
%! try
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
%! f_plot = false;
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
%! if (f_plot)
%! figure("visible", "off");
%! hold on;
%! plot(res.bearings.zi(:, 1) + 0.5 * param.L_x1, res.bearings.columns.p(:, 1), "-x;p(z);r");
%! plot(analy.x1, analy.p_wc, "-;Fowell;k");
%! xlabel("z [m]");
%! ylabel("p [Pa]");
%! grid on;
%! grid minor on;
%! title("axial pressure distribution");
%! figure("visible", "off");
%! hold on;
%! plot(res.bearings.zi(:, 1) + 0.5 * param.L_x1, res.bearings.columns.h(:, 1), "-+;h(z);r");
%! plot([0, param.L_x1], [param.h_max, param.h_min], "-;h_r(z);k");
%! xlabel("z [m]");
%! ylabel("h [m]");
%! grid on;
%! grid minor on;
%! title("radial clearance versus axial position");
%! endif
%! tol_q = 1e-3;
%! tol_p = 6e-2;
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
