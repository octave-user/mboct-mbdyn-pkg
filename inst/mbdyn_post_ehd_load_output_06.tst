## mbdyn_post_ehd_load_output.tst:06
%!test
%! try
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
%!   if (f_plot)
%!   figure("visible", "off");
%!   hold on;
%!   plot(res.bearings.zi(:, 1) + 0.5 * B, res.bearings.columns.p(:, 1), "-;p(z);r");
%!   plot(z, p_ref, "-;p_r;k");
%!   xlabel("z [m]");
%!   ylabel("p [Pa]");
%!   grid on;
%!   grid minor on;
%!   title("axial pressure distribution");
%!   figure("visible", "off");
%!   hold on;
%!   plot(res.bearings.zi(:, 1) + 0.5 * B, res.bearings.columns.h(:, 1), "-;h(z);r");
%!   plot(z, h0 + dh_dz * z, "-;h_r(z);k");
%!   xlabel("z [m]");
%!   ylabel("h [m]");
%!   grid on;
%!   grid minor on;
%!   title("radial clearance versus time");
%!   endif
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
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
