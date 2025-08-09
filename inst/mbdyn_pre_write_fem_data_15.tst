## mbdyn_pre_write_fem_data.tst:15
%!test
%! try
%! ## TEST 15
%! pkg load mboct-fem-pkg;
%! options.f_plot = false;
%! if (options.f_plot)
%!   close all;
%! endif
%! function R2j = norm_ref_frame(R2)
%!   e1 = R2(:, 1);
%!   e3 = [0; 0; 1];
%!   e2 = cross(e3, e1);
%!   e1 = cross(e2, e3);
%!   R2j = [e1 / norm(e1), e2 / norm(e2), e3 / norm(e3)];
%! endfunction
%! function [omega, r] = rotordynamics_test_case(param, options, SI_unit)
%!   cms_opt.algorithm = "shift-invert";
%!   cms_opt.refine_max_iter = int32(3);
%!   cms_opt.element.name = "elem_id_rotor";
%!   cms_opt.nodes.modal.name = "node_id_rotor";
%!   cms_opt.nodes.interfaces(1).name = "node_id_bearing1";
%!   cms_opt.nodes.interfaces(2).name = "node_id_bearing2";
%!   cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!   cms_opt.verbose = options.verbose;
%!   enable_filenames = [options.f_enable_modal, options.f_enable_beam, options.f_enable_solid];
%!   filename = "";
%! %unwind_protect
%!     filename = tempname();
%!     if (ispc())
%!       filename(filename == "\\") = "/";
%!     endif
%!     mbdyn_filename_suffix = {"cms", "beam", "solid"};
%!     for i=1:numel(mbdyn_filename_suffix)
%!       mbdyn_filenames{i} = [filename, mbdyn_filename_suffix{i}, ".mbdyn"];
%!     endfor
%!     fd = -1;
%!     if (options.f_enable_modal)
%!       unwind_protect
%!         fd = fopen(mbdyn_filenames{1}, "w");
%!         if (fd == -1)
%!           error("failed to open file \"%s\"", mbdyn_filenames{1});
%!         endif
%!         fputs(fd, "include: \"${MBDYN_ROTOR_DYN_CMS_PARAM_FILE}\";\n");
%!         fputs(fd, "set: integer ref_id_ground = 1001;\n");
%!         fputs(fd, "set: integer ref_id_shaft = 1002;\n");
%!         fputs(fd, "set: integer ref_id_rotor = 1003;\n");
%!         fputs(fd, "set: integer ref_id_bearing1 = 1004;\n");
%!         fputs(fd, "set: integer ref_id_bearing2 = 1005;\n");
%!         fputs(fd, "set: integer node_id_rotor = 2002;\n");
%!         fputs(fd, "set: integer node_id_bearing1 = 2003;\n");
%!         fputs(fd, "set: integer node_id_bearing2 = 2004;\n");
%!         fputs(fd, "set: integer body_id_unbalance = 3000;\n");
%!         fputs(fd, "set: integer elem_id_rotor = 3005;\n");
%!         fputs(fd, "set: integer joint_id_bearing1 = 3006;\n");
%!         fputs(fd, "set: integer joint_id_bearing2 = 3007;\n");
%!         fputs(fd, "set: integer joint_id_drive = 3008;\n");
%!         fputs(fd, "set: integer elem_id_inertia = 3010;\n");
%!         fputs(fd, "set: integer drive_id_rotor_speed = 5001;\n");
%!         fputs(fd, "set: integer drive_id_time_step = 5002;\n");
%!         fputs(fd, "set: real initial_time = 0.;\n");
%!         fputs(fd, "set: real final_time = 2. * pi * n / (abs(omega1 + omega0) / 2.);\n");
%!         fputs(fd, "begin: data;\n");
%!         fputs(fd, "        problem: initial value;\n");
%!         fputs(fd, "end: data;\n");
%!         fputs(fd, "begin: initial value;\n");
%!         fputs(fd, "        initial time: initial_time;\n");
%!         fputs(fd, "        final time: final_time;\n");
%!         fputs(fd, "        time step: 2. * pi / (360. * abs(omega0));\n");
%!         fputs(fd, "        strategy: change, postponed, drive_id_time_step;\n");
%!         fputs(fd, "        method: ms, 0.6;\n");
%!         fputs(fd, "        tolerance: 1e-6, test, norm, 1e-6, test, norm;\n");
%!         fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!         fputs(fd, "        max iterations: 100;\n");
%!         fputs(fd, "        derivatives max iterations: 10;\n");
%!         fputs(fd, "        derivatives coefficient: auto;\n");
%!         fputs(fd, "        output: iterations;\n");
%!         fputs(fd, "        linear solver: naive, colamd, scale, row max column max, always, max iterations, 100;\n");
%!         fputs(fd, "        nonlinear solver: nox, modified, 10,\n");
%!         fputs(fd, "             keep jacobian matrix,\n");
%!         fputs(fd, "             inner iterations before assembly, 6,\n");
%!         fputs(fd, "             jacobian operator, newton krylov,\n");
%!         fputs(fd, "             solver, line search based,\n");
%!         fputs(fd, "             forcing term, type 2,\n");
%!         fputs(fd, "             direction, newton,\n");
%!         fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!         fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!         fputs(fd, "             linear solver, gmres,\n");
%!         fputs(fd, "             linear solver max iterations, 12,\n");
%!         fputs(fd, "             krylov subspace size, 12;\n");
%!         fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!         fputs(fd, "        threads: assembly, 1;\n");
%!         fputs(fd, "end: initial value;\n");
%!         fputs(fd, "begin: control data;\n");
%!         fputs(fd, "       default output: none, structural nodes;\n");
%!         fputs(fd, "       default orientation: euler123;\n");
%!         fputs(fd, "       output precision: 16;\n");
%!         fputs(fd, "       max iterations: 0;\n");
%!         fputs(fd, "        structural nodes:\n");
%!         fputs(fd, "                +1           # modal\n");
%!         fputs(fd, "                +1           # interface 1\n");
%!         fputs(fd, "                +1              # interface 2\n");
%!         fputs(fd, "        ;\n");
%!         fputs(fd, "        joints:\n");
%!         fputs(fd, "                +1           # modal\n");
%!         fputs(fd, "                +1           # bearing1\n");
%!         fputs(fd, "                +1           # bearing2\n");
%!         fputs(fd, "                +1              # drive\n");
%!         fputs(fd, "        ;\n");
%!         fputs(fd, "        rigid bodies: 1;\n");
%!         fputs(fd, "        gravity;\n");
%!         fputs(fd, "end: control data;\n");
%!         fputs(fd, "reference: ref_id_ground,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, eye,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, null;\n");
%!         fputs(fd, "reference: ref_id_shaft,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         fputs(fd, "        reference, ref_id_ground, eye,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         fputs(fd, "        reference, ref_id_ground, 0., 0., omega0;\n");
%!         fputs(fd, "reference: ref_id_rotor,\n");
%!         fputs(fd, "        reference, ref_id_shaft, 0., 0., o,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing1,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   -0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing2,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "begin: nodes;\n");
%!         fputs(fd, "        structural: node_id_rotor, modal,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null;\n");
%!         fputs(fd, "        structural: node_id_bearing1, static,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null;\n");
%!         fputs(fd, "        structural: node_id_bearing2, static,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null;\n");
%!         fputs(fd, "end: nodes;\n");
%!         fputs(fd, "begin: elements;\n");
%!         fputs(fd, "       drive caller: drive_id_rotor_speed, string, \"(omega1 - omega0) / (final_time - initial_time) * Time + omega0\", output, yes;\n");
%!         fputs(fd, "       drive caller: drive_id_time_step, string, \"2. * pi / (max(1., 360. * abs(model::drive(drive_id_rotor_speed, Time))))\";\n");
%!         fputs(fd, "       body: body_id_unbalance,\n");
%!         fputs(fd, "             node_id_rotor,\n");
%!         fputs(fd, "                dm,\n");
%!         fputs(fd, "                reference, ref_id_rotor, dr, 0., 0.,\n");
%!         fputs(fd, "                diag, 0., 0., 0.;\n");
%!         fputs(fd, "        joint: joint_id_bearing1, total pin joint,\n");
%!         fputs(fd, "                node_id_bearing1,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, active,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "        joint: joint_id_bearing2, total pin joint,\n");
%!         fputs(fd, "               node_id_bearing2,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "               position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, inactive,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "     joint: joint_id_drive, angular velocity,\n");
%!         fputs(fd, "             # node label\n");
%!         fputs(fd, "             node_id_rotor, \n");
%!         fputs(fd, "             # direction\n");
%!         fputs(fd, "             0.,0.,1.,\n");
%!         fputs(fd, "             # angular velocity\n");
%!         fputs(fd, "             reference, drive_id_rotor_speed;\n");
%!         fputs(fd, "        include: \"${MBDYN_ROTOR_DYN_CMS_ELEM_FILE}\";\n");
%!         fputs(fd, "        inertia: elem_id_inertia,\n");
%!         fputs(fd, "                 position, reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                 orientation, reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                 body, all, joint, all, loadable, all,\n");
%!         fputs(fd, "                 output, both;\n");
%!         fputs(fd, "        gravity: uniform, 0., 0., -1., g;\n");
%!         fputs(fd, "end: elements;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!     endif
%!     fd = -1;
%!     if (options.f_enable_beam)
%!       unwind_protect
%!         fd = fopen(mbdyn_filenames{2}, "w");
%!         if (fd == -1)
%!           error("failed to open file \"%s\"", mbdyn_filenames{2});
%!         endif
%!         fputs(fd, "include: \"${MBDYN_ROTOR_DYN_CMS_PARAM_FILE}\";\n");
%!         fputs(fd, "set: integer ref_id_ground = 1001;\n");
%!         fputs(fd, "set: integer ref_id_shaft = 1002;\n");
%!         fputs(fd, "set: integer ref_id_rotor = 1003;\n");
%!         fputs(fd, "set: integer ref_id_bearing1 = 1004;\n");
%!         fputs(fd, "set: integer ref_id_bearing2 = 1005;\n");
%!         fputs(fd, "set: integer ref_id_shaft_section = 1006;\n");
%!         fputs(fd, "set: integer ref_id_bearing1_center = 1007;\n");
%!         fputs(fd, "set: integer ref_id_bearing2_center = 1008;\n");
%!         fputs(fd, "set: integer node_id_rotor = 2001;\n");
%!         fputs(fd, "set: integer node_id_bearing1 = 2002;\n");
%!         fputs(fd, "set: integer node_id_bearing1_center = 2003;\n");
%!         fputs(fd, "set: integer node_id_bearing2 = 2004;\n");
%!         fputs(fd, "set: integer node_id_bearing2_center = 2005;\n");
%!         fputs(fd, "set: integer body_id_rotor = 3001;\n");
%!         fputs(fd, "set: integer body_id_bearing1 = 3002;\n");
%!         fputs(fd, "set: integer body_id_bearing1_center = 3003;\n");
%!         fputs(fd, "set: integer body_id_bearing2 = 3004;\n");
%!         fputs(fd, "set: integer body_id_bearing2_center = 3005;\n");
%!         fputs(fd, "set: integer beam_id_shaft1 = 4001;\n");
%!         fputs(fd, "set: integer beam_id_shaft2 = 4002;\n");
%!         fputs(fd, "set: integer joint_id_bearing1 = 3006;\n");
%!         fputs(fd, "set: integer joint_id_bearing2 = 3007;\n");
%!         fputs(fd, "set: integer joint_id_drive = 3008;\n");
%!         fputs(fd, "set: integer elem_id_inertia = 4001;\n");
%!         fputs(fd, "set: integer drive_id_rotor_speed = 5001;\n");
%!         fputs(fd, "set: integer drive_id_time_step = 5002;\n");
%!         fputs(fd, "set: real l1 = 0.5 * l + o - 0.5 * w;\n");
%!         fputs(fd, "set: real l2 = 0.5 * l - o - 0.5 * w;\n");
%!         fputs(fd, "set: real rotor_m = rho * D^2 * pi / 4. * w;\n");
%!         fputs(fd, "set: real rotor_Jx = rotor_m * ((0.5 * D)^2 + w^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real rotor_Jy = rotor_Jx;\n");
%!         fputs(fd, "set: real rotor_Jz = rotor_m * (0.5 * D)^2 / 2.;\n");
%!         fputs(fd, "set: real shaft_rotor1_m = rho * d^2 * pi / 4. * (l1 / 4.);\n");
%!         fputs(fd, "set: real shaft_rotor1_Jx = shaft_rotor1_m * ((0.5 * d)^2 + (l1 / 4.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_rotor1_Jy = shaft_rotor1_Jx;\n");
%!         fputs(fd, "set: real shaft_rotor1_Jz = shaft_rotor1_m * (0.5 * d)^2 / 2.;\n");
%!         fputs(fd, "set: real shaft_rotor2_m = rho * d^2 * pi / 4. * (l2 / 4.);\n");
%!         fputs(fd, "set: real shaft_rotor2_Jx = shaft_rotor2_m * ((0.5 * d)^2 + (l2 / 4.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_rotor2_Jy = shaft_rotor2_Jx;\n");
%!         fputs(fd, "set: real shaft_rotor2_Jz = shaft_rotor2_m * (0.5 * d)^2 / 2.;\n");
%!         fputs(fd, "set: real shaft_bearing1_m = rho * d^2 * pi / 4. * (l1 / 4.);\n");
%!         fputs(fd, "set: real shaft_bearing1_Jx = shaft_bearing1_m * ((0.5 * d)^2 + (l1 / 4.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_bearing1_Jy = shaft_bearing1_Jx;\n");
%!         fputs(fd, "set: real shaft_bearing1_Jz = shaft_bearing1_m * (0.5 * d)^2 / 2.;\n");
%!         fputs(fd, "set: real shaft_bearing1_center_m = rho * d^2 * pi / 4. * (l1 / 2.);\n");
%!         fputs(fd, "set: real shaft_bearing1_center_Jx = shaft_bearing1_center_m * ((0.5 * d)^2 + (l1 / 2.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_bearing1_center_Jy = shaft_bearing1_center_Jx;\n");
%!         fputs(fd, "set: real shaft_bearing1_center_Jz = shaft_bearing1_center_m * ((0.5 * d)^2) / 2.;\n");
%!         fputs(fd, "set: real shaft_bearing2_m = rho * d^2 * pi / 4. * (l2 / 4.);\n");
%!         fputs(fd, "set: real shaft_bearing2_Jx = shaft_bearing2_m * ((0.5 * d)^2 + (l2 / 4.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_bearing2_Jy = shaft_bearing2_Jx;\n");
%!         fputs(fd, "set: real shaft_bearing2_Jz = shaft_bearing2_m * (0.5 * d)^2 / 2.;\n");
%!         fputs(fd, "set: real shaft_bearing2_center_m = rho * d^2 * pi / 4. * (l2 / 2.);\n");
%!         fputs(fd, "set: real shaft_bearing2_center_Jx = shaft_bearing2_center_m * ((0.5 * d)^2 + (l2 / 2.)^2 / 3.) / 4.;\n");
%!         fputs(fd, "set: real shaft_bearing2_center_Jy = shaft_bearing2_center_Jx;\n");
%!         fputs(fd, "set: real shaft_bearing2_center_Jz = shaft_bearing2_center_m * ((0.5 * d)^2) / 2.;\n");
%!         fputs(fd, "set: real shaft_A = d^2 * pi / 4;\n");
%!         fputs(fd, "set: real shaft_As = 9. / 10. * shaft_A;\n");
%!         fputs(fd, "set: real shaft_Iy = d^4 * pi / 64.;\n");
%!         fputs(fd, "set: real shaft_Iz = shaft_Iy;\n");
%!         fputs(fd, "set: real shaft_Ip = shaft_Iy + shaft_Iz;\n");
%!         fputs(fd, "set: real shaft_It = shaft_Ip;\n");
%!         fputs(fd, "set: real shaft_E = E;\n");
%!         fputs(fd, "set: real shaft_nu = nu;\n");
%!         fputs(fd, "set: real shaft_G = shaft_E / (2 * (1 + shaft_nu));\n");
%!         fputs(fd, "set: real shaft_rho = rho;\n");
%!         fputs(fd, "set: real shaft_damping_ratio = beta;\n");
%!         fputs(fd, "set: real initial_time = 0.;\n");
%!         fputs(fd, "set: real final_time = 2. * pi * n / (abs(omega1 + omega0) / 2.);\n");
%!         fputs(fd, "begin: data;\n");
%!         fputs(fd, "        problem: initial value;\n");
%!         fputs(fd, "end: data;\n");
%!         fputs(fd, "begin: initial value;\n");
%!         fputs(fd, "        initial time: initial_time;\n");
%!         fputs(fd, "        final time: final_time;\n");
%!         fputs(fd, "        time step: 2. * pi / (360. * abs(omega0));\n");
%!         fputs(fd, "        strategy: change, postponed, drive_id_time_step;\n");
%!         fputs(fd, "        method: ms, 0.6;\n");
%!         fputs(fd, "        tolerance: 1e-6, test, norm, 1e-6, test, norm;\n");
%!         fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!         fputs(fd, "        max iterations: 100;\n");
%!         fputs(fd, "        derivatives max iterations: 10;\n");
%!         fputs(fd, "        derivatives coefficient: auto;\n");
%!         fputs(fd, "        output: iterations;\n");
%!         fputs(fd, "        linear solver: naive, colamd, scale, row max column max, always, max iterations, 100;\n");
%!         fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!         fputs(fd, "        threads: assembly, 1;\n");
%!         fputs(fd, "        nonlinear solver: nox, modified, 10,\n");
%!         fputs(fd, "             keep jacobian matrix,\n");
%!         fputs(fd, "             inner iterations before assembly, 6,\n");
%!         fputs(fd, "             jacobian operator, newton krylov,\n");
%!         fputs(fd, "             solver, line search based,\n");
%!         fputs(fd, "             forcing term, type 2,\n");
%!         fputs(fd, "             direction, newton,\n");
%!         fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!         fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!         fputs(fd, "             linear solver, gmres,\n");
%!         fputs(fd, "             linear solver max iterations, 12,\n");
%!         fputs(fd, "             krylov subspace size, 12;\n");
%!         fputs(fd, "end: initial value;\n");
%!         fputs(fd, "begin: control data;\n");
%!         fputs(fd, "       use automatic differentiation;\n");
%!         fputs(fd, "       default output: none, structural nodes;\n");
%!         fputs(fd, "       default orientation: euler123;\n");
%!         fputs(fd, "       output precision: 16;\n");
%!         fputs(fd, "       max iterations: 0;\n");
%!         fputs(fd, "        structural nodes: 5;\n");
%!         fputs(fd, "        joints: 3;\n");
%!         fputs(fd, "        beams: 2;\n");
%!         fputs(fd, "        rigid bodies: 5;\n");
%!         if (options.f_rbk)
%!           fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!           fputs(fd, "        component,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           array, 2,\n");
%!           fputs(fd, "             mult, time, const, (omega1 - omega0) / (final_time - initial_time),\n");
%!           fputs(fd, "             const, omega0,\n");
%!           fputs(fd, "        angular acceleration,\n");
%!           fputs(fd, "        component,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           const, (omega1 - omega0) / (final_time - initial_time);\n");
%!         endif
%!         fputs(fd, "          gravity;\n");
%!         fputs(fd, "end: control data;\n");
%!         fputs(fd, "reference: ref_id_ground,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, eye,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, null;\n");
%!         fputs(fd, "reference: ref_id_shaft,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         fputs(fd, "        reference, ref_id_ground, eye,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         if (options.f_rbk)
%!           fputs(fd, "        reference, ref_id_ground, null;\n");
%!         else
%!           fputs(fd, "        reference, ref_id_ground, 0., 0., omega0;\n");
%!         endif
%!         fputs(fd, "reference: ref_id_rotor,\n");
%!         fputs(fd, "        reference, ref_id_shaft, 0., 0., o,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_shaft_section,\n");
%!         fputs(fd, "           reference, ref_id_shaft, null,\n");
%!         fputs(fd, "           reference, ref_id_shaft, 1, 0., 0., 1.,\n");
%!         fputs(fd, "                                    2, 0., 1., 0.,\n");
%!         fputs(fd, "           reference, ref_id_shaft, null,\n");
%!         fputs(fd, "           reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing1,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   -0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing1_center,\n");
%!         fputs(fd, "        reference, ref_id_bearing1,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.5 * l1,\n");
%!         fputs(fd, "        reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "        reference, ref_id_rotor, null,\n");
%!         fputs(fd, "        reference, ref_id_rotor, null;\n");
%!         fputs(fd, "reference: ref_id_bearing2,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing2_center,\n");
%!         fputs(fd, "        reference, ref_id_bearing2,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   -0.5 * l2,\n");
%!         fputs(fd, "        reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "        reference, ref_id_rotor, null,\n");
%!         fputs(fd, "        reference, ref_id_rotor, null;\n");
%!         fputs(fd, "begin: nodes;\n");
%!         fputs(fd, "        structural: node_id_rotor, dynamic,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                reference, ref_id_rotor, null;\n");
%!         fputs(fd, "        structural: node_id_bearing1, dynamic,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1, null;\n");
%!         fputs(fd, "        structural: node_id_bearing1_center, dynamic,\n");
%!         fputs(fd, "                reference, ref_id_bearing1_center, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1_center, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing1_center, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing1_center, null;\n");
%!         fputs(fd, "        structural: node_id_bearing2, dynamic,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2, null;\n");
%!         fputs(fd, "        structural: node_id_bearing2_center, dynamic,\n");
%!         fputs(fd, "                reference, ref_id_bearing2_center, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2_center, eye,\n");
%!         fputs(fd, "                reference, ref_id_bearing2_center, null,\n");
%!         fputs(fd, "                reference, ref_id_bearing2_center, null;\n");
%!         fputs(fd, "end: nodes;\n");
%!         fputs(fd, "begin: elements;\n");
%!         fputs(fd, "       drive caller: drive_id_rotor_speed, string, \"(omega1 - omega0) / (final_time - initial_time) * Time + omega0\", output, yes;\n");
%!         fputs(fd, "       drive caller: drive_id_time_step, string, \"2. * pi / (max(1., 360. * abs(model::drive(drive_id_rotor_speed, Time))))\";\n");
%!         fputs(fd, "       body: body_id_rotor,\n");
%!         fputs(fd, "             node_id_rotor,\n");
%!         fputs(fd, "             condense, 4,\n");
%!         fputs(fd, "                rotor_m,\n");
%!         fputs(fd, "                  reference, node, null,\n");
%!         fputs(fd, "                  diag, rotor_Jx, rotor_Jy, rotor_Jz,\n");
%!         fputs(fd, "                shaft_rotor1_m,\n");
%!         fputs(fd, "                  reference, node, 0., 0., -w / 2. - l1 / 8.,\n");
%!         fputs(fd, "                  diag, shaft_rotor1_Jx, shaft_rotor1_Jy, shaft_rotor1_Jz,\n");
%!         fputs(fd, "                shaft_rotor2_m,\n");
%!         fputs(fd, "                  reference, node, 0., 0., w / 2. + l2 / 8.,\n");
%!         fputs(fd, "                  diag, shaft_rotor2_Jx, shaft_rotor2_Jy, shaft_rotor2_Jz,                \n");
%!         fputs(fd, "                dm,\n");
%!         fputs(fd, "                  reference, node, dr, 0., 0.,\n");
%!         fputs(fd, "                  diag, 0., 0., 0.;\n");
%!         fputs(fd, "        body: body_id_bearing1,\n");
%!         fputs(fd, "              node_id_bearing1,\n");
%!         fputs(fd, "              shaft_bearing1_m,\n");
%!         fputs(fd, "              reference, node, 0., 0., l1 / 8.,\n");
%!         fputs(fd, "              diag, shaft_bearing1_Jx, shaft_bearing1_Jy, shaft_bearing1_Jz;\n");
%!         fputs(fd, "        body: body_id_bearing1_center,\n");
%!         fputs(fd, "              node_id_bearing1_center,\n");
%!         fputs(fd, "              shaft_bearing1_center_m,\n");
%!         fputs(fd, "              reference, node, null,\n");
%!         fputs(fd, "              diag, shaft_bearing1_center_Jx, shaft_bearing1_center_Jy, shaft_bearing1_center_Jz;\n");
%!         fputs(fd, "        body: body_id_bearing2,\n");
%!         fputs(fd, "              node_id_bearing2,\n");
%!         fputs(fd, "              shaft_bearing2_m,\n");
%!         fputs(fd, "              reference, node, 0., 0., -l2 / 8.,\n");
%!         fputs(fd, "              diag, shaft_bearing2_Jx, shaft_bearing2_Jy, shaft_bearing2_Jz;\n");
%!         fputs(fd, "        body: body_id_bearing2_center,\n");
%!         fputs(fd, "              node_id_bearing2_center,\n");
%!         fputs(fd, "              shaft_bearing2_center_m,\n");
%!         fputs(fd, "              reference, node, null,\n");
%!         fputs(fd, "              diag, shaft_bearing2_center_Jx, shaft_bearing2_center_Jy, shaft_bearing2_center_Jz;\n");
%!         fputs(fd, "    beam3: beam_id_shaft1,\n");
%!         fputs(fd, "                # node 1\n");
%!         fputs(fd, "                node_id_bearing1, position, reference, node, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # node 2\n");
%!         fputs(fd, "                node_id_bearing1_center, position, reference, node, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # node 3,\n");
%!         fputs(fd, "                node_id_rotor, position, reference, node, 0., 0., -0.5 * w,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # orientation matrix section I\n");
%!         fputs(fd, "                reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # constitutive law section I\n");
%!         fputs(fd, "                linear viscoelastic generic,\n");
%!         fputs(fd, "                diag, shaft_E * shaft_A , shaft_G * shaft_As, shaft_G * shaft_As,\n");
%!         fputs(fd, "                      shaft_G * shaft_It, shaft_E * shaft_Iy, shaft_E * shaft_Iz,\n");
%!         fputs(fd, "                proportional, shaft_damping_ratio,\n");
%!         fputs(fd, "                # orientation matrix section II\n");
%!         fputs(fd, "                same,\n");
%!         fputs(fd, "                # constitutive law section II\n");
%!         fputs(fd, "                same;\n");
%!         fputs(fd, "    beam3: beam_id_shaft2,\n");
%!         fputs(fd, "                # node 1\n");
%!         fputs(fd, "                node_id_rotor, position, reference, node, 0., 0., 0.5 * w,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,                \n");
%!         fputs(fd, "                # node 2\n");
%!         fputs(fd, "                node_id_bearing2_center, position, reference, node, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # node 3,\n");
%!         fputs(fd, "                node_id_bearing2, position, reference, node, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # orientation matrix section I\n");
%!         fputs(fd, "                reference, ref_id_shaft_section, eye,\n");
%!         fputs(fd, "                # constitutive law section I\n");
%!         fputs(fd, "                linear viscoelastic generic,\n");
%!         fputs(fd, "                diag, shaft_E * shaft_A , shaft_G * shaft_As, shaft_G * shaft_As,\n");
%!         fputs(fd, "                      shaft_G * shaft_It, shaft_E * shaft_Iy, shaft_E * shaft_Iz,\n");
%!         fputs(fd, "                proportional, shaft_damping_ratio,\n");
%!         fputs(fd, "                # orientation matrix section II\n");
%!         fputs(fd, "                same,\n");
%!         fputs(fd, "                # constitutive law section II\n");
%!         fputs(fd, "                same;\n");
%!         fputs(fd, "        joint: joint_id_bearing1, total pin joint,\n");
%!         fputs(fd, "                node_id_bearing1,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, active,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "        joint: joint_id_bearing2, total pin joint,\n");
%!         fputs(fd, "               node_id_bearing2,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "               position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, inactive,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "     joint: joint_id_drive, angular velocity,\n");
%!         fputs(fd, "             # node label\n");
%!         fputs(fd, "             node_id_rotor, \n");
%!         fputs(fd, "             # direction\n");
%!         fputs(fd, "             0.,0.,1.,\n");
%!         fputs(fd, "             # angular velocity\n");
%!         if (options.f_rbk)
%!           fputs(fd, "           null;\n");
%!         else
%!           fputs(fd, "           reference, drive_id_rotor_speed;\n");
%!         endif
%!         fputs(fd, "        inertia: elem_id_inertia,\n");
%!         fputs(fd, "                 position, reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                 orientation, reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                 body, all, joint, all, loadable, all,\n");
%!         fputs(fd, "                 output, both;\n");
%!         fputs(fd, "        gravity: uniform, 0., 0., -1., g;\n");
%!         fputs(fd, "end: elements;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!     endif
%!     fd = -1;
%!     if (options.f_enable_modal || options.f_enable_solid)
%!       geometry_file = [filename, ".geo"];
%!       unwind_protect
%!         [fd, msg] = fopen(geometry_file, "w");
%!         if (fd == -1)
%!           error("failed to open file \"%s\"", geometry_file);
%!         endif
%!         fn = fieldnames(param);
%!         for i=1:length(fn)
%!           fprintf(fd, "%s = %g;\n", fn{i}, getfield(param, fn{i}));
%!         endfor
%!         fputs(fd, "SetFactory(\"Built-in\");\n");
%!         fputs(fd, "Point(1) = {0, 0, -0.5 * l};\n");
%!         fputs(fd, "Point(2) = {0.5 * d, 0, -0.5 * l};\n");
%!         fputs(fd, "Point(3) = {0.5 * d, 0, -0.5 * w + o};\n");
%!         fputs(fd, "Point(4) = {0.5 * D, 0, -0.5 * w + o};\n");
%!         fputs(fd, "Point(5) = {0.5 * D, 0, 0.5 * w + o};\n");
%!         fputs(fd, "Point(6) = {0.5 * d, 0, 0.5 * w + o};\n");
%!         fputs(fd, "Point(7) = {0.5 * d, 0, 0.5 * l};\n");
%!         fputs(fd, "Point(8) = {0, 0, 0.5 * l};\n");
%!         fputs(fd, "Point(9) = {0, 0, 0.5 * w + o};\n");
%!         fputs(fd, "Point(10) = {0, 0, -0.5 * w + o};\n");
%!         fputs(fd, "Line(1) = {1, 2};\n");
%!         fputs(fd, "Line(2) = {2, 3};\n");
%!         fputs(fd, "Line(3) = {3, 4};\n");
%!         fputs(fd, "Line(4) = {4, 5};\n");
%!         fputs(fd, "Line(5) = {5, 6};\n");
%!         fputs(fd, "Line(6) = {6, 7};\n");
%!         fputs(fd, "Line(7) = {7, 8};\n");
%!         fputs(fd, "Line(8) = {8, 9};\n");
%!         fputs(fd, "Line(9) = {9, 10};\n");
%!         fputs(fd, "Line(10) = {10, 1};\n");
%!         fputs(fd, "Line(11) = {10,3};\n");
%!         fputs(fd, "Line(12) = {3,6};\n");
%!         fputs(fd, "Line(13) = {6,9};\n");
%!         fputs(fd, "Transfinite Curve(1) = Max(1, Round(0.5 * d / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(2) = Max(1, Round((0.5 * l - 0.5 * w + o) / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(3) = Max(1, Round((0.5 * (D - d)) / h2)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(4) = Max(1, Round(w / h2)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(5) = Max(1, Round((0.5 * (D - d)) / h2)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(6) = Max(1, Round((0.5 * l - 0.5 * w - o) / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(7) = Max(1, Round(0.5 * d / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(8) = Max(1, Round((0.5 * l - 0.5 * w - o) / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(9) = Max(1, Round(w / h2)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(10) = Max(1, Round((0.5 * l - 0.5 * w + o) / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(11) = Max(1, Round(0.5 * d / h1)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(12) = Max(1, Round(w / h2)) + 1;\n");
%!         fputs(fd, "Transfinite Curve(13) = Max(1, Round(0.5 * d / h1)) + 1;\n");
%!         fputs(fd, "Line Loop(1) = {1, 2, -11, 10};\n");
%!         fputs(fd, "Line Loop(2) = {11, 12, 13, 9};\n");
%!         fputs(fd, "Line Loop(3) = {3, 4, 5, -12};\n");
%!         fputs(fd, "Line Loop(4) = {-13, 6, 7, 8};\n");
%!         fputs(fd, "Plane Surface(1) = {1};\n");
%!         fputs(fd, "Plane Surface(2) = {2};\n");
%!         fputs(fd, "Plane Surface(3) = {3};\n");
%!         fputs(fd, "Plane Surface(4) = {4};\n");
%!         fputs(fd, "Transfinite Surface(1) = {1, 2, 3, 10};\n");
%!         fputs(fd, "Transfinite Surface(2) = {10, 3, 6, 9};\n");
%!         fputs(fd, "Transfinite Surface(3) = {3, 4, 5, 6};\n");
%!         fputs(fd, "Transfinite Surface(4) = {9, 6, 7, 8};\n");
%!         fputs(fd, "vol11[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{1}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol21[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol11[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol31[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol21[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol41[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol31[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "Recombine Surface{1, vol11[0]};\n");
%!         fputs(fd, "Recombine Surface{vol11[0], vol21[0]};\n");
%!         fputs(fd, "Recombine Surface{vol21[0], vol31[0]};\n");
%!         fputs(fd, "Recombine Surface{vol31[0], vol41[0]};\n");
%!         fputs(fd, "vol12[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{2}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol22[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol12[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol32[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol22[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol42[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol32[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "Recombine Surface{2, vol12[0]};\n");
%!         fputs(fd, "Recombine Surface{vol12[0], vol22[0]};\n");
%!         fputs(fd, "Recombine Surface{vol22[0], vol32[0]};\n");
%!         fputs(fd, "Recombine Surface{vol32[0], vol42[0]};\n");
%!         fputs(fd, "vol13[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{3}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol23[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol13[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol33[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol23[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol43[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol33[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "Recombine Surface{3, vol13[0]};\n");
%!         fputs(fd, "Recombine Surface{vol13[0], vol23[0]};\n");
%!         fputs(fd, "Recombine Surface{vol23[0], vol33[0]};\n");
%!         fputs(fd, "Recombine Surface{vol33[0], vol43[0]};\n");
%!         fputs(fd, "vol14[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{4}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol24[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol14[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol34[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol24[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "vol44[] = Extrude{{0, 0, 1},{0, 0, 0}, Pi/2}{ Surface{vol34[0]}; Layers{Round(d * Pi / (4. * h1))}; Recombine; };\n");
%!         fputs(fd, "Recombine Surface{4, vol14[0]};\n");
%!         fputs(fd, "Recombine Surface{vol14[0], vol24[0]};\n");
%!         fputs(fd, "Recombine Surface{vol24[0], vol34[0]};\n");
%!         fputs(fd, "Recombine Surface{vol34[0], vol44[0]};\n");
%!         fputs(fd, "Coherence;\n");
%!         switch (options.elem_type)
%!           case "iso8"
%!             fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!             fputs(fd, "Mesh.ElementOrder = 1;");
%!           case "iso20"
%!             fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!             fputs(fd, "Mesh.ElementOrder = 2;");
%!         endswitch
%!         fputs(fd, "Mesh 3;\n");
%!         fputs(fd, "Coherence Mesh;\n");
%!         fputs(fd, "Physical Surface(\"bearing1\", 1) = {21, 38, 55, 72};\n");
%!         fputs(fd, "Physical Surface(\"bearing2\", 2) = {260, 243, 277, 294};\n");
%!         fputs(fd, "Physical Surface(\"rotor\", 3) = {202, 180, 224, 158};\n");
%!         fputs(fd, "Physical Volume(\"shaft\", 1) = {1, 2, 4, 3, 14, 13, 15, 16, 5, 6, 8, 7};\n");
%!         fputs(fd, "Physical Volume(\"disc\", 2) = {9, 10, 12, 11};\n");
%!         fputs(fd, "Mesh.Format = 1;\n");
%!         fprintf(fd, "Save \"%s.msh\";\n", filename);
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!       ##spawn_wait(spawn("gmsh", {geometry_file})); return;
%!       fprintf(stderr, "meshing ...\n");
%!       pid = spawn("gmsh", {"-format", "msh2", "-0",geometry_file});
%!       status = spawn_wait(pid);
%!       if (status ~= 0)
%!         warning("gmsh failed with status %d", status);
%!       endif
%!       fprintf(stderr, "loading mesh ...\n");
%!       switch (options.elem_type)
%!         case "iso20"
%!           opt_mesh.elem_type = {"iso20", "penta15", "quad8", "tria6h"};
%!         case "iso8"
%!           opt_mesh.elem_type = {"iso8", "penta6", "iso4", "tria3"};
%!       endswitch
%!       mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh",opt_mesh));
%!       cms_opt.modes.number = int32(options.number_of_modes);
%!       cms_opt.nodes.modal.number = int32(rows(mesh.nodes) + 1);
%!       cms_opt.nodes.interfaces(1).number = int32(rows(mesh.nodes) + 2);
%!       cms_opt.nodes.interfaces(2).number = int32(rows(mesh.nodes) + 3);
%!       cms_opt.invariants = true;
%!       mesh.nodes(cms_opt.nodes.modal.number, :) = [0, 0, param.o, 0, 0, 0];
%!       mesh.nodes([cms_opt.nodes.interfaces.number], :) = [0, 0, -0.5 * param.l, 0, 0, 0;
%!                                                           0, 0,  0.5 * param.l, 0, 0, 0];
%!       switch (options.elem_type)
%!         case "iso20"
%!           mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 1, cms_opt.nodes.interfaces(1).number, "tria6h");
%!           mesh.elements.rbe3(3) = fem_pre_mesh_rbe3_from_surf(mesh, 2, cms_opt.nodes.interfaces(2).number, "tria6h");
%!           mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 3, cms_opt.nodes.modal.number, "quad8");
%!         case "iso8"
%!           mesh.elements.rbe3(2) = fem_pre_mesh_rbe3_from_surf(mesh, 1, cms_opt.nodes.interfaces(1).number, "iso4");
%!           mesh.elements.rbe3(3) = fem_pre_mesh_rbe3_from_surf(mesh, 2, cms_opt.nodes.interfaces(2).number, "iso4");
%!           mesh.elements.rbe3(1) = fem_pre_mesh_rbe3_from_surf(mesh, 3, cms_opt.nodes.modal.number, "iso4");
%!       endswitch
%!       load_case_dof.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
%!       switch (options.elem_type)
%!         case "iso20"
%!           mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!           mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!           mesh.materials.iso20([mesh.groups.iso20(find([mesh.groups.iso20.id] == 1)).elements]) = 1;
%!           mesh.materials.iso20([mesh.groups.iso20(find([mesh.groups.iso20.id] == 2)).elements]) = 2;
%!           mesh.materials.penta15([mesh.groups.penta15(find([mesh.groups.penta15.id] == 1)).elements]) = 1;
%!           mesh.materials.penta15([mesh.groups.penta15(find([mesh.groups.penta15.id] == 2)).elements]) = 2;
%!         case "iso8"
%!           mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!           mesh.materials.iso8([mesh.groups.iso8(find([mesh.groups.iso8.id] == 1)).elements]) = 1;
%!           mesh.materials.iso8([mesh.groups.iso8(find([mesh.groups.iso8.id] == 2)).elements]) = 2;
%!       endswitch
%!       mesh.material_data(1).rho = param.rho;
%!       mesh.material_data(1).E = param.E;
%!       mesh.material_data(1).nu = param.nu;
%!       mesh.material_data(2).rho = param.rho;
%!       mesh.material_data(2).E = 100 * param.E;
%!       mesh.material_data(2).nu = param.nu;
%!     endif
%!     if (options.f_enable_modal)
%!       fprintf(stderr, "building cms element ...\n");
%!       [mesh_cms, mat_ass, dof_map, sol_eig, cms_opt] = fem_cms_create2(mesh, load_case_dof, cms_opt);
%!       mesh_post_pro_file = sprintf("%s_post.msh", filename);
%!       fem_post_mesh_export(mesh_post_pro_file, mesh_cms);
%!       for j=1:size(sol_eig.def, 3)
%!         eig_post_pro_file_mode{j} = sprintf("%s_eig_def_%03d.msh", filename, j);
%!         fem_post_sol_step_export(eig_post_pro_file_mode{j}, sol_eig, j, j, sol_eig.f(j), options.scale_def / max(norm(sol_eig.def(:, 1:3, j), "rows")));
%!       endfor
%!       eig_post_pro_file = sprintf("%s_modes_post.geo", filename);
%!       fd = -1;
%!       unwind_protect
%!         [fd, msg] = fopen(eig_post_pro_file, "w");
%!         if (fd == -1)
%!           error("failed to open file \"%s\"", eig_post_pro_file);
%!         endif
%!         fprintf(fd, "Merge \"%s\";\n", mesh_post_pro_file);
%!         for j=1:numel(eig_post_pro_file_mode)
%!           fprintf(fd, "Merge \"%s\";\n", eig_post_pro_file_mode{j});
%!         endfor
%!         fputs(fd, "View.Type = 1;\n");
%!         fputs(fd, "View.VectorType = 5;\n");
%!         fputs(fd, "View.Visible = 1;\n");
%!         fputs(fd, "View.DisplacementFactor = 1;\n");
%!         fputs(fd, "View.ShowTime = 6;\n");
%!         fputs(fd, "View.ShowElement = 1;\n");
%!         fputs(fd, "View.IntervalsType = 3;\n");
%!         fputs(fd, "View.NbIso = 20;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!       mat_ass.Dred = param.alpha * mat_ass.Mred + param.beta * mat_ass.Kred;
%!       fem_cms_export([filename, "_cms"], mesh_cms, dof_map, mat_ass, cms_opt);
%!     endif
%!     if (options.f_enable_solid)
%!       opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!       opt_mbd_mesh.struct_nodes.type([cms_opt.nodes.modal.number]) = MBDYN_NODE_TYPE_DYNAMIC_STRUCT;
%!       opt_mbd_mesh.struct_nodes.type([cms_opt.nodes.interfaces.number]) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!       opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_shaft";
%!       solid_csl_file = [filename, "_solid.csl"];
%!       solid_nodes_file = [filename, "_solid.nod"];
%!       solid_elem_file = [filename, "_solid.elm"];
%!       load_case = struct();
%!       assert_simple(mesh.nodes(cms_opt.nodes.modal.number, 1:3), [0, 0, param.o]);
%!       assert_simple(mesh.nodes(cms_opt.nodes.interfaces(1).number, 1:3), [0, 0, -0.5 * param.l]);
%!       assert_simple(mesh.nodes(cms_opt.nodes.interfaces(2).number, 1:3), [0, 0, 0.5 * param.l]);
%!       opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, solid_csl_file, opt_mbd_mesh);
%!       opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, solid_nodes_file, opt_mbd_mesh);
%!       opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, solid_elem_file, opt_mbd_mesh);
%!       fd = -1;
%!       unwind_protect
%!         fd = fopen(mbdyn_filenames{3}, "w");
%!         if (fd == -1)
%!           error("failed to open file \"%s\"", mbdyn_filenames{2});
%!         endif
%!         fputs(fd, "include: \"${MBDYN_ROTOR_DYN_CMS_PARAM_FILE}\";\n");
%!         fputs(fd, "set: integer ref_id_ground = 1001;\n");
%!         fputs(fd, "set: integer ref_id_shaft = 1002;\n");
%!         fputs(fd, "set: integer ref_id_rotor = 1003;\n");
%!         fputs(fd, "set: integer ref_id_bearing1 = 1004;\n");
%!         fputs(fd, "set: integer ref_id_bearing2 = 1005;\n");
%!         fprintf(fd, "set: integer node_id_rotor = %d;\n", cms_opt.nodes.modal.number);
%!         fprintf(fd, "set: integer node_id_bearing1 = %d;\n", cms_opt.nodes.interfaces(1).number);
%!         fprintf(fd, "set: integer node_id_bearing2 = %d;\n", cms_opt.nodes.interfaces(2).number);
%!         fputs(fd, "set: integer body_id_unbalance = 2001;\n");
%!         fputs(fd, "set: integer joint_id_bearing1 = 3006;\n");
%!         fputs(fd, "set: integer joint_id_bearing2 = 3007;\n");
%!         fputs(fd, "set: integer joint_id_drive = 3008;\n");
%!         fputs(fd, "set: integer elem_id_inertia = 4001;\n");
%!         fputs(fd, "set: integer drive_id_rotor_speed = 5001;\n");
%!         fputs(fd, "set: integer drive_id_time_step = 5002;\n");
%!         fputs(fd, "set: real initial_time = 0.;\n");
%!         fputs(fd, "set: real final_time = 2. * pi * n / (abs(omega1 + omega0) / 2.);\n");
%!         fputs(fd, "begin: data;\n");
%!         fputs(fd, "        problem: initial value;\n");
%!         fputs(fd, "end: data;\n");
%!         fputs(fd, "begin: initial value;\n");
%!         fputs(fd, "        initial time: initial_time;\n");
%!         fputs(fd, "        final time: final_time;\n");
%!         fputs(fd, "        time step: 2. * pi / (360. * abs(omega0));\n");
%!         fputs(fd, "        strategy: change, postponed, drive_id_time_step;\n");
%!         fputs(fd, "        method: ms, 0.6;\n");
%!         fputs(fd, "        tolerance: 1e-5, test, sepnorm, 1e-5, test, norm;\n");
%!         fputs(fd, "        derivatives tolerance: 1e-3, 1e-3;\n");
%!         fputs(fd, "        max iterations: 100;\n");
%!         fputs(fd, "        derivatives max iterations: 10;\n");
%!         fputs(fd, "        derivatives coefficient: auto;\n");
%!         fputs(fd, "        output: iterations, cpu time, solver condition number, stat, yes;\n");
%!         fputs(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!         fputs(fd, "        enforce constraint equations: constraint violations;\n");
%!         fprintf(fd, "      threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!         fprintf(fd, "      threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!         fputs(fd, "        nonlinear solver: nox, modified, 25,\n");
%!         fputs(fd, "             keep jacobian matrix,\n");
%!         fputs(fd, "             inner iterations before assembly, 12,\n");
%!         fputs(fd, "             use preconditioner as solver, no,\n");
%!         fputs(fd, "             jacobian operator, newton krylov,\n");
%!         fputs(fd, "             solver, line search based,\n");
%!         fputs(fd, "             forcing term, type 2,\n");
%!         fputs(fd, "             forcing term min tolerance, 1e-10,\n");
%!         fputs(fd, "             forcing term max tolerance, 1e-8,\n");
%!         fputs(fd, "             direction, newton,\n");
%!         fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!         fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!         fputs(fd, "             linear solver, gmres,\n");
%!         fputs(fd, "             linear solver max iterations, 24,\n");
%!         fputs(fd, "             krylov subspace size, 24;\n");
%!         fputs(fd, "end: initial value;\n");
%!         fputs(fd, "begin: control data;\n");
%!         fputs(fd, "       use automatic differentiation;\n");
%!         fputs(fd, "       default output: none, structural nodes;\n");
%!         fputs(fd, "       default orientation: euler123;\n");
%!         fputs(fd, "       output precision: 16;\n");
%!         fputs(fd, "       max iterations: 0;\n");
%!         fprintf(fd, "        structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!         fprintf(fd, "        solids: %d;\n", opt_mbd_mesh.solids.number);
%!         fprintf(fd, "        joints: 3 + %d;\n", opt_mbd_mesh.joints.number);
%!         fputs(fd,   "        rigid bodies: 1;\n");
%!         if (options.f_rbk)
%!           fputs(fd, "        rigid body kinematics: drive, angular velocity,\n");
%!           fputs(fd, "        component,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           array, 2,\n");
%!           fputs(fd, "             mult, time, const, (omega1 - omega0) / (final_time - initial_time),\n");
%!           fputs(fd, "             const, omega0,\n");
%!           fputs(fd, "        angular acceleration,\n");
%!           fputs(fd, "        component,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           null,\n");
%!           fputs(fd, "           const, (omega1 - omega0) / (final_time - initial_time);\n");
%!         endif
%!         fputs(fd, "          gravity;\n");
%!         fputs(fd, "end: control data;\n");
%!         fputs(fd, "reference: ref_id_ground,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, eye,\n");
%!         fputs(fd, "        reference, global, null,\n");
%!         fputs(fd, "        reference, global, null;\n");
%!         fputs(fd, "reference: ref_id_shaft,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         fputs(fd, "        reference, ref_id_ground, eye,\n");
%!         fputs(fd, "        reference, ref_id_ground, null,\n");
%!         if (options.f_rbk)
%!           fputs(fd, "        reference, ref_id_ground, null;\n");
%!         else
%!           fputs(fd, "        reference, ref_id_ground, 0., 0., omega0;\n");
%!         endif
%!         fputs(fd, "reference: ref_id_rotor,\n");
%!         fputs(fd, "        reference, ref_id_shaft, 0., 0., o,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing1,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   -0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_bearing2,\n");
%!         fputs(fd, "        reference, ref_id_shaft,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.,\n");
%!         fputs(fd, "                   0.5 * l,\n");
%!         fputs(fd, "        reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        reference, ref_id_shaft, null;\n");
%!         fputs(fd, "begin: nodes;\n");
%!         fprintf(fd, "include: \"%s\";\n", solid_nodes_file);
%!         fputs(fd, "end: nodes;\n");
%!         fprintf(fd, "include: \"%s\";\n", solid_csl_file);
%!         fputs(fd, "begin: elements;\n");
%!         fputs(fd, "       drive caller: drive_id_rotor_speed, string, \"(omega1 - omega0) / (final_time - initial_time) * Time + omega0\", output, yes;\n");
%!         fputs(fd, "       drive caller: drive_id_time_step, string, \"2. * pi / (max(1., 360. * abs(model::drive(drive_id_rotor_speed, Time))))\";\n");
%!         fprintf(fd, "      include: \"%s\";\n", solid_elem_file);
%!         fputs(fd, "       body: body_id_unbalance,\n");
%!         fputs(fd, "             node_id_rotor,\n");
%!         fputs(fd, "                dm,\n");
%!         fputs(fd, "                  reference, ref_id_rotor, dr, 0., 0.,\n");
%!         fputs(fd, "                  diag, 0., 0., 0.;\n");
%!         fputs(fd, "        joint: joint_id_bearing1, total pin joint,\n");
%!         fputs(fd, "                node_id_bearing1,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing1, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing1, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, active,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "        joint: joint_id_bearing2, total pin joint,\n");
%!         fputs(fd, "               node_id_bearing2,\n");
%!         fputs(fd, "                        position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "                        position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "                        rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position, reference, ref_id_bearing2, null,\n");
%!         fputs(fd, "               position orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               rotation orientation, reference, ref_id_bearing2, eye,\n");
%!         fputs(fd, "               position constraint,\n");
%!         fputs(fd, "                        active, active, inactive,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "               orientation constraint,\n");
%!         fputs(fd, "                        inactive, inactive, inactive,\n");
%!         fputs(fd, "                        null;\n");
%!         fputs(fd, "     joint: joint_id_drive, angular velocity,\n");
%!         fputs(fd, "             # node label\n");
%!         fputs(fd, "             node_id_rotor, \n");
%!         fputs(fd, "             # direction\n");
%!         fputs(fd, "             0.,0.,1.,\n");
%!         fputs(fd, "             # angular velocity\n");
%!         if (~options.f_rbk)
%!           fputs(fd, "           reference, drive_id_rotor_speed;\n");
%!         else
%!           fputs(fd, "           null;\n");
%!         endif
%!         fputs(fd, "        inertia: elem_id_inertia,\n");
%!         fputs(fd, "                 position, reference, ref_id_rotor, null,\n");
%!         fputs(fd, "                 orientation, reference, ref_id_rotor, eye,\n");
%!         fputs(fd, "                 body, all, solid, all, loadable, all,\n");
%!         fputs(fd, "                 output, both;\n");
%!         fputs(fd, "        gravity: uniform, 0., 0., -1., g;\n");
%!         fputs(fd, "end: elements;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!     endif
%!     for i=1:numel(mbdyn_filenames)
%!       if (~enable_filenames(i))
%!         continue;
%!       endif
%!       options_mbd(i).output_file = sprintf("%s_%d", filename, i);
%!       if (~options.verbose)
%!         options_mbd(i).logfile = sprintf("%s_%d.stdout", filename, i);
%!       endif
%!       options_mbd(i).f_run_mbdyn2easyanim = false;
%!       param_file = sprintf("%s_%d.set", filename, i);
%!       putenv("MBDYN_ROTOR_DYN_CMS_ELEM_FILE", [filename, "_cms.elm"]);
%!       putenv("MBDYN_ROTOR_DYN_CMS_PARAM_FILE", param_file);
%!       mbdyn_pre_write_param_file(param_file, param);
%!       mbdyn_solver_run(mbdyn_filenames{i}, options_mbd(i));
%!       res(i).log_dat = mbdyn_post_load_log(options_mbd(i).output_file);
%!       [res(i).t, res(i).trajectory, res(i).deformation, res(i).velocity, res(i).acceleration, res(i).node_id] = mbdyn_post_load_output_struct(options_mbd(i).output_file);
%!       res(i).log_dat.vars = mbdyn_post_id_to_index(res(i), res(i).log_dat.vars);
%!       [res(i).drive_id, res(i).drive_data] = mbdyn_post_load_output_drv(options_mbd(i).output_file);
%!     endfor
%!     r = omega = cell(1, numel(mbdyn_filenames));
%!     if (options.f_plot)
%!       figure("visible", "off");
%!       hold on;
%!     endif
%!     for i=1:numel(mbdyn_filenames)
%!       if (~enable_filenames(i))
%!         continue;
%!       endif
%!       omega{i} = res(i).drive_data{find(res(i).drive_id == res(i).log_dat.vars.drive_id_rotor_speed)};
%!       r{i} = norm(res(i).trajectory{res(i).log_dat.vars.node_idx_rotor}(:, 1:2), "rows");
%!       if (options.f_plot)
%!         plot(omega{i} * 30 / pi * (1 / SI_unit.second), 1e3 * r{i} * SI_unit.meter, sprintf("-;%s;%d", printable_title(mbdyn_filename_suffix{i}), i));
%!       endif
%!     endfor
%!     if (options.f_plot)
%!       xlabel("n [rpm]");
%!       ylabel("r [mm]");
%!       grid on;
%!       grid minor on;
%!       title("resonance curve center of disk versus speed - magnitude");
%!       figure_list();
%!     endif
%!     tol = 0.5e-2;
%!     for i=1:numel(mbdyn_filenames)
%!       if (~enable_filenames(i))
%!         continue;
%!       endif
%!       for j=1:numel(mbdyn_filenames)
%!         if (~enable_filenames(j))
%!           continue;
%!         endif
%!         assert_simple(interp1(omega{i}, r{i}, omega{j}, "pchip", "extrap"), r{j}, tol * max(abs(r{j})));
%!       endfor
%!     endfor
%! %unwind_protect_cleanup
%!     if (numel(filename))
%!       fn = dir([filename, "*"]);
%!       for i=1:numel(fn)
%!         if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!           warning("failed to remove file \"%s\"", fn(i).name);
%!         endif
%!       endfor
%!     endif
%! %end_unwind_protect
%! endfunction
%!
%! ## Define the unit system
%! SI_unit.meter = 1e-3;
%! SI_unit.second = 1e-3;
%! SI_unit.kilogram = 1e-3;
%! SI_unit.newton = SI_unit.kilogram * SI_unit.meter / SI_unit.second^2;
%! SI_unit.pascal = SI_unit.newton / SI_unit.meter^2;
%! param.alpha = 0 / (1 / SI_unit.second);
%! param.beta = 1e-7 / (SI_unit.second);
%! param.h1 = 5e-3 / SI_unit.meter; ## fine mesh size
%! param.h2 = 5e-3 / SI_unit.meter; ## fine mesh size
%! param.l = 350e-3 / SI_unit.meter; ## bearing distance
%! param.d = 10e-3 / SI_unit.meter; ## shaft diameter
%! param.D = 150e-3 / SI_unit.meter; ##disk diameter
%! param.w = 15e-3 / SI_unit.meter; ## disk width
%! param.o = 75e-3 / SI_unit.meter; ## disk offset
%! param.ecg = 1e-3 * param.D;
%! param.dm = 1e-6 / SI_unit.kilogram;
%! param.E = 210000e6 / SI_unit.pascal;
%! param.nu = 0.3;
%! param.rho = 7850 / (SI_unit.kilogram / SI_unit.meter^3);
%! param.g = 9.81 / (SI_unit.meter / SI_unit.second^2);
%! m1 = param.D^2 * pi / 4 * param.w * param.rho;
%! param.dr = param.ecg * (m1 + param.dm) / param.dm;
%! param.omega0 = 1000 * pi / 30 / (1 / SI_unit.second);
%! param.omega1 = 10000 * pi / 30 / (1 / SI_unit.second);
%! param.n = 10;
%! options.number_of_modes = int32(10);
%! options.scale_def = 10e-3;
%! options.geo_tol = sqrt(eps);
%! options.code.use_package = false;
%! options.f_run_mbdyn = [true, true];
%! options.verbose = false;
%! options.elem_type = "iso20";
%! options.f_rbk = false;
%! options.f_enable_beam = true;
%! options.f_enable_modal = true;
%! options.f_enable_solid = false; ## disabled because of long execution time
%! [omega1, r1] = rotordynamics_test_case(param, options, SI_unit);
%! options.f_rbk = true;
%! options.f_enable_modal = false;
%! [omega2, r2] = rotordynamics_test_case(param, options, SI_unit);
%! assert_simple(r2{2}, r1{2}, 1e-3 * max(abs(r1{2})));
%! options.elem_type = "iso20";
%! options.f_enable_beam = false; ## disabled because of coarse mesh and linear elements
%! options.f_enable_solid = true;
%! options.f_enable_modal = true;
%! param.h1 = 10e-3 / SI_unit.meter; ## coarse mesh size
%! param.h2 = 100e-3 / SI_unit.meter; ## coarse mesh size
%! [omega3, r3] = rotordynamics_test_case(param, options, SI_unit);
%! options.f_rbk = false;
%! options.f_enable_modal = false;
%! [omega4, r4] = rotordynamics_test_case(param, options, SI_unit);
%! assert_simple(r4{3}, r3{3}, 1e-3 * max(abs(r3{3})));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
