## mbdyn_post_load_output_struct.tst:04
%!test
%! try
%! ## TEST4
%! ## TRACTA JOINT
%! methods = {"impliciteuler", "ms2", "ms3", "ms4", "ss2", "ss3", "ss4", "hope", "Bathe", "msstc3", "msstc4", "msstc5", "mssth3", "mssth4", "mssth5", "DIRK33", "DIRK43", "DIRK54", "hybrid,ms"};
%! stages  = [              1,     1,     1,     1,     1,     1,     1,      1,       2,        3,        4,        5,        3,        4,        5,        3,        4,        5,           1];
%! rhoinf = 0.2;
%! do_plot = false;
%! autodiff = [false, true];
%! rel_error = nan(numel(methods), numel(autodiff));
%! if (do_plot)
%!   close all;
%! endif
%! for idx_autodiff=1:numel(autodiff)
%! for idx_method=1:numel(methods)
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_struct_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " ## input parameters:\n");
%!     fprintf(fd, "set: integer stages = %d;\n",  stages(idx_method));
%!     fputs(fd, " set: real N = 360./4.;                   ## number of time steps per revolution\n");
%!     fputs(fd, " set: real omega = 40000 * pi / 30.;  ## input speed [rad/s9\n");
%!     fputs(fd, " set: real Phiz = 80. * pi / 180.;    ## prescribed angle between input shaft and output shaft [rad]\n");
%!     fputs(fd, " set: real t1 = 2. * pi / abs(omega); ## ramp up time to raise the angle from zero to Phiz [s]\n");
%!     fputs(fd, " set: real t2 = 3. * t1;             ## final time [s]\n");
%!     fputs(fd, " set: real dt = t1 / N * stages;      ## time step [s]\n");
%!     fputs(fd, " set: real R1 = 100e-3;               ## position of markers for postprocessing [m]\n");
%!     fputs(fd, " set: real R2 = 100e-3;\n");
%!     fputs(fd, " set: real R3 = 100e-3;\n");
%!     fputs(fd, " set: real R4 = 100e-3;\n");
%!     fputs(fd, " set: real l = 1000e-3;               ## length of the shafts [m]\n");
%!     fputs(fd, " set: real a = 200e-3;                ## distance between axes of revolute joints and bisector plane [m]\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: t1 + t2;\n");
%!     fputs(fd, "         time step: dt;\n");
%!     fputs(fd, "         max iterations: 100;\n");
%!     fputs(fd, "         tolerance: 1e-6, test, norm, 1e-6, test, norm;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-6, 1e-6;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     method = methods{idx_method};
%!     switch (method)
%!     case {"impliciteuler", "DIRK33", "DIRK43", "DIRK54"}
%!     otherwise
%!       method = sprintf("%s, %g", method, rhoinf);
%!     endswitch
%!     fprintf(fd, "         method: %s;\n", method);
%!     fputs(fd, "         linear solver:naive, colamd;\n");
%!     fputs(fd, "         nonlinear solver: nox, modified, 100,  keep jacobian matrix, inner iterations before assembly, 6, jacobian operator, newton krylov, forcing term, type2;\n");
%!     fputs(fd, "         threads: disable;\n");
%!     fputs(fd, "         output: iterations;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "       output meter: closest next, 0, forever, const, t1 / 36.;\n");
%!     fputs(fd, "       model: static;\n");
%!     if (autodiff(idx_autodiff))
%!       fputs(fd, "       use automatic differentiation;\n");
%!     endif
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       tolerance: 1e-8;\n");
%!     fputs(fd, "       print: equation description;\n");
%!     fputs(fd, "       structural nodes: 18;\n");
%!     fputs(fd, "       joints: 8;\n");
%!     fputs(fd, "       max iterations: 1000;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " set: integer ref_id_shaft_input = 1001;\n");
%!     fputs(fd, " set: integer ref_id_joint = 1002;\n");
%!     fputs(fd, " set: integer ref_id_hinge1 = 1003;\n");
%!     fputs(fd, " set: integer ref_id_hinge2 = 1004;\n");
%!     fputs(fd, " set: integer ref_id_hinge3 = 1005;\n");
%!     fputs(fd, " set: integer ref_id_hinge4 = 1006;\n");
%!     fputs(fd, " set: integer ref_id_hinge5 = 1007;\n");
%!     fputs(fd, " set: integer ref_id_hinge6 = 1008;\n");
%!     fputs(fd, " set: integer ref_id_hinge7 = 1009;\n");
%!     fputs(fd, " set: integer ref_id_shaft_output = 1010;\n");
%!     fputs(fd, " set: integer node_id_shaft_input = 2001;\n");
%!     fputs(fd, " set: integer node_id_shaft_intermediate1 = 2002;\n");
%!     fputs(fd, " set: integer node_id_shaft_intermediate2 = 2003;\n");
%!     fputs(fd, " set: integer node_id_shaft_output = 2004;\n");
%!     fputs(fd, " set: integer node_id_housing1 = 2005;\n");
%!     fputs(fd, " set: integer node_id_housing2 = 2006;\n");
%!     fputs(fd, " set: integer joint_id_shaft_input = 3001;\n");
%!     fputs(fd, " set: integer joint_id_hinge1 = 3002;\n");
%!     fputs(fd, " set: integer joint_id_hinge2 = 3003;\n");
%!     fputs(fd, " set: integer joint_id_hinge3 = 3004;\n");
%!     fputs(fd, " set: integer joint_id_def_hinge1 = 3005;\n");
%!     fputs(fd, " set: integer joint_id_def_hinge2 = 3006;\n");
%!     fputs(fd, " set: integer joint_id_def_hinge3 = 3007;\n");
%!     fputs(fd, " set: integer joint_id_spherical1 = 3008;\n");
%!     fputs(fd, " set: integer joint_id_spherical2 = 3009;\n");
%!     fputs(fd, " set: integer joint_id_spherical3 = 3010;\n");
%!     fputs(fd, " set: integer joint_id_shaft_output = 3011;\n");
%!     fputs(fd, " set: integer joint_id_hinge5 = 3013;\n");
%!     fputs(fd, " set: integer joint_id_hinge6 = 3014;\n");
%!     fputs(fd, " set: integer joint_id_hinge7 = 3015;\n");
%!     fputs(fd, " set: integer joint_id_hinge7_def = 3016;\n");
%!     fputs(fd, " reference: ref_id_shaft_input,\n");
%!     fputs(fd, "            position, reference, global, null,\n");
%!     fputs(fd, "            orientation, reference, global, eye,\n");
%!     fputs(fd, "            velocity, reference, global, null,\n");
%!     fputs(fd, "            angular velocity, reference, global, omega, 0., 0.;\n");
%!     fputs(fd, " reference: ref_id_joint,\n");
%!     fputs(fd, "            position, reference, ref_id_shaft_input, l, 0., 0.,\n");
%!     fputs(fd, "            orientation, reference, ref_id_shaft_input, eye,\n");
%!     fputs(fd, "            velocity, reference, ref_id_shaft_input, null,\n");
%!     fputs(fd, "            angular velocity, reference, ref_id_shaft_input, null;\n");
%!     fputs(fd, " reference: ref_id_hinge1,\n");
%!     fputs(fd, "            position, reference, ref_id_joint, -a, 0., 0.,\n");
%!     fputs(fd, "            orientation, reference, ref_id_joint, 3, 0., 0., 1.,\n");
%!     fputs(fd, "                                                  1, 1., 0., 0.,\n");
%!     fputs(fd, "            velocity, reference, ref_id_joint, null,\n");
%!     fputs(fd, "            angular velocity, reference, ref_id_joint, null;\n");
%!     fputs(fd, " reference: ref_id_hinge2,\n");
%!     fputs(fd, "            position, reference, ref_id_joint, null,\n");
%!     fputs(fd, "            orientation, reference, ref_id_joint, 3, 0., 1., 0.,\n");
%!     fputs(fd, "                                                  1, 1., 0., 0.,\n");
%!     fputs(fd, "            velocity, reference, ref_id_joint, null,\n");
%!     fputs(fd, "            angular velocity, reference, ref_id_joint, null;\n");
%!     fputs(fd, " reference: ref_id_hinge3,\n");
%!     fputs(fd, "            position, reference, ref_id_joint, a, 0., 0.,\n");
%!     fputs(fd, "            orientation, reference, ref_id_joint, 3, 0., 0., 1.,\n");
%!     fputs(fd, "                                                  1, 1., 0., 0.,\n");
%!     fputs(fd, "            velocity, reference, ref_id_joint, null,\n");
%!     fputs(fd, "            angular velocity, reference, ref_id_joint, null;\n");
%!     fputs(fd, " reference: ref_id_shaft_output,\n");
%!     fputs(fd, "            position, reference, ref_id_joint, l, 0., 0.,\n");
%!     fputs(fd, "            orientation, reference, ref_id_joint, eye,\n");
%!     fputs(fd, "            velocity, reference, ref_id_joint, null,\n");
%!     fputs(fd, "            angular velocity, reference, ref_id_joint, null;\n");
%!     fputs(fd, " reference: ref_id_hinge5,\n");
%!     fputs(fd, "            position, reference, ref_id_shaft_output, null,\n");
%!     fputs(fd, "            orientation, reference, ref_id_shaft_output, 3, 1., 0., 0.,\n");
%!     fputs(fd, "                                                         1, 0., 0., 1.,\n");
%!     fputs(fd, "            velocity, reference, ref_id_shaft_output, null,\n");
%!     fputs(fd, "            angular velocity, reference, ref_id_shaft_output, null;\n");
%!     fputs(fd, " reference: ref_id_hinge6,\n");
%!     fputs(fd, "            position, reference, ref_id_shaft_input, null,\n");
%!     fputs(fd, "            orientation, reference, ref_id_shaft_input, 3, 1., 0., 0.,\n");
%!     fputs(fd, "                                                        1, 0., 0., 1.,\n");
%!     fputs(fd, "            velocity, reference, ref_id_shaft_input, null,\n");
%!     fputs(fd, "            angular velocity, reference, ref_id_shaft_input, null;\n");
%!     fputs(fd, " reference: ref_id_hinge7,\n");
%!     fputs(fd, "            position, reference, ref_id_joint, null,\n");
%!     fputs(fd, "            orientation, reference, ref_id_joint, 3, 0., 1., 0.,\n");
%!     fputs(fd, "                                                  1, 1., 0., 0.,\n");
%!     fputs(fd, "            velocity, reference, ref_id_joint, null,\n");
%!     fputs(fd, "            angular velocity, reference, ref_id_joint, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, " structural: node_id_shaft_input, dynamic,\n");
%!     fputs(fd, "             position, reference, ref_id_shaft_input, null,\n");
%!     fputs(fd, "             orientation, reference, ref_id_shaft_input, eye,\n");
%!     fputs(fd, "             velocity, reference, ref_id_shaft_input, null,\n");
%!     fputs(fd, "             angular velocity, reference, ref_id_shaft_input, null;\n");
%!     fputs(fd, " structural: node_id_shaft_intermediate1, dynamic,\n");
%!     fputs(fd, "             position, reference, ref_id_joint, null,\n");
%!     fputs(fd, "             orientation, reference, ref_id_joint, eye,\n");
%!     fputs(fd, "             velocity, reference, ref_id_joint, null,\n");
%!     fputs(fd, "             angular velocity, reference, ref_id_joint, null;\n");
%!     fputs(fd, " structural: node_id_shaft_intermediate2, dynamic,\n");
%!     fputs(fd, "             position, reference, ref_id_joint, null,\n");
%!     fputs(fd, "             orientation, reference, ref_id_joint, eye,\n");
%!     fputs(fd, "             velocity, reference, ref_id_joint, null,\n");
%!     fputs(fd, "             angular velocity, reference, ref_id_joint, null;\n");
%!     fputs(fd, " structural: node_id_shaft_output, dynamic,\n");
%!     fputs(fd, "             position, reference, ref_id_shaft_output, null,\n");
%!     fputs(fd, "             orientation, reference, ref_id_shaft_output, eye,\n");
%!     fputs(fd, "             velocity, reference, ref_id_shaft_output, null,\n");
%!     fputs(fd, "             angular velocity, reference, ref_id_shaft_output, null;\n");
%!     fputs(fd, " structural: node_id_housing1, dynamic,\n");
%!     fputs(fd, "             position, reference, ref_id_shaft_input, null,\n");
%!     fputs(fd, "             orientation, reference, ref_id_shaft_input, eye,\n");
%!     fputs(fd, "             velocity, reference, ref_id_shaft_input, null,\n");
%!     fputs(fd, "             angular velocity, reference, global, null;\n");
%!     fputs(fd, " structural: node_id_housing2, dynamic,\n");
%!     fputs(fd, "             position, reference, ref_id_shaft_output, null,\n");
%!     fputs(fd, "             orientation, reference, ref_id_shaft_output, eye,\n");
%!     fputs(fd, "             velocity, reference, ref_id_shaft_output, null,\n");
%!     fputs(fd, "             angular velocity, reference, global, null;\n");
%!     fputs(fd, " structural: 1, dummy, node_id_shaft_input, offset,\n");
%!     fputs(fd, "             reference, ref_id_shaft_input, 0., R1, 0.,\n");
%!     fputs(fd, "             reference, ref_id_shaft_input, eye;\n");
%!     fputs(fd, " structural: 2, dummy, node_id_shaft_input, offset,\n");
%!     fputs(fd, "             reference, ref_id_shaft_input, 0., -R1, 0.,\n");
%!     fputs(fd, "             reference, ref_id_shaft_input, eye;\n");
%!     fputs(fd, " structural: 3, dummy, node_id_shaft_intermediate1, offset,\n");
%!     fputs(fd, "             reference, ref_id_joint, 0., 0., R2,\n");
%!     fputs(fd, "             reference, ref_id_joint, eye;\n");
%!     fputs(fd, " structural: 4, dummy, node_id_shaft_intermediate1, offset,\n");
%!     fputs(fd, "             reference, ref_id_joint, 0., 0., -R2,\n");
%!     fputs(fd, "             reference, ref_id_joint, eye;\n");
%!     fputs(fd, " structural: 5, dummy, node_id_shaft_intermediate1, offset,\n");
%!     fputs(fd, "             reference, ref_id_joint, 0., R2, 0.,\n");
%!     fputs(fd, "             reference, ref_id_joint, eye;\n");
%!     fputs(fd, " structural: 6, dummy, node_id_shaft_intermediate1, offset,\n");
%!     fputs(fd, "             reference, ref_id_joint, 0., -R2, 0.,\n");
%!     fputs(fd, "             reference, ref_id_joint, eye;\n");
%!     fputs(fd, " structural: 7, dummy, node_id_shaft_intermediate2, offset,\n");
%!     fputs(fd, "             reference, ref_id_joint, 0., 0., R3,\n");
%!     fputs(fd, "             reference, ref_id_joint, eye;\n");
%!     fputs(fd, " structural: 8, dummy, node_id_shaft_intermediate2, offset,\n");
%!     fputs(fd, "             reference, ref_id_joint, 0., 0., -R3,\n");
%!     fputs(fd, "             reference, ref_id_joint, eye;\n");
%!     fputs(fd, " structural: 9, dummy, node_id_shaft_intermediate2, offset,\n");
%!     fputs(fd, "             reference, ref_id_joint, 0., R3, 0.,\n");
%!     fputs(fd, "             reference, ref_id_joint, eye;\n");
%!     fputs(fd, " structural: 10, dummy, node_id_shaft_intermediate2, offset,\n");
%!     fputs(fd, "             reference, ref_id_joint, 0., -R3, 0.,\n");
%!     fputs(fd, "             reference, ref_id_joint, eye;\n");
%!     fputs(fd, " structural: 11, dummy, node_id_shaft_output, offset,\n");
%!     fputs(fd, "             reference, ref_id_shaft_output, 0., R4, 0.,\n");
%!     fputs(fd, "             reference, ref_id_shaft_output, eye;\n");
%!     fputs(fd, " structural: 12, dummy, node_id_shaft_output, offset,\n");
%!     fputs(fd, "             reference, ref_id_shaft_output, 0., -R4, 0.,\n");
%!     fputs(fd, "             reference, ref_id_shaft_output, eye;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         joint: joint_id_shaft_input, total pin joint,\n");
%!     fputs(fd, "                 node_id_shaft_input,\n");
%!     fputs(fd, "                 position, reference, ref_id_shaft_input, null,\n");
%!     fputs(fd, "                 position orientation, reference, ref_id_shaft_input, eye,\n");
%!     fputs(fd, "                 rotation orientation, reference, ref_id_shaft_input, eye,\n");
%!     fputs(fd, "                 position, reference, ref_id_shaft_input, null,\n");
%!     fputs(fd, "                 position orientation, reference, ref_id_shaft_input, eye,\n");
%!     fputs(fd, "                 rotation orientation, reference, ref_id_shaft_input, eye,\n");
%!     fputs(fd, "                 position constraint,\n");
%!     fputs(fd, "                 active, active, active, null,\n");
%!     fputs(fd, "                 orientation constraint,\n");
%!     fputs(fd, "                 angular velocity,\n");
%!     fputs(fd, "                 active,\n");
%!     fputs(fd, "                 active,\n");
%!     fputs(fd, "                 component, omega, 0., 0.;\n");
%!     fputs(fd, "         joint: joint_id_hinge1, revolute hinge,\n");
%!     fputs(fd, "                 node_id_shaft_input,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge1, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_hinge1, eye,\n");
%!     fputs(fd, "                 node_id_shaft_intermediate1,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge1, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_hinge1, eye;\n");
%!     fputs(fd, "         joint: joint_id_hinge2, total joint,\n");
%!     fputs(fd, "                 node_id_shaft_intermediate1,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge2, null,\n");
%!     fputs(fd, "                 position orientation, reference, ref_id_hinge2, eye,\n");
%!     fputs(fd, "                 rotation orientation, reference, ref_id_hinge2, eye,\n");
%!     fputs(fd, "                 node_id_shaft_intermediate2,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge2, null,\n");
%!     fputs(fd, "                 position orientation, reference, ref_id_hinge2, eye,\n");
%!     fputs(fd, "                 rotation orientation, reference, ref_id_hinge2, eye,\n");
%!     fputs(fd, "                 position constraint,\n");
%!     fputs(fd, "                 inactive, inactive, active, null,\n");
%!     fputs(fd, "                 orientation constraint,\n");
%!     fputs(fd, "                 active, active, inactive,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         joint: joint_id_hinge3, revolute hinge,\n");
%!     fputs(fd, "                 node_id_shaft_intermediate2,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge3, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_hinge3, eye,\n");
%!     fputs(fd, "                 node_id_shaft_output,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge3, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_hinge3, eye;\n");
%!     fputs(fd, "         joint: joint_id_hinge5, revolute hinge,\n");
%!     fputs(fd, "                 node_id_shaft_output,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge5, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_hinge5, eye,\n");
%!     fputs(fd, "                 node_id_housing2,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge5, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_hinge5, eye;\n");
%!     fputs(fd, "         joint: joint_id_hinge6, revolute hinge,\n");
%!     fputs(fd, "                 node_id_shaft_input,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge6, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_hinge6, eye,\n");
%!     fputs(fd, "                 node_id_housing1,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge6, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_hinge6, eye;\n");
%!     fputs(fd, "         joint: joint_id_shaft_output, total pin joint,\n");
%!     fputs(fd, "                 node_id_housing1,\n");
%!     fputs(fd, "                 position, reference, ref_id_shaft_input, null,\n");
%!     fputs(fd, "                 position orientation, reference, ref_id_shaft_input, eye,\n");
%!     fputs(fd, "                 rotation orientation, reference, ref_id_shaft_input, eye,\n");
%!     fputs(fd, "                 position, reference, ref_id_shaft_input, null,\n");
%!     fputs(fd, "                 position orientation, reference, ref_id_shaft_input, eye,\n");
%!     fputs(fd, "                 rotation orientation, reference, ref_id_shaft_input, eye,\n");
%!     fputs(fd, "                 position constraint,\n");
%!     fputs(fd, "                 inactive, inactive, inactive, null,\n");
%!     fputs(fd, "                 orientation constraint,\n");
%!     fputs(fd, "                 active,\n");
%!     fputs(fd, "                 inactive,\n");
%!     fputs(fd, "                 inactive,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         joint: joint_id_hinge7, total joint,\n");
%!     fputs(fd, "                 node_id_housing1,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge7, null,\n");
%!     fputs(fd, "                 position orientation, reference, ref_id_hinge7, eye,\n");
%!     fputs(fd, "                 rotation orientation, reference, ref_id_hinge7, eye,\n");
%!     fputs(fd, "                 node_id_housing2,\n");
%!     fputs(fd, "                 position, reference, ref_id_hinge7, null,\n");
%!     fputs(fd, "                 position orientation, reference, ref_id_hinge7, eye,\n");
%!     fputs(fd, "                 rotation orientation, reference, ref_id_hinge7, eye,\n");
%!     fputs(fd, "                 position constraint,\n");
%!     fputs(fd, "                 active, active, active, null,\n");
%!     fputs(fd, "                 orientation constraint, active, active, active, component,\n");
%!     fputs(fd, "                         const, 0.,\n");
%!     fputs(fd, "                         const, 0.,\n");
%!     fputs(fd, "                         string, \"((1 - cos(pi/2.*(Time / t1))^2) * (Time <= t1) + (Time > t1)) * Phiz\";\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [res.t, res.trajectory, res.deformation,res.velocity,res.acceleration,res.node_id, res.force_id, res.force_node_id, res.orientation_description]=mbdyn_post_load_output_struct (options.output_file);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   node_id = int32([2001, 2002, 2003, 2004]);
%!   idx_node = zeros(size(node_id), "int32");
%!   for i=1:numel(node_id)
%!     idx_node(i) = find(res.node_id == node_id(i));
%!   endfor
%!   W = cell(1, numel(idx_node));
%!   for i=1:numel(idx_node)
%!     R{i} = euler123_to_rotation_matrix(res.trajectory{idx_node(i)}(:,4:6).');
%!     W{i} = res.velocity{idx_node(i)}(:, 4:6).';
%!     Wrel{i} = zeros(size(W{i}));
%!     for j=1:columns(W{i})
%!       Wrel{i}(:, j) = R{i}(:, :, j).' * W{i}(:,j);
%!     endfor
%!   end
%!   if (do_plot)
%!   figure("visible", "off");
%!   hold on;
%!   for i=1:numel(idx_node)
%!     plot(res.t, W{i}(1, :) * 30 / pi, sprintf("-;W%dx;", i));
%!   endfor
%!   xlabel("t [s]");
%!   ylabel("n[rpm]");
%!   title("tracta joint omega1");
%!   grid on;
%!   grid minor on;
%!   figure("visible", "off");
%!   hold on;
%!   for i=1:numel(idx_node)
%!     plot(res.t, Wrel{i}(1, :) * 30 / pi, sprintf("-;Wrel%dx;", i));
%!   endfor
%!   xlabel("t [s]");
%!   ylabel("n[rpm]");
%!   title("tracta joint omega1 relative frame");
%!   grid on;
%!   grid minor on;
%!   endif
%!   t1 = log_dat.vars.t1;
%!   omega = log_dat.vars.omega;
%!   idx_t = find(res.t > 2 * t1);
%!   rel_err = -1;
%!   for i=[1,4]
%!     rel_err = max([rel_err, max(norm(Wrel{i}(:, idx_t) - [omega; zeros(2, 1)], "cols")) / abs(omega)]);
%!   endfor
%!   rel_error(idx_method, idx_autodiff) = rel_err;
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! endfor
%! endfor
%! if (do_plot)
%!   figure_list();
%! endif
%! for idx_method=1:numel(methods)
%!   fprintf(stderr, "%s: ", methods{idx_method});
%!   for idx_autodiff=1:numel(autodiff)
%!     fprintf(stderr, "%8.3e ", rel_error(idx_method, idx_autodiff));
%!   endfor
%!   fputs(stderr, "\n");
%! endfor
%! tol = 1e-7;
%! assert_simple(all(rel_error(:) < tol));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
