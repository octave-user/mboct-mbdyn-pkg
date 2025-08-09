## mbdyn_post_load_output_struct.tst:04
%!test
%! try
%! ## TEST 7
%! ## CARDANO JOINT
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! autodiff = [false, true];
%! methods = {"impliciteuler", "ms2", "ms3", "ms4", "ss2", "ss3", "ss4", "hope", "Bathe", "msstc3", "msstc4", "msstc5", "mssth3", "mssth4", "mssth5", "DIRK33", "DIRK43", "DIRK54", "hybrid,ms"};
%! stages  = [              1,     1,     1,     1,     1,     1,     1,      1,       2,        3,        4,        5,        3,        4,        5,        3,        4,        5,           1];
%! rhoinf = 0.2;
%! rel_error = nan(numel(methods), numel(autodiff));
%! fd = -1;
%! for idx_autodiff=1:numel(autodiff)
%! for idx_method=1:numel(methods)
%! %unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_struct_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fprintf(fd, "set: integer stages = %d;\n",  stages(idx_method));
%!       fputs(fd, " set: real N = 360. / stages;                  ## number of time steps\n");
%!       fputs(fd, " set: real omega = 10000 * pi / 30.; ## input speed\n");
%!       fputs(fd, " set: real oy = 20e-3;                   ## misalignment of shaft position in y-direction\n");
%!       fputs(fd, " set: real oz = 100e-3;                ## misalignment of shaft position z-direction\n");
%!       fputs(fd, " set: real l = 1000e-3;              ## length of shaft\n");
%!       fputs(fd, " set: real Phiy = 0. * pi / 180.;    ## misalignment of shaft orientation y-axis\n");
%!       fputs(fd, " set: real Phiz = 0. * pi / 180.;    ## misalignment of shaft orientation z-axis\n");
%!       fputs(fd, " set: real t1 = 2. * pi / abs(omega);\n");
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "         problem: initial value; # the default\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "         initial time: 0;\n");
%!       fputs(fd, "         final time: t1;\n");
%!       fputs(fd, "         time step: t1 / N;\n");
%!       fputs(fd, "         max iterations: 100;\n");
%!       fputs(fd, "         tolerance: 1e-8, test, norm, 1e-8, test, norm;\n");
%!       fputs(fd, "         derivatives tolerance: 1e-5, 1e-5;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         output: iterations, solver condition number, stat, yes;\n");
%!       method = methods{idx_method};
%!       switch (method)
%!       case {"impliciteuler", "DIRK33", "DIRK43", "DIRK54"}
%!       otherwise
%!         method = sprintf("%s, %g", method, rhoinf);
%!       endswitch
%!       fprintf(fd, "         method: %s;\n", method);
%!       fputs(fd, "         linear solver: umfpack, grad, scale, iterative, always;\n");
%!       fputs(fd, "         nonlinear solver: nox, modified, 100,  keep jacobian matrix, inner iterations before assembly, 8,jacobian operator, newton krylov, forcing term, type2, forcing term max tolerance, 1e-6, forcing term min tolerance, 1e-10;\n");
%!       fputs(fd, "         threads: disable;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "       output meter: closest next, 0, forever, const, t1 / 36.;\n");
%!       if (autodiff(idx_autodiff))
%!         fputs(fd, "       use automatic differentiation;\n");
%!       endif
%!       fputs(fd, "       output precision: 16;\n");
%!       fputs(fd, "       epsilon: 1;\n");
%!       fputs(fd, "       tolerance: 1e-8;\n");
%!       fputs(fd, "       print: equation description;\n");
%!       fputs(fd, "       structural nodes: 3;\n");
%!       fputs(fd, "       joints: 4;\n");
%!       fputs(fd, "       max iterations: 1000;\n");
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " set: integer ref_id_motor = 1001;\n");
%!       fputs(fd, " set: integer ref_id_shaft_input = 1002;\n");
%!       fputs(fd, " set: integer ref_id_shaft_output_static = 1003;\n");
%!       fputs(fd, " set: integer ref_id_shaft_output = 1004;\n");
%!       fputs(fd, " set: integer ref_id_shaft_axis_static = 1005;\n");
%!       fputs(fd, " set: integer ref_id_shaft_axis_dynamic = 1006;\n");
%!       fputs(fd, " set: integer node_id_shaft_input = 2001;\n");
%!       fputs(fd, " set: integer node_id_shaft_axis = 2002;\n");
%!       fputs(fd, " set: integer node_id_shaft_output = 2003;\n");
%!       fputs(fd, " set: integer joint_id_shaft_input = 3001;\n");
%!       fputs(fd, " set: integer joint_id_shaft_output = 3002;\n");
%!       fputs(fd, " set: integer joint_id_cardano_input = 3003;\n");
%!       fputs(fd, " set: integer joint_id_cardano_output = 3004;\n");
%!       fputs(fd, " reference: ref_id_motor,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            orientation, reference, global, eye,\n");
%!       fputs(fd, "            velocity, reference, global, null,\n");
%!       fputs(fd, "            angular velocity, reference, global, null;\n");
%!       fputs(fd, " reference: ref_id_shaft_input,\n");
%!       fputs(fd, "            position, reference, ref_id_motor, null,\n");
%!       fputs(fd, "            orientation, reference, ref_id_motor, eye,\n");
%!       fputs(fd, "            velocity, reference, ref_id_motor, null,\n");
%!       fputs(fd, "            angular velocity, reference, ref_id_motor, omega, 0., 0.;\n");
%!       fputs(fd, " reference: ref_id_shaft_output_static,\n");
%!       fputs(fd, "            position, reference, ref_id_motor, l, oy, oz,\n");
%!       fputs(fd, "            orientation, reference, ref_id_motor, euler123, 0., Phiy, Phiz,\n");
%!       fputs(fd, "            velocity, reference, ref_id_motor, null,\n");
%!       fputs(fd, "            angular velocity, reference, ref_id_motor, null;\n");
%!       fputs(fd, " reference: ref_id_shaft_output,\n");
%!       fputs(fd, "            position, reference, ref_id_shaft_output_static, null,\n");
%!       fputs(fd, "            orientation, reference, ref_id_shaft_output_static, eye,\n");
%!       fputs(fd, "            velocity, reference, ref_id_shaft_output_static, null,\n");
%!       fputs(fd, "            angular velocity, reference, ref_id_shaft_output_static, omega, 0., 0.;\n");
%!       fputs(fd, " reference: ref_id_shaft_axis_static,\n");
%!       fputs(fd, "            position, reference, ref_id_motor, null,\n");
%!       fputs(fd, "            orientation, reference, ref_id_motor, 1, l, oy, oz,\n");
%!       fputs(fd, "                                                  3, 0., 0., 1.,\n");
%!       fputs(fd, "            velocity, reference, ref_id_motor, null,\n");
%!       fputs(fd, "            angular velocity, reference, ref_id_motor, null;\n");
%!       fputs(fd, " reference: ref_id_shaft_axis_dynamic,\n");
%!       fputs(fd, "            position, reference, ref_id_shaft_axis_static, null,\n");
%!       fputs(fd, "            orientation, reference, ref_id_shaft_axis_static, eye,\n");
%!       fputs(fd, "            velocity, reference, ref_id_shaft_axis_static, null,\n");
%!       fputs(fd, "            angular velocity, reference, ref_id_shaft_axis_static, omega, 0., 0.;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fputs(fd, " structural: node_id_shaft_input, static,\n");
%!       fputs(fd, "             position, reference, ref_id_shaft_input, null,\n");
%!       fputs(fd, "             orientation, reference, ref_id_shaft_input, eye,\n");
%!       fputs(fd, "             velocity, reference, ref_id_shaft_input, null,\n");
%!       fputs(fd, "             angular velocity, reference, ref_id_shaft_input, null,\n");
%!       fputs(fd, "             assembly, 1e10, 1e-10, no;\n");
%!       fputs(fd, " structural: node_id_shaft_axis, static,\n");
%!       fputs(fd, "             position, reference, ref_id_shaft_axis_dynamic, null,\n");
%!       fputs(fd, "             orientation, reference, ref_id_shaft_axis_dynamic, eye,\n");
%!       fputs(fd, "             velocity, reference, ref_id_shaft_axis_dynamic, null,\n");
%!       fputs(fd, "             angular velocity, reference, ref_id_shaft_axis_dynamic, null,\n");
%!       fputs(fd, "             assembly, 1e-6, 1e-6, no;\n");
%!       fputs(fd, " structural: node_id_shaft_output, static,\n");
%!       fputs(fd, "             position, reference, ref_id_shaft_output, null,\n");
%!       fputs(fd, "             orientation, reference, ref_id_shaft_output, eye,\n");
%!       fputs(fd, "             velocity, reference, ref_id_shaft_output, null,\n");
%!       fputs(fd, "             angular velocity, reference, ref_id_shaft_output, null,\n");
%!       fputs(fd, "             assembly, 1e-3, 1e-6,no;\n");
%!       fputs(fd, " end: nodes;\n");
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "    joint: joint_id_shaft_input, total pin joint,\n");
%!       fputs(fd, "            node_id_shaft_input, \n");
%!       fputs(fd, "            position, reference, ref_id_shaft_input, null,\n");
%!       fputs(fd, "            position orientation, reference, ref_id_shaft_input, eye,\n");
%!       fputs(fd, "                 rotation orientation, reference, ref_id_shaft_input, eye,\n");
%!       fputs(fd, "            position, reference, ref_id_shaft_input, null,\n");
%!       fputs(fd, "                 position orientation, reference, ref_id_shaft_input, eye,\n");
%!       fputs(fd, "                 rotation orientation, reference, ref_id_shaft_input, eye,\n");
%!       fputs(fd, "                 position constraint,\n");
%!       fputs(fd, "                 active, active, active, null,\n");
%!       fputs(fd, "                 orientation constraint,\n");
%!       fputs(fd, "                 angular velocity,\n");
%!       fputs(fd, "                 active,\n");
%!       fputs(fd, "                 active,\n");
%!       fputs(fd, "                 component, omega, 0., 0.;\n");
%!       fputs(fd, "    joint: joint_id_shaft_output, total pin joint,\n");
%!       fputs(fd, "            node_id_shaft_output, \n");
%!       fputs(fd, "            position, reference, ref_id_shaft_output, null,\n");
%!       fputs(fd, "            position orientation, reference, ref_id_shaft_output, eye,\n");
%!       fputs(fd, "                 rotation orientation, reference, ref_id_shaft_output, eye,\n");
%!       fputs(fd, "            position, reference, ref_id_shaft_output, null,\n");
%!       fputs(fd, "                 position orientation, reference, ref_id_shaft_output, eye,\n");
%!       fputs(fd, "                 rotation orientation, reference, ref_id_shaft_output, eye,\n");
%!       fputs(fd, "                 position constraint,\n");
%!       fputs(fd, "                 inactive, active, active, null,\n");
%!       fputs(fd, "                 orientation constraint,\n");
%!       fputs(fd, "                 inactive,\n");
%!       fputs(fd, "                 active,\n");
%!       fputs(fd, "                 active,\n");
%!       fputs(fd, "                 null;\n");
%!       fputs(fd, "         joint: joint_id_cardano_input, Cardano hinge,\n");
%!       fputs(fd, "                node_id_shaft_input,\n");
%!       fputs(fd, "                position, reference, ref_id_shaft_input, null,\n");
%!       fputs(fd, "                orientation, reference, ref_id_shaft_input, eye,\n");
%!       fputs(fd, "                node_id_shaft_axis, reference, ref_id_shaft_input, null,\n");
%!       fputs(fd, "                orientation, reference, ref_id_shaft_axis_dynamic, eye;\n");
%!       fputs(fd, "         joint: joint_id_cardano_output, Cardano hinge,\n");
%!       fputs(fd, "                node_id_shaft_output,\n");
%!       fputs(fd, "                position, reference, ref_id_shaft_output, null,\n");
%!       fputs(fd, "                orientation, reference, ref_id_shaft_output, eye,\n");
%!       fputs(fd, "                node_id_shaft_axis, reference, ref_id_shaft_output, null,\n");
%!       fputs(fd, "                orientation, reference, ref_id_shaft_axis_dynamic, eye;\n");
%!       fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%! # options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [res.t, res.trajectory, res.deformation,res.velocity,res.acceleration,res.node_id, res.force_id, res.force_node_id, res.orientation_description]=mbdyn_post_load_output_struct (options.output_file);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   node_id = int32([2001, 2002, 2003]);
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
%!     plot(res.t, Wrel{i}(1, :) * 30 / pi, sprintf("-;W%dx;", i));
%!   endfor
%!   xlabel("t [s]");
%!   ylabel("n[rpm]");
%!   title("cardano joint omega1");
%!   grid on;
%!   grid minor on;
%!   endif
%!   omega = log_dat.vars.omega;
%!   idx_t = res.t > 0.5 * log_dat.vars.t1;
%!   rel_err = -1;
%!   for i=[1,3]
%!     rel_err = max([rel_err, max(norm(Wrel{i}(:, idx_t) - [omega; zeros(2, 1)], "cols")) / abs(omega)]);
%!   endfor
%!   rel_error(idx_method, idx_autodiff) = rel_err;
%! %unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! %end_unwind_protect
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
%! tol = 1e-5;
%! assert_simple(all(rel_error(:) < tol));

%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
