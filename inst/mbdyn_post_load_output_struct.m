## Copyright (C) 2014(-2020) Reinhard <octave-user@a1.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} [@var{t}, @var{trajectory}, @var{deformation}, @var{velocity}, @var{acceleration}, @var{node_id}, @var{force}, @var{force_id}, @var{force_node_id}, @var{orientation_description}] = mbdyn_post_load_output_struct(@var{mbdyn_filename}, @var{filter_node_id}, @var{filter_force_id})
##
## Loads data from MBDyn .mov, .out and .frc files ("<@var{mbdyn_filename}.mov>", "<@var{mbdyn_filename}.out>" and "<@var{mbdyn_filename}.frc>").
##
## @var{filter_node_id} @dots{} Load only the data of the corresponding nodes if present.
## If empty the data of all nodes is loaded.
##
## @var{t} @dots{} Simulation time.
##
## @var{trajectory} @dots{} Matrix of position and orientation of structural nodes.
##
## @var{deformation} @dots{} The difference between actual position/orientation and initial position/orientation.
##
## @var{velocity} @dots{} Matrix of velocity and angular velocity of structural nodes.
##
## @var{acceleration} @dots{} Matrix of acceleration and angular acceleration.
##
## @var{node_id} @dots{} Structural node numbers.
##
## @var{force} @dots{} Structural force values.
##
## @var{force_id} @dots{} Structural force element number.
##
## @var{force_node_id} @dots{} Structural node number.
##
## @var{orientation_description} @dots{} Cell array of strings. One of ("euler123", "euler321", "euler313", "phi")
##
## @end deftypefn

function [t, trajectory, deformation, velocity, acceleration, node_id, force,  force_id, force_node_id, orientation_description, t_trc, OutputFlag] = mbdyn_post_load_output_struct(mbdyn_filename, filter_node_id, filter_force_id, auto_resize_rows)

  if (nargin < 1 || nargout > 12)
    print_usage();
  endif

  if (nargin < 2)
    filter_node_id = [];
  endif

  if (nargin < 3)
    filter_force_id = [];
  endif

  if (nargin < 4)
    auto_resize_rows = false;
  endif

  [t_trc, TStep, NIter, ResErr, SolErr, SolConv, OutputFlag] = mbdyn_post_load_output_out(mbdyn_post_output_filename(mbdyn_filename), 1024, auto_resize_rows);

  t = t_trc(find(OutputFlag));

  if (nargout >= 2)
    [node_id, trajectory, velocity, acceleration, orientation_description, deformation] = mbdyn_post_load_output_mov(mbdyn_filename, filter_node_id, numel(t), 6, auto_resize_rows);
  endif

  if (nargout >= 7)
    if (2 == exist(mbdyn_post_output_filename(mbdyn_filename, ".frc"), "file"))
      [force_id, force_node_id, force] = mbdyn_post_load_output_frc(mbdyn_filename, filter_force_id, numel(t), auto_resize_rows);
    else
      force_id = [];
      force_node_id = [];
      force = {};
    endif
  endif
endfunction

%!test
%! ## TEST1
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_struct_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-1;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         nonlinear solver: nox, modified, 10,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 12,\n");
%!     fputs(fd, "             krylov subspace size, 12;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, null, 1., 0., 0, 100.;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., 9.81;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation, velocity, acceleration, node_id, force,  force_id, force_node_id, orientation_description, t_trc, output_flag] = mbdyn_post_load_output_struct(fname, [1], [1], false);
%!   g = -9.81;
%!   F = 100;
%!   m = 1;
%!   assert_simple(all(output_flag));
%!   assert_simple(node_id, int32(1));
%!   assert_simple(t, (0:0.1:1.1).', 1e-5);
%!   assert_simple(trajectory{1}, ...
%!          [0.5 * F / m * t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g * t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert_simple(velocity{1}, ...
%!         [t * F / m, ...
%!          zeros(numel(t), 1), g * t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert_simple(acceleration{1}, ...
%!         [F / m * ones(numel(t), 1), ...
%!          zeros(numel(t), 1), ...
%!          repmat(g, numel(t), 1), ...
%!          zeros(numel(t), 3)], ...
%!          1e-5);
%!  assert_simple(force{1}(:, 1:3), [repmat(F, numel(t), 1), zeros(numel(t), 2)]);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST2
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_struct_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-1;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         nonlinear solver: nox, modified, 10,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 12,\n");
%!     fputs(fd, "             krylov subspace size, 12;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, null, 1., 0., 0, 100.;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., 9.81;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation, velocity, acceleration, node_id, force,  force_id, force_node_id, orientation_description, t_trc, output_flag] = mbdyn_post_load_output_struct(fname, [1], [1], false);
%!   g = -9.81;
%!   F = 100;
%!   m = 1;
%!   assert_simple(all(output_flag));
%!   assert_simple(node_id, int32(1));
%!   assert_simple(t, (0:0.1:1.1).', 1e-5);
%!   assert_simple(trajectory{1}, ...
%!          [0.5 * F / m * t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g * t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert_simple(velocity{1}, ...
%!         [t * F / m, ...
%!          zeros(numel(t), 1), g * t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert_simple(acceleration{1}, ...
%!         [F / m * ones(numel(t), 1), ...
%!          zeros(numel(t), 1), ...
%!          repmat(g, numel(t), 1), ...
%!          zeros(numel(t), 3)], ...
%!          1e-5);
%!  assert_simple(force{1}(:, 1:3), [repmat(F, numel(t), 1), zeros(numel(t), 2)]);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST3
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_struct_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-1;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, null, 1., 0., 0, 100.;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., 9.81;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation, velocity, acceleration, node_id, force,  force_id, force_node_id, orientation_description, t_trc, output_flag] = mbdyn_post_load_output_struct(fname, [1], [1], false);
%!   g = -9.81;
%!   F = 100;
%!   m = 1;
%!   assert_simple(all(output_flag));
%!   assert_simple(node_id, int32(1));
%!   assert_simple(t, (0:0.1:1.1).', 1e-5);
%!   assert_simple(trajectory{1}, ...
%!          [0.5 * F / m * t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g * t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert_simple(velocity{1}, ...
%!         [t * F / m, ...
%!          zeros(numel(t), 1), g * t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert_simple(acceleration{1}, ...
%!         [F / m * ones(numel(t), 1), ...
%!          zeros(numel(t), 1), ...
%!          repmat(g, numel(t), 1), ...
%!          zeros(numel(t), 3)], ...
%!          1e-5);
%!  assert_simple(force{1}(:, 1:3), [repmat(F, numel(t), 1), zeros(numel(t), 2)]);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST4
%! ## TRACTA JOINT
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_struct_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " ## input parameters:\n");
%!     fputs(fd, " set: real N = 360;                   ## number of time steps per revolution\n");
%!     fputs(fd, " set: real omega = 40000 * pi / 30.;  ## input speed [rad/s9\n");
%!     fputs(fd, " set: real Phiz = 80. * pi / 180.;    ## prescribed angle between input shaft and output shaft [rad]\n");
%!     fputs(fd, " set: real t1 = 2. * pi / abs(omega); ## ramp up time to raise the angle from zero to Phiz [s]\n");
%!     fputs(fd, " set: real t2 = 10. * t1;             ## final time [s]\n");
%!     fputs(fd, " set: real dt = t1 / N;               ## time step [s]\n");
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
%!     fputs(fd, "         method: msstc5, 0.6;\n");
%!     fputs(fd, "         linear solver:naive, colamd;\n");
%!     fputs(fd, "         nonlinear solver: nox, modified, 100,  keep jacobian matrix, inner iterations before assembly, 6, jacobian operator, newton krylov, forcing term, type2;\n");
%!     fputs(fd, "         threads: disable;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "       model: static;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
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
%!   options.logfile = [fname, ".stdout"];
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
%!   figure_list();
%!   t1 = log_dat.vars.t1;
%!   omega = log_dat.vars.omega;
%!   idx_t = find(res.t > t1);
%!   tol = 1e-6;
%!   for i=[1,4]
%!     assert_simple(max(norm(Wrel{i}(:, idx_t) - [omega; zeros(2, 1)], "cols")) < tol * abs(omega));
%!   endfor
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_struct_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-1;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");

%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, null, 1., 0., 0, 100.;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., 9.81;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation, velocity, acceleration, node_id, force,  force_id, force_node_id, orientation_description, t_trc, output_flag] = mbdyn_post_load_output_struct(fname, [1], [1], false);
%!   g = -9.81;
%!   F = 100;
%!   m = 1;
%!   assert_simple(all(output_flag));
%!   assert_simple(node_id, int32(1));
%!   assert_simple(t, (0:0.1:1.1).', 1e-5);
%!   assert_simple(trajectory{1}, ...
%!          [0.5 * F / m * t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g * t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert_simple(velocity{1}, ...
%!         [t * F / m, ...
%!          zeros(numel(t), 1), g * t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert_simple(acceleration{1}, ...
%!         [F / m * ones(numel(t), 1), ...
%!          zeros(numel(t), 1), ...
%!          repmat(g, numel(t), 1), ...
%!          zeros(numel(t), 3)], ...
%!          1e-5);
%!  assert_simple(force{1}(:, 1:3), [repmat(F, numel(t), 1), zeros(numel(t), 2)]);
%!  figure("visible", "off");
%!  subplot(3, 1, 1);
%!  hold on;
%!  for i=1:3
%!    plot(t, trajectory{1}(:, i), sprintf("-;x%d;%d", i, i));
%!  endfor
%!  xlabel("t [s]");
%!  ylabel("x [m]");
%!  grid on;
%!  grid minor on;
%!  title("trajectory");
%!  subplot(3, 1, 2);
%!  hold on;
%!  for i=1:3
%!    plot(t, velocity{1}(:, i), sprintf("-;v%d;%d", i, i));
%!  endfor
%!  xlabel("t [s]");
%!  ylabel("v [m/s]");
%!  grid on;
%!  grid minor on;
%!  title("velocity");
%!  subplot(3, 1, 3);
%!  hold on;
%!  for i=1:3
%!    plot(t, acceleration{1}(:, i), sprintf("-;a%d;%d", i, i));
%!  endfor
%!  xlabel("t [s]");
%!  ylabel("a [m/s^2]");
%!  grid on;
%!  grid minor on;
%!  title("acceleration");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
