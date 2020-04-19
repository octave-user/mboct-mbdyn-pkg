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
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_struct_XXXXXX"));
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
%!   assert(all(output_flag));
%!   assert(node_id, int32(1));
%!   assert(t, (0:0.1:1.1).', 1e-5);
%!   assert(trajectory{1}, ...
%!          [0.5 * F / m * t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g * t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert(velocity{1}, ...
%!         [t * F / m, ...
%!          zeros(numel(t), 1), g * t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert(acceleration{1}, ...
%!         [F / m * ones(numel(t), 1), ...
%!          zeros(numel(t), 1), ...
%!          repmat(g, numel(t), 1), ...
%!          zeros(numel(t), 3)], ...
%!          1e-5);
%!  assert(force{1}(:, 1:3), [repmat(F, numel(t), 1), zeros(numel(t), 2)]);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_struct_XXXXXX"));
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
%!   assert(all(output_flag));
%!   assert(node_id, int32(1));
%!   assert(t, (0:0.1:1.1).', 1e-5);
%!   assert(trajectory{1}, ...
%!          [0.5 * F / m * t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g * t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert(velocity{1}, ...
%!         [t * F / m, ...
%!          zeros(numel(t), 1), g * t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert(acceleration{1}, ...
%!         [F / m * ones(numel(t), 1), ...
%!          zeros(numel(t), 1), ...
%!          repmat(g, numel(t), 1), ...
%!          zeros(numel(t), 3)], ...
%!          1e-5);
%!  assert(force{1}(:, 1:3), [repmat(F, numel(t), 1), zeros(numel(t), 2)]);
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
