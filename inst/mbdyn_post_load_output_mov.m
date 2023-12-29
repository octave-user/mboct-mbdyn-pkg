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
## @deftypefn {Function File} [@var{node_id}, @var{trajectory}, @var{velocity}, @var{acceleration}, @var{orientation_description}, @var{deformation}]=mbdyn_post_load_output_mov (@var{mbdyn_filename}, @var{filter_node_id}, @var{append_rows})
##
## Loads the MBDyn .mov file named "<@var{mbdyn_filename}.mov>".
##
## @var{filter_node_id} @dots{} Load only the data of the corresponding nodes if present.
## If empty, the data of all nodes is loaded.
##
## @var{append_rows} @dots{} This parameter is used in order to optimize the memory allocation and should be set to the expected number of rows per node.
##
## @var{append_columns} @dots{} This parameter is used in order to optimize the memory allocation and should be set to the expected number of columns.
##
## @var{node_id} @dots{} Node number of structural nodes.
##
## @var{trajectory} @dots{} Matrix of position and orientation of structural nodes.
##
## @var{velocity} @dots{} Matrix of velocity and angular velocity of structural nodes.
##
## @var{acceleration} @dots{} Matrix of acceleration and angular acceleration.
##
## @var{orientation_description} @dots{} Cell array of strings. One of ("euler123", "euler321", "euler313", "phi").
##
## @var{deformation} @dots{} The difference between actual position/orientation and initial position/orientation.
##
## @end deftypefn

function [node_id, trajectory, velocity, acceleration, orientation_description, deformation] = mbdyn_post_load_output_mov(mbdyn_filename, filter_node_id, append_rows, append_columns, auto_resize_rows)
  if (nargin < 1 || nargout > 6)
    print_usage();
  endif

  if (nargin < 2)
    filter_node_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  if (nargin < 4)
    append_columns = 6;
  endif

  if (nargin < 5)
    auto_resize_rows = true;
  endif

  if (~ischar(mbdyn_filename))
    error("mbdyn_filename must be a string!");
  endif

  if (~isscalar(append_rows))
    error("append_rows must be a scalar!");
  endif

  mov_filename = mbdyn_post_output_filename(mbdyn_filename, ".mov");
  log_filename = mbdyn_post_output_filename(mbdyn_filename, ".log");

  column_count = 0;

  if (nargout >= 2)
    column_count = 6;
  endif

  if (nargout >= 3)
    column_count = 12;
  endif

  nodes = mbdyn_post_load_log_node(log_filename);

  node_labels = [nodes.label];

  [node_id, data] = mbdyn_post_load_output(mov_filename, column_count, filter_node_id, append_rows, append_columns, 1, auto_resize_rows);

  node_label_idx = zeros(1,length(node_id));

  for i=1:length(node_id)
    node_label_idx_i = find(node_id(i) == node_labels);

    if (length(node_label_idx_i) == 1)
      node_label_idx(i) = node_label_idx_i;
    else
      error("node_id %d not found in file \"%s\"!", node_id(i), log_filename);
    endif
  endfor

  orientation_description = { nodes(node_label_idx).orientation_description };

  if (length(data) == 0)
    if (nargout >= 2)
      trajectory = {};
    endif

    if (nargout >=  3)
      velocity = {};
    endif

    if (nargout >= 4)
      acceleration = {};
    endif

    if (nargout >= 6)
      deformation = {};
    endif
  else
    for i=1:length(data)
      if (nargout >= 2)
        trajectory{i} = data{i}(:, 1:6);

        if (nargout >=  3)
          velocity{i} = data{i}(:, 7:12);
        endif

        if (nargout >= 4)
          if (columns(data{i}) >= 18)
            acceleration{i} = data{i}(:, 13:18);
          else
            acceleration{i} = zeros(rows(data{i}), 0);
          endif
        endif

        switch (orientation_description{i})
          case {"euler123", "euler321", "euler313"}
            trajectory{i}(:,4:6) *= pi / 180; % euler angles are output in degrees
          case "phi"
            ## rotation vectors are already in radians
        endswitch
      endif

      if (nargout >= 6)
        node_idx_i = find(node_id(i) == [nodes.label]);

        if (length(node_idx_i) ~= 1)
          error("node_id(%d)=%d must appear exactly one times in the .log file!",i,node_id(i));
        endif

        deformation{i} = trajectory{i} - repmat([nodes(node_idx_i).X0.', ...
                                                 nodes(node_idx_i).Phi0.'], ...
                                                rows(trajectory{i}),1);

      endif
    endfor
  endif
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_mov_XXXXXX"));
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
%!   [t, dt, niter, reserr, solerr, solconv, output_flag] = mbdyn_post_load_output_out(fname, 1024, false);
%!   [node_id, trajectory, velocity, acceleration, orientation_description, deformation] = mbdyn_post_load_output_mov(options.output_file, [1], numel(t), 6, false);
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
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_mov_XXXXXX"));
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
%!   [t, dt, niter, reserr, solerr, solconv, output_flag] = mbdyn_post_load_output_out(fname, 1024, false);
%!   [node_id, trajectory, velocity, acceleration, orientation_description, deformation] = mbdyn_post_load_output_mov(options.output_file, [1], numel(t), 6, false);
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
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_mov_XXXXXX"));
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
%!   [t, dt, niter, reserr, solerr, solconv, output_flag] = mbdyn_post_load_output_out(fname, 1024, false);
%!   [node_id, trajectory, velocity, acceleration, orientation_description, deformation] = mbdyn_post_load_output_mov(options.output_file, [1], numel(t), 6, false);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_mov_XXXXXX"));
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
%!   [t, dt, niter, reserr, solerr, solconv, output_flag] = mbdyn_post_load_output_out(fname, 1024, false);
%!   [node_id, trajectory, velocity, acceleration, orientation_description, deformation] = mbdyn_post_load_output_mov(options.output_file, [1], numel(t), 6, false);
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
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
