## Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{force_id}, @var{force_node_id}, @var{force}, @var{arm}] = mbdyn_post_load_output_frc(@var{mbdyn_filename}, @var{filter_force_id}, @var{append_rows})
##
## Loads the MBDyn .frc file named "<@var{mbdyn_filename}.frc>".
##
## @var{filter_node_id} @dots{} Load only the data of the corresponding nodes if present.
## If empty, the data of all nodes is loaded.
##
## @var{append_rows} @dots{} Hint for memory reallocation.
##
## @var{force_id} @dots{} Structural force element number
##
## @var{force_node_id} @dots{} Structural node number
##
## @var{force} @dots{} Structural force values
##
## @var{arm} @dots{} The arm of the force, in the global frame (i.e. referred to point @{0, 0, 0@} and oriented as the global frame).
##
## @end deftypefn

function [force_id, force_node_id, force, arm] = mbdyn_post_load_output_frc(mbdyn_filename, filter_force_id, append_rows, auto_resize_rows)
  if (nargin < 1 || nargout > 4)
    print_usage();
  endif

  if (nargin < 2)
    filter_force_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  if (nargin < 4)
    auto_resize_rows = true;
  endif

  frc_filename = mbdyn_post_output_filename(mbdyn_filename, ".frc");
  
  column_count = 0;

  if (nargout >= 2)
    ++column_count;
  endif

  if (nargout >= 3)
    column_count += 3;
  endif

  if (nargout >= 4)
    column_count += 3;
  endif

  [force_id, data] = mbdyn_post_load_output(frc_filename, column_count, filter_force_id, append_rows, 4, 1, auto_resize_rows);

  if (numel(data) > 0)
    for i=1:numel(data)
      if (nargout >= 2)
        if (columns(data{i}) >= 8)
          force_node_id{i} = data{i}(:, [1, 8]);
        else
          force_node_id{i} = data{i}(:, 1);
        endif
      endif

      if (nargout >= 3)
        if (columns(data{i}) >= 11)
          force{i} = data{i}(:, [2:4, 9:11]);
        elseif (columns(data{i}) >= 4)
          force{i} = data{i}(:, 2:4);          
        elseif (columns(data{i}) >= 2)
          force{i} = data{i}(:, 2:end);
        else
          force{i} = [];
        endif
      endif

      if (nargout >= 4)
        if (columns(data{i}) >= 14)
          arm{i} = data{i}(:, [5:7, 12:14]);
        elseif (columns(data{i}) >= 7)
          arm{i} = data{i}(:, 5:7);
        else
          arm{i} = [];
        endif
      endif
    endfor
  else
    force = {};
    force_node_id = {};
    arm = {};
  endif
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_frc_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 200;\n");
%!     fputs(fd, " set: real F3 = 300;\n");
%!     fputs(fd, " set: real m1 = 1;\n");
%!     fputs(fd, " set: real m2 = 2;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
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
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 3;\n");
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "                 euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "                 euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., const, F2;\n");
%!     fputs(fd, "         force: 3, absolute internal, 1, position, 30.1, 30.2, 30.3, 2, position, 40.1, 40.2, 40.3, 0., 0., 1., const, F3;\n");
%!     fputs(fd, "         joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "         joint: 2, clamp, 2, node, node;\n");
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
%!   t = mbdyn_post_load_output_out(options.output_file);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [force_id, force_node_id, force, arm] = mbdyn_post_load_output_frc(options.output_file, 1:3, numel(t));
%!   forces = mbdyn_post_load_log_force(options.output_file);
%!   X1 = log_dat.nodes(1).X0;
%!   X2 = log_dat.nodes(2).X0;
%!   R1 = log_dat.nodes(1).R0;
%!   R2 = log_dat.nodes(2).R0;
%!   F1 = log_dat.vars.F1;
%!   F2 = log_dat.vars.F2;
%!   F3 = log_dat.vars.F3;
%!   l1 = forces(1).arm1;
%!   l2 = forces(2).arm1;
%!   l3_1 = forces(3).arm1;
%!   l3_2 = forces(3).arm2;
%!   assert_simple(force_id, int32([1,2,3]));
%!   assert_simple(force_node_id{1}, repmat(1, numel(t), 1));
%!   assert_simple(force_node_id{2}, repmat(2, numel(t), 1));
%!   assert_simple(force_node_id{3}, repmat([1, 2], numel(t), 1));
%!   assert_simple(force{1}, repmat([F1, 0, 0], numel(t), 1));
%!   assert_simple(force{2}, repmat([0, F2, 0], numel(t), 1));
%!   assert_simple(force{3}, repmat([0, 0, F3, 0, 0, -F3], numel(t), 1));
%!   tol = eps^0.3;
%!   assert_simple(arm{1}, repmat((R1 * l1 + X1).', numel(t), 1), tol);
%!   assert_simple(arm{2}, repmat((R2 * l2 + X2).', numel(t), 1), tol);
%!   assert_simple(arm{3}, repmat([(R1 * l3_1 + X1).', (R2 * l3_2 + X2).'], numel(t), 1), tol);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_frc_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 200;\n");
%!     fputs(fd, " set: real F3 = 300;\n");
%!     fputs(fd, " set: real m1 = 1;\n");
%!     fputs(fd, " set: real m2 = 2;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
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
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 3;\n");
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "                 euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "                 euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., const, F2;\n");
%!     fputs(fd, "         force: 3, absolute internal, 1, position, 30.1, 30.2, 30.3, 2, position, 40.1, 40.2, 40.3, 0., 0., 1., const, F3;\n");
%!     fputs(fd, "         joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "         joint: 2, clamp, 2, node, node;\n");
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
%!   t = mbdyn_post_load_output_out(options.output_file);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [force_id, force_node_id, force, arm] = mbdyn_post_load_output_frc(options.output_file, 1:3, numel(t));
%!   forces = mbdyn_post_load_log_force(options.output_file);
%!   X1 = log_dat.nodes(1).X0;
%!   X2 = log_dat.nodes(2).X0;
%!   R1 = log_dat.nodes(1).R0;
%!   R2 = log_dat.nodes(2).R0;
%!   F1 = log_dat.vars.F1;
%!   F2 = log_dat.vars.F2;
%!   F3 = log_dat.vars.F3;
%!   l1 = forces(1).arm1;
%!   l2 = forces(2).arm1;
%!   l3_1 = forces(3).arm1;
%!   l3_2 = forces(3).arm2;
%!   assert_simple(force_id, int32([1,2,3]));
%!   assert_simple(force_node_id{1}, repmat(1, numel(t), 1));
%!   assert_simple(force_node_id{2}, repmat(2, numel(t), 1));
%!   assert_simple(force_node_id{3}, repmat([1, 2], numel(t), 1));
%!   assert_simple(force{1}, repmat([F1, 0, 0], numel(t), 1));
%!   assert_simple(force{2}, repmat([0, F2, 0], numel(t), 1));
%!   assert_simple(force{3}, repmat([0, 0, F3, 0, 0, -F3], numel(t), 1));
%!   tol = eps^0.3;
%!   assert_simple(arm{1}, repmat((R1 * l1 + X1).', numel(t), 1), tol);
%!   assert_simple(arm{2}, repmat((R2 * l2 + X2).', numel(t), 1), tol);
%!   assert_simple(arm{3}, repmat([(R1 * l3_1 + X1).', (R2 * l3_2 + X2).'], numel(t), 1), tol);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_frc_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 200;\n");
%!     fputs(fd, " set: real F3 = 300;\n");
%!     fputs(fd, " set: real m1 = 1;\n");
%!     fputs(fd, " set: real m2 = 2;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
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
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 3;\n");
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "                 euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "                 euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., const, F2;\n");
%!     fputs(fd, "         force: 3, absolute internal, 1, position, 30.1, 30.2, 30.3, 2, position, 40.1, 40.2, 40.3, 0., 0., 1., const, F3;\n");
%!     fputs(fd, "         joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "         joint: 2, clamp, 2, node, node;\n");
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
%!   t = mbdyn_post_load_output_out(options.output_file);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [force_id, force_node_id, force, arm] = mbdyn_post_load_output_frc(options.output_file, 1:3, numel(t));
%!   forces = mbdyn_post_load_log_force(options.output_file);
%!   X1 = log_dat.nodes(1).X0;
%!   X2 = log_dat.nodes(2).X0;
%!   R1 = log_dat.nodes(1).R0;
%!   R2 = log_dat.nodes(2).R0;
%!   F1 = log_dat.vars.F1;
%!   F2 = log_dat.vars.F2;
%!   F3 = log_dat.vars.F3;
%!   l1 = forces(1).arm1;
%!   l2 = forces(2).arm1;
%!   l3_1 = forces(3).arm1;
%!   l3_2 = forces(3).arm2;
%!   assert_simple(force_id, int32([1,2,3]));
%!   assert_simple(force_node_id{1}, repmat(1, numel(t), 1));
%!   assert_simple(force_node_id{2}, repmat(2, numel(t), 1));
%!   assert_simple(force_node_id{3}, repmat([1, 2], numel(t), 1));
%!   assert_simple(force{1}, repmat([F1, 0, 0], numel(t), 1));
%!   assert_simple(force{2}, repmat([0, F2, 0], numel(t), 1));
%!   assert_simple(force{3}, repmat([0, 0, F3, 0, 0, -F3], numel(t), 1));
%!   tol = eps^0.3;
%!   assert_simple(arm{1}, repmat((R1 * l1 + X1).', numel(t), 1), tol);
%!   assert_simple(arm{2}, repmat((R2 * l2 + X2).', numel(t), 1), tol);
%!   assert_simple(arm{3}, repmat([(R1 * l3_1 + X1).', (R2 * l3_2 + X2).'], numel(t), 1), tol);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_frc_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 200;\n");
%!     fputs(fd, " set: real F3 = 300;\n");
%!     fputs(fd, " set: real m1 = 1;\n");
%!     fputs(fd, " set: real m2 = 2;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
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
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 3;\n");
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "                 euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "                 euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., const, F2;\n");
%!     fputs(fd, "         force: 3, absolute internal, 1, position, 30.1, 30.2, 30.3, 2, position, 40.1, 40.2, 40.3, 0., 0., 1., const, F3;\n");
%!     fputs(fd, "         joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "         joint: 2, clamp, 2, node, node;\n");
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
%!   t = mbdyn_post_load_output_out(options.output_file);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [force_id, force_node_id, force, arm] = mbdyn_post_load_output_frc(options.output_file, 1:3, numel(t));
%!   forces = mbdyn_post_load_log_force(options.output_file);
%!   X1 = log_dat.nodes(1).X0;
%!   X2 = log_dat.nodes(2).X0;
%!   R1 = log_dat.nodes(1).R0;
%!   R2 = log_dat.nodes(2).R0;
%!   F1 = log_dat.vars.F1;
%!   F2 = log_dat.vars.F2;
%!   F3 = log_dat.vars.F3;
%!   l1 = forces(1).arm1;
%!   l2 = forces(2).arm1;
%!   l3_1 = forces(3).arm1;
%!   l3_2 = forces(3).arm2;
%!   assert_simple(force_id, int32([1,2,3]));
%!   assert_simple(force_node_id{1}, repmat(1, numel(t), 1));
%!   assert_simple(force_node_id{2}, repmat(2, numel(t), 1));
%!   assert_simple(force_node_id{3}, repmat([1, 2], numel(t), 1));
%!   assert_simple(force{1}, repmat([F1, 0, 0], numel(t), 1));
%!   assert_simple(force{2}, repmat([0, F2, 0], numel(t), 1));
%!   assert_simple(force{3}, repmat([0, 0, F3, 0, 0, -F3], numel(t), 1));
%!   tol = eps^0.3;
%!   assert_simple(arm{1}, repmat((R1 * l1 + X1).', numel(t), 1), tol);
%!   assert_simple(arm{2}, repmat((R2 * l2 + X2).', numel(t), 1), tol);
%!   assert_simple(arm{3}, repmat([(R1 * l3_1 + X1).', (R2 * l3_2 + X2).'], numel(t), 1), tol);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
