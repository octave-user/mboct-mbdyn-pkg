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
## @deftypefn {Function File} [@var{joint_id}, @var{local_reaction}, @var{global_reaction}] = mbdyn_post_load_output_jnt(@var{mbdyn_filename}, @var{filter_joint_id}, @var{append_rows})
##
## Loads the MBDyn .jnt file named "<@var{mbdyn_filename}.jnt>".
##
## @var{joint_id} @dots{} The label of joint element.
##
## @var{local_reaction} @dots{} The three components of the reaction force and reaction couple in a local reference frame.
##
## @var{global_reaction} @dots{} The three components of the reaction force and reaction couple in the global frame.
##
## @end deftypefn

function [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(mbdyn_filename, filter_joint_id, append_rows, auto_resize_rows)
  if (nargin < 1 || nargin > 4 || nargout > 3)
    print_usage();
  endif

  if (nargin < 2)
    filter_joint_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  if (nargin < 4)
    auto_resize_rows = true;
  endif

  jnt_filename = mbdyn_post_output_filename(mbdyn_filename, ".jnt");
  
  column_count = 0;
  
  if (nargout >= 2)
    column_count = 6;
  endif
  
  if (nargout >= 3)
    column_count = 12;
  endif
  
  [joint_id, data] = mbdyn_post_load_output(jnt_filename, column_count, filter_joint_id, append_rows, 12, 1, auto_resize_rows);
  
  if (length(joint_id) == 0)
    local_reaction = {};
    global_reaction = {};
  else
    for i=1:length(data)
      if (nargout >= 2)
	local_reaction{i} = data{i}(:, 1:6);
      endif
      if (nargout >= 3)
	global_reaction{i} = data{i}(:, 7:12);
      endif
    endfor
  endif
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_jnt_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 200;\n");
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
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static,\n");
%!     fputs(fd, "                 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "                 euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 2, static,\n");
%!     fputs(fd, "                 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "                 euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., const, F2;\n");
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
%!   forces = mbdyn_post_load_log_force(options.output_file);
%!   [joint_id, Fl, Fg] = mbdyn_post_load_output_jnt(options.output_file, 1:2, numel(t));
%!   X1 = log_dat.nodes(1).X0;
%!   X2 = log_dat.nodes(2).X0;
%!   R1 = log_dat.nodes(1).R0;
%!   R2 = log_dat.nodes(2).R0;
%!   F1 = repmat([log_dat.vars.F1, 0, 0], numel(t), 1);
%!   F2 = repmat([0, log_dat.vars.F2, 0], numel(t), 1);
%!   l1 = forces(1).arm1;
%!   l2 = forces(2).arm1;
%!   M1 = cross(repmat((R1 * l1).', numel(t), 1), F1);
%!   M2 = cross(repmat((R2 * l2).', numel(t), 1), F2);
%!   tol = eps^0.3;
%!   assert_simple(Fg{1}(:, 1:3), F1, tol * norm(F1));
%!   assert_simple(Fg{1}(:, 4:6), M1, tol * norm(M1));
%!   assert_simple(Fl{1}(:, 1:3), F1 * R1, tol * norm(F1));
%!   assert_simple(Fl{1}(:, 4:6), M1 * R1, tol * norm(M1));
%!   assert_simple(Fg{2}(:, 1:3), F2, tol * norm(F2));
%!   assert_simple(Fg{2}(:, 4:6), M2, tol * norm(M2));
%!   assert_simple(Fl{2}(:, 1:3), F2 * R2, tol * norm(F2));
%!   assert_simple(Fl{2}(:, 4:6), M2 * R2, tol * norm(M2));
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_jnt_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 200;\n");
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
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static,\n");
%!     fputs(fd, "                 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "                 euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 2, static,\n");
%!     fputs(fd, "                 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "                 euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., const, F2;\n");
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
%!   forces = mbdyn_post_load_log_force(options.output_file);
%!   [joint_id, Fl, Fg] = mbdyn_post_load_output_jnt(options.output_file, 1:2, numel(t));
%!   X1 = log_dat.nodes(1).X0;
%!   X2 = log_dat.nodes(2).X0;
%!   R1 = log_dat.nodes(1).R0;
%!   R2 = log_dat.nodes(2).R0;
%!   F1 = repmat([log_dat.vars.F1, 0, 0], numel(t), 1);
%!   F2 = repmat([0, log_dat.vars.F2, 0], numel(t), 1);
%!   l1 = forces(1).arm1;
%!   l2 = forces(2).arm1;
%!   M1 = cross(repmat((R1 * l1).', numel(t), 1), F1);
%!   M2 = cross(repmat((R2 * l2).', numel(t), 1), F2);
%!   tol = eps^0.3;
%!   assert_simple(Fg{1}(:, 1:3), F1, tol * norm(F1));
%!   assert_simple(Fg{1}(:, 4:6), M1, tol * norm(M1));
%!   assert_simple(Fl{1}(:, 1:3), F1 * R1, tol * norm(F1));
%!   assert_simple(Fl{1}(:, 4:6), M1 * R1, tol * norm(M1));
%!   assert_simple(Fg{2}(:, 1:3), F2, tol * norm(F2));
%!   assert_simple(Fg{2}(:, 4:6), M2, tol * norm(M2));
%!   assert_simple(Fl{2}(:, 1:3), F2 * R2, tol * norm(F2));
%!   assert_simple(Fl{2}(:, 4:6), M2 * R2, tol * norm(M2));
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_jnt_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 200;\n");
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
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static,\n");
%!     fputs(fd, "                 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "                 euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 2, static,\n");
%!     fputs(fd, "                 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "                 euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., const, F2;\n");
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
%!   forces = mbdyn_post_load_log_force(options.output_file);
%!   [joint_id, Fl, Fg] = mbdyn_post_load_output_jnt(options.output_file, 1:2, numel(t));
%!   X1 = log_dat.nodes(1).X0;
%!   X2 = log_dat.nodes(2).X0;
%!   R1 = log_dat.nodes(1).R0;
%!   R2 = log_dat.nodes(2).R0;
%!   F1 = repmat([log_dat.vars.F1, 0, 0], numel(t), 1);
%!   F2 = repmat([0, log_dat.vars.F2, 0], numel(t), 1);
%!   l1 = forces(1).arm1;
%!   l2 = forces(2).arm1;
%!   M1 = cross(repmat((R1 * l1).', numel(t), 1), F1);
%!   M2 = cross(repmat((R2 * l2).', numel(t), 1), F2);
%!   tol = eps^0.3;
%!   assert_simple(Fg{1}(:, 1:3), F1, tol * norm(F1));
%!   assert_simple(Fg{1}(:, 4:6), M1, tol * norm(M1));
%!   assert_simple(Fl{1}(:, 1:3), F1 * R1, tol * norm(F1));
%!   assert_simple(Fl{1}(:, 4:6), M1 * R1, tol * norm(M1));
%!   assert_simple(Fg{2}(:, 1:3), F2, tol * norm(F2));
%!   assert_simple(Fg{2}(:, 4:6), M2, tol * norm(M2));
%!   assert_simple(Fl{2}(:, 1:3), F2 * R2, tol * norm(F2));
%!   assert_simple(Fl{2}(:, 4:6), M2 * R2, tol * norm(M2));
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_jnt_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 200;\n");
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
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static,\n");
%!     fputs(fd, "                 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "                 euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 2, static,\n");
%!     fputs(fd, "                 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "                 euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., const, F2;\n");
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
%!   forces = mbdyn_post_load_log_force(options.output_file);
%!   [joint_id, Fl, Fg] = mbdyn_post_load_output_jnt(options.output_file, 1:2, numel(t));
%!   X1 = log_dat.nodes(1).X0;
%!   X2 = log_dat.nodes(2).X0;
%!   R1 = log_dat.nodes(1).R0;
%!   R2 = log_dat.nodes(2).R0;
%!   F1 = repmat([log_dat.vars.F1, 0, 0], numel(t), 1);
%!   F2 = repmat([0, log_dat.vars.F2, 0], numel(t), 1);
%!   l1 = forces(1).arm1;
%!   l2 = forces(2).arm1;
%!   M1 = cross(repmat((R1 * l1).', numel(t), 1), F1);
%!   M2 = cross(repmat((R2 * l2).', numel(t), 1), F2);
%!   tol = eps^0.3;
%!   assert_simple(Fg{1}(:, 1:3), F1, tol * norm(F1));
%!   assert_simple(Fg{1}(:, 4:6), M1, tol * norm(M1));
%!   assert_simple(Fl{1}(:, 1:3), F1 * R1, tol * norm(F1));
%!   assert_simple(Fl{1}(:, 4:6), M1 * R1, tol * norm(M1));
%!   assert_simple(Fg{2}(:, 1:3), F2, tol * norm(F2));
%!   assert_simple(Fg{2}(:, 4:6), M2, tol * norm(M2));
%!   assert_simple(Fl{2}(:, 1:3), F2 * R2, tol * norm(F2));
%!   assert_simple(Fl{2}(:, 4:6), M2 * R2, tol * norm(M2));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
