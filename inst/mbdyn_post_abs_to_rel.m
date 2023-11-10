## Copyright (C) 2011(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} mbdyn_post_abs_to_rel(@var{ref_node}, @var{abs_mov_file}, @var{rel_mov_file}, @var{ref_only}, @var{input_mode}, @var{output_mode})
##
## Converts absolute node positions an velocities stored in file "<@var{abs_mov_file}>" to relative positions and velocities with respect to node <@var{ref_node}>.
## The result will be stored in file <@var{rel_mov_file}>.
##
## @end deftypefn

function mbdyn_post_abs_to_rel(ref_node, abs_mov_file, rel_mov_file, ref_only, input_mode, output_mode)
  if (nargin < 3 || nargin > 6)
    print_usage();
  endif

  if (nargin < 4)
    ref_only = 0;
  endif

  if (nargin < 5)
    input_mode = "euler123";
  endif

  if (nargin < 6)
    output_mode = "euler123";
  endif

  if (~isscalar(ref_node))
    print_usage();
  endif

  if (~ischar(abs_mov_file))
    print_usage();
  endif

  if (~ischar(rel_mov_file))
    print_usage();
  endif

  mbdyn_path_init();
  
  command = sprintf("exec awk -f abs2rel.awk -v RefNode=%d -v RefOnly=%d -v InputMode=%s -v OutputMode=%s '%s.mov' > '%s.mov'", ...
		    ref_node, ...
		    ref_only, ...
		    input_mode, ...
		    output_mode, ...
		    mbdyn_convert_path(abs_mov_file), ...
		    mbdyn_convert_path(rel_mov_file));

  rc = shell(command);
  
  if (0 ~= rc)
    error("command \"%s\" failed with status %d", command, rc);
  endif
endfunction
  
%!test
%! ## TEST1
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
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
%!     fputs(fd, "         final time: 0.1;\n");
%!     fputs(fd, "         time step: 1e-5;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-3;\n");
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
%!     fputs(fd, "             output: iterations;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");
%!     fputs(fd, "     output meter: closest next, 0, forever, const, 1e-2;\n");
%!     fputs(fd, "     output precision: 16;\n");
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " reference: 1,\n");
%!     fputs(fd, "     position, reference, global, 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "     orientation, reference, global, euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "     velocity, reference, global, 1., 2., 3.,\n");
%!     fputs(fd, "     angular velocity, reference, global, 0.1, 0.2, 0.3;\n");
%!     fputs(fd, " reference: 2,\n");
%!     fputs(fd, "     position, reference, 1, 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "     orientation, reference, 1, euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "     velocity, reference, 1, null,\n");
%!     fputs(fd, "     angular velocity, reference, 1, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, eye,");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 reference, 2, eye,");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., mult, time, const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., mult, time, const, F2;\n");
%!     fputs(fd, "         joint: 1, total joint,\n");
%!     fputs(fd, "                   1,\n");
%!     fputs(fd, "                    position, reference, 1, null,\n");
%!     fputs(fd, "                    position orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    rotation orientation, reference, 1, eye,\n");
%!     fputs(fd, "                   2,\n");
%!     fputs(fd, "                    position, reference, 1, null,\n");
%!     fputs(fd, "                    position orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    rotation orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    position constraint, active, active, active, null,\n");
%!     fputs(fd, "                    orientation constraint, active, active, active, null;\n");
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
%!   mbdyn_post_abs_to_rel(2, options.output_file, [options.output_file, "_rel"], false);
%!   copyfile([options.output_file, ".log"], [options.output_file, "_rel.log"]);
%!   copyfile([options.output_file, ".out"], [options.output_file, "_rel.out"]);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory_abs, deformation_abs, velocity_abs, acceleration_abs, node_id] = mbdyn_post_load_output_struct(options.output_file, 1:2);
%!   [t, trajectory_rel, deformation_rel, velocity_rel, acceleration_rel, node_id] = mbdyn_post_load_output_struct([options.output_file, "_rel"], 1:2);
%!   tol = eps^0.4;
%!   for i=1:numel(node_id)
%!     assert(trajectory_rel{i}, repmat(trajectory_rel{i}(1,:), numel(t), 1), tol * max(max(abs(trajectory_abs{i}))));
%!     assert(velocity_rel{i}, zeros(numel(t), 6), tol * max(max(abs(velocity_abs{i}))));
%!     ## FIXME: bug in abs2rel.awk?
%!     ## assert(acceleration_rel{i}, zeros(numel(t), 6), tol * max(max(abs(acceleration_abs{i}))));
%!   endfor
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
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
%!     fputs(fd, "         final time: 0.1;\n");
%!     fputs(fd, "         time step: 1e-5;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-3;\n");
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
%!     fputs(fd, "             output: iterations;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     output meter: closest next, 0, forever, const, 1e-2;\n");
%!     fputs(fd, "     output precision: 16;\n");
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " reference: 1,\n");
%!     fputs(fd, "     position, reference, global, 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "     orientation, reference, global, euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "     velocity, reference, global, 1., 2., 3.,\n");
%!     fputs(fd, "     angular velocity, reference, global, 0.1, 0.2, 0.3;\n");
%!     fputs(fd, " reference: 2,\n");
%!     fputs(fd, "     position, reference, 1, 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "     orientation, reference, 1, euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "     velocity, reference, 1, null,\n");
%!     fputs(fd, "     angular velocity, reference, 1, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, eye,");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 reference, 2, eye,");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., mult, time, const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., mult, time, const, F2;\n");
%!     fputs(fd, "         joint: 1, total joint,\n");
%!     fputs(fd, "                   1,\n");
%!     fputs(fd, "                    position, reference, 1, null,\n");
%!     fputs(fd, "                    position orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    rotation orientation, reference, 1, eye,\n");
%!     fputs(fd, "                   2,\n");
%!     fputs(fd, "                    position, reference, 1, null,\n");
%!     fputs(fd, "                    position orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    rotation orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    position constraint, active, active, active, null,\n");
%!     fputs(fd, "                    orientation constraint, active, active, active, null;\n");
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
%!   mbdyn_post_abs_to_rel(2, options.output_file, [options.output_file, "_rel"], false);
%!   copyfile([options.output_file, ".log"], [options.output_file, "_rel.log"]);
%!   copyfile([options.output_file, ".out"], [options.output_file, "_rel.out"]);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory_abs, deformation_abs, velocity_abs, acceleration_abs, node_id] = mbdyn_post_load_output_struct(options.output_file, 1:2);
%!   [t, trajectory_rel, deformation_rel, velocity_rel, acceleration_rel, node_id] = mbdyn_post_load_output_struct([options.output_file, "_rel"], 1:2);
%!   tol = eps^0.4;
%!   for i=1:numel(node_id)
%!     assert(trajectory_rel{i}, repmat(trajectory_rel{i}(1,:), numel(t), 1), tol * max(max(abs(trajectory_abs{i}))));
%!     assert(velocity_rel{i}, zeros(numel(t), 6), tol * max(max(abs(velocity_abs{i}))));
%!     ## FIXME: bug in abs2rel.awk?
%!     ## assert(acceleration_rel{i}, zeros(numel(t), 6), tol * max(max(abs(acceleration_abs{i}))));
%!   endfor
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
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
%!     fputs(fd, "         final time: 0.1;\n");
%!     fputs(fd, "         time step: 1e-5;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, "         output: iterations;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     output meter: closest next, 0, forever, const, 1e-2;\n");
%!     fputs(fd, "     output precision: 16;\n");
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " reference: 1,\n");
%!     fputs(fd, "     position, reference, global, 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "     orientation, reference, global, euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "     velocity, reference, global, 1., 2., 3.,\n");
%!     fputs(fd, "     angular velocity, reference, global, 0.1, 0.2, 0.3;\n");
%!     fputs(fd, " reference: 2,\n");
%!     fputs(fd, "     position, reference, 1, 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "     orientation, reference, 1, euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "     velocity, reference, 1, null,\n");
%!     fputs(fd, "     angular velocity, reference, 1, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, eye,");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 reference, 2, eye,");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., mult, time, const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., mult, time, const, F2;\n");
%!     fputs(fd, "         joint: 1, total joint,\n");
%!     fputs(fd, "                   1,\n");
%!     fputs(fd, "                    position, reference, 1, null,\n");
%!     fputs(fd, "                    position orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    rotation orientation, reference, 1, eye,\n");
%!     fputs(fd, "                   2,\n");
%!     fputs(fd, "                    position, reference, 1, null,\n");
%!     fputs(fd, "                    position orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    rotation orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    position constraint, active, active, active, null,\n");
%!     fputs(fd, "                    orientation constraint, active, active, active, null;\n");
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
%!   mbdyn_post_abs_to_rel(2, options.output_file, [options.output_file, "_rel"], false);
%!   copyfile([options.output_file, ".log"], [options.output_file, "_rel.log"]);
%!   copyfile([options.output_file, ".out"], [options.output_file, "_rel.out"]);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory_abs, deformation_abs, velocity_abs, acceleration_abs, node_id] = mbdyn_post_load_output_struct(options.output_file, 1:2);
%!   [t, trajectory_rel, deformation_rel, velocity_rel, acceleration_rel, node_id] = mbdyn_post_load_output_struct([options.output_file, "_rel"], 1:2);
%!   tol = eps^0.4;
%!   for i=1:numel(node_id)
%!     assert(trajectory_rel{i}, repmat(trajectory_rel{i}(1,:), numel(t), 1), tol * max(max(abs(trajectory_abs{i}))));
%!     assert(velocity_rel{i}, zeros(numel(t), 6), tol * max(max(abs(velocity_abs{i}))));
%!     ## FIXME: bug in abs2rel.awk?
%!     ## assert(acceleration_rel{i}, zeros(numel(t), 6), tol * max(max(abs(acceleration_abs{i}))));
%!   endfor
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST4
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!  fputs(fd, "set: real r = 10.;\n");
%!  fputs(fd, "set: real m = 1.;\n");
%!  fputs(fd, "set: real gx = 0;\n");
%!  fputs(fd, "set: real gy = -9.81;\n");
%!  fputs(fd, "set: real gz = 0.;\n");
%!  fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!  fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!  fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!  fputs(fd, "set: real initial_time = 0;\n");
%!  fputs(fd, "set: real final_time = 0.1;\n");
%!  fputs(fd, "set: integer N = 720;\n");
%!  fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!  fputs(fd, "begin: data;\n");
%!  fputs(fd, "	problem: initial value; # the default\n");
%!  fputs(fd, "end: data;\n");
%!  fputs(fd, "begin: initial value;\n");
%!  fputs(fd, "    threads: assembly, 1;\n");
%!  fputs(fd, "	initial time: initial_time;\n");
%!  fputs(fd, "	final time: final_time;\n");
%!  fputs(fd, "	time step: time_step;\n");
%!  fputs(fd, "    max iterations: 100;\n");
%!  fputs(fd, "    method: ms, 0.;\n");
%!  fputs(fd, "    nonlinear solver: nox, modified, 1000,\n");
%!  fputs(fd, "        keep jacobian matrix,\n");
%!  fputs(fd, "        inner iterations before assembly, 6,\n");
%!  fputs(fd, "        jacobian operator, newton krylov,\n");
%!  fputs(fd, "        solver, line search based,\n");
%!  fputs(fd, "        forcing term, type 2,\n");
%!  fputs(fd, "        verbose, yes,\n");
%!  fputs(fd, "        print convergence info, no,\n");
%!  fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!  fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!  fputs(fd, "        linear solver, gmres,\n");
%!  fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!  fputs(fd, "        linear solver max iterations, 10,\n");
%!  fputs(fd, "        krylov subspace size, 10,\n");
%!  fputs(fd, "        minimum step, 1e-12,\n");
%!  fputs(fd, "        use preconditioner as solver, no;\n");
%!  fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!  fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!  fputs(fd, "        derivatives max iterations: 20;\n");
%!  fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!  fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!  fputs(fd, "end: initial value;\n");
%!  fputs(fd, "begin: control data;\n");
%!  fputs(fd, "       tolerance: 1e-10;\n");
%!  fputs(fd, "       max iterations: 100;\n");
%!  fputs(fd, "	structural nodes: 1;\n");
%!  fputs(fd, "	rigid bodies: 1;\n");
%!  fputs(fd, "	joints: 1;\n");
%!  fputs(fd, "    forces: 1;\n");
%!  fputs(fd, "    gravity;\n");
%!  fputs(fd, "end: control data;\n");
%!  fputs(fd, "reference: 1,\n");
%!  fputs(fd, "    position, reference, global, 0.05, 0.02, 0.03,\n");
%!  fputs(fd, "    orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!  fputs(fd, "    velocity, reference, global, null,\n");
%!  fputs(fd, "    angular velocity, reference, global, null;\n");
%!  fputs(fd, "reference: 2,\n");
%!  fputs(fd, "   position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!  fputs(fd, "   orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!  fputs(fd, "   velocity, reference, 1, null,\n");
%!  fputs(fd, "   angular velocity, reference, 1, 0., 0., omega0;\n");
%!  fputs(fd, "begin: nodes;\n");
%!  fputs(fd, "	structural: 1, dynamic,\n");
%!  fputs(fd, "		position, reference, 2, 0.07, 0.08, 0.02,\n");
%!  fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!  fputs(fd, "		velocity, reference, 2, null,\n");
%!  fputs(fd, "		angular velocity, reference, 2, null;\n");
%!  fputs(fd, "end: nodes;\n");
%!  fputs(fd, "begin: elements;\n");
%!  fputs(fd, "	body: 1, 1, \n");
%!  fputs(fd, "		# mass\n");
%!  fputs(fd, "		m,\n");
%!  fputs(fd, "		# center of mass\n");
%!  fputs(fd, "		null,\n");
%!  fputs(fd, "		# inertia matrix\n");
%!  fputs(fd, "		diag, .1, 1., 100.;\n");
%!  fputs(fd, "	joint: 1, total pin joint,\n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		# relative position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "		# relative orientation matrix\n");
%!  fputs(fd, "		position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!  fputs(fd, "		# absolute position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "                position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, eye,\n");
%!  fputs(fd, "                position constraint,\n");
%!  fputs(fd, "                active, active, active, null,\n");
%!  fputs(fd, "                orientation constraint,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                inactive,\n");
%!  fputs(fd, "                null;\n");
%!  fputs(fd, "	couple: 1, absolute, \n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		position, reference, node, null, \n");
%!  fputs(fd, "		component, 100., \n");
%!  fputs(fd, "			   100., \n");
%!  fputs(fd, "			   1000.;\n");
%!  fputs(fd, "    gravity: uniform, component,\n");
%!  fputs(fd, "        gx, gy, gz;\n");
%!  fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST5
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!  fputs(fd, "set: real r = 10.;\n");
%!  fputs(fd, "set: real m = 1.;\n");
%!  fputs(fd, "set: real gx = 0;\n");
%!  fputs(fd, "set: real gy = -9.81;\n");
%!  fputs(fd, "set: real gz = 0.;\n");
%!  fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!  fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!  fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!  fputs(fd, "set: real initial_time = 0;\n");
%!  fputs(fd, "set: real final_time = 0.1;\n");
%!  fputs(fd, "set: integer N = 720;\n");
%!  fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!  fputs(fd, "begin: data;\n");
%!  fputs(fd, "	problem: initial value; # the default\n");
%!  fputs(fd, "end: data;\n");
%!  fputs(fd, "begin: initial value;\n");
%!  fputs(fd, "    threads: assembly, 1;\n");
%!  fputs(fd, "	initial time: initial_time;\n");
%!  fputs(fd, "	final time: final_time;\n");
%!  fputs(fd, "	time step: time_step;\n");
%!  fputs(fd, "    max iterations: 100;\n");
%!  fputs(fd, "    linear solver: naive, colamd;\n");
%!  fputs(fd, "    method: ms, 0.;\n");
%!  fputs(fd, "    nonlinear solver: nox, modified, 1000,\n");
%!  fputs(fd, "        keep jacobian matrix,\n");
%!  fputs(fd, "        inner iterations before assembly, 6,\n");
%!  fputs(fd, "        jacobian operator, newton krylov,\n");
%!  fputs(fd, "        solver, line search based,\n");
%!  fputs(fd, "        forcing term, type 2,\n");
%!  fputs(fd, "        verbose, yes,\n");
%!  fputs(fd, "        print convergence info, no,\n");
%!  fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!  fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!  fputs(fd, "        linear solver, gmres,\n");
%!  fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!  fputs(fd, "        linear solver max iterations, 10,\n");
%!  fputs(fd, "        krylov subspace size, 10,\n");
%!  fputs(fd, "        minimum step, 1e-12,\n");
%!  fputs(fd, "        use preconditioner as solver, no;\n");
%!  fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!  fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!  fputs(fd, "        derivatives max iterations: 20;\n");
%!  fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!  fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!  fputs(fd, "end: initial value;\n");
%!  fputs(fd, "begin: control data;\n");
%!  fputs(fd, "       tolerance: 1e-10;\n");
%!  fputs(fd, "       max iterations: 100;\n");
%!  fputs(fd, "	structural nodes: 1;\n");
%!  fputs(fd, "	rigid bodies: 1;\n");
%!  fputs(fd, "	joints: 1;\n");
%!  fputs(fd, "    forces: 1;\n");
%!  fputs(fd, "    gravity;\n");
%!  fputs(fd, "end: control data;\n");
%!  fputs(fd, "reference: 1,\n");
%!  fputs(fd, "    position, reference, global, 0.05, 0.02, 0.03,\n");
%!  fputs(fd, "    orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!  fputs(fd, "    velocity, reference, global, null,\n");
%!  fputs(fd, "    angular velocity, reference, global, null;\n");
%!  fputs(fd, "reference: 2,\n");
%!  fputs(fd, "   position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!  fputs(fd, "   orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!  fputs(fd, "   velocity, reference, 1, null,\n");
%!  fputs(fd, "   angular velocity, reference, 1, 0., 0., omega0;\n");
%!  fputs(fd, "begin: nodes;\n");
%!  fputs(fd, "	structural: 1, dynamic,\n");
%!  fputs(fd, "		position, reference, 2, 0.07, 0.08, 0.02,\n");
%!  fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!  fputs(fd, "		velocity, reference, 2, null,\n");
%!  fputs(fd, "		angular velocity, reference, 2, null;\n");
%!  fputs(fd, "end: nodes;\n");
%!  fputs(fd, "begin: elements;\n");
%!  fputs(fd, "	body: 1, 1, \n");
%!  fputs(fd, "		# mass\n");
%!  fputs(fd, "		m,\n");
%!  fputs(fd, "		# center of mass\n");
%!  fputs(fd, "		null,\n");
%!  fputs(fd, "		# inertia matrix\n");
%!  fputs(fd, "		diag, .1, 1., 100.;\n");
%!  fputs(fd, "	joint: 1, total pin joint,\n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		# relative position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "		# relative orientation matrix\n");
%!  fputs(fd, "		position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!  fputs(fd, "		# absolute position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "                position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, eye,                \n");
%!  fputs(fd, "                position constraint,\n");
%!  fputs(fd, "                active, active, active, null,\n");
%!  fputs(fd, "                orientation constraint,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                inactive,\n");
%!  fputs(fd, "                null;\n");
%!  fputs(fd, "	couple: 1, absolute, \n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		position, reference, node, null, \n");
%!  fputs(fd, "		component, 100., \n");
%!  fputs(fd, "			   100., \n");
%!  fputs(fd, "			   1000.;\n");
%!  fputs(fd, "    gravity: uniform, component,\n");
%!  fputs(fd, "        gx, gy, gz;\n");
%!  fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST6
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!  fputs(fd, "set: real r = 10.;\n");
%!  fputs(fd, "set: real m = 1.;\n");
%!  fputs(fd, "set: real gx = 0;\n");
%!  fputs(fd, "set: real gy = -9.81;\n");
%!  fputs(fd, "set: real gz = 0.;\n");
%!  fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!  fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!  fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!  fputs(fd, "set: real initial_time = 0;\n");
%!  fputs(fd, "set: real final_time = 0.1;\n");
%!  fputs(fd, "set: integer N = 720;\n");
%!  fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!  fputs(fd, "begin: data;\n");
%!  fputs(fd, "	problem: initial value; # the default\n");
%!  fputs(fd, "end: data;\n");
%!  fputs(fd, "begin: initial value;\n");
%!  fputs(fd, "    threads: assembly, 1;\n");
%!  fputs(fd, "	initial time: initial_time;\n");
%!  fputs(fd, "	final time: final_time;\n");
%!  fputs(fd, "	time step: time_step;\n");
%!  fputs(fd, "    max iterations: 100;\n");
%!  fputs(fd, "    method: ms, 0.;\n");
%!  fputs(fd, "    nonlinear solver: nox, modified, 1000,\n");
%!  fputs(fd, "        keep jacobian matrix,\n");
%!  fputs(fd, "        inner iterations before assembly, 6,\n");
%!  fputs(fd, "        jacobian operator, newton,\n");
%!  fputs(fd, "        solver, line search based,\n");
%!  fputs(fd, "        forcing term, type 2,\n");
%!  fputs(fd, "        verbose, yes,\n");
%!  fputs(fd, "        print convergence info, no,\n");
%!  fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!  fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!  fputs(fd, "        linear solver, gmres,\n");
%!  fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!  fputs(fd, "        linear solver max iterations, 10,\n");
%!  fputs(fd, "        krylov subspace size, 10,\n");
%!  fputs(fd, "        minimum step, 1e-12,\n");
%!  fputs(fd, "        use preconditioner as solver, no;\n");
%!  fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!  fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!  fputs(fd, "        derivatives max iterations: 20;\n");
%!  fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!  fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!  fputs(fd, "end: initial value;\n");
%!  fputs(fd, "begin: control data;\n");
%!  fputs(fd, "       tolerance: 1e-10;\n");
%!  fputs(fd, "       max iterations: 100;\n");
%!  fputs(fd, "	structural nodes: 1;\n");
%!  fputs(fd, "	rigid bodies: 1;\n");
%!  fputs(fd, "	joints: 1;\n");
%!  fputs(fd, "    forces: 1;\n");
%!  fputs(fd, "    gravity;\n");
%!  fputs(fd, "end: control data;\n");
%!  fputs(fd, "reference: 1,\n");
%!  fputs(fd, "    position, reference, global, 0.05, 0.02, 0.03,\n");
%!  fputs(fd, "    orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!  fputs(fd, "    velocity, reference, global, null,\n");
%!  fputs(fd, "    angular velocity, reference, global, null;\n");
%!  fputs(fd, "reference: 2,\n");
%!  fputs(fd, "   position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!  fputs(fd, "   orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!  fputs(fd, "   velocity, reference, 1, null,\n");
%!  fputs(fd, "   angular velocity, reference, 1, 0., 0., omega0;\n");
%!  fputs(fd, "begin: nodes;\n");
%!  fputs(fd, "	structural: 1, dynamic,\n");
%!  fputs(fd, "		position, reference, 2, 0.07, 0.08, 0.02,\n");
%!  fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!  fputs(fd, "		velocity, reference, 2, null,\n");
%!  fputs(fd, "		angular velocity, reference, 2, null;\n");
%!  fputs(fd, "end: nodes;\n");
%!  fputs(fd, "begin: elements;\n");
%!  fputs(fd, "	body: 1, 1, \n");
%!  fputs(fd, "		# mass\n");
%!  fputs(fd, "		m,\n");
%!  fputs(fd, "		# center of mass\n");
%!  fputs(fd, "		null,\n");
%!  fputs(fd, "		# inertia matrix\n");
%!  fputs(fd, "		diag, .1, 1., 100.;\n");
%!  fputs(fd, "	joint: 1, total pin joint,\n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		# relative position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "		# relative orientation matrix\n");
%!  fputs(fd, "		position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!  fputs(fd, "		# absolute position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "                position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, eye,                \n");
%!  fputs(fd, "                position constraint,\n");
%!  fputs(fd, "                active, active, active, null,\n");
%!  fputs(fd, "                orientation constraint,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                inactive,\n");
%!  fputs(fd, "                null;\n");
%!  fputs(fd, "	couple: 1, absolute, \n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		position, reference, node, null, \n");
%!  fputs(fd, "		component, 100., \n");
%!  fputs(fd, "			   100., \n");
%!  fputs(fd, "			   1000.;\n");
%!  fputs(fd, "    gravity: uniform, component,\n");
%!  fputs(fd, "        gx, gy, gz;\n");
%!  fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST7
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!  fputs(fd, "set: real r = 10.;\n");
%!  fputs(fd, "set: real m = 1.;\n");
%!  fputs(fd, "set: real gx = 0;\n");
%!  fputs(fd, "set: real gy = -9.81;\n");
%!  fputs(fd, "set: real gz = 0.;\n");
%!  fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!  fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!  fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!  fputs(fd, "set: real initial_time = 0;\n");
%!  fputs(fd, "set: real final_time = 0.1;\n");
%!  fputs(fd, "set: integer N = 720;\n");
%!  fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!  fputs(fd, "begin: data;\n");
%!  fputs(fd, "	problem: initial value; # the default\n");
%!  fputs(fd, "end: data;\n");
%!  fputs(fd, "begin: initial value;\n");
%!  fputs(fd, "    threads: assembly, 1;\n");
%!  fputs(fd, "	initial time: initial_time;\n");
%!  fputs(fd, "	final time: final_time;\n");
%!  fputs(fd, "	time step: time_step;\n");
%!  fputs(fd, "    max iterations: 100;\n");
%!  fputs(fd, "    method: ms, 0.;\n");
%!  fputs(fd, "    nonlinear solver: nox, modified, 1000,\n");
%!  fputs(fd, "        keep jacobian matrix,\n");
%!  fputs(fd, "        inner iterations before assembly, 6,\n");
%!  fputs(fd, "        jacobian operator, newton,\n");
%!  fputs(fd, "        solver, line search based,\n");
%!  fputs(fd, "        line search method, more thuente,\n");
%!  fputs(fd, "        sufficient decrease condition, armijo goldstein,\n");
%!  fputs(fd, "        forcing term, type 2,\n");
%!  fputs(fd, "        direction, newton,\n");
%!  fputs(fd, "        verbose, yes,\n");
%!  fputs(fd, "        print convergence info, no,\n");
%!  fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!  fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!  fputs(fd, "        linear solver, gmres,\n");
%!  fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!  fputs(fd, "        linear solver max iterations, 10,\n");
%!  fputs(fd, "        krylov subspace size, 10,\n");
%!  fputs(fd, "        minimum step, 1e-12,\n");
%!  fputs(fd, "        use preconditioner as solver, no;\n");
%!  fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!  fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!  fputs(fd, "        derivatives max iterations: 20;\n");
%!  fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!  fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!  fputs(fd, "end: initial value;\n");
%!  fputs(fd, "begin: control data;\n");
%!  fputs(fd, "       tolerance: 1e-10;\n");
%!  fputs(fd, "       max iterations: 100;\n");
%!  fputs(fd, "	structural nodes: 1;\n");
%!  fputs(fd, "	rigid bodies: 1;\n");
%!  fputs(fd, "	joints: 1;\n");
%!  fputs(fd, "    forces: 1;\n");
%!  fputs(fd, "    gravity;\n");
%!  fputs(fd, "end: control data;\n");
%!  fputs(fd, "reference: 1,\n");
%!  fputs(fd, "    position, reference, global, 0.05, 0.02, 0.03,\n");
%!  fputs(fd, "    orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!  fputs(fd, "    velocity, reference, global, null,\n");
%!  fputs(fd, "    angular velocity, reference, global, null;\n");
%!  fputs(fd, "reference: 2,\n");
%!  fputs(fd, "   position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!  fputs(fd, "   orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!  fputs(fd, "   velocity, reference, 1, null,\n");
%!  fputs(fd, "   angular velocity, reference, 1, 0., 0., omega0;\n");
%!  fputs(fd, "begin: nodes;\n");
%!  fputs(fd, "	structural: 1, dynamic,\n");
%!  fputs(fd, "		position, reference, 2, 0.07, 0.08, 0.02,\n");
%!  fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!  fputs(fd, "		velocity, reference, 2, null,\n");
%!  fputs(fd, "		angular velocity, reference, 2, null;\n");
%!  fputs(fd, "end: nodes;\n");
%!  fputs(fd, "begin: elements;\n");
%!  fputs(fd, "	body: 1, 1, \n");
%!  fputs(fd, "		# mass\n");
%!  fputs(fd, "		m,\n");
%!  fputs(fd, "		# center of mass\n");
%!  fputs(fd, "		null,\n");
%!  fputs(fd, "		# inertia matrix\n");
%!  fputs(fd, "		diag, .1, 1., 100.;\n");
%!  fputs(fd, "	joint: 1, total pin joint,\n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		# relative position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "		# relative orientation matrix\n");
%!  fputs(fd, "		position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!  fputs(fd, "		# absolute position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "                position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, eye,                \n");
%!  fputs(fd, "                position constraint,\n");
%!  fputs(fd, "                active, active, active, null,\n");
%!  fputs(fd, "                orientation constraint,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                inactive,\n");
%!  fputs(fd, "                null;\n");
%!  fputs(fd, "	couple: 1, absolute, \n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		position, reference, node, null, \n");
%!  fputs(fd, "		component, 100., \n");
%!  fputs(fd, "			   100., \n");
%!  fputs(fd, "			   1000.;\n");
%!  fputs(fd, "    gravity: uniform, component,\n");
%!  fputs(fd, "        gx, gy, gz;\n");
%!  fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST8
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!  fputs(fd, "set: real r = 10.;\n");
%!  fputs(fd, "set: real m = 1.;\n");
%!  fputs(fd, "set: real gx = 0;\n");
%!  fputs(fd, "set: real gy = -9.81;\n");
%!  fputs(fd, "set: real gz = 0.;\n");
%!  fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!  fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!  fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!  fputs(fd, "set: real initial_time = 0;\n");
%!  fputs(fd, "set: real final_time = 0.1;\n");
%!  fputs(fd, "set: integer N = 720;\n");
%!  fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!  fputs(fd, "begin: data;\n");
%!  fputs(fd, "	problem: initial value; # the default\n");
%!  fputs(fd, "end: data;\n");
%!  fputs(fd, "begin: initial value;\n");
%!  fputs(fd, "    threads: assembly, 1;\n");
%!  fputs(fd, "	initial time: initial_time;\n");
%!  fputs(fd, "	final time: final_time;\n");
%!  fputs(fd, "	time step: time_step;\n");
%!  fputs(fd, "    max iterations: 100;\n");
%!  fputs(fd, "    method: ms, 0.;\n");
%!  fputs(fd, "    nonlinear solver: nox, modified, 1000,\n");
%!  fputs(fd, "        keep jacobian matrix,\n");
%!  fputs(fd, "        inner iterations before assembly, 6,\n");
%!  fputs(fd, "        jacobian operator, newton,\n");
%!  fputs(fd, "        solver, line search based,\n");
%!  fputs(fd, "        forcing term, type 2,\n");
%!  fputs(fd, "        direction, broyden,\n");
%!  fputs(fd, "        verbose, yes,\n");
%!  fputs(fd, "        print convergence info, no,\n");
%!  fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!  fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!  fputs(fd, "        linear solver, gmres,\n");
%!  fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!  fputs(fd, "        linear solver max iterations, 10,\n");
%!  fputs(fd, "        krylov subspace size, 10,\n");
%!  fputs(fd, "        minimum step, 1e-12,\n");
%!  fputs(fd, "        use preconditioner as solver, no;\n");
%!  fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!  fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!  fputs(fd, "        derivatives max iterations: 20;\n");
%!  fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!  fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!  fputs(fd, "end: initial value;\n");
%!  fputs(fd, "begin: control data;\n");
%!  fputs(fd, "       tolerance: 1e-10;\n");
%!  fputs(fd, "       max iterations: 100;\n");
%!  fputs(fd, "	structural nodes: 1;\n");
%!  fputs(fd, "	rigid bodies: 1;\n");
%!  fputs(fd, "	joints: 1;\n");
%!  fputs(fd, "    forces: 1;\n");
%!  fputs(fd, "    gravity;\n");
%!  fputs(fd, "end: control data;\n");
%!  fputs(fd, "reference: 1,\n");
%!  fputs(fd, "    position, reference, global, 0.05, 0.02, 0.03,\n");
%!  fputs(fd, "    orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!  fputs(fd, "    velocity, reference, global, null,\n");
%!  fputs(fd, "    angular velocity, reference, global, null;\n");
%!  fputs(fd, "reference: 2,\n");
%!  fputs(fd, "   position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!  fputs(fd, "   orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!  fputs(fd, "   velocity, reference, 1, null,\n");
%!  fputs(fd, "   angular velocity, reference, 1, 0., 0., omega0;\n");
%!  fputs(fd, "begin: nodes;\n");
%!  fputs(fd, "	structural: 1, dynamic,\n");
%!  fputs(fd, "		position, reference, 2, 0.07, 0.08, 0.02,\n");
%!  fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!  fputs(fd, "		velocity, reference, 2, null,\n");
%!  fputs(fd, "		angular velocity, reference, 2, null;\n");
%!  fputs(fd, "end: nodes;\n");
%!  fputs(fd, "begin: elements;\n");
%!  fputs(fd, "	body: 1, 1, \n");
%!  fputs(fd, "		# mass\n");
%!  fputs(fd, "		m,\n");
%!  fputs(fd, "		# center of mass\n");
%!  fputs(fd, "		null,\n");
%!  fputs(fd, "		# inertia matrix\n");
%!  fputs(fd, "		diag, .1, 1., 100.;\n");
%!  fputs(fd, "	joint: 1, total pin joint,\n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		# relative position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "		# relative orientation matrix\n");
%!  fputs(fd, "		position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!  fputs(fd, "		# absolute position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "                position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, eye,                \n");
%!  fputs(fd, "                position constraint,\n");
%!  fputs(fd, "                active, active, active, null,\n");
%!  fputs(fd, "                orientation constraint,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                inactive,\n");
%!  fputs(fd, "                null;\n");
%!  fputs(fd, "	couple: 1, absolute, \n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		position, reference, node, null, \n");
%!  fputs(fd, "		component, 100., \n");
%!  fputs(fd, "			   100., \n");
%!  fputs(fd, "			   1000.;\n");
%!  fputs(fd, "    gravity: uniform, component,\n");
%!  fputs(fd, "        gx, gy, gz;\n");
%!  fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST9
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!  fputs(fd, "set: real r = 10.;\n");
%!  fputs(fd, "set: real m = 1.;\n");
%!  fputs(fd, "set: real gx = 0;\n");
%!  fputs(fd, "set: real gy = -9.81;\n");
%!  fputs(fd, "set: real gz = 0.;\n");
%!  fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!  fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!  fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!  fputs(fd, "set: real initial_time = 0;\n");
%!  fputs(fd, "set: real final_time = 0.1;\n");
%!  fputs(fd, "set: integer N = 720;\n");
%!  fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!  fputs(fd, "begin: data;\n");
%!  fputs(fd, "	problem: initial value; # the default\n");
%!  fputs(fd, "end: data;\n");
%!  fputs(fd, "begin: initial value;\n");
%!  fputs(fd, "    threads: assembly, 1;\n");
%!  fputs(fd, "	initial time: initial_time;\n");
%!  fputs(fd, "	final time: final_time;\n");
%!  fputs(fd, "	time step: time_step;\n");
%!  fputs(fd, "    max iterations: 100;\n");
%!  fputs(fd, "    method: ms, 0.;\n");
%!  fputs(fd, "    nonlinear solver: nox, modified, 1000,\n");
%!  fputs(fd, "        keep jacobian matrix,\n");
%!  fputs(fd, "        inner iterations before assembly, 6,\n");
%!  fputs(fd, "        jacobian operator, newton,\n");
%!  fputs(fd, "        solver, inexact trust region based,\n");
%!  fputs(fd, "        forcing term, type 2,\n");
%!  fputs(fd, "        direction, newton,\n");
%!  fputs(fd, "        verbose, yes,\n");
%!  fputs(fd, "        print convergence info, no,\n");
%!  fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!  fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!  fputs(fd, "        linear solver, gmres,\n");
%!  fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!  fputs(fd, "        linear solver max iterations, 10,\n");
%!  fputs(fd, "        krylov subspace size, 10,\n");
%!  fputs(fd, "        minimum step, 1e-12,\n");
%!  fputs(fd, "        use preconditioner as solver, no;\n");
%!  fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!  fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!  fputs(fd, "        derivatives max iterations: 20;\n");
%!  fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!  fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!  fputs(fd, "end: initial value;\n");
%!  fputs(fd, "begin: control data;\n");
%!  fputs(fd, "       tolerance: 1e-10;\n");
%!  fputs(fd, "       max iterations: 100;\n");
%!  fputs(fd, "	structural nodes: 1;\n");
%!  fputs(fd, "	rigid bodies: 1;\n");
%!  fputs(fd, "	joints: 1;\n");
%!  fputs(fd, "    forces: 1;\n");
%!  fputs(fd, "    gravity;\n");
%!  fputs(fd, "end: control data;\n");
%!  fputs(fd, "reference: 1,\n");
%!  fputs(fd, "    position, reference, global, 0.05, 0.02, 0.03,\n");
%!  fputs(fd, "    orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!  fputs(fd, "    velocity, reference, global, null,\n");
%!  fputs(fd, "    angular velocity, reference, global, null;\n");
%!  fputs(fd, "reference: 2,\n");
%!  fputs(fd, "   position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!  fputs(fd, "   orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!  fputs(fd, "   velocity, reference, 1, null,\n");
%!  fputs(fd, "   angular velocity, reference, 1, 0., 0., omega0;\n");
%!  fputs(fd, "begin: nodes;\n");
%!  fputs(fd, "	structural: 1, dynamic,\n");
%!  fputs(fd, "		position, reference, 2, 0.07, 0.08, 0.02,\n");
%!  fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!  fputs(fd, "		velocity, reference, 2, null,\n");
%!  fputs(fd, "		angular velocity, reference, 2, null;\n");
%!  fputs(fd, "end: nodes;\n");
%!  fputs(fd, "begin: elements;\n");
%!  fputs(fd, "	body: 1, 1, \n");
%!  fputs(fd, "		# mass\n");
%!  fputs(fd, "		m,\n");
%!  fputs(fd, "		# center of mass\n");
%!  fputs(fd, "		null,\n");
%!  fputs(fd, "		# inertia matrix\n");
%!  fputs(fd, "		diag, .1, 1., 100.;\n");
%!  fputs(fd, "	joint: 1, total pin joint,\n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		# relative position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "		# relative orientation matrix\n");
%!  fputs(fd, "		position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!  fputs(fd, "		# absolute position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "                position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, eye,                \n");
%!  fputs(fd, "                position constraint,\n");
%!  fputs(fd, "                active, active, active, null,\n");
%!  fputs(fd, "                orientation constraint,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                inactive,\n");
%!  fputs(fd, "                null;\n");
%!  fputs(fd, "	couple: 1, absolute, \n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		position, reference, node, null, \n");
%!  fputs(fd, "		component, 100., \n");
%!  fputs(fd, "			   100., \n");
%!  fputs(fd, "			   1000.;\n");
%!  fputs(fd, "    gravity: uniform, component,\n");
%!  fputs(fd, "        gx, gy, gz;\n");
%!  fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST10
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!  fputs(fd, "set: real r = 10.;\n");
%!  fputs(fd, "set: real m = 1.;\n");
%!  fputs(fd, "set: real gx = 0;\n");
%!  fputs(fd, "set: real gy = -9.81;\n");
%!  fputs(fd, "set: real gz = 0.;\n");
%!  fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!  fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!  fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!  fputs(fd, "set: real initial_time = 0;\n");
%!  fputs(fd, "set: real final_time = 0.1;\n");
%!  fputs(fd, "set: integer N = 720;\n");
%!  fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!  fputs(fd, "begin: data;\n");
%!  fputs(fd, "	problem: initial value; # the default\n");
%!  fputs(fd, "end: data;\n");
%!  fputs(fd, "begin: initial value;\n");
%!  fputs(fd, "    threads: assembly, 1;\n");
%!  fputs(fd, "	initial time: initial_time;\n");
%!  fputs(fd, "	final time: final_time;\n");
%!  fputs(fd, "	time step: time_step;\n");
%!  fputs(fd, "    max iterations: 100;\n");
%!  fputs(fd, "    method: ms, 0.;\n");
%!  fputs(fd, "    nonlinear solver: nox, modified, 1000,\n");
%!  fputs(fd, "        keep jacobian matrix,\n");
%!  fputs(fd, "        inner iterations before assembly, 6,\n");
%!  fputs(fd, "        jacobian operator, newton,\n");
%!  fputs(fd, "        solver, tensor based,\n");
%!  fputs(fd, "        forcing term, type 2,\n");
%!  fputs(fd, "        direction, newton,\n");
%!  fputs(fd, "        verbose, yes,\n");
%!  fputs(fd, "        print convergence info, no,\n");
%!  fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!  fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!  fputs(fd, "        linear solver, gmres,\n");
%!  fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!  fputs(fd, "        linear solver max iterations, 10,\n");
%!  fputs(fd, "        krylov subspace size, 10,\n");
%!  fputs(fd, "        minimum step, 1e-12,\n");
%!  fputs(fd, "        use preconditioner as solver, no;\n");
%!  fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!  fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!  fputs(fd, "        derivatives max iterations: 20;\n");
%!  fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!  fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!  fputs(fd, "end: initial value;\n");
%!  fputs(fd, "begin: control data;\n");
%!  fputs(fd, "       tolerance: 1e-10;\n");
%!  fputs(fd, "       max iterations: 100;\n");
%!  fputs(fd, "	structural nodes: 1;\n");
%!  fputs(fd, "	rigid bodies: 1;\n");
%!  fputs(fd, "	joints: 1;\n");
%!  fputs(fd, "    forces: 1;\n");
%!  fputs(fd, "    gravity;\n");
%!  fputs(fd, "end: control data;\n");
%!  fputs(fd, "reference: 1,\n");
%!  fputs(fd, "    position, reference, global, 0.05, 0.02, 0.03,\n");
%!  fputs(fd, "    orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!  fputs(fd, "    velocity, reference, global, null,\n");
%!  fputs(fd, "    angular velocity, reference, global, null;\n");
%!  fputs(fd, "reference: 2,\n");
%!  fputs(fd, "   position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!  fputs(fd, "   orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!  fputs(fd, "   velocity, reference, 1, null,\n");
%!  fputs(fd, "   angular velocity, reference, 1, 0., 0., omega0;\n");
%!  fputs(fd, "begin: nodes;\n");
%!  fputs(fd, "	structural: 1, dynamic,\n");
%!  fputs(fd, "		position, reference, 2, 0.07, 0.08, 0.02,\n");
%!  fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!  fputs(fd, "		velocity, reference, 2, null,\n");
%!  fputs(fd, "		angular velocity, reference, 2, null;\n");
%!  fputs(fd, "end: nodes;\n");
%!  fputs(fd, "begin: elements;\n");
%!  fputs(fd, "	body: 1, 1, \n");
%!  fputs(fd, "		# mass\n");
%!  fputs(fd, "		m,\n");
%!  fputs(fd, "		# center of mass\n");
%!  fputs(fd, "		null,\n");
%!  fputs(fd, "		# inertia matrix\n");
%!  fputs(fd, "		diag, .1, 1., 100.;\n");
%!  fputs(fd, "	joint: 1, total pin joint,\n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		# relative position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "		# relative orientation matrix\n");
%!  fputs(fd, "		position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!  fputs(fd, "		# absolute position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "                position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, eye,                \n");
%!  fputs(fd, "                position constraint,\n");
%!  fputs(fd, "                active, active, active, null,\n");
%!  fputs(fd, "                orientation constraint,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                inactive,\n");
%!  fputs(fd, "                null;\n");
%!  fputs(fd, "	couple: 1, absolute, \n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		position, reference, node, null, \n");
%!  fputs(fd, "		component, 100., \n");
%!  fputs(fd, "			   100., \n");
%!  fputs(fd, "			   1000.;\n");
%!  fputs(fd, "    gravity: uniform, component,\n");
%!  fputs(fd, "        gx, gy, gz;\n");
%!  fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST11
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!  fputs(fd, "set: real r = 10.;\n");
%!  fputs(fd, "set: real m = 1.;\n");
%!  fputs(fd, "set: real gx = 0;\n");
%!  fputs(fd, "set: real gy = -9.81;\n");
%!  fputs(fd, "set: real gz = 0.;\n");
%!  fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!  fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!  fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!  fputs(fd, "set: real initial_time = 0;\n");
%!  fputs(fd, "set: real final_time = 0.1;\n");
%!  fputs(fd, "set: integer N = 720;\n");
%!  fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!  fputs(fd, "begin: data;\n");
%!  fputs(fd, "	problem: initial value; # the default\n");
%!  fputs(fd, "end: data;\n");
%!  fputs(fd, "begin: initial value;\n");
%!  fputs(fd, "    threads: assembly, 1;\n");
%!  fputs(fd, "	initial time: initial_time;\n");
%!  fputs(fd, "	final time: final_time;\n");
%!  fputs(fd, "	time step: time_step;\n");
%!  fputs(fd, "    max iterations: 100;\n");
%!  fputs(fd, "    method: ms, 0.;\n");
%!  fputs(fd, "    nonlinear solver: nox, modified, 1000,\n");
%!  fputs(fd, "        keep jacobian matrix,\n");
%!  fputs(fd, "        inner iterations before assembly, 6,\n");
%!  fputs(fd, "        jacobian operator, newton,\n");
%!  fputs(fd, "        solver, trust region based,\n");
%!  fputs(fd, "        forcing term, type 2,\n");
%!  fputs(fd, "        direction, newton,\n");
%!  fputs(fd, "        verbose, yes,\n");
%!  fputs(fd, "        print convergence info, no,\n");
%!  fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!  fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!  fputs(fd, "        linear solver, gmres,\n");
%!  fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!  fputs(fd, "        linear solver max iterations, 10,\n");
%!  fputs(fd, "        krylov subspace size, 10,\n");
%!  fputs(fd, "        minimum step, 1e-12,\n");
%!  fputs(fd, "        use preconditioner as solver, no;\n");
%!  fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!  fputs(fd, "        derivatives tolerance: 1e-4, 1e-4;\n");
%!  fputs(fd, "        derivatives max iterations: 3000;\n");
%!  fputs(fd, "        derivatives coefficient: 1e-4;\n");
%!  fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!  fputs(fd, "end: initial value;\n");
%!  fputs(fd, "begin: control data;\n");
%!  fputs(fd, "       tolerance: 1e-10;\n");
%!  fputs(fd, "       max iterations: 100;\n");
%!  fputs(fd, "	structural nodes: 1;\n");
%!  fputs(fd, "	rigid bodies: 1;\n");
%!  fputs(fd, "	joints: 1;\n");
%!  fputs(fd, "    forces: 1;\n");
%!  fputs(fd, "    gravity;\n");
%!  fputs(fd, "end: control data;\n");
%!  fputs(fd, "reference: 1,\n");
%!  fputs(fd, "    position, reference, global, 0.05, 0.02, 0.03,\n");
%!  fputs(fd, "    orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!  fputs(fd, "    velocity, reference, global, null,\n");
%!  fputs(fd, "    angular velocity, reference, global, null;\n");
%!  fputs(fd, "reference: 2,\n");
%!  fputs(fd, "   position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!  fputs(fd, "   orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!  fputs(fd, "   velocity, reference, 1, null,\n");
%!  fputs(fd, "   angular velocity, reference, 1, 0., 0., omega0;\n");
%!  fputs(fd, "begin: nodes;\n");
%!  fputs(fd, "	structural: 1, dynamic,\n");
%!  fputs(fd, "		position, reference, 2, 0.07, 0.08, 0.02,\n");
%!  fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!  fputs(fd, "		velocity, reference, 2, null,\n");
%!  fputs(fd, "		angular velocity, reference, 2, null;\n");
%!  fputs(fd, "end: nodes;\n");
%!  fputs(fd, "begin: elements;\n");
%!  fputs(fd, "	body: 1, 1, \n");
%!  fputs(fd, "		# mass\n");
%!  fputs(fd, "		m,\n");
%!  fputs(fd, "		# center of mass\n");
%!  fputs(fd, "		null,\n");
%!  fputs(fd, "		# inertia matrix\n");
%!  fputs(fd, "		diag, .1, 1., 100.;\n");
%!  fputs(fd, "	joint: 1, total pin joint,\n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		# relative position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "		# relative orientation matrix\n");
%!  fputs(fd, "		position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!  fputs(fd, "		# absolute position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "                position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, eye,                \n");
%!  fputs(fd, "                position constraint,\n");
%!  fputs(fd, "                active, active, active, null,\n");
%!  fputs(fd, "                orientation constraint,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                inactive,\n");
%!  fputs(fd, "                null;\n");
%!  fputs(fd, "	couple: 1, absolute, \n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		position, reference, node, null, \n");
%!  fputs(fd, "		component, 100., \n");
%!  fputs(fd, "			   100., \n");
%!  fputs(fd, "			   1000.;\n");
%!  fputs(fd, "    gravity: uniform, component,\n");
%!  fputs(fd, "        gx, gy, gz;\n");
%!  fputs(fd, "end: elements;\n");
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
%!   # options.mbdyn_command = "valgrind mbdyn";
%!   mbdyn_solver_run(fname, options);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST12
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!  fputs(fd, "set: real r = 10.;\n");
%!  fputs(fd, "set: real m = 1.;\n");
%!  fputs(fd, "set: real gx = 0;\n");
%!  fputs(fd, "set: real gy = -9.81;\n");
%!  fputs(fd, "set: real gz = 0.;\n");
%!  fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!  fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!  fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!  fputs(fd, "set: real initial_time = 0;\n");
%!  fputs(fd, "set: real final_time = 0.1;\n");
%!  fputs(fd, "set: integer N = 720;\n");
%!  fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!  fputs(fd, "begin: data;\n");
%!  fputs(fd, "	problem: initial value; # the default\n");
%!  fputs(fd, "end: data;\n");
%!  fputs(fd, "begin: initial value;\n");
%!  fputs(fd, "    threads: assembly, 1;\n");
%!  fputs(fd, "	initial time: initial_time;\n");
%!  fputs(fd, "	final time: final_time;\n");
%!  fputs(fd, "	time step: time_step;\n");
%!  fputs(fd, "    max iterations: 100;\n");
%!  fputs(fd, "    linear solver: naive, colamd;\n");
%!  fputs(fd, "    method: ms, 0.;\n");
%!  fputs(fd, "    nonlinear solver: nox, modified, 1000,\n");
%!  fputs(fd, "        keep jacobian matrix,\n");
%!  fputs(fd, "        inner iterations before assembly, 6,\n");
%!  fputs(fd, "        jacobian operator, newton krylov,\n");
%!  fputs(fd, "        solver, line search based,\n");
%!  fputs(fd, "        forcing term, type 2,\n");
%!  fputs(fd, "        verbose, yes,\n");
%!  fputs(fd, "        print convergence info, no,\n");
%!  fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!  fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!  fputs(fd, "        linear solver, gmres,\n");
%!  fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!  fputs(fd, "        linear solver max iterations, 10,\n");
%!  fputs(fd, "        krylov subspace size, 10,\n");
%!  fputs(fd, "        minimum step, 1e-12,\n");
%!  fputs(fd, "        use preconditioner as solver, yes;\n");
%!  fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!  fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!  fputs(fd, "        derivatives max iterations: 20;\n");
%!  fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!  fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!  fputs(fd, "end: initial value;\n");
%!  fputs(fd, "begin: control data;\n");
%!  fputs(fd, "       tolerance: 1e-10;\n");
%!  fputs(fd, "       max iterations: 100;\n");
%!  fputs(fd, "	structural nodes: 1;\n");
%!  fputs(fd, "	rigid bodies: 1;\n");
%!  fputs(fd, "	joints: 1;\n");
%!  fputs(fd, "    forces: 1;\n");
%!  fputs(fd, "    gravity;\n");
%!  fputs(fd, "end: control data;\n");
%!  fputs(fd, "reference: 1,\n");
%!  fputs(fd, "    position, reference, global, 0.05, 0.02, 0.03,\n");
%!  fputs(fd, "    orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!  fputs(fd, "    velocity, reference, global, null,\n");
%!  fputs(fd, "    angular velocity, reference, global, null;\n");
%!  fputs(fd, "reference: 2,\n");
%!  fputs(fd, "   position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!  fputs(fd, "   orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!  fputs(fd, "   velocity, reference, 1, null,\n");
%!  fputs(fd, "   angular velocity, reference, 1, 0., 0., omega0;\n");
%!  fputs(fd, "begin: nodes;\n");
%!  fputs(fd, "	structural: 1, dynamic,\n");
%!  fputs(fd, "		position, reference, 2, 0.07, 0.08, 0.02,\n");
%!  fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!  fputs(fd, "		velocity, reference, 2, null,\n");
%!  fputs(fd, "		angular velocity, reference, 2, null;\n");
%!  fputs(fd, "end: nodes;\n");
%!  fputs(fd, "begin: elements;\n");
%!  fputs(fd, "	body: 1, 1, \n");
%!  fputs(fd, "		# mass\n");
%!  fputs(fd, "		m,\n");
%!  fputs(fd, "		# center of mass\n");
%!  fputs(fd, "		null,\n");
%!  fputs(fd, "		# inertia matrix\n");
%!  fputs(fd, "		diag, .1, 1., 100.;\n");
%!  fputs(fd, "	joint: 1, total pin joint,\n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		# relative position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "		# relative orientation matrix\n");
%!  fputs(fd, "		position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!  fputs(fd, "		# absolute position\n");
%!  fputs(fd, "		position, reference, 2, null,\n");
%!  fputs(fd, "                position orientation, reference, 2, eye,\n");
%!  fputs(fd, "                rotation orientation, reference, 2, eye,\n");
%!  fputs(fd, "                position constraint,\n");
%!  fputs(fd, "                active, active, active, null,\n");
%!  fputs(fd, "                orientation constraint,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                active,\n");
%!  fputs(fd, "                inactive,\n");
%!  fputs(fd, "                null;\n");
%!  fputs(fd, "	couple: 1, absolute, \n");
%!  fputs(fd, "		1, \n");
%!  fputs(fd, "		position, reference, node, null, \n");
%!  fputs(fd, "		component, 100., \n");
%!  fputs(fd, "			   100., \n");
%!  fputs(fd, "			   1000.;\n");
%!  fputs(fd, "    gravity: uniform, component,\n");
%!  fputs(fd, "        gx, gy, gz;\n");
%!  fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST13
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, "set: real r = 10.;\n");
%!     fputs(fd, "set: real m = 1.;\n");
%!     fputs(fd, "set: real gx = 0;\n");
%!     fputs(fd, "set: real gy = -9.81;\n");
%!     fputs(fd, "set: real gz = 0.;\n");
%!     fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!     fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!     fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!     fputs(fd, "set: real initial_time = 0;\n");
%!     fputs(fd, "set: real final_time = 0.1;\n");
%!     fputs(fd, "set: integer N = 720;\n");
%!     fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "    threads: assembly, 1;\n");
%!     fputs(fd, "	   initial time: initial_time;\n");
%!     fputs(fd, "	   final time: final_time;\n");
%!     fputs(fd, "	   time step: time_step;\n");
%!     fputs(fd, "    max iterations: 100;\n");
%!     fputs(fd, "    linear solver: umfpack, grad;\n");
%!     fputs(fd, "    method: ms, 0.;\n");
%!     fputs(fd, "    nonlinear solver: nox,\n");
%!     fputs(fd, "        jacobian operator, newton,\n");
%!     fputs(fd, "        solver, line search based,\n");
%!     fputs(fd, "        forcing term, type 2,\n");
%!     fputs(fd, "        verbose, yes,\n");
%!     fputs(fd, "        print convergence info, no,\n");
%!     fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!     fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!     fputs(fd, "        linear solver, gmres,\n");
%!     fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!     fputs(fd, "        linear solver max iterations, 10,\n");
%!     fputs(fd, "        krylov subspace size, 10,\n");
%!     fputs(fd, "        minimum step, 1e-12,\n");
%!     fputs(fd, "        use preconditioner as solver, yes;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       tolerance: 1e-10;\n");
%!     fputs(fd, "       max iterations: 100;\n");
%!     fputs(fd, "	      structural nodes: 1;\n");
%!     fputs(fd, "	      rigid bodies: 1;\n");
%!     fputs(fd, "	      joints: 1;\n");
%!     fputs(fd, "       forces: 1;\n");
%!     fputs(fd, "       gravity;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: 1,\n");
%!     fputs(fd, "       position, reference, global, 0.05, 0.02, 0.03,\n");
%!     fputs(fd, "       orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!     fputs(fd, "       velocity, reference, global, null,\n");
%!     fputs(fd, "       angular velocity, reference, global, null;\n");
%!     fputs(fd, "reference: 2,\n");
%!     fputs(fd, "       position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!     fputs(fd, "       orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!     fputs(fd, "       velocity, reference, 1, null,\n");
%!     fputs(fd, "       angular velocity, reference, 1, 0., 0., omega0;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "	      structural: 1, dynamic,\n");
%!     fputs(fd, "		       position, reference, 2, 0.07, 0.08, 0.02,\n");
%!     fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!     fputs(fd, "		       velocity, reference, 2, null,\n");
%!     fputs(fd, "		       angular velocity, reference, 2, null;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "	 body: 1, 1, \n");
%!     fputs(fd, "	          # mass\n");
%!     fputs(fd, "		  m,\n");
%!     fputs(fd, "		  # center of mass\n");
%!     fputs(fd, "		  null,\n");
%!     fputs(fd, "		  # inertia matrix\n");
%!     fputs(fd, "		  diag, .1, 1., 100.;\n");
%!     fputs(fd, "	joint: 1, total pin joint,\n");
%!     fputs(fd, "		1, \n");
%!     fputs(fd, "		# relative position\n");
%!     fputs(fd, "		position, reference, 2, null,\n");
%!     fputs(fd, "		# relative orientation matrix\n");
%!     fputs(fd, "		position orientation, reference, 2, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!     fputs(fd, "		# absolute position\n");
%!     fputs(fd, "		position, reference, 2, null,\n");
%!     fputs(fd, "                position orientation, reference, 2, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, 2, eye,\n");
%!     fputs(fd, "                position constraint,\n");
%!     fputs(fd, "                active, active, active, null,\n");
%!     fputs(fd, "                orientation constraint,\n");
%!     fputs(fd, "                active,\n");
%!     fputs(fd, "                active,\n");
%!     fputs(fd, "                inactive,\n");
%!     fputs(fd, "                null;\n");
%!     fputs(fd, "	couple: 1, absolute, \n");
%!     fputs(fd, "		1, \n");
%!     fputs(fd, "		position, reference, node, null, \n");
%!     fputs(fd, "		component, 100., \n");
%!     fputs(fd, "			   100., \n");
%!     fputs(fd, "			   1000.;\n");
%!     fputs(fd, " gravity: uniform, component,\n");
%!     fputs(fd, "        gx, gy, gz;\n");
%!     fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST14
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, "set: real r = 10.;\n");
%!     fputs(fd, "set: real m = 1.;\n");
%!     fputs(fd, "set: real gx = 0;\n");
%!     fputs(fd, "set: real gy = -9.81;\n");
%!     fputs(fd, "set: real gz = 0.;\n");
%!     fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!     fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!     fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!     fputs(fd, "set: real initial_time = 0;\n");
%!     fputs(fd, "set: real final_time = 0.1;\n");
%!     fputs(fd, "set: integer N = 720;\n");
%!     fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "    threads: assembly, 1;\n");
%!     fputs(fd, "	   initial time: initial_time;\n");
%!     fputs(fd, "	   final time: final_time;\n");
%!     fputs(fd, "	   time step: time_step;\n");
%!     fputs(fd, "    max iterations: 100;\n");
%!     fputs(fd, "    linear solver: umfpack, cc;\n");
%!     fputs(fd, "    method: ms, 0.;\n");
%!     fputs(fd, "    nonlinear solver: nox,\n");
%!     fputs(fd, "        jacobian operator, newton,\n");
%!     fputs(fd, "        solver, line search based,\n");
%!     fputs(fd, "        forcing term, type 2,\n");
%!     fputs(fd, "        verbose, yes,\n");
%!     fputs(fd, "        print convergence info, no,\n");
%!     fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!     fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!     fputs(fd, "        linear solver, gmres,\n");
%!     fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!     fputs(fd, "        linear solver max iterations, 10,\n");
%!     fputs(fd, "        krylov subspace size, 10,\n");
%!     fputs(fd, "        minimum step, 1e-12,\n");
%!     fputs(fd, "        use preconditioner as solver, yes;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "       tolerance: 1e-10;\n");
%!     fputs(fd, "       max iterations: 100;\n");
%!     fputs(fd, "	      structural nodes: 1;\n");
%!     fputs(fd, "	      rigid bodies: 1;\n");
%!     fputs(fd, "	      joints: 1;\n");
%!     fputs(fd, "       forces: 1;\n");
%!     fputs(fd, "       gravity;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: 1,\n");
%!     fputs(fd, "       position, reference, global, 0.05, 0.02, 0.03,\n");
%!     fputs(fd, "       orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!     fputs(fd, "       velocity, reference, global, null,\n");
%!     fputs(fd, "       angular velocity, reference, global, null;\n");
%!     fputs(fd, "reference: 2,\n");
%!     fputs(fd, "       position, reference, 1, r * cos(Phi0), r * sin(Phi0), 0.,\n");
%!     fputs(fd, "       orientation, reference, 1, euler123, 20. * pi / 180., 16 * pi / 180., 45. * pi / 180.,\n");
%!     fputs(fd, "       velocity, reference, 1, null,\n");
%!     fputs(fd, "       angular velocity, reference, 1, 0., 0., omega0;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "	      structural: 1, dynamic,\n");
%!     fputs(fd, "		       position, reference, 2, 0.07, 0.08, 0.02,\n");
%!     fputs(fd, "                orientation, reference, 2, euler123, 10. * pi / 180., 20. * pi / 180., 90. * pi / 180., \n");
%!     fputs(fd, "		       velocity, reference, 2, null,\n");
%!     fputs(fd, "		       angular velocity, reference, 2, null;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "	 body: 1, 1, \n");
%!     fputs(fd, "	          # mass\n");
%!     fputs(fd, "		  m,\n");
%!     fputs(fd, "		  # center of mass\n");
%!     fputs(fd, "		  null,\n");
%!     fputs(fd, "		  # inertia matrix\n");
%!     fputs(fd, "		  diag, .1, 1., 100.;\n");
%!     fputs(fd, "	joint: 1, total pin joint,\n");
%!     fputs(fd, "		1, \n");
%!     fputs(fd, "		# relative position\n");
%!     fputs(fd, "		position, reference, 2, null,\n");
%!     fputs(fd, "		# relative orientation matrix\n");
%!     fputs(fd, "		position orientation, reference, 2, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, 2, euler123, 0.1*pi/180., 0.2*pi/180., 0.3*pi/180.,\n");
%!     fputs(fd, "		# absolute position\n");
%!     fputs(fd, "		position, reference, 2, null,\n");
%!     fputs(fd, "                position orientation, reference, 2, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, 2, eye,\n");
%!     fputs(fd, "                position constraint,\n");
%!     fputs(fd, "                active, active, active, null,\n");
%!     fputs(fd, "                orientation constraint,\n");
%!     fputs(fd, "                active,\n");
%!     fputs(fd, "                active,\n");
%!     fputs(fd, "                inactive,\n");
%!     fputs(fd, "                null;\n");
%!     fputs(fd, "	couple: 1, absolute, \n");
%!     fputs(fd, "		1, \n");
%!     fputs(fd, "		position, reference, node, null, \n");
%!     fputs(fd, "		component, 100., \n");
%!     fputs(fd, "			   100., \n");
%!     fputs(fd, "			   1000.;\n");
%!     fputs(fd, " gravity: uniform, component,\n");
%!     fputs(fd, "        gx, gy, gz;\n");
%!     fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST15
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, "set: real r = 10.;\n");
%!     fputs(fd, "set: real m = 1.;\n");
%!     fputs(fd, "set: real gx = 0;\n");
%!     fputs(fd, "set: real gy = -9.81;\n");
%!     fputs(fd, "set: real gz = 0.;\n");
%!     fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!     fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!     fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!     fputs(fd, "set: real initial_time = 0;\n");
%!     fputs(fd, "set: real final_time = 0.1;\n");
%!     fputs(fd, "set: integer N = 720;\n");
%!     fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "    threads: assembly, 1;\n");
%!     fputs(fd, "	   initial time: initial_time;\n");
%!     fputs(fd, "	   final time: final_time;\n");
%!     fputs(fd, "	   time step: time_step;\n");
%!     fputs(fd, "    max iterations: 100;\n");
%!     fputs(fd, "    linear solver: aztecoo, tolerance, 1e-8, max iterations, 3, preconditioner, klu;\n");
%!     fputs(fd, "    method: ms, 0.;\n");
%!     fputs(fd, "    nonlinear solver: nox,\n");
%!     fputs(fd, "        jacobian operator, newton,\n");
%!     fputs(fd, "        solver, line search based,\n");
%!     fputs(fd, "        forcing term, type 2,\n");
%!     fputs(fd, "        verbose, yes,\n");
%!     fputs(fd, "        print convergence info, no,\n");
%!     fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!     fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!     fputs(fd, "        linear solver, gmres,\n");
%!     fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!     fputs(fd, "        linear solver max iterations, 10,\n");
%!     fputs(fd, "        krylov subspace size, 10,\n");
%!     fputs(fd, "        minimum step, 1e-12,\n");
%!     fputs(fd, "        use preconditioner as solver, yes;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "       tolerance: 1e-10;\n");
%!     fputs(fd, "       max iterations: 100;\n");
%!     fputs(fd, "	      structural nodes: 1;\n");
%!     fputs(fd, "	      rigid bodies: 1;\n");
%!     fputs(fd, "	      joints: 1;\n");
%!     fputs(fd, "       forces: 1;\n");
%!     fputs(fd, "       gravity;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: 1,\n");
%!     fputs(fd, "       position, reference, global, 0.5, 0.02, 0.03,\n");
%!     fputs(fd, "       orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!     fputs(fd, "       velocity, reference, global, null,\n");
%!     fputs(fd, "       angular velocity, reference, global, null;\n");
%!     fputs(fd, "reference: 2,\n");
%!     fputs(fd, "       position, reference, 1, null,\n");
%!     fputs(fd, "       orientation, reference, 1, eye,\n");
%!     fputs(fd, "       velocity, reference, 1, null,\n");
%!     fputs(fd, "       angular velocity, reference, 1, 0., 0., omega0;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "	      structural: 1, dynamic,\n");
%!     fputs(fd, "		       position, reference, 2, 1.2, 2.3, 4.5,\n");
%!     fputs(fd, "                orientation, reference, 2, euler123, 45. * pi / 180., 30. * pi / 180, 70. * pi / 180., \n");
%!     fputs(fd, "		       velocity, reference, 2, null,\n");
%!     fputs(fd, "		       angular velocity, reference, 2, null;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "	 body: 1, 1, \n");
%!     fputs(fd, "	          # mass\n");
%!     fputs(fd, "		  m,\n");
%!     fputs(fd, "		  # center of mass\n");
%!     fputs(fd, "		  null,\n");
%!     fputs(fd, "		  # inertia matrix\n");
%!     fputs(fd, "		  diag, .1, 1., 100.;\n");
%!     fputs(fd, "	joint: 1, total pin joint,\n");
%!     fputs(fd, "		1, \n");
%!     fputs(fd, "		# relative position\n");
%!     fputs(fd, "		position, reference, 2, null,\n");
%!     fputs(fd, "		# relative orientation matrix\n");
%!     fputs(fd, "		position orientation, reference, 2, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, 2, eye,\n");
%!     fputs(fd, "		# absolute position\n");
%!     fputs(fd, "		position, reference, 2, null,\n");
%!     fputs(fd, "                position orientation, reference, 2, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, 2, eye,\n");
%!     fputs(fd, "                position constraint,\n");
%!     fputs(fd, "                active, active, active, null,\n");
%!     fputs(fd, "                orientation constraint,\n");
%!     fputs(fd, "                active,\n");
%!     fputs(fd, "                active,\n");
%!     fputs(fd, "                inactive,\n");
%!     fputs(fd, "                null;\n");
%!     fputs(fd, "	couple: 1, absolute, \n");
%!     fputs(fd, "		1, \n");
%!     fputs(fd, "		position, reference, node, null, \n");
%!     fputs(fd, "		component, 100., \n");
%!     fputs(fd, "			   100., \n");
%!     fputs(fd, "			   1000.;\n");
%!     fputs(fd, " gravity: uniform, component,\n");
%!     fputs(fd, "        gx, gy, gz;\n");
%!     fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST16
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, "set: real r = 10.;\n");
%!     fputs(fd, "set: real m = 1.;\n");
%!     fputs(fd, "set: real gx = 0;\n");
%!     fputs(fd, "set: real gy = -9.81;\n");
%!     fputs(fd, "set: real gz = 0.;\n");
%!     fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!     fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!     fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!     fputs(fd, "set: real initial_time = 0;\n");
%!     fputs(fd, "set: real final_time = 0.1;\n");
%!     fputs(fd, "set: integer N = 720;\n");
%!     fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "    threads: assembly, 1;\n");
%!     fputs(fd, "	   initial time: initial_time;\n");
%!     fputs(fd, "	   final time: final_time;\n");
%!     fputs(fd, "	   time step: time_step;\n");
%!     fputs(fd, "    max iterations: 100;\n");
%!     fputs(fd, "    linear solver: aztecoo, tolerance, 1e-8, max iterations, 3, preconditioner, klu;\n");
%!     fputs(fd, "    method: ms, 0.;\n");
%!     fputs(fd, "    nonlinear solver: nox,\n");
%!     fputs(fd, "        jacobian operator, newton,\n");
%!     fputs(fd, "        solver, line search based,\n");
%!     fputs(fd, "        forcing term, type 2,\n");
%!     fputs(fd, "        verbose, yes,\n");
%!     fputs(fd, "        print convergence info, no,\n");
%!     fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!     fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!     fputs(fd, "        linear solver, gmres,\n");
%!     fputs(fd, "        linear solver tolerance, 1e-10,\n");
%!     fputs(fd, "        linear solver max iterations, 10,\n");
%!     fputs(fd, "        krylov subspace size, 10,\n");
%!     fputs(fd, "        minimum step, 1e-12,\n");
%!     fputs(fd, "        use preconditioner as solver, yes;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       tolerance: 1e-10;\n");
%!     fputs(fd, "       max iterations: 100;\n");
%!     fputs(fd, "	      structural nodes: 1;\n");
%!     fputs(fd, "	      rigid bodies: 1;\n");
%!     fputs(fd, "	      joints: 1;\n");
%!     fputs(fd, "       forces: 1;\n");
%!     fputs(fd, "       gravity;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: 1,\n");
%!     fputs(fd, "       position, reference, global, 0.5, 0.02, 0.03,\n");
%!     fputs(fd, "       orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!     fputs(fd, "       velocity, reference, global, null,\n");
%!     fputs(fd, "       angular velocity, reference, global, null;\n");
%!     fputs(fd, "reference: 2,\n");
%!     fputs(fd, "       position, reference, 1, null,\n");
%!     fputs(fd, "       orientation, reference, 1, eye,\n");
%!     fputs(fd, "       velocity, reference, 1, null,\n");
%!     fputs(fd, "       angular velocity, reference, 1, 0., 0., omega0;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "	      structural: 1, dynamic,\n");
%!     fputs(fd, "		       position, reference, 2, 1.2, 2.3, 4.5,\n");
%!     fputs(fd, "                orientation, reference, 2, euler123, 45. * pi / 180., 30. * pi / 180, 70. * pi / 180., \n");
%!     fputs(fd, "		       velocity, reference, 2, null,\n");
%!     fputs(fd, "		       angular velocity, reference, 2, null;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "	 body: 1, 1, \n");
%!     fputs(fd, "	          # mass\n");
%!     fputs(fd, "		  m,\n");
%!     fputs(fd, "		  # center of mass\n");
%!     fputs(fd, "		  null,\n");
%!     fputs(fd, "		  # inertia matrix\n");
%!     fputs(fd, "		  diag, .1, 1., 100.;\n");
%!     fputs(fd, "	joint: 1, total pin joint,\n");
%!     fputs(fd, "		1, \n");
%!     fputs(fd, "		# relative position\n");
%!     fputs(fd, "		position, reference, 2, null,\n");
%!     fputs(fd, "		# relative orientation matrix\n");
%!     fputs(fd, "		position orientation, reference, 2, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, 2, eye,\n");
%!     fputs(fd, "		# absolute position\n");
%!     fputs(fd, "		position, reference, 2, null,\n");
%!     fputs(fd, "                position orientation, reference, 2, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, 2, eye,\n");
%!     fputs(fd, "                position constraint,\n");
%!     fputs(fd, "                active, active, active, null,\n");
%!     fputs(fd, "                orientation constraint,\n");
%!     fputs(fd, "                active,\n");
%!     fputs(fd, "                active,\n");
%!     fputs(fd, "                inactive,\n");
%!     fputs(fd, "                null;\n");
%!     fputs(fd, "	couple: 1, absolute, \n");
%!     fputs(fd, "		1, \n");
%!     fputs(fd, "		position, reference, node, null, \n");
%!     fputs(fd, "		component, 100., \n");
%!     fputs(fd, "			   100., \n");
%!     fputs(fd, "			   1000.;\n");
%!     fputs(fd, " gravity: uniform, component,\n");
%!     fputs(fd, "        gx, gy, gz;\n");
%!     fputs(fd, "end: elements;\n");
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
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST17
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, "set: real r = 10.;\n");
%!     fputs(fd, "set: real m = 1.;\n");
%!     fputs(fd, "set: real gx = 0;\n");
%!     fputs(fd, "set: real gy = -9.81;\n");
%!     fputs(fd, "set: real gz = 0.;\n");
%!     fputs(fd, "set: real Phi0 = -90 * pi / 180;\n");
%!     fputs(fd, "set: real omega0 = 1300 * pi / 180;\n");
%!     fputs(fd, "set: real delta_Phi = 2*pi/360.;\n");
%!     fputs(fd, "set: real initial_time = 0;\n");
%!     fputs(fd, "set: real final_time = 0.1;\n");
%!     fputs(fd, "set: integer N = 720;\n");
%!     fputs(fd, "set: real time_step = 2 * pi / (N * max(1,abs(omega0)));\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "    threads: assembly, 1;\n");
%!     fputs(fd, "	   initial time: initial_time;\n");
%!     fputs(fd, "	   final time: final_time;\n");
%!     fputs(fd, "	   time step: time_step;\n");
%!     fputs(fd, "    max iterations: 100;\n");
%!     fputs(fd, "    linear solver: aztecoo, tolerance, 1e-6, max iterations, 3, preconditioner, klu;\n");
%!     fputs(fd, "    method: ms, 0.;\n");
%!     fputs(fd, "    nonlinear solver: nox,\n");
%!     fputs(fd, "        jacobian operator, newton,\n");
%!     fputs(fd, "        solver, line search based,\n");
%!     fputs(fd, "        forcing term, type 2,\n");
%!     fputs(fd, "        forcing term min tolerance, 1e-8,\n");
%!     fputs(fd, "        forcing term max tolerance, 1e-3,\n");
%!     fputs(fd, "        verbose, yes,\n");
%!     fputs(fd, "        print convergence info, no,\n");
%!     fputs(fd, "        weighted rms absolute tolerance, 1e-3,\n");
%!     fputs(fd, "        weighted rms relative tolerance, 1e-2,\n");
%!     fputs(fd, "        linear solver, gmres,\n");
%!     fputs(fd, "        linear solver max iterations, 30,\n");
%!     fputs(fd, "        krylov subspace size, 30,\n");
%!     fputs(fd, "        minimum step, 1e-12,\n");
%!     fputs(fd, "        use preconditioner as solver, no;\n");
%!     fputs(fd, "        tolerance: 1e-6, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-3, 0;\n");
%!     fputs(fd, "        derivatives max iterations: 20;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-6;\n");
%!     fputs(fd, "        output: iterations, CPU time, solver condition number, stat, yes;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       tolerance: 1e-10;\n");
%!     fputs(fd, "       max iterations: 100;\n");
%!     fputs(fd, "	      structural nodes: 1;\n");
%!     fputs(fd, "	      rigid bodies: 1;\n");
%!     fputs(fd, "	      joints: 1;\n");
%!     fputs(fd, "       forces: 1;\n");
%!     fputs(fd, "       gravity;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: 1,\n");
%!     fputs(fd, "       position, reference, global, 0.5, 0.02, 0.03,\n");
%!     fputs(fd, "       orientation, reference, global, euler123, 30. * pi / 180., 40. * pi / 180., 70. * pi / 180.,\n");
%!     fputs(fd, "       velocity, reference, global, null,\n");
%!     fputs(fd, "       angular velocity, reference, global, null;\n");
%!     fputs(fd, "reference: 2,\n");
%!     fputs(fd, "       position, reference, 1, null,\n");
%!     fputs(fd, "       orientation, reference, 1, eye,\n");
%!     fputs(fd, "       velocity, reference, 1, null,\n");
%!     fputs(fd, "       angular velocity, reference, 1, 0., 0., omega0;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "	      structural: 1, dynamic,\n");
%!     fputs(fd, "		       position, reference, 2, 1.2, 2.3, 4.5,\n");
%!     fputs(fd, "                orientation, reference, 2, euler123, 45. * pi / 180., 30. * pi / 180, 70. * pi / 180., \n");
%!     fputs(fd, "		       velocity, reference, 2, null,\n");
%!     fputs(fd, "		       angular velocity, reference, 2, null;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "	 body: 1, 1, \n");
%!     fputs(fd, "	          # mass\n");
%!     fputs(fd, "		  m,\n");
%!     fputs(fd, "		  # center of mass\n");
%!     fputs(fd, "		  null,\n");
%!     fputs(fd, "		  # inertia matrix\n");
%!     fputs(fd, "		  diag, .1, 1., 100.;\n");
%!     fputs(fd, "	joint: 1, total pin joint,\n");
%!     fputs(fd, "		1, \n");
%!     fputs(fd, "		# relative position\n");
%!     fputs(fd, "		position, reference, 2, null,\n");
%!     fputs(fd, "		# relative orientation matrix\n");
%!     fputs(fd, "		position orientation, reference, 2, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, 2, eye,\n");
%!     fputs(fd, "		# absolute position\n");
%!     fputs(fd, "		position, reference, 2, null,\n");
%!     fputs(fd, "                position orientation, reference, 2, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, 2, eye,\n");
%!     fputs(fd, "                position constraint,\n");
%!     fputs(fd, "                active, active, active, null,\n");
%!     fputs(fd, "                orientation constraint,\n");
%!     fputs(fd, "                active,\n");
%!     fputs(fd, "                active,\n");
%!     fputs(fd, "                inactive,\n");
%!     fputs(fd, "                null;\n");
%!     fputs(fd, "	couple: 1, absolute, \n");
%!     fputs(fd, "		1, \n");
%!     fputs(fd, "		position, reference, node, null, \n");
%!     fputs(fd, "		component, 100., \n");
%!     fputs(fd, "			   100., \n");
%!     fputs(fd, "			   100000.;\n");
%!     fputs(fd, " gravity: uniform, component,\n");
%!     fputs(fd, "        gx, gy, gz;\n");
%!     fputs(fd, "end: elements;\n");
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
%!   options.mbdyn_command = "mbdyn";
%!   mbdyn_solver_run(fname, options);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory, deformation, velocity, acceleration, node_id] = mbdyn_post_load_output_struct(options.output_file);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_abs_to_rel_XXXXXX"));
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
%!     fputs(fd, "         final time: 0.1;\n");
%!     fputs(fd, "         time step: 1e-5;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "         output: iterations;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");
%!     fputs(fd, "     output meter: closest next, 0, forever, const, 1e-2;\n");
%!     fputs(fd, "     output precision: 16;\n");
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " reference: 1,\n");
%!     fputs(fd, "     position, reference, global, 1.1, 1.2, 1.3,\n");
%!     fputs(fd, "     orientation, reference, global, euler123, -0.7 * pi, 0.5 * pi, -1.2 * pi,\n");
%!     fputs(fd, "     velocity, reference, global, 1., 2., 3.,\n");
%!     fputs(fd, "     angular velocity, reference, global, 0.1, 0.2, 0.3;\n");
%!     fputs(fd, " reference: 2,\n");
%!     fputs(fd, "     position, reference, 1, 2.1, 2.2, 2.3,\n");
%!     fputs(fd, "     orientation, reference, 1, euler123, 0.1 * pi, 0.2 * pi, 0.3 * pi,\n");
%!     fputs(fd, "     velocity, reference, 1, null,\n");
%!     fputs(fd, "     angular velocity, reference, 1, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, eye,");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 reference, 2, eye,");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 reference, 2, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 10.1, 10.2, 10.3, 1., 0., 0., mult, time, const, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 20.1, 20.2, 20.3, 0., 1., 0., mult, time, const, F2;\n");
%!     fputs(fd, "         joint: 1, total joint,\n");
%!     fputs(fd, "                   1,\n");
%!     fputs(fd, "                    position, reference, 1, null,\n");
%!     fputs(fd, "                    position orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    rotation orientation, reference, 1, eye,\n");
%!     fputs(fd, "                   2,\n");
%!     fputs(fd, "                    position, reference, 1, null,\n");
%!     fputs(fd, "                    position orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    rotation orientation, reference, 1, eye,\n");
%!     fputs(fd, "                    position constraint, active, active, active, null,\n");
%!     fputs(fd, "                    orientation constraint, active, active, active, null;\n");
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
%!   mbdyn_post_abs_to_rel(2, options.output_file, [options.output_file, "_rel"], false);
%!   copyfile([options.output_file, ".log"], [options.output_file, "_rel.log"]);
%!   copyfile([options.output_file, ".out"], [options.output_file, "_rel.out"]);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [t, trajectory_abs, deformation_abs, velocity_abs, acceleration_abs, node_id] = mbdyn_post_load_output_struct(options.output_file, 1:2);
%!   [t, trajectory_rel, deformation_rel, velocity_rel, acceleration_rel, node_id] = mbdyn_post_load_output_struct([options.output_file, "_rel"], 1:2);
%!   tol = eps^0.4;
%!   for i=1:numel(node_id)
%!     assert(trajectory_rel{i}, repmat(trajectory_rel{i}(1,:), numel(t), 1), tol * max(max(abs(trajectory_abs{i}))));
%!     assert(velocity_rel{i}, zeros(numel(t), 6), tol * max(max(abs(velocity_abs{i}))));
%!   endfor
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
