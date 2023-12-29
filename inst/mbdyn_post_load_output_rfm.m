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
## @deftypefn {Function File} [ @var{ref_id}, @var{position}, @var{orientation}, @var{velocity}, @var{angular_velocity} ] = mbdyn_post_load_output_rfm(@var{mbdyn_filename},@var{filter_ref_id})
##
## Loads the MBDyn .rfm file named "<@var{mbdyn_filename}.rfm>".
##
## @var{ref_id} @dots{} Vector of reference frame numbers used in the MBDyn input file.
##
## @var{position}(@var{i}) @dots{} Position of the reference frame <@var{ref_id}(@var{i})>
##
## @var{orientation}(@var{i}) @dots{} Orientation of the reference frame <@var{ref_id}(@var{i})>
##
## @var{velocity}(@var{i}) @dots{} Velocity of the reference frame <@var{ref_id}(@var{i})>
##
## @var{angular_velocity}(@var{i}) @dots{} Angular velocity of the reference frame <@var{ref_id}(@var{i})>
##
## @var{filter_ref_id} @dots{} Load only the data of the corresponding reference frames if present.
## If empty, the data of all reference frames is loaded.
##
## @end deftypefn

function [ref_id, position, orientation, velocity, angular_velocity] = mbdyn_post_load_output_rfm(mbdyn_filename, filter_ref_id)
  if (nargin < 1 || nargin > 2 || nargout > 5)
    print_usage();
  endif

  if (nargin < 2)
    filter_ref_id = [];
  endif

  if (~ischar(mbdyn_filename))
    error("mbdyn_filename must be a string!");
  endif

  rfm_filename = mbdyn_post_output_filename(mbdyn_filename, ".rfm");

  def_col_size = 0;

  if (nargout >= 2)
    def_col_size = 3;
  endif

  if (nargout >= 3)
    def_col_size = 6;
  endif

  if (nargout >= 4)
    def_col_size = 9;
  endif

  if (nargout >= 5)
    def_col_size = 12;
  endif

  [ref_id, data] = mbdyn_post_load_output(rfm_filename, def_col_size, filter_ref_id);

  for i=1:length(data)
    if (nargout >= 2)
      position{i} = data{i}(:,1:3);
    endif

    if (nargout >= 3)
      ## FIXME: we have to consider the orientation description of each reference frame!
      orientation{i} = data{i}(:, 4:6) * pi / 180;
    endif

    if (nargout >= 4)
      velocity{i} = data{i}(:,7:9);
    endif

    if (nargout >= 5)
      angular_velocity{i} = data{i}(:,10:12);
    endif
  endfor
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_rfm_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1e-4;\n");
%!     fputs(fd, "         time step: 1e-6;\n");
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
%!     fputs(fd, "     default output: reference frames;\n");
%!     fputs(fd, "     default orientation: euler123;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " reference: 1, position, reference, global, 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "               orientation, reference, global, euler123, 0.11, 0.22, 0.33,\n");
%!     fputs(fd, "               velocity, reference, global, 10.1, 20.2, 30.3,\n");
%!     fputs(fd, "               angular velocity, reference, global, 100.1, 200.2, 300.3;\n");
%!     fputs(fd, " reference: 2, position, reference, global, 10.1, 20.2, 30.3,\n");
%!     fputs(fd, "               orientation, reference, global, euler123, 1.11e-2, 2.22e-2, 3.33e-2,\n");
%!     fputs(fd, "               velocity, reference, global, 100.1, 200.2, 300.3,\n");
%!     fputs(fd, "               angular velocity, reference, global, 2100.1, 2200.2, 2300.3;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, eye,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
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
%!   [ref_id, X, Phi, v, omega] = mbdyn_post_load_output_rfm(options.output_file, 1:2);
%!   tol = eps^0.4;
%!   assert_simple(ref_id, int32([1,2]));
%!   assert_simple(X{1}, [1.1, 2.2, 3.3], tol);
%!   assert_simple(Phi{1}, [0.11, 0.22, 0.33], 2 * pi * tol);
%!   assert_simple(v{1}, [10.1, 20.2, 30.3], tol);
%!   assert_simple(omega{1}, [100.1, 200.2, 300.3], tol);
%!   assert_simple(X{2}, [10.1, 20.2, 30.3], tol);
%!   assert_simple(Phi{2}, [1.11e-2, 2.22e-2, 3.33e-2], 2 * pi * tol);
%!   assert_simple(v{2}, [100.1, 200.2, 300.3], tol);
%!   assert_simple(omega{2}, [2100.1, 2200.2, 2300.3], tol);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_rfm_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1e-4;\n");
%!     fputs(fd, "         time step: 1e-6;\n");
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
%!     fputs(fd, "     default output: reference frames;\n");
%!     fputs(fd, "     default orientation: euler123;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " reference: 1, position, reference, global, 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "               orientation, reference, global, euler123, 0.11, 0.22, 0.33,\n");
%!     fputs(fd, "               velocity, reference, global, 10.1, 20.2, 30.3,\n");
%!     fputs(fd, "               angular velocity, reference, global, 100.1, 200.2, 300.3;\n");
%!     fputs(fd, " reference: 2, position, reference, global, 10.1, 20.2, 30.3,\n");
%!     fputs(fd, "               orientation, reference, global, euler123, 1.11e-2, 2.22e-2, 3.33e-2,\n");
%!     fputs(fd, "               velocity, reference, global, 100.1, 200.2, 300.3,\n");
%!     fputs(fd, "               angular velocity, reference, global, 2100.1, 2200.2, 2300.3;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, eye,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
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
%!   [ref_id, X, Phi, v, omega] = mbdyn_post_load_output_rfm(options.output_file, 1:2);
%!   tol = eps^0.4;
%!   assert_simple(ref_id, int32([1,2]));
%!   assert_simple(X{1}, [1.1, 2.2, 3.3], tol);
%!   assert_simple(Phi{1}, [0.11, 0.22, 0.33], 2 * pi * tol);
%!   assert_simple(v{1}, [10.1, 20.2, 30.3], tol);
%!   assert_simple(omega{1}, [100.1, 200.2, 300.3], tol);
%!   assert_simple(X{2}, [10.1, 20.2, 30.3], tol);
%!   assert_simple(Phi{2}, [1.11e-2, 2.22e-2, 3.33e-2], 2 * pi * tol);
%!   assert_simple(v{2}, [100.1, 200.2, 300.3], tol);
%!   assert_simple(omega{2}, [2100.1, 2200.2, 2300.3], tol);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_rfm_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1e-4;\n");
%!     fputs(fd, "         time step: 1e-6;\n");
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
%!     fputs(fd, "     default output: reference frames;\n");
%!     fputs(fd, "     default orientation: euler123;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " reference: 1, position, reference, global, 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "               orientation, reference, global, euler123, 0.11, 0.22, 0.33,\n");
%!     fputs(fd, "               velocity, reference, global, 10.1, 20.2, 30.3,\n");
%!     fputs(fd, "               angular velocity, reference, global, 100.1, 200.2, 300.3;\n");
%!     fputs(fd, " reference: 2, position, reference, global, 10.1, 20.2, 30.3,\n");
%!     fputs(fd, "               orientation, reference, global, euler123, 1.11e-2, 2.22e-2, 3.33e-2,\n");
%!     fputs(fd, "               velocity, reference, global, 100.1, 200.2, 300.3,\n");
%!     fputs(fd, "               angular velocity, reference, global, 2100.1, 2200.2, 2300.3;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, eye,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
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
%!   [ref_id, X, Phi, v, omega] = mbdyn_post_load_output_rfm(options.output_file, 1:2);
%!   tol = eps^0.4;
%!   assert_simple(ref_id, int32([1,2]));
%!   assert_simple(X{1}, [1.1, 2.2, 3.3], tol);
%!   assert_simple(Phi{1}, [0.11, 0.22, 0.33], 2 * pi * tol);
%!   assert_simple(v{1}, [10.1, 20.2, 30.3], tol);
%!   assert_simple(omega{1}, [100.1, 200.2, 300.3], tol);
%!   assert_simple(X{2}, [10.1, 20.2, 30.3], tol);
%!   assert_simple(Phi{2}, [1.11e-2, 2.22e-2, 3.33e-2], 2 * pi * tol);
%!   assert_simple(v{2}, [100.1, 200.2, 300.3], tol);
%!   assert_simple(omega{2}, [2100.1, 2200.2, 2300.3], tol);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_rfm_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1e-4;\n");
%!     fputs(fd, "         time step: 1e-6;\n");
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
%!     fputs(fd, "     default output: reference frames;\n");
%!     fputs(fd, "     default orientation: euler123;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " reference: 1, position, reference, global, 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "               orientation, reference, global, euler123, 0.11, 0.22, 0.33,\n");
%!     fputs(fd, "               velocity, reference, global, 10.1, 20.2, 30.3,\n");
%!     fputs(fd, "               angular velocity, reference, global, 100.1, 200.2, 300.3;\n");
%!     fputs(fd, " reference: 2, position, reference, global, 10.1, 20.2, 30.3,\n");
%!     fputs(fd, "               orientation, reference, global, euler123, 1.11e-2, 2.22e-2, 3.33e-2,\n");
%!     fputs(fd, "               velocity, reference, global, 100.1, 200.2, 300.3,\n");
%!     fputs(fd, "               angular velocity, reference, global, 2100.1, 2200.2, 2300.3;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, eye,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
%!     fputs(fd, "                 reference, 1, null,\n");
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
%!   [ref_id, X, Phi, v, omega] = mbdyn_post_load_output_rfm(options.output_file, 1:2);
%!   tol = eps^0.4;
%!   assert_simple(ref_id, int32([1,2]));
%!   assert_simple(X{1}, [1.1, 2.2, 3.3], tol);
%!   assert_simple(Phi{1}, [0.11, 0.22, 0.33], 2 * pi * tol);
%!   assert_simple(v{1}, [10.1, 20.2, 30.3], tol);
%!   assert_simple(omega{1}, [100.1, 200.2, 300.3], tol);
%!   assert_simple(X{2}, [10.1, 20.2, 30.3], tol);
%!   assert_simple(Phi{2}, [1.11e-2, 2.22e-2, 3.33e-2], 2 * pi * tol);
%!   assert_simple(v{2}, [100.1, 200.2, 300.3], tol);
%!   assert_simple(omega{2}, [2100.1, 2200.2, 2300.3], tol);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
