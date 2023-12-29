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
## @deftypefn {Function File} [@var{inertia}] = mbdyn_post_inertia_compute(@var{body_groups}, @var{bodies}, @var{nodes})
## @deftypefnx {} [@var{inertia}] = mbdyn_post_inertia_compute(@var{body_groups}, @var{mbdyn_filename})
##
## Computes the mass, center of gravity and momentum of inertia of groups of bodies.
##
## @var{body_groups}.name @dots{} Name of the group of bodies
##
## @var{body_groups}.labels @dots{} Array of labels of all bodies which belong to this group
##
## @var{bodies} @dots{} Data structure returned from mbdyn_post_load_log_body
##
## @var{nodes} @dots{} Data structure returned from mbdyn_post_load_log
##
## @var{mbdyn_filename} @dots{} Name of mbdyn output files
##
## @end deftypefn

function [inertia] = mbdyn_post_inertia_compute(body_groups = "all", varargin)
  if (nargin < 2 || nargin > 3 || nargout > 1)
    print_usage();
  endif

  inertia = struct();

  if (~isstruct(body_groups) && ~ischar(body_groups))
    error("body_groups must be a cell array or a string!");
  endif

  if (length(varargin) == 1 && ischar(varargin{1}))
    mbdyn_filename = varargin{1};
    bodies = mbdyn_post_load_log_body(mbdyn_filename);
    nodes = mbdyn_post_load_log_node(mbdyn_filename);
  elseif (length(varargin) == 2 && isstruct(varargin{1}) && isstruct(varargin{2}))
    bodies = varargin{1};
    nodes = varargin{2};
  else
    print_usage();
    return;
  endif

  N = length(bodies);

  if (ischar(body_groups))
    switch(body_groups)
      case 'all'
        body_groups = struct();
        body_groups(1).labels =  [ bodies.label ];
        body_groups(1).name = 'all';
      case 'every'
        body_groups = struct();
        for i=1:length(bodies)
          body_groups(i).labels = [bodies(i).label];
          body_groups(i).name = sprintf('body: %d',bodies(i).label);
        endfor
      otherwise
        error("invalid argument: body_groups=\"%s\"",body_groups);
    endswitch
  endif

  for i=1:length(body_groups)

    inertia(i).name = body_groups(i).name;
    inertia(i).dm = 0;
    inertia(i).Xgc = zeros(1,3);
    inertia(i).J = zeros(3,3);
    Mgc = zeros(3,1);

    for j=1:length(body_groups(i).labels)
      idx_body = find([bodies.label] == body_groups(i).labels(j));

      if (length(idx_body) ~= 1)
        error("bodies.label is not unique");
      endif

      if (length(idx_body) == 0)
        error("body %d not found!",body_groups(i).labels(j));
      endif

      idx_node = find([nodes.label] == bodies(idx_body).node);

      if (length(idx_body) ~= 1)
        error("body label not found or not unique!");
      endif

      if (length(idx_node)  ~= 1)
        error("node label not found or not unique!");
      endif

      m_j = bodies(idx_body).dm;
      J_j = bodies(idx_body).J;
      R_j = nodes(idx_node).R0;
      dXgc_j = R_j * bodies(idx_body).Xgc;
      Xgc_j = dXgc_j + nodes(idx_node).X0;

      Mgc +=  Xgc_j * m_j;
      inertia(i).dm += m_j;
      inertia(i).J += R_j * J_j * R_j.' + (skew(dXgc_j) * skew(dXgc_j) - skew(Xgc_j) * skew(Xgc_j)) * m_j;
    endfor
    inertia(i).Xgc = Mgc / inertia(i).dm;
    inertia(i).J += skew(inertia(i).Xgc) * skew(inertia(i).Xgc) * inertia(i).dm;
    inertia(i).bodies = sort(body_groups(i).labels);
  endfor
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_inertia_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real m1 = 123;\n");
%!     fputs(fd, " set: real m2 = 0.456;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1e-3;\n");
%!     fputs(fd, "         time step: 1e-4;\n");
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
%!     fputs(fd, "     output precision: 16;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 0.5, 0.6, 0.7,\n");
%!     fputs(fd, "                 euler123, 0.3 * pi, 0.7 * pi, -1.4 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "                 euler123, pi / 4., -pi / 3., pi / 2.,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, reference, global, null, diag, 0.456, 0.789, 0.123, inertial, reference, global, eye;\n");
%!     fputs(fd, "         body: 2, 2, m2, reference, global, null, diag, 45, 78, 123, inertial, reference, global, eye;\n");
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
%!   bodies = mbdyn_post_load_log_body(fname);
%!   log_dat = mbdyn_post_load_log(fname);
%!   body_groups(1).labels = [1, 2];
%!   body_groups(1).name = "body 1 + body 2";
%!   inertia = mbdyn_post_inertia_compute(body_groups, bodies, log_dat.nodes);
%!   mbdyn_post_inertia_print(inertia, [fname, "inertia.dat"]);
%!   tol = 1e-5;
%!   assert_simple(inertia(1).dm, 123.456, tol * 123.456);
%!   assert_simple(inertia(1).Xgc, zeros(3, 1), tol);
%!   assert_simple(inertia(1).J, diag([45.456, 78.789, 123.123]), tol * 123.123);
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
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_inertia_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real m1 = 123;\n");
%!     fputs(fd, " set: real m2 = 0.456;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1e-3;\n");
%!     fputs(fd, "         time step: 1e-4;\n");
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
%!     fputs(fd, "     output precision: 16;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 0.5, 0.6, 0.7,\n");
%!     fputs(fd, "                 euler123, 0.3 * pi, 0.7 * pi, -1.4 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "                 euler123, pi / 4., -pi / 3., pi / 2.,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, reference, global, null, diag, 0.456, 0.789, 0.123, inertial, reference, global, eye;\n");
%!     fputs(fd, "         body: 2, 2, m2, reference, global, null, diag, 45, 78, 123, inertial, reference, global, eye;\n");
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
%!   bodies = mbdyn_post_load_log_body(fname);
%!   log_dat = mbdyn_post_load_log(fname);
%!   body_groups(1).labels = [1, 2];
%!   body_groups(1).name = "body 1 + body 2";
%!   inertia = mbdyn_post_inertia_compute(body_groups, bodies, log_dat.nodes);
%!   mbdyn_post_inertia_print(inertia, [fname, "inertia.dat"]);
%!   tol = 1e-5;
%!   assert_simple(inertia(1).dm, 123.456, tol * 123.456);
%!   assert_simple(inertia(1).Xgc, zeros(3, 1), tol);
%!   assert_simple(inertia(1).J, diag([45.456, 78.789, 123.123]), tol * 123.123);
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
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_inertia_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real m1 = 123;\n");
%!     fputs(fd, " set: real m2 = 0.456;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1e-3;\n");
%!     fputs(fd, "         time step: 1e-4;\n");
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
%!     fputs(fd, "     output precision: 16;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 0.5, 0.6, 0.7,\n");
%!     fputs(fd, "                 euler123, 0.3 * pi, 0.7 * pi, -1.4 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "                 euler123, pi / 4., -pi / 3., pi / 2.,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, reference, global, null, diag, 0.456, 0.789, 0.123, inertial, reference, global, eye;\n");
%!     fputs(fd, "         body: 2, 2, m2, reference, global, null, diag, 45, 78, 123, inertial, reference, global, eye;\n");
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
%!   bodies = mbdyn_post_load_log_body(fname);
%!   log_dat = mbdyn_post_load_log(fname);
%!   body_groups(1).labels = [1, 2];
%!   body_groups(1).name = "body 1 + body 2";
%!   inertia = mbdyn_post_inertia_compute(body_groups, bodies, log_dat.nodes);
%!   mbdyn_post_inertia_print(inertia, [fname, "inertia.dat"]);
%!   tol = 1e-5;
%!   assert_simple(inertia(1).dm, 123.456, tol * 123.456);
%!   assert_simple(inertia(1).Xgc, zeros(3, 1), tol);
%!   assert_simple(inertia(1).J, diag([45.456, 78.789, 123.123]), tol * 123.123);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_inertia_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real m1 = 123;\n");
%!     fputs(fd, " set: real m2 = 0.456;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1e-3;\n");
%!     fputs(fd, "         time step: 1e-4;\n");
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

%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     output precision: 16;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 0.5, 0.6, 0.7,\n");
%!     fputs(fd, "                 euler123, 0.3 * pi, 0.7 * pi, -1.4 * pi,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "                 euler123, pi / 4., -pi / 3., pi / 2.,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, reference, global, null, diag, 0.456, 0.789, 0.123, inertial, reference, global, eye;\n");
%!     fputs(fd, "         body: 2, 2, m2, reference, global, null, diag, 45, 78, 123, inertial, reference, global, eye;\n");
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
%!   bodies = mbdyn_post_load_log_body(fname);
%!   log_dat = mbdyn_post_load_log(fname);
%!   body_groups(1).labels = [1, 2];
%!   body_groups(1).name = "body 1 + body 2";
%!   inertia = mbdyn_post_inertia_compute(body_groups, bodies, log_dat.nodes);
%!   mbdyn_post_inertia_print(inertia, [fname, "inertia.dat"]);
%!   tol = 1e-5;
%!   assert_simple(inertia(1).dm, 123.456, tol * 123.456);
%!   assert_simple(inertia(1).Xgc, zeros(3, 1), tol);
%!   assert_simple(inertia(1).J, diag([45.456, 78.789, 123.123]), tol * 123.123);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
