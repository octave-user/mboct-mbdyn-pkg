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
## @deftypefn {Function File} [@var{node_id}, @var{beta}, @var{gamma}, @var{beta_dot}, @var{gamma_dot}] = mbdyn_post_load_output_ine(@var{mbdyn_filename}, @var{filter_node_id})
##
## Loads the MBDyn .ine file named "<@var{mbdyn_filename}.ine>".
##
## @var{node_id} @dots{} The node identifier
##
## @var{beta} @dots{} Three components of the momentum in the absolute reference frame
##
## @var{gamma} @dots{} The three components of the momenta moment in the absolute reference frame,
## with respect to the coordinates of the node, thus to a moving frame
##
## @var{beta_dot} @dots{} The three components of the derivative of the momentum
##
## @var{gamma_dot} @dots{} The three components of the derivative of the momentum moment
##
## @end deftypefn

function [node_id, beta, gamma, beta_dot, gamma_dot] = mbdyn_post_load_output_ine(mbdyn_filename, filter_node_id, append_rows)
  if (nargin < 1 || nargin > 3 || nargout > 5)
    print_usage();
  endif

  if (nargin < 2)
    filter_node_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif
  
  ine_filename = mbdyn_post_output_filename(mbdyn_filename, ".ine");

  column_count = (nargout - 1) * 3;

  [node_id, data] = mbdyn_post_load_output(ine_filename, column_count, filter_node_id, append_rows);

  for i=1:length(data)
    if (nargout >= 2)
      beta{i} = data{i}(:, 1:3);
    endif
    
    if (nargout >= 3)
      gamma{i} = data{i}(:, 4:6);
    endif
    
    if (nargout >= 4)
      beta_dot{i} = data{i}(:, 7:9);
    endif
    
    if (nargout >= 5)
      gamma_dot{i} = data{i}(:, 10:12);
    endif
  endfor
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_ine_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real m = 1.234;\n");
%!     fputs(fd, " set: real J11 = 0.011;\n");
%!     fputs(fd, " set: real J22 = 0.022;\n");
%!     fputs(fd, " set: real J33 = 0.033;\n");
%!     fputs(fd, " set: real omega1 = 1.123;\n");
%!     fputs(fd, " set: real omega2 = 2.234;\n");
%!     fputs(fd, " set: real omega3 = 3.456;\n");
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
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 omega1, omega2, omega3,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m, null, diag, J11, J22, J33;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
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
%!   [t] = mbdyn_post_load_output_out(fname, 1024, false);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   g = log_dat.vars.g;
%!   m = log_dat.vars.m;
%!   J11 = log_dat.vars.J11;
%!   J22 = log_dat.vars.J22;
%!   J33 = log_dat.vars.J33;
%!   J = diag([J11, J22, J33]);
%!   omega = [log_dat.vars.omega1; log_dat.vars.omega2; log_dat.vars.omega3];
%!   vref = [0; 0; -g] * t.';
%!   betaref = m * vref;
%!   gammaref = J * omega;
%!   betadotref = [0; 0; -m * g];
%!   gammadotref = zeros(3, 1);
%!   [node_id, beta, gamma, beta_dot, gamma_dot] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta, gamma, beta_dot] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta, gamma] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   assert_simple(node_id, int32(1));
%!   tol = eps^0.4;
%!   assert_simple(beta{1}, betaref.', tol * max(max(abs(betaref))));
%!   assert_simple(gamma{1}, repmat(gammaref.', numel(t), 1), tol * norm(gammaref));
%!   assert_simple(beta_dot{1}, repmat(betadotref.', numel(t), 1), tol * norm(betadotref));
%!   assert_simple(gamma_dot{1}, repmat(gammadotref.', numel(t), 1), tol * max([1,norm(gammadotref)]));
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_ine_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real m = 1.234;\n");
%!     fputs(fd, " set: real J11 = 0.011;\n");
%!     fputs(fd, " set: real J22 = 0.022;\n");
%!     fputs(fd, " set: real J33 = 0.033;\n");
%!     fputs(fd, " set: real omega1 = 1.123;\n");
%!     fputs(fd, " set: real omega2 = 2.234;\n");
%!     fputs(fd, " set: real omega3 = 3.456;\n");
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
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 omega1, omega2, omega3,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m, null, diag, J11, J22, J33;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
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
%!   [t] = mbdyn_post_load_output_out(fname, 1024, false);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   g = log_dat.vars.g;
%!   m = log_dat.vars.m;
%!   J11 = log_dat.vars.J11;
%!   J22 = log_dat.vars.J22;
%!   J33 = log_dat.vars.J33;
%!   J = diag([J11, J22, J33]);
%!   omega = [log_dat.vars.omega1; log_dat.vars.omega2; log_dat.vars.omega3];
%!   vref = [0; 0; -g] * t.';
%!   betaref = m * vref;
%!   gammaref = J * omega;
%!   betadotref = [0; 0; -m * g];
%!   gammadotref = zeros(3, 1);
%!   [node_id, beta, gamma, beta_dot, gamma_dot] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta, gamma, beta_dot] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta, gamma] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   assert_simple(node_id, int32(1));
%!   tol = eps^0.4;
%!   assert_simple(beta{1}, betaref.', tol * max(max(abs(betaref))));
%!   assert_simple(gamma{1}, repmat(gammaref.', numel(t), 1), tol * norm(gammaref));
%!   assert_simple(beta_dot{1}, repmat(betadotref.', numel(t), 1), tol * norm(betadotref));
%!   assert_simple(gamma_dot{1}, repmat(gammadotref.', numel(t), 1), tol * max([1,norm(gammadotref)]));
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_ine_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real m = 1.234;\n");
%!     fputs(fd, " set: real J11 = 0.011;\n");
%!     fputs(fd, " set: real J22 = 0.022;\n");
%!     fputs(fd, " set: real J33 = 0.033;\n");
%!     fputs(fd, " set: real omega1 = 1.123;\n");
%!     fputs(fd, " set: real omega2 = 2.234;\n");
%!     fputs(fd, " set: real omega3 = 3.456;\n");
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
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 omega1, omega2, omega3,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m, null, diag, J11, J22, J33;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
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
%!   [t] = mbdyn_post_load_output_out(fname, 1024, false);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   g = log_dat.vars.g;
%!   m = log_dat.vars.m;
%!   J11 = log_dat.vars.J11;
%!   J22 = log_dat.vars.J22;
%!   J33 = log_dat.vars.J33;
%!   J = diag([J11, J22, J33]);
%!   omega = [log_dat.vars.omega1; log_dat.vars.omega2; log_dat.vars.omega3];
%!   vref = [0; 0; -g] * t.';
%!   betaref = m * vref;
%!   gammaref = J * omega;
%!   betadotref = [0; 0; -m * g];
%!   gammadotref = zeros(3, 1);
%!   [node_id, beta, gamma, beta_dot, gamma_dot] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta, gamma, beta_dot] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta, gamma] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   assert_simple(node_id, int32(1));
%!   tol = eps^0.4;
%!   assert_simple(beta{1}, betaref.', tol * max(max(abs(betaref))));
%!   assert_simple(gamma{1}, repmat(gammaref.', numel(t), 1), tol * norm(gammaref));
%!   assert_simple(beta_dot{1}, repmat(betadotref.', numel(t), 1), tol * norm(betadotref));
%!   assert_simple(gamma_dot{1}, repmat(gammadotref.', numel(t), 1), tol * max([1,norm(gammadotref)]));
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_ine_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real m = 1.234;\n");
%!     fputs(fd, " set: real J11 = 0.011;\n");
%!     fputs(fd, " set: real J22 = 0.022;\n");
%!     fputs(fd, " set: real J33 = 0.033;\n");
%!     fputs(fd, " set: real omega1 = 1.123;\n");
%!     fputs(fd, " set: real omega2 = 2.234;\n");
%!     fputs(fd, " set: real omega3 = 3.456;\n");
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
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 omega1, omega2, omega3,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m, null, diag, J11, J22, J33;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
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
%!   [t] = mbdyn_post_load_output_out(fname, 1024, false);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   g = log_dat.vars.g;
%!   m = log_dat.vars.m;
%!   J11 = log_dat.vars.J11;
%!   J22 = log_dat.vars.J22;
%!   J33 = log_dat.vars.J33;
%!   J = diag([J11, J22, J33]);
%!   omega = [log_dat.vars.omega1; log_dat.vars.omega2; log_dat.vars.omega3];
%!   vref = [0; 0; -g] * t.';
%!   betaref = m * vref;
%!   gammaref = J * omega;
%!   betadotref = [0; 0; -m * g];
%!   gammadotref = zeros(3, 1);
%!   [node_id, beta, gamma, beta_dot, gamma_dot] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta, gamma, beta_dot] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta, gamma] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id, beta] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   [node_id] = mbdyn_post_load_output_ine(options.output_file, 1, numel(t));
%!   assert_simple(node_id, int32(1));
%!   tol = eps^0.4;
%!   assert_simple(beta{1}, betaref.', tol * max(max(abs(betaref))));
%!   assert_simple(gamma{1}, repmat(gammaref.', numel(t), 1), tol * norm(gammaref));
%!   assert_simple(beta_dot{1}, repmat(betadotref.', numel(t), 1), tol * norm(betadotref));
%!   assert_simple(gamma_dot{1}, repmat(gammadotref.', numel(t), 1), tol * max([1,norm(gammadotref)]));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
