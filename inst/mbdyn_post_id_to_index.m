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
## @deftypefn {Function File} @var{vars} = mbdyn_post_id_to_index(@var{res},@var{vars_in})
##
## Convert node and element id numbers to index values (e.g. vars.node_id_* => vars.node_idx_*).
## 
## @var{res} @dots{} Data structure with id numbers to be converted.
##
## @var{res}.node_id @dots{} Vector of node id's.
##
## @var{res}.force_id @dots{} Vector of force id's.
##
## @var{res}.torque_id @dots{} Vector of torque id's.
##
## @var{res}.joint_id @dots{} Vector of joint id's.
##
## @var{res}.elem_id @dots{} Vector of element id's.
##
## @var{vars_in} @dots{} Data structure returned from mbdyn_post_load_log_vars.
##
## @var{vars} @dots{} Data structure with the same information like <@var{vars_in}>
## plus the requested indices (e.g. node_id_xxx -> node_idx_xxx).
##
## @end deftypefn

function vars = mbdyn_post_id_to_index(res, vars_in)
  if (nargin ~= 2 || nargout > 1)
    print_usage();
  endif
  
  vars = vars_in;

  id_names = {"node",       "node";
              "abs_node",   "abs_node";
              "elec_node",  "elec_node";
              "prm_node",   "prm_node";
              "hyd_node",   "hyd_node";
              "therm_node", "therm_node";
              "force",      "force";
              "force",      "torque";
              "joint",      "joint";
              "elem",       "elem";
              "genel",      "genel";
              "drive",      "drive";
              "trace",      "trace"};

  for i=1:rows(id_names)
    if (isfield(res, [id_names{i, 1}, "_id"]))
      vars = mbdyn_post_node_id_to_node_index(getfield(res, [id_names{i, 1}, "_id"]), vars, [id_names{i, 2}, "_id_"], [id_names{i, 2}, "_idx_"], vars);
    endif
  endfor
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_mov_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer node_id_1 = 1001;\n");
%!     fputs(fd, " set: integer body_id_1 = 2001;\n");
%!     fputs(fd, " set: integer force_id_1 = 3001;\n");
%!     fputs(fd, " set: real m1 = 2.5;\n");
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real J11 = 1.;\n");
%!     fputs(fd, " set: real J22 = 1.;\n");
%!     fputs(fd, " set: real J33 = 1.;\n");
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
%!     fputs(fd, "         structural: node_id_1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: body_id_1, node_id_1, m1, null, diag, J11, J22, J33;\n");
%!     fputs(fd, "         force: force_id_1, absolute, node_id_1, position, null, 1., 0., 0, F1;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = true;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation, velocity, acceleration, res.node_id, force, res.force_id, force_node_id] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   log_dat.vars = mbdyn_post_id_to_index(res, log_dat.vars);
%!   g = -log_dat.vars.g;
%!   F = force{log_dat.vars.force_idx_1}(:, 1);
%!   m = log_dat.vars.m1;
%!   assert_simple(res.node_id, int32(log_dat.vars.node_id_1));
%!   assert_simple(t, (0:0.1:1.1).', 1e-5);
%!   assert_simple(trajectory{log_dat.vars.node_idx_1}, ...
%!          [0.5 * F / m .* t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g .* t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert_simple(velocity{log_dat.vars.node_idx_1}, ...
%!         [t .* F / m, ...
%!          zeros(numel(t), 1), g .* t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert_simple(acceleration{log_dat.vars.node_idx_1}, ...
%!         [F / m, ...
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_mov_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer node_id_1 = 1001;\n");
%!     fputs(fd, " set: integer body_id_1 = 2001;\n");
%!     fputs(fd, " set: integer force_id_1 = 3001;\n");
%!     fputs(fd, " set: real m1 = 2.5;\n");
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real J11 = 1.;\n");
%!     fputs(fd, " set: real J22 = 1.;\n");
%!     fputs(fd, " set: real J33 = 1.;\n");
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
%!     fputs(fd, "         structural: node_id_1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: body_id_1, node_id_1, m1, null, diag, J11, J22, J33;\n");
%!     fputs(fd, "         force: force_id_1, absolute, node_id_1, position, null, 1., 0., 0, F1;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = true;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation, velocity, acceleration, res.node_id, force, res.force_id, force_node_id] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   log_dat.vars = mbdyn_post_id_to_index(res, log_dat.vars);
%!   g = -log_dat.vars.g;
%!   F = force{log_dat.vars.force_idx_1}(:, 1);
%!   m = log_dat.vars.m1;
%!   assert_simple(res.node_id, int32(log_dat.vars.node_id_1));
%!   assert_simple(t, (0:0.1:1.1).', 1e-5);
%!   assert_simple(trajectory{log_dat.vars.node_idx_1}, ...
%!          [0.5 * F / m .* t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g .* t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert_simple(velocity{log_dat.vars.node_idx_1}, ...
%!         [t .* F / m, ...
%!          zeros(numel(t), 1), g .* t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert_simple(acceleration{log_dat.vars.node_idx_1}, ...
%!         [F / m, ...
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_mov_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer node_id_1 = 1001;\n");
%!     fputs(fd, " set: integer body_id_1 = 2001;\n");
%!     fputs(fd, " set: integer force_id_1 = 3001;\n");
%!     fputs(fd, " set: real m1 = 2.5;\n");
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real J11 = 1.;\n");
%!     fputs(fd, " set: real J22 = 1.;\n");
%!     fputs(fd, " set: real J33 = 1.;\n");
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
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: node_id_1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: body_id_1, node_id_1, m1, null, diag, J11, J22, J33;\n");
%!     fputs(fd, "         force: force_id_1, absolute, node_id_1, position, null, 1., 0., 0, F1;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = true;
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation, velocity, acceleration, res.node_id, force, res.force_id, force_node_id] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   log_dat.vars = mbdyn_post_id_to_index(res, log_dat.vars);
%!   g = -log_dat.vars.g;
%!   F = force{log_dat.vars.force_idx_1}(:, 1);
%!   m = log_dat.vars.m1;
%!   assert_simple(res.node_id, int32(log_dat.vars.node_id_1));
%!   assert_simple(t, (0:0.1:1.1).', 1e-5);
%!   assert_simple(trajectory{log_dat.vars.node_idx_1}, ...
%!          [0.5 * F / m .* t.^2, ...
%!           zeros(numel(t), 1), ...
%!           0.5 * g .* t.^2, ...
%!           zeros(numel(t), 3)], 1e-5);
%!  assert_simple(velocity{log_dat.vars.node_idx_1}, ...
%!         [t .* F / m, ...
%!          zeros(numel(t), 1), g .* t, ...
%!          zeros(numel(t), 3)], 1e-5);
%!  assert_simple(acceleration{log_dat.vars.node_idx_1}, ...
%!         [F / m, ...
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
