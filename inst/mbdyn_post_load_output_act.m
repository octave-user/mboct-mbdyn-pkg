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
## @deftypefn {Function File} [@var{beam_id}, @var{F_I}, @var{M_I}, @var{F_II}, @var{M_II}] = mbdyn_post_load_output_act(@var{mbdyn_filename})
## @deftypefnx {} [@dots{}] = mbdyn_post_load_output_act(@var{mbdyn_filename}, @var{filter_beam_id})
## @deftypefnx {} [@dots{}] = mbdyn_post_load_output_act(@var{mbdyn_filename}, @var{filter_beam_id}, @var{append_rows})
##
## Loads internal forces and moments of beam2 and beam3 elements from MBDyn output file "<@var{mbdyn_filename}>.act".
##
## @var{filter_beam_id} @dots{} Array of beam labels to be loaded. If this array is empty, all beams are loaded.
##
## @var{append_rows} @dots{} Hint for memory reallocation.
##
## @var{beam_id} @dots{} The label of the beam.
##
## @var{F_I} @dots{} The three components of the force at the first evaluation point,
## oriented according to the reference frame of that beam section.
##
## @var{M_I} @dots{} The three components of the moment at the first evaluation point,
## oriented according to the reference frame of that beam section.
##
## The three-node beam element generates six additional columns:
##
## @var{F_II} @dots{} The three components of the force at the second evaluation point,
## oriented according to the reference frame of that beam section;
##
## @var{M_II} @dots{} The three components of the moment at the second evaluation point,
## oriented according to the reference frame of that beam section.
##
## @end deftypefn

function [beam_id, F_I, M_I, F_II, M_II] = mbdyn_post_load_output_act(mbdyn_filename, filter_beam_id, append_rows)
  if (nargin < 1 || nargin > 3 || nargout > 5)
    print_usage();
  endif

  if (nargin < 2)
    filter_beam_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  act_filename = mbdyn_post_output_filename(mbdyn_filename, ".act");

  column_count = (nargout - 1) * 3;

  [beam_id, data] = mbdyn_post_load_output(act_filename, column_count, filter_beam_id, append_rows);

  for i=1:length(data)
    if (nargout >= 2)
      F_I{i} = data{i}(:,1:3);
    endif
    if (nargout >= 3)
      M_I{i} = data{i}(:,4:6);
    endif
    if (nargout >= 4)
      F_II{i} = data{i}(:,7:9);
    endif
    if (nargout >= 5)
      M_II{i} = data{i}(:,10:12);
    endif
  endfor
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_act_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real L = 0.5;\n");
%!     fputs(fd, " set: real F1 = 10;\n");
%!     fputs(fd, " set: real D = 50e-3;\n");
%!     fputs(fd, " set: real A = D^2 * pi / 4.;\n");
%!     fputs(fd, " set: real Ay = 9. / 10. * A;\n");
%!     fputs(fd, " set: real Az = Ay;\n");
%!     fputs(fd, " set: real Iy = D^4 * pi / 64.;\n");
%!     fputs(fd, " set: real Iz = Iy;\n");
%!     fputs(fd, " set: real It = Iy + Iz;\n");
%!     fputs(fd, " set: real E = 210000e6;\n");
%!     fputs(fd, " set: real G = 81500e6;\n");
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
%!     fputs(fd, "         tolerance: 1.e-5;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 3;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     beams: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 2, static,\n");
%!     fputs(fd, "                 0.5 * L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 3, static,\n");
%!     fputs(fd, "                 L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "         beam3: 1, 1, null, 2, null, 3, null, eye,\n");
%!     fputs(fd, "                linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz,\n");
%!     fputs(fd, "                same, same;\n");
%!     fputs(fd, "         force: 1, absolute, 3, position, null, 0., 0., -1, mult, time, F1;\n");
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
%!   log_dat = mbdyn_post_load_log(fname);
%!   [beam_id, F_I, M_I, F_II, M_II] = mbdyn_post_load_output_act(fname, 1, numel(t));
%!   tol = eps^0.4;
%!   assert(F_I{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert(F_II{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert(M_I{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 + 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
%!   assert(M_II{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 - 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_act_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real L = 0.5;\n");
%!     fputs(fd, " set: real F1 = 10;\n");
%!     fputs(fd, " set: real D = 50e-3;\n");
%!     fputs(fd, " set: real A = D^2 * pi / 4.;\n");
%!     fputs(fd, " set: real Ay = 9. / 10. * A;\n");
%!     fputs(fd, " set: real Az = Ay;\n");
%!     fputs(fd, " set: real Iy = D^4 * pi / 64.;\n");
%!     fputs(fd, " set: real Iz = Iy;\n");
%!     fputs(fd, " set: real It = Iy + Iz;\n");
%!     fputs(fd, " set: real E = 210000e6;\n");
%!     fputs(fd, " set: real G = 81500e6;\n");
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
%!     fputs(fd, "         tolerance: 1.e-5;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 3;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     beams: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 2, static,\n");
%!     fputs(fd, "                 0.5 * L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, "         structural: 3, static,\n");
%!     fputs(fd, "                 L, 0., 0.,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "         beam3: 1, 1, null, 2, null, 3, null, eye,\n");
%!     fputs(fd, "                linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz,\n");
%!     fputs(fd, "                same, same;\n");
%!     fputs(fd, "         force: 1, absolute, 3, position, null, 0., 0., -1, mult, time, F1;\n");
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
%!   log_dat = mbdyn_post_load_log(fname);
%!   [beam_id, F_I, M_I, F_II, M_II] = mbdyn_post_load_output_act(fname, 1, numel(t));
%!   tol = eps^0.4;
%!   assert(F_I{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert(F_II{1}(:, 3), -log_dat.vars.F1 * t, tol * abs(log_dat.vars.F1));
%!   assert(M_I{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 + 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
%!   assert(M_II{1}(:, 2), log_dat.vars.F1 * t * 0.5 * (1 - 1/sqrt(3)) * log_dat.vars.L, tol * abs(log_dat.vars.F1 * log_dat.vars.L));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
