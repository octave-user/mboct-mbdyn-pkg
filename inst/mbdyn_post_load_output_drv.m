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
## @deftypefn {Function File} [@var{drive_id}, @var{drive_value}] = mbdyn_post_load_output_drv(@var{mbdyn_filename}, @var{filter_drive_id}, @var{append_rows})
##
## Loads output from drive callers from MBDyn output file "<@var{mbdyn_filename}.drv>".
##
## @var{filter_drive_id} @dots{} Load only drives with drive numbers in <@var{filter_drive_id}>. If empty, everything is loaded.
##
## @var{drive_id} @dots{} Drive identifier of the loaded drives.
##
## @var{append_rows} @dots{} The expected number of rows for which memory should be allocated.
##
## @end deftypefn

function [drive_id, drive_value] = mbdyn_post_load_output_drv(mbdyn_filename, filter_drive_id, append_rows, auto_resize_rows)
  if (nargin < 1 || nargout > 2)
    print_usage();
  endif

  if (nargin < 2)
    filter_drive_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  if (nargin < 4)
    auto_resize_rows = true;
  endif

  drv_filename = mbdyn_post_output_filename(mbdyn_filename, ".drv");

  column_count = 1;

  [drive_id, drive_value] = mbdyn_post_load_output(drv_filename, column_count, filter_drive_id, append_rows, 2, 1, auto_resize_rows);
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_drv_XXXXXX"));
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
%!     fputs(fd, "         drive caller: 1, node, 1, structural, string, \"X[1]\", direct, output, value, yes;\n");
%!     fputs(fd, "         drive caller: 2, node, 1, structural, string, \"X[3]\", direct, output, value, yes;\n");
%!     fputs(fd, "         drive caller: 3, node, 1, structural, string, \"XP[1]\", direct, output, value, yes;\n");
%!     fputs(fd, "         drive caller: 4, node, 1, structural, string, \"XP[3]\", direct, output, value, yes;\n");
%!     fputs(fd, "         drive caller: 5, node, 1, structural, string, \"XPP[1]\", direct, output, value, yes;\n");
%!     fputs(fd, "         drive caller: 6, node, 1, structural, string, \"XPP[3]\", direct, output, value, yes;\n");
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
%!   t = mbdyn_post_load_output_out(fname, 1024, false);
%!   [drive_id, drive_value] = mbdyn_post_load_output_drv(fname, 1:6, numel(t));
%!   g = -9.81;
%!   F = 100;
%!   m = 1;
%!   tol = 1e-12;
%!   assert(drive_value{1}, 0.5 * F / m * t.^2, tol);
%!   assert(drive_value{2}, 0.5 * g * t.^2, tol);
%!   assert(drive_value{3}, t * F / m, tol);
%!   assert(drive_value{4}, g * t, tol);
%!   assert(drive_value{5}, repmat(F / m, numel(t), 1), tol);
%!   assert(drive_value{6}, repmat(g, numel(t), 1), tol);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_drv_XXXXXX"));
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
%!     fputs(fd, "         drive caller: 1, node, 1, structural, string, \"X[1]\", direct, output, value, yes;\n");
%!     fputs(fd, "         drive caller: 2, node, 1, structural, string, \"X[3]\", direct, output, value, yes;\n");
%!     fputs(fd, "         drive caller: 3, node, 1, structural, string, \"XP[1]\", direct, output, value, yes;\n");
%!     fputs(fd, "         drive caller: 4, node, 1, structural, string, \"XP[3]\", direct, output, value, yes;\n");
%!     fputs(fd, "         drive caller: 5, node, 1, structural, string, \"XPP[1]\", direct, output, value, yes;\n");
%!     fputs(fd, "         drive caller: 6, node, 1, structural, string, \"XPP[3]\", direct, output, value, yes;\n");
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
%!   t = mbdyn_post_load_output_out(fname, 1024, false);
%!   [drive_id, drive_value] = mbdyn_post_load_output_drv(fname, 1:6, numel(t));
%!   g = -9.81;
%!   F = 100;
%!   m = 1;
%!   tol = 1e-12;
%!   assert(drive_value{1}, 0.5 * F / m * t.^2, tol);
%!   assert(drive_value{2}, 0.5 * g * t.^2, tol);
%!   assert(drive_value{3}, t * F / m, tol);
%!   assert(drive_value{4}, g * t, tol);
%!   assert(drive_value{5}, repmat(F / m, numel(t), 1), tol);
%!   assert(drive_value{6}, repmat(g, numel(t), 1), tol);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
