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
## @deftypefn {Function File} [@var{node_id}, @var{y}, @var{ydot}] = mbdyn_post_load_output_abs(@var{output_file})
##
## Load output data from abstract nodes from MBDyn output file @var{output_file}.
##
## @var{node_id} @dots{} Array of node numbers.
##
## @var{y} @dots{} Cell array of abstract node values.
##
## @var{ydot} @dots{} Cell array of abstract node time derivatives.
##
## @end deftypefn

function [node_id, y, ydot] = mbdyn_post_load_output_abs(abs_filename, filter_node_id, append_rows)
  if (nargin < 1 || nargin > 3 || nargout > 3)
    print_usage();
  endif

  if (nargin < 2)
    filter_node_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif
  
  abs_filename = mbdyn_post_output_filename(abs_filename, ".abs");

  [node_id, data] = mbdyn_post_load_output(abs_filename, 2, filter_node_id, append_rows);
  
  for i=1:numel(data)
    y{i} = data{i}(:,1);
    
    if (nargout >= 3)
      if (columns(data{i}) >= 2)
        ydot{i} = data{i}(:, 2);
      else
        ydot{i} = [];
      endif
    endif
  endfor
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_abs_XXXXXX"));
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
%!     fputs(fd, "     abstract nodes: 2;\n");
%!     fputs(fd, "     genels: 2;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         abstract: 1, algebraic;\n");
%!     fputs(fd, "         abstract: 2, algebraic;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         genel: 1, clamp, 1, abstract, element, 1, body, string,\"E\", direct;\n");
%!     fputs(fd, "         genel: 2, clamp, 2, abstract, element, 1, body, string,\"V\", direct;\n");
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
%!   [t] = mbdyn_post_load_output_out(fname, 1024, false);
%!   g = -9.81;
%!   m = 1;
%!   v = g * t;
%!   s = 0.5 * g * t.^2;
%!   E = 0.5 * m * v.^2;
%!   V = -m * g * s;
%!   [node_id, abs_value, abs_deriv] = mbdyn_post_load_output_abs(options.output_file, [1, 2], numel(t));
%!   assert(abs_value{1}, E, 1e-4);
%!   assert(abs_value{2}, V, 1e-4);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_abs_XXXXXX"));
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
%!     fputs(fd, "     abstract nodes: 2;\n");
%!     fputs(fd, "     genels: 2;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         abstract: 1, algebraic;\n");
%!     fputs(fd, "         abstract: 2, algebraic;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         genel: 1, clamp, 1, abstract, element, 1, body, string,\"E\", direct;\n");
%!     fputs(fd, "         genel: 2, clamp, 2, abstract, element, 1, body, string,\"V\", direct;\n");
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
%!   [t] = mbdyn_post_load_output_out(fname, 1024, false);
%!   g = -9.81;
%!   m = 1;
%!   v = g * t;
%!   s = 0.5 * g * t.^2;
%!   E = 0.5 * m * v.^2;
%!   V = -m * g * s;
%!   [node_id, abs_value, abs_deriv] = mbdyn_post_load_output_abs(options.output_file, [1, 2], numel(t));
%!   assert(abs_value{1}, E, 1e-4);
%!   assert(abs_value{2}, V, 1e-4);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
