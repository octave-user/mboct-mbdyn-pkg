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
## @deftypefn {Function File} [@var{forces}] = mbdyn_post_load_log_force(@var{mbdyn_filename})
##
## Parses the MBDyn log-file "<@var{mbdyn_filename}>.log" and returns information about the structural forces.
##
## @var{forces}.label @dots{} The label of the structural force or couple.
##
## @var{forces}.node1 @dots{} The label of the first structural node the force is applied to.
##
## @var{forces}.arm1 @dots{}  The distance between the force and the first structural node.
##
## @var{forces}.node2 @dots{} If the force is an internal structural force node2 is the label of the second node the force is applied to.
##
## @var{forces}.arm2 @dots{}  If the force is an internal structural force arm2 is the distance between the force and the second node.
##
## @var{forces}.internal @dots{} Equal to one if the force is an internal structural force.
##
## @var{forces}.absolute @dots{} Equal to one if the force is a absolute force. Equal to zero if the force is a follower force.
##
## @var{force}.couple @dots{} Equal to zero if it is a force. Equal to one if it is a couple.
##
## @end deftypefn

function [forces] = mbdyn_post_load_log_force(mbdyn_filename)
  if (nargin < 1)
    print_usage();
  endif

  log_filename = mbdyn_post_output_filename(mbdyn_filename, ".log");

  fid = -1;

  unwind_protect
    [fid, msg] = fopen(log_filename, "rt");

    if (fid == -1)
      error("could not open file \"%s\": %s", log_filename, msg);
    endif

    i = 0;
    line_no = 0;

    forces = struct()([]);

    while (1)
      line = fgets(fid);

      if (~ischar(line) && line == -1)
        break;
      endif

      ++line_no;

      tag_end = find(line == ':');

      if (length(tag_end) >= 1)
        tag_end = tag_end(1);
        tag = line(1:tag_end);
        data = line(length(tag)+1:end);
        switch (tag)
          case {"structural internal follower force:",
                "structural internal follower couple:",
                "structural internal absolute force:",
                "structural internal absolute couple:",
                "structural follower force:",
                "structural follower couple:",
                "structural absolute force:",
                "structural absolute couple:"}

            ++i;

            [label, node1, arm1x, arm1y, arm1z, node2, arm2x, arm2y, arm2z, count] = sscanf(data,"%d %d %g %g %g %d %g %g %g", "C");

            if (count ~= 5 && count ~= 9)
              error("parse error in file \"%s\": line %d",log_filename,line_no);
            endif

            if (count >= 5)
              forces(i).label = int32(label);
              forces(i).node1 = int32(node1);
              forces(i).arm1 = [ arm1x; arm1y; arm1z ];
            endif

            if (count >= 9)
              forces(i).node2 = int32(node2);
              forces(i).arm2 = [ arm2x; arm2y; arm2z ];
            endif

            forces(i).internal = logical(length(strfind(tag, "internal")) > 0);
            forces(i).absolute = logical(length(strfind(tag, "absolute")) > 0);
            forces(i).couple   = logical(length(strfind(tag, "couple")) > 0);
        endswitch
      endif
    endwhile
  unwind_protect_cleanup
    if (fid ~= -1)
      fclose(fid);
    endif
  end_unwind_protect

  if (length(forces) > 0)
    [force_labels,force_idx] = sort([forces.label]);
    forces = forces(force_idx);
  endif
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_log_force_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 20;\n");
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
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1., 2., 3.,\n");
%!     fputs(fd, "                 euler123, 0.1, 0.2, 0.3,\n");
%!     fputs(fd, "                 10., 20., 30.,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 1.1, 1.2, 1.3, 1., 0., 0, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 2.1, 2.2, 2.3, 0., 1., 0, F2;\n");
%!     fputs(fd, "         force: 3, follower internal, 1, position, 3.1, 3.2, 3.3, 2, position, 4.1, 4.2, 4.3, 0., 1., 0, F2;\n");
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
%!   forces = mbdyn_post_load_log_force(fname);
%!   assert_simple(forces(1).label, int32(1));
%!   assert_simple(forces(1).node1, int32(1));
%!   assert_simple(forces(1).arm1, [1.1; 1.2; 1.3]);
%!   assert_simple(forces(2).label, int32(2));
%!   assert_simple(forces(2).node1, int32(2));
%!   assert_simple(forces(2).arm1, [2.1; 2.2; 2.3]);
%!   assert_simple(forces(3).label, int32(3));
%!   assert_simple(forces(3).node1, int32(1));
%!   assert_simple(forces(3).node2, int32(2));
%!   assert_simple(forces(3).arm1, [3.1; 3.2; 3.3]);
%!   assert_simple(forces(3).arm2, [4.1; 4.2; 4.3]);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_log_force_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 20;\n");
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
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1., 2., 3.,\n");
%!     fputs(fd, "                 euler123, 0.1, 0.2, 0.3,\n");
%!     fputs(fd, "                 10., 20., 30.,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 1.1, 1.2, 1.3, 1., 0., 0, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 2.1, 2.2, 2.3, 0., 1., 0, F2;\n");
%!     fputs(fd, "         force: 3, follower internal, 1, position, 3.1, 3.2, 3.3, 2, position, 4.1, 4.2, 4.3, 0., 1., 0, F2;\n");
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
%!   forces = mbdyn_post_load_log_force(fname);
%!   assert_simple(forces(1).label, int32(1));
%!   assert_simple(forces(1).node1, int32(1));
%!   assert_simple(forces(1).arm1, [1.1; 1.2; 1.3]);
%!   assert_simple(forces(2).label, int32(2));
%!   assert_simple(forces(2).node1, int32(2));
%!   assert_simple(forces(2).arm1, [2.1; 2.2; 2.3]);
%!   assert_simple(forces(3).label, int32(3));
%!   assert_simple(forces(3).node1, int32(1));
%!   assert_simple(forces(3).node2, int32(2));
%!   assert_simple(forces(3).arm1, [3.1; 3.2; 3.3]);
%!   assert_simple(forces(3).arm2, [4.1; 4.2; 4.3]);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_log_force_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 20;\n");
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
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1., 2., 3.,\n");
%!     fputs(fd, "                 euler123, 0.1, 0.2, 0.3,\n");
%!     fputs(fd, "                 10., 20., 30.,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 1.1, 1.2, 1.3, 1., 0., 0, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 2.1, 2.2, 2.3, 0., 1., 0, F2;\n");
%!     fputs(fd, "         force: 3, follower internal, 1, position, 3.1, 3.2, 3.3, 2, position, 4.1, 4.2, 4.3, 0., 1., 0, F2;\n");
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
%!   forces = mbdyn_post_load_log_force(fname);
%!   assert_simple(forces(1).label, int32(1));
%!   assert_simple(forces(1).node1, int32(1));
%!   assert_simple(forces(1).arm1, [1.1; 1.2; 1.3]);
%!   assert_simple(forces(2).label, int32(2));
%!   assert_simple(forces(2).node1, int32(2));
%!   assert_simple(forces(2).arm1, [2.1; 2.2; 2.3]);
%!   assert_simple(forces(3).label, int32(3));
%!   assert_simple(forces(3).node1, int32(1));
%!   assert_simple(forces(3).node2, int32(2));
%!   assert_simple(forces(3).arm1, [3.1; 3.2; 3.3]);
%!   assert_simple(forces(3).arm2, [4.1; 4.2; 4.3]);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_log_force_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 20;\n");
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
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1., 2., 3.,\n");
%!     fputs(fd, "                 euler123, 0.1, 0.2, 0.3,\n");
%!     fputs(fd, "                 10., 20., 30.,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, 1.1, 1.2, 1.3, 1., 0., 0, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, 2.1, 2.2, 2.3, 0., 1., 0, F2;\n");
%!     fputs(fd, "         force: 3, follower internal, 1, position, 3.1, 3.2, 3.3, 2, position, 4.1, 4.2, 4.3, 0., 1., 0, F2;\n");
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
%!   forces = mbdyn_post_load_log_force(fname);
%!   assert_simple(forces(1).label, int32(1));
%!   assert_simple(forces(1).node1, int32(1));
%!   assert_simple(forces(1).arm1, [1.1; 1.2; 1.3]);
%!   assert_simple(forces(2).label, int32(2));
%!   assert_simple(forces(2).node1, int32(2));
%!   assert_simple(forces(2).arm1, [2.1; 2.2; 2.3]);
%!   assert_simple(forces(3).label, int32(3));
%!   assert_simple(forces(3).node1, int32(1));
%!   assert_simple(forces(3).node2, int32(2));
%!   assert_simple(forces(3).arm1, [3.1; 3.2; 3.3]);
%!   assert_simple(forces(3).arm2, [4.1; 4.2; 4.3]);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
