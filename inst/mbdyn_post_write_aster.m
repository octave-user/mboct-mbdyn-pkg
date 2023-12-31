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
## @deftypefn {Function File} mbdyn_post_write_aster(@var{data}, @var{mesh_file})
##
## Converts a MBDyn model to Code_Aster .mail format.
##
## @var{data} @dots{} The result obtained from mbdyn_post_load_log.
##
## @var{mesh_file} @dots{} The name of the Code_Aster .mail file.
##
## @end deftypefn

function mbdyn_post_write_aster(data, mesh_file)
  if (nargin < 1 || nargin > 2)
    print_usage();
  endif

  if (nargin < 2)
    mesh_file = stdout;
  endif

  owns_fd = false;
  fd = -1;

  unwind_protect
    if (ischar(mesh_file))
      owns_fd = true;

      [fd, msg] = fopen(mesh_file, "wt");

      if (fd == -1)
        error("failed to open file \"%s\": %s", mesh_file, msg);
      endif
    else
      fd = mesh_file;

      if (~isscalar(fd))
        error("mesh_file must be a open file descriptor or a file name");
      endif
    endif


    fprintf(fd, "%% --------------------------------------------------------------------------------\n");
    fprintf(fd, " TITRE\n");
    fprintf(fd, "%s\n", ctime(time()));
    fprintf(fd, " FINSF\n");

    if (length(data.nodes) > 0)
      fprintf(fd, " %%\n");
      fprintf(fd, " COOR_3D\n");

      for i=1:length(data.nodes)
        fprintf(fd, " N%-8d ", data.nodes(i).label);

        for j=1:3
          fprintf(fd, "%.14e ", data.nodes(i).X0(j));
        endfor

        fprintf(fd, "\n");
      endfor

      fprintf(fd, " FINSF\n");
    endif

    if (length(data.beams2) > 0)
      fprintf(fd, " %%\n");

      fprintf(fd, " SEG2\n");

      for i=1:length(data.beams2)
        fprintf(fd, " M%-8d ", data.beams2(i).label);

        for j=1:2
          fprintf(fd, "N%-8d ", data.beams2(i).nodes(j).label);
        endfor

        fprintf(fd, "\n");
      endfor

      fprintf(fd, " FINSF\n");
    endif

    if (length(data.beams3) > 0)
      fprintf(fd, " %%\n");

      fprintf(fd, " SEG3\n");

      for i=1:length(data.beams3)
        fprintf(fd, " M%-8d ", data.beams3(i).label);

        for j=1:3
          fprintf(fd, "N%-8d ", data.beams3(i).nodes([1,3,2])(j).label);
        endfor

        fprintf(fd, "\n");
      endfor

      fprintf(fd, " FINSF\n");
    endif

    fprintf(fd, " FIN\n");
  unwind_protect_cleanup
    if (fd ~= -1 && owns_fd)
      fclose(fd);
    endif
  end_unwind_protect
endfunction

%!test
%! f_print_output = false;
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_write_aster_XXXXXX"));
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
%!   log_dat = mbdyn_post_load_log(fname);
%!   aster_file = [options.output_file, ".mail"];
%!   mbdyn_post_write_aster(log_dat, aster_file);
%!   if (f_print_output)
%!     spawn_wait(spawn("cat", {aster_file}));
%!   endif
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
%! f_print_output = false;
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_write_aster_XXXXXX"));
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
%!   log_dat = mbdyn_post_load_log(fname);
%!   aster_file = [options.output_file, ".mail"];
%!   mbdyn_post_write_aster(log_dat, aster_file);
%!   if (f_print_output)
%!     spawn_wait(spawn("cat", {aster_file}));
%!   endif
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
%! f_print_output = false;
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_write_aster_XXXXXX"));
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
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
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
%!   log_dat = mbdyn_post_load_log(fname);
%!   aster_file = [options.output_file, ".mail"];
%!   mbdyn_post_write_aster(log_dat, aster_file);
%!   if (f_print_output)
%!     spawn_wait(spawn("cat", {aster_file}));
%!   endif
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
%! f_print_output = true;
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_write_aster_XXXXXX"));
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
%!     fputs(fd, "     use automatic differentiation;\n");

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
%!   log_dat = mbdyn_post_load_log(fname);
%!   aster_file = [options.output_file, ".mail"];
%!   mbdyn_post_write_aster(log_dat, aster_file);
%!   if (f_print_output)
%!     spawn_wait(spawn("cat", {aster_file}));
%!   endif
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
