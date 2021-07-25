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
## @deftypefn {Function File} mbdyn_pre_write_fem_data(@var{fem_filename}, @var{Mred}, @var{Dred}, @var{Sred}, @var{Tred}, @var{X0}, @var{ured0}, @var{dured_dt0}, @var{diagM}, @var{m}, @var{Xgc}, @var{Jgc}, @var{node_list}, @var{Inv3}, @var{Inv4}, @var{Inv5}, @var{Inv8}, @var{Inv9})
##
## Creates a FEM input file for MBDyn's model element.
##
## @var{fem_filename} @dots{} output filename
##
## @var{Mred} @dots{} reduced order mass matrix
##
## @var{Sred} @dots{} reduced order stiffness matrix
##
## @var{Tred} @dots{} matrix of mode shapes at exposed nodes
##
## @var{X0} @dots{} initial position of exposed nodes
##
## @var{ured0} @dots{} initial modal displacements (optional, default zero)
##
## @var{dured_dt0} @dots{} initial modal velocities (optional, default zero)
##
## @var{diagM} @dots{} diagonal of the lumped mass matrix at exposed nodes
##
## @var{m} @dots{} total mass of the FEM model
##
## @var{Xgc} @dots{} center of gravity of the FEM model
##
## @var{Jgc} @dots{} momentum of inertia of the FEM model at center of gravity
##
## @var{node_list} @dots{} list of node numbers for exposed nodes
##
## @var{Inv3} @dots{} Static coupling between rigid body and FEM node displacements
##
## @var{Inv4} @dots{} Static coupling between rigid body rotations and FEM node displacements
##
## @var{Inv5} @dots{} Static coupling between FEM node displacements
##
## @var{Inv8} @dots{} Linear coupling terms
##
## @var{Inv9} @dots{} Quadratic coupling terms
##
## References:
## @nospell{MBDyn Theory and Developer's Manual Version 1.7.1}
## @nospell{Pierangelo Masarati}
## @nospell{Dipartimento di Ingegneria Aerospaziale}
## @nospell{Politecnico di Milano}
## @nospell{January 9, 2017}
## @end deftypefn

function mbdyn_pre_write_fem_data(fem_filename, Mred, Dred, Sred, Tred, X0, ured0, dured_dt0, diagM, m, Xgc, Jgc, node_list, Inv3, Inv4, Inv5, Inv8, Inv9)
  if (nargin < 5)
    print_usage();
  endif

  if (~ischar(fem_filename))
    error("fem_filename must be a string");
  endif

  if (~ismatrix(Mred) || ~issquare(Mred) || ~isreal(Mred))
    error("Mred must be square real matrix");
  endif

  N = columns(Mred);

  if (~ismatrix(Sred) || ~isreal(Sred) || rows(Sred) ~= N || columns(Sred) ~= N)
    error("Sred must be a square real matrix and must have the same dimensions like Mred");
  endif

  if (~ismatrix(Dred) || ~isreal(Dred) || numel(Dred) > 0 && (rows(Dred) ~= N || columns(Dred) ~= N))
    error("Dred must be a square real matrix and must have the same dimensions like Mred");
  endif

  if (~ismatrix(Tred) || ~isreal(Tred) || columns(Tred) ~= N)
    error("Tred must be a square real matrix with the same numer of columns like Mred");
  endif

  if (~ismatrix(X0) || ~isreal(X0) || rows(X0) ~= 3 || columns(X0) * 6 ~= rows(Tred))
    error("X0 must be a real 3 x rows(Tred) / 6 matrix");
  endif

  if (nargin >= 6)
    if (~isvector(ured0) || ~isreal(ured0) || columns(ured0) ~= 1 || rows(ured0) ~= N)
      error("ured0 must be real vector with the same number of rows like Mred");
    endif
  endif

  if (nargin >= 7)
    if (~isvector(dured_dt0) || ~isreal(dured_dt0) || columns(dured_dt0) ~= 1 || rows(dured_dt0) ~= N)
      error("dured_dt0 must be real vector with the same number of rows like Mred");
    endif
  endif

  if (nargin >= 8 && numel(diagM) ~= 0)
    if (~isvector(diagM) || ~isreal(diagM) || numel(diagM) ~= columns(X0) * 6)
      error("diagM must be a real vector with the same length as the number of columns of X0");
    endif
  endif

  if (nargin >= 14 && numel(Inv3) ~= 0)
    if (~ismatrix(Inv3) || ~isreal(Inv3) || rows(Inv3) ~= 3 || columns(Inv3) ~= N)
      error("Inv3 must be a real 3 x N matrix");
    endif
  endif

  if (nargin >= 15 && numel(Inv4) ~= 0)
    if (~ismatrix(Inv4) || ~isreal(Inv4) || rows(Inv4) ~= 3 || columns(Inv4) ~= N)
      error("Inv4 must be a real 3 x N matrix");
    endif
  endif

  if (nargin >= 16 && numel(Inv5) ~= 0)
    if (~isreal(Inv5) || rows(Inv5) ~= 3 || columns(Inv5) ~= N || size(Inv5, 3) ~= N)
      error("Inv5 must be a real 3 x N x N matrix");
    endif
  endif

  if (nargin >= 17 && numel(Inv8) ~= 0)
    if (~isreal(Inv8) || rows(Inv8) ~= 3 || columns(Inv8) ~= 3 || size(Inv8, 3) ~= N)
      error("Inv8 must be a real 3 x 3 x N matrix");
    endif
  endif

  if (nargin >= 18 && numel(Inv9) ~= 0)
    if (~isreal(Inv9) || rows(Inv9) ~= 3 || columns(Inv9) ~= 3 || size(Inv9, 3) ~= N  || size(Inv9, 4) ~= N)
      error("Inv9 must be a real 3 x 3 x N x N matrix");
    endif
  endif

  fd = -1;
  
  unwind_protect
    [fd, msg] = fopen(fem_filename, "wt");

    if (fd == -1)
      error("could not open file \"%s\": %s", fem_filename, msg);
    endif

    REVISION = "REV0";
    NODE = columns(X0);
    NORMAL = columns(Tred);
    ATTACHED = 0;
    CONSTRAINT = 0;
    REJECTED = 0;

    if (nargin < 13 || ~numel(node_list))
      node_list = 1:NODE;
    endif

    fprintf(fd, "** MBDyn MODAL DATA FILE\n");
    fprintf(fd, "** NODE SET \"ALL\"\n");

    fprintf(fd, "\n** RECORD GROUP 1, HEADER\n");
    fprintf(fd, "**   REVISION,  NODE,  NORMAL, ATTACHMENT, CONSTRAINT, REJECTED MODES.\n");
    fprintf(fd, "      %s             %d       %d          %d                 %d                 %d\n", REVISION, NODE, NORMAL, ATTACHED, CONSTRAINT, REJECTED);
    fprintf(fd, "**\n");

    fprintf(fd, "\n** RECORD GROUP 2, FINITE ELEMENT NODE LIST\n");
    fprintf(fd, "%d ", node_list);
    fprintf(fd, "\n**\n");

    if (nargin >= 6)
      fprintf(fd, "\n** RECORD GROUP 3, INITIAL MODAL DISPLACEMENTS\n");
      fprintf(fd, "%.16e ", ured0);
      fprintf(fd, "\n**\n");
    endif

    if (nargin >= 7)
      fprintf(fd, "\n** RECORD GROUP 4, INITIAL MODALVELOCITIES\n");
      fprintf(fd, "%.16e ", dured_dt0);
      fprintf(fd, "\n**\n");
    endif

    for i=1:rows(X0)
      fprintf(fd, "\n** RECORD GROUP %d, NODAL %s COORDINATES\n", i + 4, {"X", "Y", "Z"}{i});
      fprintf(fd, "%.16e\n", X0(i, :));
      fprintf(fd, "**\n");
    endfor

    fprintf(fd, "** RECORD GROUP 8, MODE SHAPES\n");

    for i=1:columns(Tred)
      fprintf(fd, "**    NORMAL MODE SHAPE #  %d\n", i);
      fprintf(fd, "%.16e %.16e %.16e %.16e %.16e %.16e\n", Tred(:, i));
    endfor

    fprintf(fd, "**\n");
    fprintf(fd, "** RECORD GROUP 9, MODAL MASS MATRIX\n");

    for i=1:rows(Mred)
      fprintf(fd, "%.16e ", Mred(i, :));
      fprintf(fd, "\n");
    endfor

    fprintf(fd, "**\n");
    fprintf(fd, "** RECORD GROUP 10, MODAL STIFFNESS MATRIX\n");

    for i=1:rows(Sred)
      fprintf(fd, "%.16e ", Sred(i, :));
      fprintf(fd, "\n");
    endfor

    fprintf(fd, "**\n");

    if (nargin >= 8 && numel(diagM) > 0)
      fprintf(fd, "** RECORD GROUP 11, DIAGONAL OF LUMPED MASS MATRIX\n");
      fprintf(fd, "%.16e %.16e %.16e %.16e %.16e %.16e\n", diagM);
      fprintf(fd, "**\n");
    endif

    if (nargin >= 11 && numel(Jgc) > 0)
      fprintf(fd, "** RECORD GROUP 12, RIGID BODY INERTIA MATRIX\n");
      fprintf(fd, "%.16e\n", m);
      fprintf(fd, "%.16e %.16e %.16e\n", Xgc);
      fprintf(fd, "%.16e %.16e %.16e\n", Jgc.');
      fprintf(fd, "**\n");
    endif

    if (numel(Dred) > 0)
      fprintf(fd, "**\n");
      fprintf(fd, "** RECORD GROUP 13, MODAL DAMPING MATRIX\n");

      for i=1:rows(Dred)
        fprintf(fd, "%.16e ", Dred(i, :));
        fprintf(fd, "\n");
      endfor

      fprintf(fd, "**\n");
    endif

    if (nargin >= 14 && numel(Inv3) > 0)
      fprintf(fd, "**\n");
      fprintf(fd, "** RECORD GROUP 14, INVARIANT 3\n");

      for i=1:rows(Inv3)
        fprintf(fd, "%.16e ", Inv3(i,:));
        fprintf(fd, "\n");
      endfor

      fprintf(fd, "**\n");
    endif

    if (nargin >= 15 && numel(Inv4) > 0)
      fprintf(fd, "**\n");
      fprintf(fd, "** RECORD GROUP 15, INVARIANT 4\n");

      for i=1:rows(Inv4)
        fprintf(fd, "%.16e ", Inv4(i, :));
        fprintf(fd, "\n");
      endfor

      fprintf(fd,"**\n");
    endif

    if (nargin >= 17 && numel(Inv8) > 0)
      fprintf(fd, "**\n");
      fprintf(fd, "** RECORD GROUP 16, INVARIANT 8\n");

      for i=1:rows(Inv8)
        for j=1:size(Inv8, 3)
          fprintf(fd, "%.16e ", Inv8(i, :,  j));
        endfor
        fprintf(fd, "\n");
      endfor

      fprintf(fd, "**\n");
    endif

    if (nargin >= 16 && numel(Inv5) > 0)
      fprintf(fd, "**\n");
      fprintf(fd, "** RECORD GROUP 17, INVARIANT 5\n");

      for i=1:rows(Inv5)
        for j=1:size(Inv5, 3)
          fprintf(fd, "%.16e ", Inv5(i, :, j));
        endfor
        fprintf(fd, "\n");
      endfor

      fprintf(fd, "**\n");
    endif

    if (nargin >= 16 && numel(Inv9) > 0)
      fprintf(fd, "**\n");
      fprintf(fd, "** RECORD GROUP 18, INVARIANT 9\n");

      for i=1:rows(Inv9)
        for j=1:size(Inv9, 3)
          for k=1:size(Inv9, 4)
            fprintf(fd, "%.16e ", Inv9(i, :, j, k));
          endfor
        endfor
        fprintf(fd, "\n");
      endfor

      fprintf(fd, "**\n");
    endif
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction

%!test
%! fd = -1;
%! omega1 = [  0,  0,  0,  0,   0];
%! omega2 = [0.5, 50,  5, 10, 100];
%! omega3 = [0.2, 10, 20, 50, 500];
%! m = 0.1;
%! s1 = 0;
%! s2 = 0;
%! s3 = 0.05;
%! J1 = 5e-3;
%! J2 = 5e-3;
%! J3 = 10e-3;
%! h3 = 0.075;
%! nsteps = 14400;
%! nrev = 1;
%! n0 = 0.1;
%! tolres = 1e-9;
%! silent = true;
%! for i=1:numel(omega1)
%!   unwind_protect
%!     unwind_protect
%!       [fd, fname1] = mkstemp(fullfile(tempdir(), "mbdyn_pre_write_fem_data1_XXXXXX"));
%!       if (fd == -1)
%!         error("failed to open temporary file");
%!       endif
%!       fprintf(fd, " set: real omega1 = %e;\n", omega1(i));
%!       fprintf(fd, " set: real omega2 = %e;\n", omega2(i));
%!       fprintf(fd, " set: real omega3 = %e;\n", omega3(i));
%!       fputs(fd, " set: real omega = sqrt(omega1^2 + omega2^2 + omega3^2);\n");
%!       fprintf(fd, " set: real m = %e;\n", m);
%!       fprintf(fd, " set: real s1 = %e;\n", s1);
%!       fprintf(fd, " set: real s2 = %e;\n", s2);
%!       fprintf(fd, " set: real s3 = %e;\n", s3);
%!       fprintf(fd, " set: real J1 = %e;\n", J1);
%!       fprintf(fd, " set: real J2 = %e;\n", J2);
%!       fprintf(fd, " set: real J3 = %e;\n", J3);
%!       fprintf(fd, " set: real h3 = %e;\n", h3);
%!       fprintf(fd, " set: real nsteps = %e;\n", nsteps);
%!       fprintf(fd, " set: real nrev = %e;\n", nrev);
%!       fprintf(fd, " set: real n0 = %e;\n", n0);
%!       fprintf(fd, " set: real tolres = %e;\n", tolres);
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "         problem: initial value;\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "         initial time: 0;\n");
%!       fputs(fd, "         final time: 2. * pi * nrev / abs(omega);\n");
%!       fputs(fd, "         time step: 2. * pi / (nsteps * abs(omega));\n");
%!       fputs(fd, "         linear solver: naive, colamd, scale, row max column max, always;\n");
%!       fputs(fd, "         method: ms, 0.6;\n");
%!       fputs(fd, "         max iterations: 10;\n");
%!       fputs(fd, "         tolerance: tolres;\n");
%!       fputs(fd, "         threads: assembly, 1;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: auto;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "     max iterations: 0;\n");
%!       fputs(fd, "     structural nodes: 1;\n");
%!       fputs(fd, "     joints: 1;\n");
%!       fputs(fd, "     rigid bodies: 1;\n");
%!       fputs(fd, "     gravity;\n");
%!       fputs(fd, "     output meter: closest next, 2 * pi * n0 / abs(omega), forever, 2 * pi / (100 * abs(omega));\n");
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fputs(fd, "         structural: 1, dynamic, 0., 0., h3, eye, reference, global, h3 * omega2, -h3 * omega1, 0., reference, global, omega1, omega2, omega3, accelerations, yes;\n");
%!       fputs(fd, " end: nodes;\n");
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "         body: 1, 1, m, reference, node, s1, s2, s3, diag, J1, J2, J3;\n");
%!       fputs(fd, "         joint: 1, total pin joint, 1,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position constraint, active, active, active, null,\n");
%!       fputs(fd, "            orientation constraint, inactive, inactive, inactive, null;\n");
%!       fputs(fd, "            gravity: uniform, 0., 0., -1., 9.81;\n");
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!     end_unwind_protect
%!     options.output_file = fname1;
%!     options.verbose = false;
%!     if (silent)
%!       options.logfile = [fname1, ".stdout"];
%!     endif
%!     mbdyn_solver_run(fname1, options);
%!     [t1, trajectory1, deformation1, velocity1, acceleration1, node_id1] = mbdyn_post_load_output_struct(fname1);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname1);
%!       files = dir([fname1, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
%!     unwind_protect
%!       [fd, fname2] = mkstemp(fullfile(tempdir(), "mbdyn_pre_write_fem_data2_XXXXXX"));
%!       if (fd == -1)
%!         error("failed to open temporary file");
%!       endif
%!       fprintf(fd, " set: real omega1 = %e;\n", omega1(i));
%!       fprintf(fd, " set: real omega2 = %e;\n", omega2(i));
%!       fprintf(fd, " set: real omega3 = %e;\n", omega3(i));
%!       fputs(fd, " set: real omega = sqrt(omega1^2 + omega2^2 + omega3^2);\n");
%!       fprintf(fd, " set: real m = %e;\n", m);
%!       fprintf(fd, " set: real s1 = %e;\n", s1);
%!       fprintf(fd, " set: real s2 = %e;\n", s2);
%!       fprintf(fd, " set: real s3 = %e;\n", s3);
%!       fprintf(fd, " set: real J1 = %e;\n", J1);
%!       fprintf(fd, " set: real J2 = %e;\n", J2);
%!       fprintf(fd, " set: real J3 = %e;\n", J3);
%!       fprintf(fd, " set: real h3 = %e;\n", h3);
%!       fprintf(fd, " set: real nsteps = %e;\n", nsteps);
%!       fprintf(fd, " set: real nrev = %e;\n", nrev);
%!       fprintf(fd, " set: real n0 = %e;\n", n0);
%!       fprintf(fd, " set: real tolres = %e;\n", tolres);
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "         problem: initial value;\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "         initial time: 0;\n");
%!       fputs(fd, "         final time: 2. * pi * nrev / abs(omega);\n");
%!       fputs(fd, "         time step: 2. * pi / (nsteps * abs(omega));\n");
%!       fputs(fd, "         linear solver: naive, colamd, scale, row max column max, always;\n");
%!       fputs(fd, "         method: ms, 0.6;\n");
%!       fputs(fd, "         max iterations: 10;\n");
%!       fputs(fd, "         tolerance: tolres;\n");
%!       fputs(fd, "         threads: assembly, 1;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: auto;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "     max iterations: 0;\n");
%!       fputs(fd, "     structural nodes: 1;\n");
%!       fputs(fd, "     joints: 2;\n");
%!       fputs(fd, "     rigid bodies: 1;\n");
%!       fputs(fd, "     gravity;\n");
%!       fputs(fd, "     output meter: closest next, 2 * pi * n0 / abs(omega), forever, 2 * pi / (100 * abs(omega));\n");
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fputs(fd, "         structural: 1, modal, 0., 0., h3, eye, reference, global, h3 * omega2, -h3 * omega1, 0., reference, global, omega1, omega2, omega3, accelerations, yes;\n");
%!       fputs(fd, " end: nodes;\n");
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "         body: 1, 1, m, reference, node, s1, s2, s3, diag, J1, J2, J3;\n");
%!       fputs(fd, "         joint: 1, total pin joint, 1,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position constraint, active, active, active, null,\n");
%!       fputs(fd, "            orientation constraint, inactive, inactive, inactive, null;\n");
%!       fputs(fd, "         joint: 2, modal, 1,\n");
%!       fputs(fd, "            1, from file,\n");
%!       fputs(fd, "            damping from file,\n");
%!       fprintf(fd, "          \"%s.fem\",\n", fname2);
%!       fputs(fd, "            origin node, 1000,\n");
%!       fputs(fd, "            0;\n");
%!       fputs(fd, "            gravity: uniform, 0., 0., -1., 9.81;\n");
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!     end_unwind_protect
%!     MRED = 1;
%!     DRED = 0;
%!     SRED = 1e10;
%!     TRED = zeros(6, 1);
%!     X0 = zeros(3, 1);
%!     URED0 = 0;
%!     DURED_DT0 = 0;
%!     DIAGM = [];
%!     M = 0;
%!     XGC = zeros(3, 1);
%!     JGC = zeros(3, 3);
%!     NODE_LIST = 1000;
%!     mbdyn_pre_write_fem_data([fname2, ".fem"], MRED, DRED, SRED, TRED, X0, URED0, DURED_DT0, DIAGM, M, XGC, JGC, NODE_LIST);
%!     options.output_file = fname2;
%!     if (silent)
%!       options.logfile = [fname2, ".stdout"];
%!     endif
%!     mbdyn_solver_run(fname2, options);
%!     [t2, trajectory2, deformation2, velocity2, acceleration2, node_id2] = mbdyn_post_load_output_struct(fname2);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname2);
%!       files = dir([fname2, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
%!   fd = -1;
%!   unwind_protect
%!     unwind_protect
%!       [fd, fname3] = mkstemp(fullfile(tempdir(), "mbdyn_pre_write_fem_data3_XXXXXX"));
%!       if (fd == -1)
%!         error("failed to open temporary file");
%!       endif
%!       fprintf(fd, " set: real omega1 = %e;\n", omega1(i));
%!       fprintf(fd, " set: real omega2 = %e;\n", omega2(i));
%!       fprintf(fd, " set: real omega3 = %e;\n", omega3(i));
%!       fputs(fd, " set: real omega = sqrt(omega1^2 + omega2^2 + omega3^2);\n");
%!       fprintf(fd, " set: real m = %e;\n", m);
%!       fprintf(fd, " set: real s1 = %e;\n", s1);
%!       fprintf(fd, " set: real s2 = %e;\n", s2);
%!       fprintf(fd, " set: real s3 = %e;\n", s3);
%!       fprintf(fd, " set: real J1 = %e;\n", J1);
%!       fprintf(fd, " set: real J2 = %e;\n", J2);
%!       fprintf(fd, " set: real J3 = %e;\n", J3);
%!       fprintf(fd, " set: real h3 = %e;\n", h3);
%!       fprintf(fd, " set: real nsteps = %e;\n", nsteps);
%!       fprintf(fd, " set: real nrev = %e;\n", nrev);
%!       fprintf(fd, " set: real n0 = %e;\n", n0);
%!       fprintf(fd, " set: real tolres = %e;\n", tolres);
%!       fputs(fd, " begin: data;\n");
%!       fputs(fd, "         problem: initial value;\n");
%!       fputs(fd, " end: data;\n");
%!       fputs(fd, " begin: initial value;\n");
%!       fputs(fd, "         initial time: 0;\n");
%!       fputs(fd, "         final time: 2. * pi * nrev / abs(omega);\n");
%!       fputs(fd, "         time step: 2. * pi / (nsteps * abs(omega));\n");
%!       fputs(fd, "         linear solver: naive, colamd, scale, row max column max, always;\n");
%!       fputs(fd, "         method: ms, 0.6;\n");
%!       fputs(fd, "         max iterations: 10;\n");
%!       fputs(fd, "         tolerance: tolres;\n");
%!       fputs(fd, "         threads: assembly, 1;\n");
%!       fputs(fd, "         derivatives max iterations: 10;\n");
%!       fputs(fd, "         derivatives coefficient: auto;\n");
%!       fputs(fd, "         #output: iterations, residual;\n");
%!       fputs(fd, "         #abort after: derivatives;\n");
%!       fputs(fd, " end: initial value;\n");
%!       fputs(fd, " begin: control data;\n");
%!       fputs(fd, "     max iterations: 0;\n");
%!       fputs(fd, "     structural nodes: 1;\n");
%!       fputs(fd, "     joints: 2;\n");
%!       fputs(fd, "     gravity;\n");
%!       fputs(fd, "     output meter: closest next, 2 * pi * n0 / abs(omega), forever, 2 * pi / (100 * abs(omega));\n");
%!       fputs(fd, " end: control data;\n");
%!       fputs(fd, " begin: nodes;\n");
%!       fputs(fd, "         structural: 1, modal, 0., 0., h3, eye, reference, global, h3 * omega2, -h3 * omega1, 0., reference, global, omega1, omega2, omega3, accelerations, yes;\n");
%!       fputs(fd, " end: nodes;\n");
%!       fputs(fd, " begin: elements;\n");
%!       fputs(fd, "         joint: 1, total pin joint, 1,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position, reference, global, null,\n");
%!       fputs(fd, "            position orientation, eye,\n");
%!       fputs(fd, "            rotation orientation, eye,\n");
%!       fputs(fd, "            position constraint, active, active, active, null,\n");
%!       fputs(fd, "            orientation constraint, inactive, inactive, inactive, null;\n");
%!       fputs(fd, "         joint: 2, modal, 1,\n");
%!       fputs(fd, "            1, from file,\n");
%!       fputs(fd, "            damping from file,\n");
%!       fprintf(fd, "          \"%s.fem\",\n", fname3);
%!       fputs(fd, "            origin node, 1000,\n");
%!       fputs(fd, "            0;\n");
%!       fputs(fd, "            gravity: uniform, 0., 0., -1., 9.81;\n");
%!       fputs(fd, " end: elements;\n");
%!     unwind_protect_cleanup
%!       if (fd ~= -1)
%!         fclose(fd);
%!       endif
%!     end_unwind_protect
%!     MRED = 1;
%!     DRED = 0;
%!     SRED = 1e10;
%!     TRED = zeros(6, 1);
%!     X0 = zeros(3, 1);
%!     URED0 = 0;
%!     DURED_DT0 = 0;
%!     DIAGM = [];
%!     M = m;
%!     XGC = [s1; s2; s3];
%!     JGC = diag([J1, J2, J3]);
%!     NODE_LIST = 1000;
%!     mbdyn_pre_write_fem_data([fname3, ".fem"], MRED, DRED, SRED, TRED, X0, URED0, DURED_DT0, DIAGM, M, XGC, JGC, NODE_LIST);
%!     options.output_file = fname3;
%!     if (silent)
%!       options.logfile = [fname3, ".stdout"];
%!     endif
%!     mbdyn_solver_run(fname3, options);
%!     [t3, trajectory3, deformation3, velocity3, acceleration3, node_id3] = mbdyn_post_load_output_struct(fname3);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       unlink(fname3);
%!       files = dir([fname3, "*"]);
%!       for j=1:numel(files)
%!         unlink(fullfile(files(j).folder, files(j).name));
%!       endfor
%!     endif
%!   end_unwind_protect
%!   tol = 1e-4;
%!   for j=1:numel(node_id2)
%!     assert(trajectory2{j}, trajectory1{j}, tol * max(max(abs(trajectory1{j}))));
%!     assert(velocity2{j}, velocity1{j}, tol * max(max(abs(velocity1{j}))));
%!     assert(acceleration2{j}, acceleration1{j}, tol * max(max(abs(acceleration1{j}))));
%!     assert(trajectory3{j}, trajectory1{j}, tol * max(max(abs(trajectory1{j}))));
%!     assert(velocity3{j}, velocity1{j}, tol * max(max(abs(velocity1{j}))));
%!     assert(acceleration3{j}, acceleration1{j}, tol * max(max(abs(acceleration1{j}))));
%!   endfor
%! endfor
