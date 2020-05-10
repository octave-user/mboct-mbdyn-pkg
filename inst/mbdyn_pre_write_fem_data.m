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
      fprintf(fd, "%e %e %e %e %e %e\n", Tred(:, i));
    endfor

    fprintf(fd, "**\n");
    fprintf(fd, "** RECORD GROUP 9, MODAL MASS MATRIX\n");

    for i=1:rows(Mred)
      fprintf(fd, "%e ", Mred(i, :));
      fprintf(fd, "\n");
    endfor

    fprintf(fd, "**\n");
    fprintf(fd, "** RECORD GROUP 10, MODAL STIFFNESS MATRIX\n");

    for i=1:rows(Sred)
      fprintf(fd, "%e ", Sred(i, :));
      fprintf(fd, "\n");
    endfor

    fprintf(fd, "**\n");

    if (nargin >= 8 && numel(diagM) > 0)
      fprintf(fd, "** RECORD GROUP 11, DIAGONAL OF LUMPED MASS MATRIX\n");
      fprintf(fd, "%e %e %e %e %e %e\n", diagM);
      fprintf(fd, "**\n");
    endif

    if (nargin >= 11 && numel(Jgc) > 0)
      fprintf(fd, "** RECORD GROUP 12, RIGID BODY INERTIA MATRIX\n");
      fprintf(fd, "%e\n", m);
      fprintf(fd, "%e %e %e\n", Xgc);
      fprintf(fd, "%e %e %e\n", Jgc.');
      fprintf(fd, "**\n");
    endif

    if (numel(Dred) > 0)
      fprintf(fd, "**\n");
      fprintf(fd, "** RECORD GROUP 13, MODAL DAMPING MATRIX\n");

      for i=1:rows(Dred)
        fprintf(fd, "%e ", Dred(i, :));
        fprintf(fd, "\n");
      endfor

      fprintf(fd, "**\n");
    endif

    if (nargin >= 14 && numel(Inv3) > 0)
      fprintf(fd, "**\n");
      fprintf(fd, "** RECORD GROUP 14, INVARIANT 3\n");

      for i=1:rows(Inv3)
        fprintf(fd, "%e ", Inv3(i,:));
        fprintf(fd, "\n");
      endfor

      fprintf(fd, "**\n");
    endif

    if (nargin >= 15 && numel(Inv4) > 0)
      fprintf(fd, "**\n");
      fprintf(fd, "** RECORD GROUP 15, INVARIANT 4\n");

      for i=1:rows(Inv4)
        fprintf(fd, "%e ", Inv4(i, :));
        fprintf(fd, "\n");
      endfor

      fprintf(fd,"**\n");
    endif

    if (nargin >= 17 && numel(Inv8) > 0)
      fprintf(fd, "**\n");
      fprintf(fd, "** RECORD GROUP 16, INVARIANT 8\n");

      for i=1:rows(Inv8)
        for j=1:size(Inv8, 3)
          fprintf(fd, "%e ", Inv8(i, :,  j));
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
          fprintf(fd, "%e ", Inv5(i, :, j));
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
            fprintf(fd, "%e ", Inv9(i, :, j, k));
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
