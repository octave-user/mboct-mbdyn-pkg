## Copyright (C) 2016 Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} mbdyn_post_print_eig(@var{modal}, @var{log_dat}, @var{output_file})
##
## Print all eigenvalues from an eigenanalysis to a file.
##
## @var{modal} @dots{} Return value of mbdyn_post_load_output_eig.
##
## @var{log_dat} @dots{} Return value of mbdyn_post_load_log.
##
## @var{output_file} @dots{} File name for results.
##
## @end deftypefn

function mbdyn_post_print_eig(modal, log_dat, varargin)
  fd = -1;
  owns_fd = false;

  if (length(varargin) >= 1)
    if (ischar(varargin{1}))
      owns_fd = true;
      [fd, msg] = fopen(varargin{1}, "wt");
      if (fd == -1)
        error("failed to open file \"%s\": %s", varargin{1}, msg);
      endif
    else
      fd = varargin{1};
    endif
  else
    fd = stdout;
  endif

  unwind_protect
    if (isfield(log_dat.dof_info, "struct_node_labels"))
      for j=1:columns(modal.VR)
        fprintf(fd, "f=%.3fHz\t", modal.f(j));
        for i=1:length(log_dat.dof_info.struct_node_labels)
          fprintf(fd, "node %d\t", log_dat.dof_info.struct_node_labels(i));
        endfor
        fprintf(fd, "\n");

        norm_VR = max(abs(modal.VR(:, j)));

        for k=1:6
          switch (k)
            case {1, 2, 3}
              fprintf(fd, "X%d\t", k);
            case {4, 5, 6}
              fprintf(fd, "g%d\t", k - 3);
          endswitch
          for i=1:length(log_dat.dof_info.struct_node_dofs)
            fprintf(fd, "%.6f\t", abs(modal.VR(log_dat.dof_info.struct_node_dofs(i) + k, j)) / norm_VR);
          endfor
          fprintf(fd, "\n");
        endfor
      endfor
    endif
  unwind_protect_cleanup
    if (fd ~= -1 && owns_fd)
      fclose(fd);
    endif
  end_unwind_protect
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     s1 = 5000;
%!     m1 = 2.5;
%!     d1 = 10;
%!     s2 = 1000;
%!     m2 = 3;
%!     d2 = 20;
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_print_eig_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fprintf(fd, "set: real s1 = %g;\n", s1);
%!     fprintf(fd, "set: real m1 = %g;\n", m1);
%!     fprintf(fd, "set: real d1 = %g;\n", d1);
%!     fprintf(fd, "set: real s2 = %g;\n", s2);
%!     fprintf(fd, "set: real m2 = %g;\n", m2);
%!     fprintf(fd, "set: real d2 = %g;\n", d2);
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
%!     fputs(fd, "         eigenanalysis: 0,\n");
%!     fputs(fd, "          suffix format, \"%02d\",\n");
%!     fputs(fd, "          output full matrices,\n");
%!     fputs(fd, "          output eigenvectors,\n");
%!     fputs(fd, "          results output precision, 16,\n");
%!     fputs(fd, "          parameter, 0.01, use lapack, balance, permute;\n");
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
%!     fputs(fd, "     genels: 4;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static, null, eye, null, null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         genel: 1, spring support, 1, structural, 1, algebraic, linear viscoelastic generic, s1, d1;\n");
%!     fputs(fd, "         genel: 2, mass, 1, structural, 1, algebraic, m1;\n");
%!     fputs(fd, "         genel: 3, spring support, 1, structural, 2, algebraic, linear viscoelastic generic, s2, d2;\n");
%!     fputs(fd, "         genel: 4, mass, 1, structural, 2, algebraic, m2;\n");
%!     fputs(fd, "         joint: 1, total pin joint, 1,\n");
%!     fputs(fd, "            position, null,\n");
%!     fputs(fd, "            position orientation, eye,\n");
%!     fputs(fd, "            rotation orientation, eye,\n");
%!     fputs(fd, "            position, null,\n");
%!     fputs(fd, "            position orientation, eye,\n");
%!     fputs(fd, "            rotation orientation, eye,\n");
%!     fputs(fd, "            position constraint, inactive, inactive, active, null,\n");
%!     fputs(fd, "            orientation constraint, active, active, active, null;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.positive_frequencies = false;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   modal = mbdyn_post_load_output_eig(fname, options, 0);
%!   log_dat = mbdyn_post_load_log(fname);
%!   omega1 = sqrt(s1 / m1 - (d1 / (2 * m1))^2);
%!   omega2 = sqrt(s2 / m2 - (d2 / (2 * m2))^2);
%!   alpha1 = -d1 / (2 * m1);
%!   alpha2 = -d2 / (2 * m2);
%!   lambda1 = alpha1 + [-1j; 1j] * omega1;
%!   lambda2 = alpha2 + [-1j; 1j] * omega2;
%!   lambda = [lambda1; lambda2];
%!   [dummy, idx] = sort(imag(lambda), "ascend");
%!   lambda = lambda(idx);
%!   assert_simple(modal.lambda, lambda, eps^0.9 * max(abs(lambda)));
%!   mbdyn_post_print_eig(modal, log_dat, [fname, "_res.txt"]);
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
%!     s1 = 5000;
%!     m1 = 2.5;
%!     d1 = 10;
%!     s2 = 1000;
%!     m2 = 3;
%!     d2 = 20;
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_print_eig_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fprintf(fd, "set: real s1 = %g;\n", s1);
%!     fprintf(fd, "set: real m1 = %g;\n", m1);
%!     fprintf(fd, "set: real d1 = %g;\n", d1);
%!     fprintf(fd, "set: real s2 = %g;\n", s2);
%!     fprintf(fd, "set: real m2 = %g;\n", m2);
%!     fprintf(fd, "set: real d2 = %g;\n", d2);
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
%!     fputs(fd, "         eigenanalysis: 0,\n");
%!     fputs(fd, "          suffix format, \"%02d\",\n");
%!     fputs(fd, "          output full matrices,\n");
%!     fputs(fd, "          output eigenvectors,\n");
%!     fputs(fd, "          results output precision, 16,\n");
%!     fputs(fd, "          parameter, 0.01, use lapack, balance, permute;\n");
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
%!     fputs(fd, "     genels: 4;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static, null, eye, null, null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         genel: 1, spring support, 1, structural, 1, algebraic, linear viscoelastic generic, s1, d1;\n");
%!     fputs(fd, "         genel: 2, mass, 1, structural, 1, algebraic, m1;\n");
%!     fputs(fd, "         genel: 3, spring support, 1, structural, 2, algebraic, linear viscoelastic generic, s2, d2;\n");
%!     fputs(fd, "         genel: 4, mass, 1, structural, 2, algebraic, m2;\n");
%!     fputs(fd, "         joint: 1, total pin joint, 1,\n");
%!     fputs(fd, "            position, null,\n");
%!     fputs(fd, "            position orientation, eye,\n");
%!     fputs(fd, "            rotation orientation, eye,\n");
%!     fputs(fd, "            position, null,\n");
%!     fputs(fd, "            position orientation, eye,\n");
%!     fputs(fd, "            rotation orientation, eye,\n");
%!     fputs(fd, "            position constraint, inactive, inactive, active, null,\n");
%!     fputs(fd, "            orientation constraint, active, active, active, null;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.positive_frequencies = false;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   modal = mbdyn_post_load_output_eig(fname, options, 0);
%!   log_dat = mbdyn_post_load_log(fname);
%!   omega1 = sqrt(s1 / m1 - (d1 / (2 * m1))^2);
%!   omega2 = sqrt(s2 / m2 - (d2 / (2 * m2))^2);
%!   alpha1 = -d1 / (2 * m1);
%!   alpha2 = -d2 / (2 * m2);
%!   lambda1 = alpha1 + [-1j; 1j] * omega1;
%!   lambda2 = alpha2 + [-1j; 1j] * omega2;
%!   lambda = [lambda1; lambda2];
%!   [dummy, idx] = sort(imag(lambda), "ascend");
%!   lambda = lambda(idx);
%!   assert_simple(modal.lambda, lambda, eps^0.9 * max(abs(lambda)));
%!   mbdyn_post_print_eig(modal, log_dat, [fname, "_res.txt"]);
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
%!     s1 = 5000;
%!     m1 = 2.5;
%!     d1 = 10;
%!     s2 = 1000;
%!     m2 = 3;
%!     d2 = 20;
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_print_eig_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fprintf(fd, "set: real s1 = %g;\n", s1);
%!     fprintf(fd, "set: real m1 = %g;\n", m1);
%!     fprintf(fd, "set: real d1 = %g;\n", d1);
%!     fprintf(fd, "set: real s2 = %g;\n", s2);
%!     fprintf(fd, "set: real m2 = %g;\n", m2);
%!     fprintf(fd, "set: real d2 = %g;\n", d2);
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
%!     fputs(fd, "         eigenanalysis: 0,\n");
%!     fputs(fd, "          suffix format, \"%02d\",\n");
%!     fputs(fd, "          output full matrices,\n");
%!     fputs(fd, "          output eigenvectors,\n");
%!     fputs(fd, "          results output precision, 16,\n");
%!     fputs(fd, "          parameter, 0.01, use lapack, balance, permute;\n");
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     genels: 4;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static, null, eye, null, null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         genel: 1, spring support, 1, structural, 1, algebraic, linear viscoelastic generic, s1, d1;\n");
%!     fputs(fd, "         genel: 2, mass, 1, structural, 1, algebraic, m1;\n");
%!     fputs(fd, "         genel: 3, spring support, 1, structural, 2, algebraic, linear viscoelastic generic, s2, d2;\n");
%!     fputs(fd, "         genel: 4, mass, 1, structural, 2, algebraic, m2;\n");
%!     fputs(fd, "         joint: 1, total pin joint, 1,\n");
%!     fputs(fd, "            position, null,\n");
%!     fputs(fd, "            position orientation, eye,\n");
%!     fputs(fd, "            rotation orientation, eye,\n");
%!     fputs(fd, "            position, null,\n");
%!     fputs(fd, "            position orientation, eye,\n");
%!     fputs(fd, "            rotation orientation, eye,\n");
%!     fputs(fd, "            position constraint, inactive, inactive, active, null,\n");
%!     fputs(fd, "            orientation constraint, active, active, active, null;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.positive_frequencies = false;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   modal = mbdyn_post_load_output_eig(fname, options, 0);
%!   log_dat = mbdyn_post_load_log(fname);
%!   omega1 = sqrt(s1 / m1 - (d1 / (2 * m1))^2);
%!   omega2 = sqrt(s2 / m2 - (d2 / (2 * m2))^2);
%!   alpha1 = -d1 / (2 * m1);
%!   alpha2 = -d2 / (2 * m2);
%!   lambda1 = alpha1 + [-1j; 1j] * omega1;
%!   lambda2 = alpha2 + [-1j; 1j] * omega2;
%!   lambda = [lambda1; lambda2];
%!   [dummy, idx] = sort(imag(lambda), "ascend");
%!   lambda = lambda(idx);
%!   assert_simple(modal.lambda, lambda, eps^0.9 * max(abs(lambda)));
%!   mbdyn_post_print_eig(modal, log_dat, [fname, "_res.txt"]);
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
%!     s1 = 5000;
%!     m1 = 2.5;
%!     d1 = 10;
%!     s2 = 1000;
%!     m2 = 3;
%!     d2 = 20;
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_print_eig_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fprintf(fd, "set: real s1 = %g;\n", s1);
%!     fprintf(fd, "set: real m1 = %g;\n", m1);
%!     fprintf(fd, "set: real d1 = %g;\n", d1);
%!     fprintf(fd, "set: real s2 = %g;\n", s2);
%!     fprintf(fd, "set: real m2 = %g;\n", m2);
%!     fprintf(fd, "set: real d2 = %g;\n", d2);
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
%!     fputs(fd, "         eigenanalysis: 0,\n");
%!     fputs(fd, "          suffix format, \"%02d\",\n");
%!     fputs(fd, "          output full matrices,\n");
%!     fputs(fd, "          output eigenvectors,\n");
%!     fputs(fd, "          results output precision, 16,\n");
%!     fputs(fd, "          parameter, 0.01, use lapack, balance, permute;\n");
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
%!     fputs(fd, "     genels: 4;\n");
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, static, null, eye, null, null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         genel: 1, spring support, 1, structural, 1, algebraic, linear viscoelastic generic, s1, d1;\n");
%!     fputs(fd, "         genel: 2, mass, 1, structural, 1, algebraic, m1;\n");
%!     fputs(fd, "         genel: 3, spring support, 1, structural, 2, algebraic, linear viscoelastic generic, s2, d2;\n");
%!     fputs(fd, "         genel: 4, mass, 1, structural, 2, algebraic, m2;\n");
%!     fputs(fd, "         joint: 1, total pin joint, 1,\n");
%!     fputs(fd, "            position, null,\n");
%!     fputs(fd, "            position orientation, eye,\n");
%!     fputs(fd, "            rotation orientation, eye,\n");
%!     fputs(fd, "            position, null,\n");
%!     fputs(fd, "            position orientation, eye,\n");
%!     fputs(fd, "            rotation orientation, eye,\n");
%!     fputs(fd, "            position constraint, inactive, inactive, active, null,\n");
%!     fputs(fd, "            orientation constraint, active, active, active, null;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.positive_frequencies = false;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   modal = mbdyn_post_load_output_eig(fname, options, 0);
%!   log_dat = mbdyn_post_load_log(fname);
%!   omega1 = sqrt(s1 / m1 - (d1 / (2 * m1))^2);
%!   omega2 = sqrt(s2 / m2 - (d2 / (2 * m2))^2);
%!   alpha1 = -d1 / (2 * m1);
%!   alpha2 = -d2 / (2 * m2);
%!   lambda1 = alpha1 + [-1j; 1j] * omega1;
%!   lambda2 = alpha2 + [-1j; 1j] * omega2;
%!   lambda = [lambda1; lambda2];
%!   [dummy, idx] = sort(imag(lambda), "ascend");
%!   lambda = lambda(idx);
%!   assert_simple(modal.lambda, lambda, eps^0.9 * max(abs(lambda)));
%!   mbdyn_post_print_eig(modal, log_dat, stdout);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
