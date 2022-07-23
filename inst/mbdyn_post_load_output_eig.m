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
## @deftypefn {Function File} @var{modal} = mbdyn_post_load_output_eig(@var{mbdyn_output_file})
## @deftypefnx {} @dots{} = mbdyn_post_load_output_eig(@var{mbdyn_output_file}, @var{options})
## @deftypefnx {} @dots{} = mbdyn_post_load_output_eig(@var{mbdyn_output_file}, @var{options}, @var{index})
##
## Loads data form an eigenanalysis from MBDyn output file "<@var{mbdyn_output_file}_%02d.m>".
##
## @var{mbdyn_output_file} @dots{} Name of the MBDyn output file without extension.
##
## @var{options}.positive_frequencies @dots{} Return only eigenvalues with imaginary values greater than zero.
##
## @var{options}.solve_qz @dots{} Solve the eigenvalue problem within Octave.
##
## @var{index} @dots{} If there are more than one eigenanalysis results (e.g. "mod_01.m", "mod_02.m", @dots{}) then load the file corresponding to the index <@var{index}>.
##
## @var{modal} @dots{} Structure that contains mode shapes, eigenvalues and node positions.
##
## @end deftypefn

function modal = mbdyn_post_load_output_eig(mbdyn_output_file, options, index)
  if (nargin < 1 || nargin > 3 || nargout > 1)
    print_usage();
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (nargin < 3)
    index = 0;
  endif
  
  if (~isfield(options, "positive_frequencies"))
    options.positive_frequencies = true;
  endif

  if (~isfield(options, "solve_qz"))
    options.solve_qz = false;
  endif
  
  [inp_dir, inp_name, inp_ext] = fileparts(mbdyn_output_file);

  mbdyn_output_file = fullfile(inp_dir, sprintf("%s_%02d.m", inp_name, index));

  [info, err, msg] = stat(mbdyn_output_file);

  if (err ~= 0)
    mbdyn_output_file = fullfile(inp_dir, cstrcat(inp_name, ".m"));
    [info, err, msg] = stat(mbdyn_output_file);
  endif

  if (err ~= 0)
    error("file not found \"%s\"", mbdyn_output_file);
  endif
  
  source(mbdyn_output_file);
  
  if (1 == exist("dCoef", "var"))
    modal.dCoef = dCoef;
  endif
  
  if (1 == exist("Aplus", "var"))
    modal.Aplus = Aplus;
  endif
  
  if (1 == exist("Aminus", "var"))
    modal.Aminus = Aminus;
  endif

  if ((options.solve_qz || ...
       1 ~= exist("alpha", "var")) && ...
      1 == exist("Aplus", "var") && ...
      1 == exist("Aminus", "var"))
    clear alpha VR VL;
    [VR, alpha2] = eig(Aminus, Aplus);
    LAMBDA2 = -diag(alpha2);
    modal.lambda = (LAMBDA2 + 1) ./ (dCoef * (LAMBDA2 - 1));
    modal.f = imag(modal.lambda) / (2 * pi);  
  elseif (1 == exist("alpha", "var"))
    LAMBDA = (alpha(:,1) + 1j * alpha(:, 2)) ./ alpha(:, 3);
    modal.lambda = 1 / dCoef * (LAMBDA - 1) ./ (LAMBDA + 1);
    modal.f = imag(modal.lambda) / (2 * pi);
  else
    modal.lambda = [];
    modal.f = [];
  endif
  
  if (length(modal.f) > 0)
    [modal.f, idx_f] = sort(modal.f);
    if (options.positive_frequencies)
      idx_gtz = find(modal.f > 0);
      modal.f = modal.f(idx_gtz);
      idx_f = idx_f(idx_gtz);
    endif
    modal.lambda = modal.lambda(idx_f);
  endif

  if (isfield(modal, "lambda"))
    delta = -real(modal.lambda);
    omegad = imag(modal.lambda);
    a0 = (delta ./ omegad).^2;
    modal.D = sqrt(a0 ./ (1 + a0));
    modal.omega0 = omegad ./ sqrt(1 - modal.D.^2);
  endif
  
  if (1 == exist("dTime", "var"))
    modal.dTime = dTime;
  endif
  
  if (1 == exist("X0", "var"))
    modal.X0 = X0;
  endif
  
  if (1 == exist("idx", "var"))
    modal.idx = idx;
  endif

  if (1 == exist("labels", "var"))
    modal.labels = labels;
  endif
  
  if (1 == exist("VR", "var"))
    modal.VR = VR(:, idx_f);
  else
    modal.VR = [];
  endif
  
  if (1 == exist("VL", "var"))
    modal.VL = VL(:, idx_f);
  else
    modal.VL = [];
  endif  
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_eig_XXXXXX"));
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
%!     fputs(fd, "     abstract nodes: 2;\n");
%!     fputs(fd, "     genels: 4;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         abstract: 1, differential;\n");
%!     fputs(fd, "         abstract: 2, differential;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         genel: 1, spring support, 1, abstract, algebraic, linear viscoelastic generic, s1, d1;\n");
%!     fputs(fd, "         genel: 2, mass, 1, abstract, algebraic, m1;\n");
%!     fputs(fd, "         genel: 3, spring support, 2, abstract, algebraic, linear viscoelastic generic, s2, d2;\n");
%!     fputs(fd, "         genel: 4, mass, 2, abstract, algebraic, m2;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   options.f_run_mbdyn2easyanim = false;
%!   options.positive_frequencies = false;
%!   mbdyn_solver_run(fname, options);
%!   modal = mbdyn_post_load_output_eig(fname, options, 0);
%!   omega1 = sqrt(s1 / m1 - (d1 / (2 * m1))^2);
%!   omega2 = sqrt(s2 / m2 - (d2 / (2 * m2))^2);
%!   alpha1 = -d1 / (2 * m1);
%!   alpha2 = -d2 / (2 * m2);
%!   lambda1 = alpha1 + [-1j; 1j] * omega1;
%!   lambda2 = alpha2 + [-1j; 1j] * omega2;
%!   lambda = [lambda1; lambda2];
%!   [dummy, idx] = sort(imag(lambda), "ascend");
%!   lambda = lambda(idx);
%!   assert(modal.lambda, lambda, eps^0.9 * max(abs(lambda)));
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_eig_XXXXXX"));
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
%!     fputs(fd, "     abstract nodes: 2;\n");
%!     fputs(fd, "     genels: 4;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         abstract: 1, differential;\n");
%!     fputs(fd, "         abstract: 2, differential;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         genel: 1, spring support, 1, abstract, algebraic, linear viscoelastic generic, s1, d1;\n");
%!     fputs(fd, "         genel: 2, mass, 1, abstract, algebraic, m1;\n");
%!     fputs(fd, "         genel: 3, spring support, 2, abstract, algebraic, linear viscoelastic generic, s2, d2;\n");
%!     fputs(fd, "         genel: 4, mass, 2, abstract, algebraic, m2;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   options.f_run_mbdyn2easyanim = false;
%!   options.positive_frequencies = false;
%!   mbdyn_solver_run(fname, options);
%!   modal = mbdyn_post_load_output_eig(fname, options, 0);
%!   omega1 = sqrt(s1 / m1 - (d1 / (2 * m1))^2);
%!   omega2 = sqrt(s2 / m2 - (d2 / (2 * m2))^2);
%!   alpha1 = -d1 / (2 * m1);
%!   alpha2 = -d2 / (2 * m2);
%!   lambda1 = alpha1 + [-1j; 1j] * omega1;
%!   lambda2 = alpha2 + [-1j; 1j] * omega2;
%!   lambda = [lambda1; lambda2];
%!   [dummy, idx] = sort(imag(lambda), "ascend");
%!   lambda = lambda(idx);
%!   assert(modal.lambda, lambda, eps^0.9 * max(abs(lambda)));
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_eig_XXXXXX"));
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
%!     fputs(fd, "     abstract nodes: 2;\n");
%!     fputs(fd, "     genels: 4;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         abstract: 1, differential;\n");
%!     fputs(fd, "         abstract: 2, differential;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         genel: 1, spring support, 1, abstract, algebraic, linear viscoelastic generic, s1, d1;\n");
%!     fputs(fd, "         genel: 2, mass, 1, abstract, algebraic, m1;\n");
%!     fputs(fd, "         genel: 3, spring support, 2, abstract, algebraic, linear viscoelastic generic, s2, d2;\n");
%!     fputs(fd, "         genel: 4, mass, 2, abstract, algebraic, m2;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   options.f_run_mbdyn2easyanim = false;
%!   options.positive_frequencies = false;
%!   mbdyn_solver_run(fname, options);
%!   modal = mbdyn_post_load_output_eig(fname, options, 0);
%!   omega1 = sqrt(s1 / m1 - (d1 / (2 * m1))^2);
%!   omega2 = sqrt(s2 / m2 - (d2 / (2 * m2))^2);
%!   alpha1 = -d1 / (2 * m1);
%!   alpha2 = -d2 / (2 * m2);
%!   lambda1 = alpha1 + [-1j; 1j] * omega1;
%!   lambda2 = alpha2 + [-1j; 1j] * omega2;
%!   lambda = [lambda1; lambda2];
%!   [dummy, idx] = sort(imag(lambda), "ascend");
%!   lambda = lambda(idx);
%!   assert(modal.lambda, lambda, eps^0.9 * max(abs(lambda)));
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_output_eig_XXXXXX"));
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
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");

%!     fputs(fd, "     abstract nodes: 2;\n");
%!     fputs(fd, "     genels: 4;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         abstract: 1, differential;\n");
%!     fputs(fd, "         abstract: 2, differential;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         genel: 1, spring support, 1, abstract, algebraic, linear viscoelastic generic, s1, d1;\n");
%!     fputs(fd, "         genel: 2, mass, 1, abstract, algebraic, m1;\n");
%!     fputs(fd, "         genel: 3, spring support, 2, abstract, algebraic, linear viscoelastic generic, s2, d2;\n");
%!     fputs(fd, "         genel: 4, mass, 2, abstract, algebraic, m2;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   options.f_run_mbdyn2easyanim = false;
%!   options.positive_frequencies = false;
%!   mbdyn_solver_run(fname, options);
%!   modal = mbdyn_post_load_output_eig(fname, options, 0);
%!   omega1 = sqrt(s1 / m1 - (d1 / (2 * m1))^2);
%!   omega2 = sqrt(s2 / m2 - (d2 / (2 * m2))^2);
%!   alpha1 = -d1 / (2 * m1);
%!   alpha2 = -d2 / (2 * m2);
%!   lambda1 = alpha1 + [-1j; 1j] * omega1;
%!   lambda2 = alpha2 + [-1j; 1j] * omega2;
%!   lambda = [lambda1; lambda2];
%!   [dummy, idx] = sort(imag(lambda), "ascend");
%!   lambda = lambda(idx);
%!   assert(modal.lambda, lambda, eps^0.9 * max(abs(lambda)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
