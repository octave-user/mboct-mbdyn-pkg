## Copyright (C) 2014(-2023) Reinhard <octave-user@a1.net>
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

  if (~isfield(options, "use_netcdf"))
    options.use_netcdf = false;
  endif

  [inp_dir, inp_name, inp_ext] = fileparts(mbdyn_output_file);

  if (options.use_netcdf)
    mbdyn_output_file = fullfile(inp_dir, [inp_name, ".nc"]);

    [info, err, msg] = stat(mbdyn_output_file);

    if (err ~= 0)
      error("file not found \"%s\"", mbdyn_output_file);
    endif

    modal = mbdyn_post_load_modal_data_nc(mbdyn_output_file, index);
  else
    mbdyn_output_file = fullfile(inp_dir, sprintf("%s_%02d.m", inp_name, index));

    [info, err, msg] = stat(mbdyn_output_file);

    if (err ~= 0)
      mbdyn_output_file = fullfile(inp_dir, cstrcat(inp_name, ".m"));
      [info, err, msg] = stat(mbdyn_output_file);
    endif

    if (err ~= 0)
      error("file not found \"%s\"", mbdyn_output_file);
    endif

    modal = mbdyn_post_load_modal_data_m(mbdyn_output_file);
  endif

  if (isfield(modal, "alpha"))
    LAMBDA = (modal.alpha(:,1) + 1j * modal.alpha(:, 2)) ./ modal.alpha(:, 3);
    modal.lambda = 1 / modal.dCoef * (LAMBDA - 1) ./ (LAMBDA + 1);
    modal.f = imag(modal.lambda) / (2 * pi);
  else
    modal.lambda = [];
    modal.f = [];
  endif

  if (~isempty(modal.f))
    [modal.f, idx_f] = sort(modal.f);
    if (options.positive_frequencies)
      idx_gtz = find(modal.f > 0);
      modal.f = modal.f(idx_gtz);
      idx_f = idx_f(idx_gtz);
    endif
    modal.lambda = modal.lambda(idx_f);
    if (isfield(modal, "VR"))
      modal.VR = modal.VR(:, idx_f);
    endif

    if (isfield(modal, "VL"))
      modal.VL = modal.VL(:, idx_f);
    endif
  endif

  if (isfield(modal, "lambda"))
    delta = -real(modal.lambda);
    omegad = imag(modal.lambda);
    a0 = (delta ./ omegad).^2;
    modal.D = sqrt(a0 ./ (1 + a0));
    modal.omega0 = omegad ./ sqrt(1 - modal.D.^2);
  endif
endfunction

function modal = mbdyn_post_load_modal_data_m(mbdyn_output_file)
  source(mbdyn_output_file);

  var_names = {"dTime", "dCoef", "Aplus", "Aminus", "VR", "VL", "alpha", "X0", "idx", "labels"};

  modal = struct();

  for i=1:numel(var_names)
    if (1 == exist(var_names{i}, "var"))
      modal = setfield(modal, var_names{i}, eval(var_names{i}));
    endif
  endfor
endfunction

function modal = mbdyn_post_load_modal_data_nc(mbdyn_output_file, index)
  pkg load netcdf;

  info = ncinfo(mbdyn_output_file);

  modal = struct();

  prefix = sprintf("eig.%d.", index);

  var_names = {[prefix, "time"],   "dTime";
               [prefix, "dCoef"],  "dCoef";
               [prefix, "Aplus"],  "Aplus";
               [prefix, "Aminus"], "Aminus";
               [prefix, "VR"],     "VR";
               [prefix, "VL"],     "VL";
               [prefix, "alpha"],  "alpha";
               "eig.idx",          "idx";
               "eig.labels",       "labels";
               "eig.joint.idx",    "joint_idx";
               "eig.joint.labels", "joint_labels"};

  for i=1:rows(var_names)
    try
      value = ncread(mbdyn_output_file, var_names{i, 1});
    catch
      switch (lasterror.message)
        case "NetCDF: Variable not found"
          continue
        otherwise
          rethrow(lasterror());
      endswitch
    end_try_catch
    modal = setfield(modal, var_names{i, 2}, value);
  endfor

  if (isfield(modal, "time"))
    modal.dTime = modal.time;
    modal = rmfield(modal, "time");
  endif

  if (isfield(modal, "alpha"))
    modal.alpha = modal.alpha.';
  endif

  if (isfield(modal, "VR"))
    modal.VR = complex(modal.VR(:, :, 1), modal.VR(:, :, 2));
  endif

  if (isfield(modal, "VL"))
    modal.VL = complex(modal.VL(:, :, 1), modal.VL(:, :, 2));
  endif

  if (isfield(modal, "Aplus"))
    Aplus_attr = ncreadatt(mbdyn_output_file, [prefix, "Aplus"], "description");

    Aplus_format = "sparse";

    if (strncmp(Aplus_attr, Aplus_format, length(Aplus_format)))
      modal.Aplus = spconvert(modal.Aplus.');
      modal.Aminus = spconvert(modal.Aminus.');
    endif
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
%!     fputs(fd, "     output results: netcdf, text;\n");
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
%!   options.use_netcdf = true;
%!   modalnc = mbdyn_post_load_output_eig(fname, options, 0);
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
%!   assert(modalnc.lambda, lambda, eps^0.9 * max(abs(lambda)))
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
%!     fputs(fd, "     output results: netcdf, text;\n");
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
%!   options.use_netcdf = true;
%!   modalnc = mbdyn_post_load_output_eig(fname, options, 0);
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
%!   assert(modalnc.lambda, lambda, eps^0.9 * max(abs(lambda)));
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
%!     fputs(fd, "     output results: netcdf, text;\n");
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
%!   options.use_netcdf = true;
%!   modalnc = mbdyn_post_load_output_eig(fname, options, 0);
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
%!   assert(modalnc.lambda, lambda, eps^0.9 * max(abs(lambda)));
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
%!     fputs(fd, "     output results: netcdf, text;\n");
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
%!   options.use_netcdf = true;
%!   modalnc = mbdyn_post_load_output_eig(fname, options, 0);
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
%!   assert(modalnc.lambda, lambda, eps^0.9 * max(abs(lambda)));
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
%! ## TEST 330
%! ## VIBRATIONS OF COMPLETE SPHERICAL SHELLS WITH IMPERFECTIONS
%! ## Thomas A. Duffey
%! ## Jason E. Pepin
%! ## Amy N. Robertson
%! ## Michael L. Steinzig
%! ## Internatial Modal Analysis Conference (IMAC-XXIII)
%! ## Orlando, Florida
%! ## January 31-February 3, 2005
%! ## Los Alamos
%! ## NATIONAL LABORATORY
%! pkg load mboct-fem-pkg;
%! do_plot = false;
%! if (do_plot)
%!   close all;
%! endif
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   unwind_protect
%!     [fd, msg] = fopen([filename, ".geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s.geo\"", filename);
%!     endif
%!     options.verbose = true;
%!     SI_unit_meter = 1;
%!     SI_unit_second = 1;
%!     SI_unit_kilogram = 1;
%!     SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!     SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!     SI_unit_rad = 1;
%!     E = 28e6 * 6895 / SI_unit_pascal;
%!     nu = 0.28;
%!     rho = 0.000751 * 4.4482 / (25.4e-3^4) / (SI_unit_kilogram / SI_unit_meter^3);
%!     R = 4.4688 * 25.4e-3 / SI_unit_meter;
%!     t = 0.0625 * 25.4e-3 / SI_unit_meter;
%!     W0 = [0; 0; 0] / (SI_unit_rad / SI_unit_second);
%!     fmin = 1. / SI_unit_second^-1;
%!     fmax = 10000 / SI_unit_second^-1;
%!     number_of_modes = 39*2;
%!     t1 = 1 / SI_unit_second;
%!     N = 1;
%!     mesh_size = 2 * R * pi / 36 / SI_unit_meter;
%!     fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fprintf(fd, "R = %g;\n", R);
%!     fprintf(fd, "t = %g;\n", t);
%!     fprintf(fd, "h = %g;\n", mesh_size);
%!     fputs(fd, "Point(1) = {0,0,-R - 0.5 * t,h};\n");
%!     fputs(fd, "Point(2) = {0,0,-R + 0.5 * t,h};\n");
%!     fputs(fd, "Line(1) = {1,2};\n");
%!     fputs(fd, "A[] = Extrude{{0,1,0},{0,0,0},Pi}{Line{1}; Layers{Ceil(R * Pi / h)};Recombine;};\n");
%!     fputs(fd, "V1[] = Extrude{{0,0,1},{0,0,0},Pi}{Surface{A[1]};Layers{Ceil(R * Pi/h)};Recombine;};\n");
%!     fputs(fd, "V2[] = Extrude{{0,0,1},{0,0,0},Pi}{Surface{V1[0]};Layers{Ceil(R * Pi/h)};Recombine;};\n");
%!     fputs(fd, "Physical Volume(\"volume\", 1) = {V1[1],V2[1]};\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   [~] = unlink([filename, ".msh"]);
%!   pid = spawn("gmsh", {"-format", "msh2", "-0", "-order", "2", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   opt_msh.elem_type = {"iso20r", "penta15"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   mesh.materials.iso20r = ones(rows(mesh.elements.iso20r), 1, "int32");
%!   mesh.materials.penta15 = ones(rows(mesh.elements.penta15), 1, "int32");
%!   mesh.material_data.rho = rho;
%!   mesh.material_data.C = fem_pre_mat_isotropic(E, nu);
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   elem_file = [filename, ".elem"];
%!   mbdyn_file = [filename, ".mbdyn"];
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_center";
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fputs(fd, "set: integer ref_id_center = 1001;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fprintf(fd, "    final time: %g;\n", t1);
%!     fprintf(fd, "    time step: %g;\n", t1 / N);
%!     fputs(fd, "    max iterations: 100;\n");
%!     fprintf(fd, "    tolerance: %g, test, norm, %g, test, norm;\n", 1e-4 / SI_unit_newton, 1e-6 / SI_unit_meter);
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 100;\n");
%!     fputs(fd, "    method: ms4, 0.6;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30, linear solver max iterations, 60, verbose, 3;\n");
%!     fputs(fd, "    threads: assembly, 4;\n");
%!     fputs(fd, "    threads: solver, 4;\n");
%!     fputs(fd, "    eigenanalysis: list, 1, 0.,\n");
%!     fputs(fd, "    output eigenvectors,\n");
%!     fputs(fd, "        output geometry,\n");
%!     fprintf(fd, "        lower frequency limit, %g, upper frequency limit, %g,\n", fmin, fmax);
%!     fprintf(fd, "    use arpack,%d,%d,0.,suffix format,\"%%02d\";\n", number_of_modes, 2 * number_of_modes + 1);
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, "     output results: netcdf, text;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " reference: ref_id_center,\n");
%!     fputs(fd, "   position, reference, global, null,\n");
%!     fputs(fd, "   orientation, reference, global, eye,\n");
%!     fputs(fd, "   velocity, reference, global, null,\n");
%!     fprintf(fd, "   angular velocity, reference, global%s;\n", sprintf(", %g", W0));
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!     options_mbd.output_file = sprintf("%s_mbd", filename);
%!     options_eig.positive_frequencies = false;
%!     if (~options.verbose)
%!       options_mbd.logfile = [options_mbd.output_file, ".stdout"];
%!     endif
%!     options_mbd.mbdyn_command = "mbdyn";
%!     info = mbdyn_solver_run(mbdyn_file, options_mbd);
%!     [mesh_sol, sol] = mbdyn_post_load_output_sol(options_mbd.output_file);
%!     modal = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig);
%!     options_eig.use_netcdf = true;
%!     modalnc = mbdyn_post_load_output_eig(options_mbd.output_file, options_eig);
%!     fref = [5078; 6005; 6378; 6729] / SI_unit_second^-1;
%!     assert(modal.f([7, 12, 19, 39] - 6), fref, 5e-4 * max(fref));
%!     assert(modalnc.f([7, 12, 19, 39] - 6), fref, 5e-4 * max(fref));
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       [~] = unlink(fullfile(fn(i).folder, fn(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
