## mbdyn_post_load_output_eig.tst:05
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_load_output_eig_XXXXXX"));
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
%!   options.use_netcdf = false;
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
%!   assert_simple(modal.lambda, lambda, eps^0.9 * max(abs(lambda)));
%!   assert_simple(modalnc.lambda, lambda, eps^0.9 * max(abs(lambda)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
