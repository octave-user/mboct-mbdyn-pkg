## mbdyn_post_frequency_response.tst:01
%!test
%! try
%! ## TEST 1
%! ##  Frequency response of a cantilever beam
%! ##  Robert Gash, Klaus Knothe
%! ##  Strukturdynamik Band 2
%! ##  Page 25 figure 2.10, equation 2.38b
%! f_print_input_file = false;
%! f_plot_response = false;
%! options.verbose = false;
%! if (f_plot_response)
%!   close("all");
%! endif
%! E = [70000e6, 210000e6, 4e6, 80000e6];
%! nu = [0.3, 0.3, 0.5, 0.4];
%! rho = [2700, 7850, 1000, 1];
%! which_eigenvalues = {"smallest magnitude", "largest imaginary part"};
%! solvers = {"lapack", "arpack"};
%! linear_solvers = {"umfpack, grad", "klu, grad", "naive, colamd"};
%! assert_simple(numel(E), numel(rho));
%! for n=1:numel(linear_solvers)
%!   for m=1:numel(solvers)
%!     for k=1:numel(which_eigenvalues)
%!       switch (solvers{m})
%!         case "lapack"
%!           switch(k)
%!             case 1
%!             otherwise
%!                     # skip this test because option is not used at all
%!               continue;
%!           endswitch
%!           switch (n)
%!             case 1
%!             otherwise
%!                     # skip this test because option is not used at all
%!               continue
%!           endswitch
%!       endswitch
%!       for i=1:numel(E)
%!         N = 30;
%!         M = 500;
%!         l = 4;
%!         X = [linspace(0, l, N);
%!              zeros(2, N)];
%!         g = 9.81;
%!         D = 10e-3;
%!         A = D^2 * pi / 4.;
%!         Ay = 9. / 10. * A;
%!         Az = Ay;
%!         Iy = D^4 * pi / 64.;
%!         Iz = Iy;
%!         It = Iy + Iz;
%!         G = E(i) / (2 * (1 + nu(i)));
%!         B = E(i) * Iy;
%!         mu = rho(i) * A;
%!         P0 = (0.3 * l) / l^3 * 3 * B;
%!         omega1 = sqrt(B / (mu * l^4));
%!         wstat = P0 * l^3 / (3 * B);
%!         laml = linspace(0, 50, 1000);
%!         flambda = @(laml) cos(laml) + 1 ./ cosh(laml);
%!         flam = flambda(laml);
%!         idx = find(sign(flam(2:end)) ~= sign(flam(1:end-1)));
%!         omega_crit = zeros(1, numel(idx));
%!         for j=1:numel(idx)
%!           omega_crit(j) = fzero(flambda, laml([idx(j), idx(j)+1]))^2 * omega1;
%!         endfor
%!         omega = linspace(1e-5, 1, M) * 100 * omega_crit(1);
%!         DELTA = sqrt(omega ./ omega1) / l;
%!         V = 3 ./ (DELTA * l).^3 .* (sin(DELTA * l) .* cosh(DELTA * l) - cos(DELTA * l) .* sinh(DELTA * l)) ./ (1 + cos(DELTA * l) .* cosh(DELTA * l));
%!         options.A = "A";
%!         options.Ip = "It";
%!         options.Iy = "Iy";
%!         options.Iz = "Iz";
%!         options.constitutive_law_number = "const_law_id_beam";
%!         options.rho = "rho";
%!         beam = mbdyn_pre_beam_compute(X, N, 20);
%!         fd = -1;
%!         unwind_protect
%!           unwind_protect
%!             [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_frequency_response_XXXXXX"));
%!             if (fd == -1)
%!               error("failed to open temporary file");
%!             endif
%!             fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!             fprintf(fd, " set: real A = %g;\n", A);
%!             fprintf(fd, " set: real Ay = %g;\n", Ay);
%!             fprintf(fd, " set: real Az = %g;\n", Az);
%!             fprintf(fd, " set: real Iy = %g;\n", Iy);
%!             fprintf(fd, " set: real Iz = %g;\n", Iz);
%!             fprintf(fd, " set: real It = %g;\n", It);
%!             fprintf(fd, " set: real E = %g;\n", E(i));
%!             fprintf(fd, " set: real G = %g;\n", G);
%!             fprintf(fd, " set: real rho = %g;\n", rho(i));
%!             fprintf(fd, " set: real fmin = %g;\n", 1e-1 * min(omega_crit) / (2 * pi));
%!             fprintf(fd, " set: real fmax = %g;\n", 1.5 * max(omega_crit) / (2 * pi));
%!             fputs(fd, " begin: data;\n");
%!             fputs(fd, "         problem: initial value;\n");
%!             fputs(fd, " end: data;\n");
%!             fputs(fd, " begin: initial value;\n");
%!             fputs(fd, "         initial time: 0;\n");
%!             fputs(fd, "         final time: 1e-2;\n");
%!             fputs(fd, "         time step: 1e-2;\n");
%!             fprintf(fd, "         linear solver: %s, scale, iterative, once, max iterations, 10;\n", linear_solvers{n});
%!             fputs(fd, "         method: ms, 0.6;\n");
%!             fputs(fd, "         max iterations: 10;\n");
%!             fputs(fd, "         tolerance: 1.e-6;\n");
%!             fputs(fd, "         threads: assembly, 1;\n");
%!             fputs(fd, "         derivatives max iterations: 10;\n");
%!             fputs(fd, "         derivatives coefficient: auto;\n");
%!             fputs(fd, "         eigenanalysis: 0,\n");
%!             fputs(fd, "          suffix format, \"%02d\",\n");
%!             fputs(fd, "          output matrices,\n");
%!             fputs(fd, "          output eigenvectors,\n");
%!             fputs(fd, "          output geometry,\n");
%!             fputs(fd, "          lower frequency limit, fmin,\n");
%!             fputs(fd, "          upper frequency limit, fmax,\n");
%!             fputs(fd, "          results output precision, 16,\n");
%!             fprintf(fd, "        mode, %s", which_eigenvalues{k});
%!             switch (solvers{m})
%!               case "arpack"
%!                 nev = 5 * numel(omega_crit);
%!                 ncv = 2 * nev + 1;
%!                 fprintf(fd, ",\n        use arpack, %d, %d, 0", nev, ncv);
%!               case "lapack"
%!                 fprintf(fd, ",\n        use lapack");
%!             endswitch
%!             fputs(fd, ";\n");
%!             fputs(fd, "         nonlinear solver: nox, modified, 10,\n");
%!             fputs(fd, "             keep jacobian matrix,\n");
%!             fputs(fd, "             inner iterations before assembly, 6,\n");
%!             fputs(fd, "             jacobian operator, newton krylov,\n");
%!             fputs(fd, "             solver, line search based,\n");
%!             fputs(fd, "             forcing term, type 2,\n");
%!             fputs(fd, "             direction, newton,\n");
%!             fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!             fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!             fputs(fd, "             linear solver, gmres,\n");
%!             fputs(fd, "             linear solver max iterations, 12,\n");
%!             fputs(fd, "             krylov subspace size, 12;\n");
%!             fputs(fd, " end: initial value;\n");
%!             fputs(fd, " begin: control data;\n");
%!             fputs(fd, "     use automatic differentiation;\n");
%!             fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!             fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!             fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!             fputs(fd, "     joints: 1;\n");
%!             fputs(fd, "     output results: netcdf, text;\n");
%!             fputs(fd, " end: control data;\n");
%!             fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!             mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!             fputs(fd, " begin: nodes;\n");
%!             mbdyn_pre_beam_write_nodes(beam, fd, options);
%!             fputs(fd, " end: nodes;\n");
%!             fputs(fd, " begin: elements;\n");
%!             mbdyn_pre_beam_write_bodies(beam, fd, options);
%!             mbdyn_pre_beam_write_beams(beam, fd, options);
%!             fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!             fputs(fd, " end: elements;\n");
%!           unwind_protect_cleanup
%!             if (fd ~= -1)
%!               fclose(fd);
%!             endif
%!           end_unwind_protect
%!           options.output_file = fname;
%!           if (f_print_input_file)
%!             spawn_wait(spawn("cat", {fname}));
%!           endif
%!           if (~options.verbose)
%! # options.logfile = [fname, ".stdout"];
%!           endif
%!           info = mbdyn_solver_run(fname, options);
%!           log_dat = mbdyn_post_load_log(fname);
%!           nc = [false, true];
%!           modal = cell(size(nc));
%!           F = cell(size(nc));
%!           for idxnc=1:numel(nc)
%!             modal{idxnc} = mbdyn_post_load_output_eig(options.output_file, struct("use_netcdf", nc(idxnc)));
%!             excitation.node_label = columns(beam.Xn);
%!             excitation.offset = [0; 0; 0];
%!             excitation.direction = [0; 0; 1];
%!             response.node_label = columns(beam.Xn);
%!             response.offset = [0; 0; 0];
%!             response.direction = [0; 0; 1];
%!             options.matrix_type = "wre";
%!             options.number_of_processors = int32(1);
%!             F{idxnc} = mbdyn_post_frequency_response(modal{idxnc}, log_dat.dof_info, excitation, response, omega, P0, options);
%!             if (f_plot_response)
%!               figure("visible", "off");
%!               subplot(2, 1, 1);
%!               hold on;
%!               semilogy(omega / (2 * pi), abs(F{idxnc}), "-;abs(F);r");
%!               semilogy(omega / (2 * pi), abs(V * wstat), "-;abs(V);k");
%!               grid minor on;
%!               xlabel("f [Hz]");
%!               ylabel("abs(F) [m]");
%!               title("frequency response magnitude");
%!               subplot(2, 1, 2);
%!               hold on;
%!               plot(omega / (2 * pi), 180 / pi * arg(F{idxnc}), "-;arg(F);r");
%!               plot(omega / (2 * pi), 180 / pi * arg(V * wstat), "-;arg(V);k");
%!               grid minor on;
%!               xlabel("f [Hz]");
%!               ylabel("arg(F) [deg]");
%!               title("frequency response phase");
%!             endif
%!           endfor
%!         unwind_protect_cleanup
%!           if (fd ~= -1)
%!             unlink(fname);
%!             files = dir([fname, "*"]);
%!             for j=1:numel(files)
%!               unlink(fullfile(files(j).folder, files(j).name));
%!             endfor
%!           endif
%!         end_unwind_protect
%!         tol = 2e-2;
%!         for idxnc=1:numel(F)
%!           assert_simple(sum(abs(F{idxnc} - V.' * wstat).^2) / sum(abs(V * wstat).^2) < tol);
%!           wyz = zeros(1, 2 * numel(omega_crit));
%!           idxyz = 0;
%!           idx_node = find(modal{idxnc}.labels == excitation.node_label);
%!           for o=1:numel(modal{idxnc}.lambda)
%!             Vi = abs(modal{idxnc}.VR(modal{idxnc}.idx(idx_node) + (1:6), o));
%!             switch (max(Vi))
%!               case {Vi(2), Vi(6), Vi(3), Vi(5)}
%!                 if (idxyz < numel(wyz))
%!                   wyz(++idxyz) = 2 * pi * modal{idxnc}.f(o); ## extract only the bending modes in the x-y and y-z plane
%!                 endif
%!             endswitch
%!           endfor
%!           erryz = max(max(abs(wyz(1:2:end-1) ./ omega_crit - 1)), max(abs(wyz(2:2:end) ./ omega_crit - 1)));
%!           if (options.verbose)
%!             fprintf(stderr, "err(%d:%d:%d:%d): %.3f%%\n", n, m, k, i, 100 * erryz);
%!           endif
%!           assert_simple(erryz < tol);
%!         endfor
%!         fn = fieldnames(modal{1});
%!         for idxfn=1:numel(fn)
%!           assert_simple(getfield(modal{1}, fn{idxfn}), getfield(modal{2}, fn{idxfn}));
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%! endfor
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
