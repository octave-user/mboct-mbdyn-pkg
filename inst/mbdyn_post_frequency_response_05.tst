## mbdyn_post_frequency_response.tst:05
%!test
%! try
%! ##  Frequency response of a cantilever beam
%! ##  Robert Gash, Klaus Knothe
%! ##  Strukturdynamik Band 2
%! ##  Page 25 figure 2.10, equation 2.38b
%! f_print_input_file = false;
%! f_plot_response = false;
%! if (f_plot_response)
%!   close("all");
%! endif
%! N = 40;
%! M = 200;
%! l = 4;
%! X = [linspace(0, l, N);
%!      zeros(2, N)];
%! g = 9.81;
%! D = 10e-3;
%! A = D^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = D^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! B = E * Iy;
%! mu = rho * A;
%! P0 = (0.3 * l) / l^3 * 3 * B;
%! omega1 = sqrt(B / (mu * l^4));
%! wstat = P0 * l^3 / (3 * B);
%! omega_crit = [3.516, 22.035, 61.697] * omega1;
%! omega = linspace(1e-5, 1, M) * 0.5 * omega_crit(end);
%! DELTA = sqrt(omega ./ omega1) / l;
%! V = 3 ./ (DELTA * l).^3 .* (sin(DELTA * l) .* cosh(DELTA * l) - cos(DELTA * l) .* sinh(DELTA * l)) ./ (1 + cos(DELTA * l) .* cosh(DELTA * l));
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_frequency_response_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1e-2;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         eigenanalysis: 0,\n");
%!     fputs(fd, "          suffix format, \"%02d\",\n");
%!     fputs(fd, "          output sparse matrices,\n");
%!     fputs(fd, "          output eigenvectors,\n");
%!     fputs(fd, "          results output precision, 16,\n");
%!     fputs(fd, "          parameter, 0.01, use arpack, 100, 200, 0;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     output results: netcdf, text;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   if (f_print_input_file)
%!     spawn_wait(spawn("cat", {fname}));
%!   endif
%!   options.verbose = false;
%!   if (~options.verbose)
%! # options.logfile = [fname, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(fname, options);
%!   log_dat = mbdyn_post_load_log(fname);
%!   nc = [true, false];
%!   modal = cell(size(nc));
%!   F = cell(size(nc));
%!   for idxnc=1:numel(nc)
%!     modal{idxnc} = mbdyn_post_load_output_eig(options.output_file, struct("use_netcdf", nc(idxnc)));
%!     excitation.node_label = columns(beam.Xn);
%!     excitation.offset = [0; 0; 0];
%!     excitation.direction = [0; 0; 1];
%!     response.node_label = columns(beam.Xn);
%!     response.offset = [0; 0; 0];
%!     response.direction = [0; 0; 1];
%!     options.matrix_type = "wre";
%!     options.solver = "pardiso";
%!     F{idxnc} = mbdyn_post_frequency_response(modal{idxnc}, log_dat.dof_info, excitation, response, omega, P0, options);
%!     if (f_plot_response)
%!       figure("visible", "off");
%!       subplot(2, 1, 1);
%!       hold on;
%!       semilogy(omega / (2 * pi), abs(F{idxnc}), "-;abs(F);r");
%!       semilogy(omega / (2 * pi), abs(V * wstat), "-;abs(V);k");
%!       grid minor on;
%!       xlabel("f [Hz]");
%!       ylabel("abs(F) [m]");
%!       title("frequency response magnitude");
%!       subplot(2, 1, 2);
%!       hold on;
%!       plot(omega / (2 * pi), 180 / pi * arg(F{idxnc}), "-;arg(F);r");
%!       plot(omega / (2 * pi), 180 / pi * arg(V * wstat), "-;arg(V);k");
%!       grid minor on;
%!       xlabel("f [Hz]");
%!       ylabel("arg(F) [deg]");
%!       title("frequency response phase");
%!     endif
%!     tol = 1e-2;
%!     assert_simple(sum(abs(F{idxnc} - V.' * wstat).^2) / sum(abs(V * wstat).^2) < tol);
%!   endfor
%!   fn = fieldnames(modal{1});
%!   for idxfn=1:numel(fn)
%!     assert_simple(getfield(modal{1}, fn{idxfn}), getfield(modal{2}, fn{idxfn}));
%!   endfor
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
