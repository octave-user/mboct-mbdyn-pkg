%!test
%! N = 4;

%! Theta2 = 88*pi/180;
%! Theta3 = -10*pi/180;

%! X = [ linspace(1,7,N);
%!       zeros(1,N);
%!       zeros(1,N) ];

%! R = [cos(Theta2)*cos(Theta3),-sin(Theta3),sin(Theta2)*cos(Theta3);
%!      cos(Theta2)*sin(Theta3),cos(Theta3),sin(Theta2)*sin(Theta3);
%!      -sin(Theta2),0,cos(Theta2)];

%! n = [ 1, 1;
%!       0, 0;
%!       0, 0 ];

%! beam = mbdyn_pre_beam_compute(R * X, N);
%! assert_simple(R, beam.Rn(:,:,1), sqrt(eps));

%!test
%! N = 100;
%! X = [ 1,2,3,4;
%!       0,3,0,0;
%!       0,4,0,0 ];
%! beam = mbdyn_pre_beam_compute(X, N, 40);
%! norm_dXn = norm(beam.Xn(:,2:end) - beam.Xn(:,1:end-1), 2, 'cols');
%! f = max(abs(1 - norm_dXn/mean(norm_dXn)));
%! assert_simple(f < 1e-2);

%!test
%! close all;
%! N = 4;

%! Theta2 = 88*pi/180;
%! Theta3 = -10*pi/180;

%! X = [ linspace(1,7,N);
%!       zeros(1,N);
%!       zeros(1,N) ];

%! R = [cos(Theta2)*cos(Theta3),-sin(Theta3),sin(Theta2)*cos(Theta3);
%!      cos(Theta2)*sin(Theta3),cos(Theta3),sin(Theta2)*sin(Theta3);
%!      -sin(Theta2),0,cos(Theta2)];

%! n = [ 1, 1;
%!       0, 0;
%!       0, 0 ];
%! beam = mbdyn_pre_beam_compute(R*X,N);
%! figure("visible", "off");
%! mbdyn_pre_beam_plot(beam,struct("s",0.001,"Rn",false,"Rg",false));
%! title("straight beam");
%! assert_simple(R,beam.Rn(:,:,1),sqrt(eps));

%!test
%! f_print_input_file = false;
%! F = 10;
%! d = 1e-3;
%! A = d^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = d^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! D = 10e-3;
%! L = 1e-3;
%! n = 1;
%! N = n * 7;
%! Phi = linspace(0, 2 * pi * n, n * 36);
%! X = [0.5 * D * cos(Phi);
%!      0.5 * D * sin(Phi);
%!      linspace(0, L, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 10);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_beam_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
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
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, "     forces: 1;\n");
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
%!     fputs(fd, "joint: 2, total pin joint,\n");
%!     fprintf(fd, " %d,\n", columns(beam.Xn));
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fputs(fd, "     position constraint, active, active, inactive, null,\n");
%!     fputs(fd, "     orientation constraint, active, active, active, null;\n");
%!     fprintf(fd, "force: 1, absolute, %d, position, reference, %d, null, 0., 0., -1., mult, time, F;\n", columns(beam.Xn), columns(beam.Xn));
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
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert_simple(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert_simple(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert_simple(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert_simple(bodies(i).node, int32(i));
%!   endfor
%!   R = G * d^4 / (8 * n * D^3);
%!   wref = -F / R;
%!   tol = 1e-2;
%!   assert_simple(deformation{end}(end, 3), wref, tol * abs(wref));
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
%! f_print_input_file = false;
%! F = 10;
%! d = 1e-3;
%! A = d^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = d^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! D = 10e-3;
%! L = 1e-3;
%! n = 1;
%! N = n * 7;
%! Phi = linspace(0, 2 * pi * n, n * 36);
%! X = [0.5 * D * cos(Phi);
%!      0.5 * D * sin(Phi);
%!      linspace(0, L, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 10);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_beam_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
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
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, "     forces: 1;\n");
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
%!     fputs(fd, "joint: 2, total pin joint,\n");
%!     fprintf(fd, " %d,\n", columns(beam.Xn));
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fputs(fd, "     position constraint, active, active, inactive, null,\n");
%!     fputs(fd, "     orientation constraint, active, active, active, null;\n");
%!     fprintf(fd, "force: 1, absolute, %d, position, reference, %d, null, 0., 0., -1., mult, time, F;\n", columns(beam.Xn), columns(beam.Xn));
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
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert_simple(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert_simple(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert_simple(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert_simple(bodies(i).node, int32(i));
%!   endfor
%!   R = G * d^4 / (8 * n * D^3);
%!   wref = -F / R;
%!   tol = 1e-2;
%!   assert_simple(deformation{end}(end, 3), wref, tol * abs(wref));
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
%! f_print_input_file = false;
%! F = 10;
%! d = 1e-3;
%! A = d^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = d^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! D = 10e-3;
%! L = 1e-3;
%! n = 1;
%! N = n * 7;
%! Phi = linspace(0, 2 * pi * n, n * 36);
%! X = [0.5 * D * cos(Phi);
%!      0.5 * D * sin(Phi);
%!      linspace(0, L, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 10);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_beam_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
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
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, "     forces: 1;\n");
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
%!     fputs(fd, "joint: 2, total pin joint,\n");
%!     fprintf(fd, " %d,\n", columns(beam.Xn));
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fputs(fd, "     position constraint, active, active, inactive, null,\n");
%!     fputs(fd, "     orientation constraint, active, active, active, null;\n");
%!     fprintf(fd, "force: 1, absolute, %d, position, reference, %d, null, 0., 0., -1., mult, time, F;\n", columns(beam.Xn), columns(beam.Xn));
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
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert_simple(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert_simple(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert_simple(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert_simple(bodies(i).node, int32(i));
%!   endfor
%!   R = G * d^4 / (8 * n * D^3);
%!   wref = -F / R;
%!   tol = 1e-2;
%!   assert_simple(deformation{end}(end, 3), wref, tol * abs(wref));
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
%! close all;
%! N = 10;
%! X = [ 1,2,3,4;
%!       0,3,0,0;
%!       0,4,0,0 ];
%! beam = mbdyn_pre_beam_compute(X,N,50);

%! figure("visible", "off");
%! plot(beam.ti,beam.si,'-;s(t);1');
%! xlabel('t []');
%! ylabel('s [m]');
%! title('curve length versus parameter');
%! grid on;
%! grid minor on;

%! norm_dXn = norm(beam.Xn(:,2:end) - beam.Xn(:,1:end-1),2,'cols');
%! f = max(abs(1-norm_dXn/mean(norm_dXn)));
%! figure("visible","off");
%! stem(norm_dXn,'o-;norm(dXn);1');
%! xlabel('node #');
%! ylabel('norm(dXn) [m]');
%! assert_simple(f < 0.18);
%! mbdyn_pre_beam_plot(beam,struct("s",0.1,"X",true,"Rn",true,"Rg",true));

%!test
%! close all;
%! r = 0.5;
%! v = 0.5 * r;
%! h = 2;
%! n = 1;
%! N = 4;
%! k = 200;
%! Phi = linspace(0,2*pi*n,k*N+1);
%! X = [ r*cos(Phi);
%!       r*sin(Phi);
%!       h * Phi./(2*pi) ];
%!
%! Xt = [ -r*sin(Phi);
%!         r*cos(Phi);
%!         repmat(h / (2 * pi),1,length(Phi)) ];
%!
%! t = 2*pi*n/N*Xt(:,[1,end]);
%! beam = mbdyn_pre_beam_compute(X,N);
%! figure("visible", "off");
%! hold on;
%! plot3(X(1,1:k:k*N+1),X(2,1:k:k*N+1),X(3,1:k:k*N+1),'x;curve;3');
%! set(plot3(X(1,:),X(2,:),X(3,:),'-;curve;3'),'linewidth',1);
%! set(gca(),'dataaspectratio',[1,1,1]);
%! xlabel('x');
%! ylabel('y');
%! zlabel('z');
%! grid on;
%! grid minor on;
%! title('helix');
%! mbdyn_pre_beam_plot(beam,struct("s",0.08,"Rn",true,"Rg",true));
%! hold on;
%! plot3(X(1,:),X(2,:),X(3,:),'--;curve;0');

%!demo
%!demo
%! f_print_input_file = false;
%! F = 8;
%! d = 1.3e-3;
%! A = d^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = d^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! D = 10e-3;
%! L = 5e-3;
%! n = 5;
%! N = n * 21;
%! Phi = linspace(0, 2 * pi * n, n * 3600);
%! X = [0.5 * D * cos(Phi);
%!      0.5 * D * sin(Phi);
%!      linspace(0, L, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 10);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_beam_compute_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
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
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");

%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, "     forces: 1;\n");
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
%!     fputs(fd, "joint: 2, total pin joint,\n");
%!     fprintf(fd, " %d,\n", columns(beam.Xn));
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fputs(fd, "     position constraint, active, active, inactive, null,\n");
%!     fputs(fd, "     orientation constraint, active, active, active, null;\n");
%!     fprintf(fd, "force: 1, absolute, %d, position, reference, %d, null, 0., 0., -1., mult, time, F;\n", columns(beam.Xn), columns(beam.Xn));
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
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options.output_file);
%!   res.log_dat = mbdyn_post_load_log(fname);
%!   res.bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(res.log_dat.nodes)
%!     assert_simple(res.log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert_simple(res.log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(res.log_dat.beams3)
%!     for j=1:3
%!       assert_simple(res.log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(res.bodies)
%!     assert_simple(res.bodies(i).node, int32(i));
%!   endfor
%!   R = G * d^4 / (8 * n * D^3);
%!   wref = -F / R;
%!   tol = 1e-2;
%!   assert_simple(res.deformation{end}(end, 3), wref, tol * abs(wref));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
