%!demo
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 20;
%! l = 1000e-3;
%! h = 500e-3;
%! g = 9.81*1000;
%! D = 20e-3;
%! A = D^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = D^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! X = [linspace(0, l, N);
%!      zeros(2, N)];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! options.first_node_number = 3;
%! options.start_node = "1";
%! options.end_node = "2";
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_plot_model_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real g = %g;\n", g);
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fprintf(fd, " set: real l = %g;\n", l);
%!     fprintf(fd, " set: real h = %g;\n", h);
%!     fprintf(fd, " set: integer N = 10000;\n");
%!     fputs(fd, "set: real t1 = N;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: t1;\n");
%!     fputs(fd, "         time step: t1 / N;\n");
%!     fputs(fd, "         linear solver: umfpack, colamd, scale, row max column max, always, max iterations, 10;\n");
%!     fputs(fd, "         method: implicit euler;\n");
%!     fputs(fd, "         max iterations: 100;\n");
%!     fputs(fd, "         tolerance: 1.e-4, 1e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         output: iterations, solver condition number, stat, yes;\n");
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
%!     fputs(fd, "       output meter: closest next, 0., forever, const, t1 / 10;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, " structural: 1, static, position, 0., 0., h, orientation, eye, velocity, null, angular velocity, null;\n");
%!     fputs(fd, " structural: 2, static, position, l, 0., h, orientation, eye, velocity, null, angular velocity, null;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, " joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., mult, time, g / t1;\n");
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
%!   [res.t, res.trajectory, res.deformation, res.velocity, res.acceleration, res.node_id] = mbdyn_post_load_output_struct(options.output_file);
%!   res.log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   [mesh, sol] = mbdyn_post_export_model(res);
%!   if (~isempty(pkg("list", "mboct-fem-pkg")))
%!     pkg load mboct-fem-pkg;
%!     opts.print_and_exit = true;
%!     opts.print_to_file = fname;
%!     opts.rotation_angle = [pi/2, 0, 0];
%!     opts.skin_only = false;
%!     opts.show_element = false;
%!     opts.animation_delay = 0;
%!     opts.scale_def = 1;
%!     fem_post_sol_external(mesh, sol, opts);
%!     fn = dir([opts.print_to_file, "*.jpg"]);
%!     for i=1:numel(fn)
%!       [img, map, alpha] = imread(fullfile(fn(i).folder, fn(i).name));
%!       figure("visible", "off");
%!       imshow(img, map);
%!       title(sprintf("deformed mesh step %d", i));
%!     endfor
%!   endif
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
