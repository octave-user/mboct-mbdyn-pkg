## mbdyn_post_deformations_scale.tst:03
%!test
%! try
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
%! %unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_deformation_scale_XXXXXX"));
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
%!   for k=1:2
%!     node_groups.X0_node_id = log_dat.nodes(1).label;
%!     switch (k)
%!     case 1
%!       node_groups.R0ref_node_id = log_dat.nodes(1).label;
%!     case 2
%!       node_groups.R0ref_node_id = [log_dat.nodes([1, ceil(end / 2), end]).label];
%!     endswitch
%!     node_groups.scale = 100;
%!     node_groups.node_id_scale = [log_dat.nodes.label];
%!     mbdyn_post_deformations_scale([fname, ".log"], [fname, ".mov"], [fname, "_scaled.mov"], node_groups);
%!     copyfile([fname, ".log"], [fname, "_scaled.log"]);
%!     copyfile([fname, ".out"], [fname, "_scaled.out"]);
%!     [scaled{k}.t, scaled{k}.trajectory, scaled{k}.deformation] = mbdyn_post_load_output_struct([fname, "_scaled"]);
%!   endfor
%!   R = G * d^4 / (8 * n * D^3);
%!   wref = -F / R;
%!   tol = 1e-2;
%!   assert_simple(deformation{end}(end, 3), wref, tol * abs(wref));
%!   for k=1:numel(scaled)
%!     assert_simple(scaled{k}.deformation{end}(end, 3), node_groups.scale * wref, tol * abs(node_groups.scale * wref));
%!   endfor
%! %unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! %end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
