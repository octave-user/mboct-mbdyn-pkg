## mbdyn_pre_beam_write_bodies.tst:03
%!test
%! try
%! f_plot_res = false;
%! fd = -1;
%! unwind_protect
%! shape = "straight x";
%! param.initial_time = 0;
%! param.final_time = 1;
%! param.E = 210000e6;
%! param.G = 81500e6;
%! d = 12e-3;
%! l = 5000e-3;
%! param.rho = 7850;
%! param.number_of_elements = 50;
%! interpolation_points = 1000;
%! [fd, input_file] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_beam_write_bodies_XXXXXX"));
%! if (fd == -1)
%!   error("failed to open temporary file");
%! endif
%! unwind_protect
%! options_pre.reference_frame = "ref_id_beam";
%! options_pre.first_reference_frame_number = "ref_id_1";
%! options_pre.first_node_number = "node_id_1";
%! options_pre.first_beam_number = "beam_id_1";
%! options_pre.first_body_number = "body_id_1";
%! options_pre.constitutive_law_number = "const_law_id_1";
%! options_pre.rho = "rho";
%! options_pre.A = "A";
%! options_pre.Iy = "Iy";
%! options_pre.Iz = "Iz";
%! options_pre.Ip = "Ip";
%! options_pre.output_flag = "output_nodes";
%! options_pre.node_type = "dynamic";
%! options_pre.reorient_cross_sections = true;
%! options_run.f_run_mbdyn = true;
%! options_run.f_run_mbdyn2easyanim = true;
%! options_run.f_runEasyAnim = false;
%! options_run.verbose = false;
%! options_run.logfile = [input_file, ".stdout"];
%! param.A = d^2 * pi / 4;
%! param.Ay = param.Az = param.A * 9 / 10;
%! param.Iy = param.Iz = d^4 * pi / 64;
%! param.It = param.Ip = param.Iy + param.Iz;
%! param.output_nodes = uint8(true);
%! param.F1 = 0.001;
%! param.x0 = 0;
%! param.y0 = 0;
%! param.z0 = 0;
%! param.Phix0 = 0;
%! param.Phiy0 = 0;
%! param.Phiz0 = 0;
%! switch( shape )
%!     case "straight x"
%!         X = [ 0, l/3, l/2, l;
%!               0, 0,   0,  0;
%!               0, 0,   0,  0 ];
%!     case "straight y"
%!         X = [ 0, 0,    0,   0;
%!               0, l/3,  l/2, l;
%!               0, 0,    0,   0 ];
%!     case "straight z"
%!         X = [ 0, 0 ,  0,   0;
%!               0, 0,   0,   0;
%!               0, l/3, l/2, l ];
%!     case "straight -x"
%!         X = [ 0, -l/3, -l/2, -l;
%!               0, 0,   0,  0;
%!               0, 0,   0,  0 ];
%!     case "straight -y"
%!         X = [ 0, 0,    0,   0;
%!               0, -l/3,  -l/2, -l;
%!               0, 0,    0,   0 ];
%!     case "straight -z"
%!         X = [ 0, 0 ,  0,   0;
%!               0, 0,   0,   0;
%!               0, -l/3, -l/2, -l ];
%!     case "circular arc xy"
%!         R = l;
%!         Phi = linspace(0,pi/2,interpolation_points);
%!         X = [ -R * cos(Phi);
%!                R * sin(Phi);
%!                zeros(1,columns(Phi)) ];
%!     otherwise
%!         error("unknown shape \"%s\"", shape);
%! endswitch
%! param.number_of_nodes = 2 * param.number_of_elements + 1;
%! switch (options_pre.node_type)
%!     case "dynamic"
%!         param.number_of_bodies = param.number_of_nodes;
%!     case "static"
%!         param.number_of_bodies = 0;
%!     otherwise
%!         error("invalid node type \"%s\"", options_pre.node_type);
%! endswitch
%! beam = mbdyn_pre_beam_compute(X, param.number_of_elements, ceil(interpolation_points/columns(X)));
%! if (f_plot_res)
%!   figure("visible", "off");
%!   mbdyn_pre_beam_plot(beam,struct("Rn",true,"Rg",false));
%!   set(plot3(X(1,:),X(2,:),X(3,:),'-;X;c'),'linewidth',3);
%! endif
%! mbdyn_pre_write_param_file(fd, param);
%! fputs(fd, "begin: data;\n");
%! fputs(fd, "    problem: initial value; # the default\n");
%! fputs(fd, "end: data;\n");
%! fputs(fd, "begin: initial value;\n");
%! fputs(fd, "    initial time: 0;\n");
%! fputs(fd, "    final time: 1;\n");
%! fputs(fd, "    time step: 1e-2;\n");
%! fputs(fd, "    max iterations: 10;\n");
%! fputs(fd, "    derivatives coefficient: 1e-8;\n");
%! fputs(fd, "    derivatives max iterations: 10;\n");
%! fputs(fd, "    derivatives tolerance: 1e-5;\n");
%! fputs(fd, "    tolerance: 1e-5;\n");
%! fputs(fd, "    eigenanalysis: initial_time,\n");
%! fputs(fd, "    suffix format, \"%02d\",\n");
%! fputs(fd, "    output sparse matrices, \n");
%! fputs(fd, "    parameter, 0.01,\n");
%! fputs(fd, "    output eigenvectors, output geometry,lower frequency limit, 0.001,upper frequency limit, 10000,\n");
%! fputs(fd, "    use arpack,12,120,0;\n");
%! fputs(fd, "    output: counter, messages;\n");
%! fputs(fd, "    threads: assembly, 1;\n");
%! fputs(fd, "    nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%! fputs(fd, "end: initial value;\n");
%! fputs(fd, "begin: control data;\n");
%! fputs(fd, "    model: static;\n");
%! fputs(fd, "    structural nodes: number_of_nodes;\n");
%! fputs(fd, "    rigid bodies: number_of_bodies;\n");
%! fputs(fd, "    beams: number_of_elements;\n");
%! fputs(fd, "    default output: none, forces;\n");
%! fputs(fd, "    forces: 1;\n");
%! fputs(fd, "    joints: 1;\n");
%! fputs(fd, "end: control data;\n");
%! fputs(fd, "set: integer node_id_1  = 1001;\n");
%! fputs(fd, "set: integer body_id_1  = 2001;\n");
%! fputs(fd, "set: integer beam_id_1  = 3001;\n");
%! fputs(fd, "set: integer ref_id_beam = 4001;\n");
%! fputs(fd, "set: integer ref_id_1   = 4002;\n");
%! fputs(fd, "set: integer force_id_1 = 5001;\n");
%! fputs(fd, "set: integer const_law_id_1 = 6001;\n");
%! fputs(fd, "set: integer joint_id_clamp = 7001;\n");
%! fputs(fd, "constitutive law: const_law_id_1,\n");
%! fputs(fd, "    6, linear viscoelastic generic,\n");
%! fputs(fd, "    diag, E * A, G * Ay, G * Az, \n");
%! fputs(fd, "          G * It, E * Iy, E * Iz,\n");
%! fputs(fd, "    proportional, 0;\n");
%! fputs(fd, "reference: ref_id_beam,\n");
%! fputs(fd, "    position, reference, global, x0, y0, z0,\n");
%! fputs(fd, "    orientation, reference, global, euler123, Phix0, Phiy0, Phiz0,\n");
%! fputs(fd, "    velocity, reference, global, null,\n");
%! fputs(fd, "    angular velocity, reference, global, null;\n");
%! mbdyn_pre_beam_write_reference_frames(beam, fd, options_pre);
%! fputs(fd, "begin: nodes;\n");
%! mbdyn_pre_beam_write_nodes(beam, fd, options_pre);
%! fputs(fd, "end: nodes;\n");
%! fputs(fd, "begin: elements;\n");
%! mbdyn_pre_beam_write_beams(beam, fd, options_pre);
%! switch (options_pre.node_type)
%!     case "dynamic"
%!         mbdyn_pre_beam_write_bodies(beam, fd, options_pre);
%!     case "static"
%!     otherwise
%!         error("invalid node type \"%s\"", options_pre.node_type);
%! endswitch
%! fputs(fd, "    joint: joint_id_clamp, clamp,\n");
%! fputs(fd, "    node_id_1, node, node;\n");
%! fputs(fd, "    force: force_id_1, absolute, node_id_1 + number_of_nodes - 1,\n");
%! fputs(fd, "        position, reference, node, null,\n");
%! fputs(fd, "        0., 1., 0.,\n");
%! fputs(fd, "        string,\"F1 * ( (Time <= final_time ) * ( Time - initial_time ) / ( final_time - initial_time ) + ( Time > final_time ) )\";\n");
%! fputs(fd, "end: elements;\n");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     fclose(fd);
%!   endif
%! end_unwind_protect
%! options_run.output_file = input_file;
%! mbdyn_solver_run(input_file, options_run);
%! vars = mbdyn_post_load_log_vars(options_run.output_file);
%! filter_force_id = vars.force_id_1;
%! filter_node_id = vars.node_id_1 + param.number_of_nodes - 1;
%! [t, trajectory, deformation, velocity, acceleration, node_id, force, force_id, force_node_id] = mbdyn_post_load_output_struct (options_run.output_file,filter_node_id, filter_force_id);
%! modal = mbdyn_post_load_output_eig(options_run.output_file);
%! options_modal.f_run_mbdyn2easyanim = true;
%! options_modal.f_run_EasyAnim = true;
%! options_modal.scale = 0.5;
%! mbdyn_post_eig_to_mov_file(options_run.output_file, [options_run.output_file, "_%d.mov"], options_modal, modal);
%! [bodies] = mbdyn_post_load_log_body(options_run.output_file);
%! [nodes, dof_info] = mbdyn_post_load_log_node(options_run.output_file);
%! beam_bodies(1).name = 'beam';
%! beam_bodies(1).labels = vars.body_id_1 + ( 1:param.number_of_bodies ) - 1;
%! [inertia] = mbdyn_post_inertia_compute(beam_bodies,bodies,nodes);
%! node_idx_node_1 = find(node_id == filter_node_id);
%! force_idx_1 = find(force_id == filter_force_id);
%! if (f_plot_res)
%!   figure("visible", "off");
%!   plot(t, deformation{node_idx_node_1}(:,1:3), '-;u;');
%!   grid on;
%!   grid minor on;
%!   xlabel('t [s]');
%!   ylabel('u [m]');
%!   title('deformation at the end node');
%! endif
%! f_a = param.F1 * l^3 / (3 * param.E * param.Iz);
%! f = deformation{node_idx_node_1}(end, 2);
%! m_a = param.rho*param.A*l;
%! r_a = d/2;
%! Jx_a = param.rho*param.Ip*l;
%! Jz_a = Jy_a = m_a * (r_a^2/4+l^2/12);
%! assert_simple(f, f_a, 1e-5 * abs(f_a));
%! assert_simple(inertia.dm, m_a, 1e-5 * abs(m_a));
%! assert_simple(inertia.J(1,1), Jx_a, 1e-2 * abs(Jx_a));
%! assert_simple(inertia.J(2,2), Jy_a, 1e-2 * abs(Jy_a));
%! assert_simple(inertia.J(3,3), Jz_a, 1e-2 * abs(Jz_a));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(input_file);
%!     files = dir([input_file, "*"]);
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
