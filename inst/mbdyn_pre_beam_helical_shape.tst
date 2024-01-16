%!test
%! close all;
%! Di = [ 10.5e-3; 3e-3; 10.5e-3 ];
%! d = 1.5e-3;
%! n = 7.5;
%! options.L = 30e-3;
%! options.D = Di + d;
%! options.z = options.L * [ 0; 0.05; 0.95; 1 ];
%! options.z_D = options.L * [ 0; 0.5; 1 ];
%! options.Phi = 2 * pi * n * [ 0; 0.2; 0.8; 1 ];
%! options.N = 3600;
%! options.interpolation = 'linear';
%! options.interpolation_D = 'spline';
%! [X, shape] = mbdyn_pre_beam_helical_shape(options);
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! plot(options.Phi, options.z, 'x;z_i(Phi_i);3');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! plot(options.z_D, options.D, 'x;D_i(z_i);3');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! daspect(ones(1,3));

%!test
%! close all;
%! Di = 10.5e-3;
%! options.d = 1.5e-3;
%! options.na = 2;
%! options.ni = 2;
%! options.L = 30e-3;
%! options.D = Di + options.d;
%! options.N = 3600;
%! options.type = "inactive coils, const pitch, const diameter";
%! number_of_elements_per_coil = 7;
%! N = ceil(number_of_elements_per_coil*options.na + 2 * options.ni);
%! interpolation_points = 1;
%! [X, shape] = mbdyn_pre_beam_helical_shape(options);
%! beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
%! mbdyn_pre_beam_plot(beam,struct("Rn",true, "Rg",false, "s",options.d));
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! daspect(ones(1,3));

%!test
%! close all;
%! Di = 10.5e-3;
%! options.d = 1.5e-3;
%! options.na = 2;
%! options.ni = 3-270/360;
%! options.L = 30e-3;
%! options.D = Di + options.d;
%! options.N = 3600;
%! options.N_begin = 10;
%! options.N_end = 10;
%! options.type = "const pitch, const diameter";
%! number_of_elements_per_coil = 7;
%! N = ceil(number_of_elements_per_coil*options.na);
%! interpolation_points = 1;
%! [X, shape] = mbdyn_pre_beam_helical_shape(options);
%! beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! nfig = figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! mbdyn_pre_beam_plot(beam,struct("Rn",true, "Rg",false, "s",options.d, "figure",nfig));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! daspect(ones(1,3));

%!test
%! close all;
%! options.d = 1.5e-3;
%! options.na = 7.25;
%! options.ni = 0.25;
%! options.L = 30e-3;
%! options.Di = [12.45e-3, 10.45e-3, 12.45e-3 ];
%! options.type = "const pitch, 3 diameters";
%! number_of_elements_per_coil = 7;
%! N = ceil(number_of_elements_per_coil*options.na);
%! interpolation_points = 1;
%! [X, shape] = mbdyn_pre_beam_helical_shape(options);
%! beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
%! mbdyn_pre_beam_plot(beam,struct("Rn",true, "Rg",false, "s",options.d));
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! plot(shape.Phi, shape.z, 'x;z_i(Phi_i);3');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! daspect(ones(1, 3));

%!test
%! close all;
%! options.d = 1.6e-3;
%! options.na = 6;
%! options.ni = 3;
%! options.L = 34.7e-3;
%! options.Di = 12.4e-3;
%! options.z_rel = [ 0
%! 0.125106383
%! 0.301702128
%! 0.513191489
%! 0.726382979
%! 0.892340426
%! 1];
%! options.interpolation = "linear";
%! options.type = "relative variable pitch, const diameter";
%! number_of_elements_per_coil = 7;
%! N = ceil(number_of_elements_per_coil*options.na);
%! interpolation_points = 1;
%! [X, shape] = mbdyn_pre_beam_helical_shape(options);
%! beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
%! mbdyn_pre_beam_plot(beam,struct("Rn",true, "Rg",false, "s",options.d));
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! % plot(shape.Phi, shape.z, 'x;z_i(Phi_i);3');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! daspect(ones(1, 3));

%!test
%! close all;
%! param.shape.d = 1.6e-3;
%! param.shape.na = 5.5;
%! param.shape.ni = 3;
%! param.shape.L = 34.48e-3;
%! param.shape.Di = [13.45e-3, 12.45e-3, 13.45e-3];
%! param.shape.z_rel_D = [0, 0.5, 1];
%! param.shape.z_rel = [ 0, 0.38, 0.62, 1 ];
%! param.shape.interpolation_z = "linear";
%! param.shape.interpolation_D = "linear";
%! param.shape.type = "relative variable pitch, variable diameter";
%! param.shape.inactive_coil_alignment = false;
%! number_of_elements_per_coil = 7;
%! N = ceil(number_of_elements_per_coil*param.shape.na);
%! interpolation_points = 1;
%! [X, shape] = mbdyn_pre_beam_helical_shape(param.shape);
%! beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
%! mbdyn_pre_beam_plot(beam,struct("Rn",true, "Rg",false, "s",param.shape.d));
%! fd = -1;
%! unwind_protect
%! [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_beam_helical_shape_XXXXXX"), true);
%! mbdyn_pre_beam_print_nodes(beam, fd);
%! unwind_protect_cleanup
%! fclose(fd);
%! assert_simple(unlink(fname), 0);
%! end_unwind_protect
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! daspect(ones(1, 3));

%!test
%! close all;
%! spring.d = 2.6e-3;
%! spring.L = 31.4e-3;
%! spring.Di = [ 13.55e-3; 13.55e-3; 13.9e-3; 13.9e-3; 13.55e-3; 13.55e-3];
%! spring.na = 4.6;
%! spring.ni = 2.7;
%! spring.z_D = [ -spring.d; 0; 6.5e-3; spring.L - 6.5e-3; spring.L; spring.L + spring.d ];
%! spring.N = int32(3600);
%! spring.interpolation_D = 'linear';
%! spring.type = "inactive coils, const pitch, variable diameter";
%! spring.ground_surface = 340 * pi / 180;  # angle of ground area at each side
%! [X, shape] = mbdyn_pre_beam_helical_shape(spring);
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! daspect(ones(1, 3));

%!demo
%!demo
%! close all;
%! spring.d = 2.6e-3;
%! spring.L = 31.4e-3;
%! spring.Di = [ 13.55e-3; 13.55e-3; 13.9e-3; 13.9e-3; 13.55e-3; 13.55e-3];
%! spring.na = 4.6;
%! spring.ni = 2.7;
%! spring.z_D = [ -spring.d; 0; 6.5e-3; spring.L - 6.5e-3; spring.L; spring.L + spring.d ];
%! spring.N = int32(3600);
%! spring.interpolation_D = 'linear';
%! spring.type = "inactive coils, const pitch, variable diameter";
%! spring.ground_surface = 340 * pi / 180;  # angle of ground area at each side
%! [X, shape] = mbdyn_pre_beam_helical_shape(spring);
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:), X(2, :), X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! daspect(ones(1, 3));
