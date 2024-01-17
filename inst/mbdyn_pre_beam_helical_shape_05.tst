## mbdyn_pre_beam_helical_shape.tst:05
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
