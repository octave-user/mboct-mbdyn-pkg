## mbdyn_pre_beam_helical_shape.tst:03
%!test
%! f_plot = false;
%! if (f_plot)
%!   close all;
%! endif
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
%! if (f_plot)
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);r');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);r');
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
%! endif