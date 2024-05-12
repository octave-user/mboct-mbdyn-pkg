## mbdyn_pre_beam_helical_shape.tst:04
%!test
%! try
%! f_plot = false;
%! if (f_plot)
%! close all;
%! endif
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
%! if (f_plot)
%! mbdyn_pre_beam_plot(beam,struct("Rn",true, "Rg",false, "s",options.d));
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);r');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! plot(shape.Phi, shape.z, 'x;z_i(Phi_i);b');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);r');
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
%! endif
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
