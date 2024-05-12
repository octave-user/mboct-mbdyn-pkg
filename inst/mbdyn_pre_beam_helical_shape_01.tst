## mbdyn_pre_beam_helical_shape.tst:01
%!test
%! try
%! f_plot = false;
%! if (f_plot)
%!   close all;
%! endif
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
%! if (f_plot)
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);r');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! plot(options.Phi, options.z, 'x;z_i(Phi_i);b');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);r');
%! plot(options.z_D, options.D, 'x;D_i(z_i);b');
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
%! endif
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
