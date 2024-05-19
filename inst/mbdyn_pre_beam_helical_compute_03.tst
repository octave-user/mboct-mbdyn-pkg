## mbdyn_pre_beam_helical_compute.tst:03
%!test
%! try
%! f_plot = false;
%! if (f_plot)
%!   close all;
%! endif
%! options.d = 1.5e-3;
%! options.na = 1.2;
%! options.ni = 3;
%! options.L = 30e-3;
%! options.Di = 10.6e-3;
%! options.N_coil = 20;
%! options.type = "const pitch, const diameter";
%! number_of_elements_per_coil = 7;
%! interpolation_points = 1;
%! [beam, X, shape] = mbdyn_pre_beam_helical_compute(options, number_of_elements_per_coil, interpolation_points);
%! if (f_plot)
%! figure("visible","off");
%! hold on;
%! plot(180 / pi * shape.Phi, shape.z,'-;z(Phi);r');
%! grid on;
%! grid minor on;
%! xlabel('Phi [deg]');
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
%! figure_list();
%! endif
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
