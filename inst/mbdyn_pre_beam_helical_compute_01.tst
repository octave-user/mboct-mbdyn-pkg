## mbdyn_pre_beam_helical_compute.tst:01
%!test
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
%! nfig = figure("visible", "off");
%! hold on;
%! mbdyn_pre_beam_plot(beam,struct("Rn",true,"Rg",false,"s",options.d,"figure",nfig));
%! plot3(X(1,:),X(2,:),X(3,:),'-;X;b');
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
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! set(gca(),"DataAspectRatio",[1,1,1]);
%! endif