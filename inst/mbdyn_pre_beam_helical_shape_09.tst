## mbdyn_pre_beam_helical_shape.tst:09
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
