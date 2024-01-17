## mbdyn_pre_beam_compute.tst:03
%!test
%! close all;
%! N = 4;
%! Theta2 = 88*pi/180;
%! Theta3 = -10*pi/180;
%! X = [ linspace(1,7,N);
%!       zeros(1,N);
%!       zeros(1,N) ];
%! R = [cos(Theta2)*cos(Theta3),-sin(Theta3),sin(Theta2)*cos(Theta3);
%!      cos(Theta2)*sin(Theta3),cos(Theta3),sin(Theta2)*sin(Theta3);
%!      -sin(Theta2),0,cos(Theta2)];
%! n = [ 1, 1;
%!       0, 0;
%!       0, 0 ];
%! beam = mbdyn_pre_beam_compute(R*X,N);
%! figure("visible", "off");
%! mbdyn_pre_beam_plot(beam,struct("s",0.001,"Rn",false,"Rg",false));
%! title("straight beam");
%! assert_simple(R,beam.Rn(:,:,1),sqrt(eps));
