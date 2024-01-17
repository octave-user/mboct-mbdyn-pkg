## mbdyn_pre_beam_compute.tst:01
%!test
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
%! beam = mbdyn_pre_beam_compute(R * X, N);
%! assert_simple(R, beam.Rn(:,:,1), sqrt(eps));
