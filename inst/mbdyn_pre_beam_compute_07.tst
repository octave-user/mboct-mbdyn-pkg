## mbdyn_pre_beam_compute.tst:07
%!test
%! f_plot = false;
%! if (f_plot)
%!   close all;
%! endif
%! N = 10;
%! X = [ 1,2,3,4;
%!       0,3,0,0;
%!       0,4,0,0 ];
%! beam = mbdyn_pre_beam_compute(X,N,50);
%! if (f_plot)
%! figure("visible", "off");
%!   plot(beam.ti,beam.si,'-;s(t);r');
%!   xlabel('t []');
%!   ylabel('s [m]');
%!   title('curve length versus parameter');
%!   grid on;
%!   grid minor on;
%! endif
%! norm_dXn = norm(beam.Xn(:,2:end) - beam.Xn(:,1:end-1),2,'cols');
%! f = max(abs(1-norm_dXn/mean(norm_dXn)));
%! if (f_plot)
%! figure("visible","off");
%! stem(norm_dXn,'o-;norm(dXn);r');
%! xlabel('node #');
%! ylabel('norm(dXn) [m]');
%! endif
%! assert_simple(f < 0.18);
%! mbdyn_pre_beam_plot(beam,struct("s",0.1,"X",true,"Rn",true,"Rg",true));
