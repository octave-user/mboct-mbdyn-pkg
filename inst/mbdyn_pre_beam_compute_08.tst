## mbdyn_pre_beam_compute.tst:08
%!test
%! try
%! f_plot = false;
%! if (f_plot)
%!   close all;
%! endif
%! r = 0.5;
%! v = 0.5 * r;
%! h = 2;
%! n = 1;
%! N = 4;
%! k = 200;
%! Phi = linspace(0,2*pi*n,k*N+1);
%! X = [ r*cos(Phi);
%!       r*sin(Phi);
%!       h * Phi./(2*pi) ];
%!
%! Xt = [ -r*sin(Phi);
%!         r*cos(Phi);
%!         repmat(h / (2 * pi),1,length(Phi)) ];
%!
%! t = 2*pi*n/N*Xt(:,[1,end]);
%! beam = mbdyn_pre_beam_compute(X,N);
%! if (f_plot)
%!   figure("visible", "off");
%!   hold on;
%!   plot3(X(1,1:k:k*N+1),X(2,1:k:k*N+1),X(3,1:k:k*N+1),'x;curve;b');
%!   set(plot3(X(1,:),X(2,:),X(3,:),'-;curve;b'),'linewidth',1);
%!   set(gca(),'dataaspectratio',[1,1,1]);
%!   xlabel('x');
%!   ylabel('y');
%!   zlabel('z');
%!   grid on;
%!   grid minor on;
%!   title('helix');
%! endif
%! mbdyn_pre_beam_plot(beam,struct("s",0.08,"Rn",true,"Rg",true));
%! if (f_plot)
%!   hold on;
%!   plot3(X(1,:),X(2,:),X(3,:),'--;curve;k');
%! endif
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
