## mbdyn_pre_beam_compute.tst:02
%!test
%! N = 100;
%! X = [ 1,2,3,4;
%!       0,3,0,0;
%!       0,4,0,0 ];
%! beam = mbdyn_pre_beam_compute(X, N, 40);
%! norm_dXn = norm(beam.Xn(:,2:end) - beam.Xn(:,1:end-1), 2, 'cols');
%! f = max(abs(1 - norm_dXn/mean(norm_dXn)));
%! assert_simple(f < 1e-2);
