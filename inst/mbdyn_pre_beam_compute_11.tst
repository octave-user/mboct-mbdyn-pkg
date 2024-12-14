## mbdyn_pre_beam_compute.tst:01
%!test
%! close all;
%! try
%! N1 = 10;
%! N2 = 10;
%! N3 = 10;
%! N = 10;
%! x1 = linspace(0, 2, N1);
%! y1 = zeros(1, N1);
%! z1 = zeros(1, N1);
%! x2 = repmat(2.5, 1, N2);
%! y2 = linspace(0.5, 2.5, N2);
%! z2 = zeros(1, N2);
%! x3 = repmat(2.5, 1, N3);
%! y3 = repmat(3, 1, N3);
%! z3 = linspace(0.5, 2.5, N3);
%! X = [x1, x2, x3;
%!      y1, y2, y3;
%!      z1, z2, z3];
%! opts.smooth_curvature = true;
%! beam = mbdyn_pre_beam_compute(X, N, 1, opts);
%! assert(max(abs(beam.cosPhi - 1)) < sqrt(eps));
%! opts.Rn = true;
%! opts.Rg = true;
%! mbdyn_pre_beam_plot(beam, opts);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch


## mbdyn_pre_beam_compute.tst:01
%!test
%! close all;
%! try
%! N1 = 10;
%! N2 = 10;
%! N3 = 10;
%! N = 5;
%! x1 = linspace(0, 2, N1);
%! y1 = 0.05 * sin(x1 / 2 * 2 * pi);
%! z1 = zeros(1, N1);
%! X = [x1;
%!      y1;
%!      z1];
%! opts.smooth_curvature = true;
%! beam = mbdyn_pre_beam_compute(X, N, 1, opts);
%! assert(max(abs(beam.cosPhi - 1)) < sqrt(eps));
%! opts.Rn = true;
%! opts.Rg = true;
%! mbdyn_pre_beam_plot(beam, opts);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch

%!test
%! close all;
%! try
%! N1 = 10;
%! N2 = 10;
%! N3 = 10;
%! N = 5;
%! x1 = linspace(0, 2, N1);
%! y1 = 0.05 * sin(x1 / 2 * 2 * pi);
%! z1 = 0.05 * cos(x1 / 2 * 4 * pi);
%! X = [x1;
%!      y1;
%!      z1];
%! opts.smooth_curvature = true;
%! beam = mbdyn_pre_beam_compute(X, N, 1, opts);
%! assert(max(abs(beam.cosPhi - 1)) < sqrt(eps));
%! opts.Rn = true;
%! opts.Rg = true;
%! mbdyn_pre_beam_plot(beam, opts);
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
