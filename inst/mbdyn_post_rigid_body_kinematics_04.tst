## mbdyn_post_rigid_body_kinematics.tst:04
%!demo
%! try
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
