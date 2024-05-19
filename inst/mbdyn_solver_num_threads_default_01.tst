## mbdyn_solver_num_threads_default.tst:01
%!test
%! try
%! save_env = getenv("MBD_NUM_THREADS");
%! unwind_protect
%!   putenv("MBD_NUM_THREADS", "4");
%!   assert_simple(mbdyn_solver_num_threads_default(), int32(4));
%! unwind_protect_cleanup
%!   putenv("MBD_NUM_THREADS", save_env);
%! end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
