%!test
%! save_env = getenv("MBD_NUM_THREADS");
%! unwind_protect
%!   putenv("MBD_NUM_THREADS", "4");
%!   assert_simple(mbdyn_solver_num_threads_default(), int32(4));
%! unwind_protect_cleanup
%!   putenv("MBD_NUM_THREADS", save_env);
%! end_unwind_protect
