## mbdyn_pre_write_param_file.tst:03
%!test
%! try
%! fd = -1;
%! %unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_write_param_file_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     param.x = -1.5;
%!     param.i = int32(2);
%!     param.b1 = true;
%!     param.b2 = false;
%!     param.s = "123 456 789";
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "   initial time: 0;\n");
%!     fputs(fd, "   final time: 1;\n");
%!     fputs(fd, "   time step: 0.1;\n");
%!     fputs(fd, "   nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "   abstract nodes: 1;\n");
%!     fputs(fd, "   genels: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, " abstract: 1, differential;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, " genel: 1, mass, 1, abstract, algebraic, 1.;\n");
%!     fputs(fd, " genel: 2, spring support, 1, abstract, algebraic, 1.;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!     opts.output_file = fname;
%!     opts.verbose = false;
%!     opts.logfile = [fname, ".stdout"];
%!     mbdyn_solver_run(fname, opts);
%!     log_dat = mbdyn_post_load_log(opts.output_file);
%!     assert_simple(log_dat.vars.x, param.x);
%!     assert_simple(log_dat.vars.i, param.i);
%!     assert_simple(log_dat.vars.b1, param.b1);
%!     assert_simple(log_dat.vars.b2, param.b2);
%!     assert_simple(log_dat.vars.s, param.s);
%! %unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!   endif
%!   files = dir([fname, ".*"]);
%!   for i=1:numel(files)
%!     unlink(fullfile(files(i).folder, files(i).name));
%!   endfor
%! %end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
