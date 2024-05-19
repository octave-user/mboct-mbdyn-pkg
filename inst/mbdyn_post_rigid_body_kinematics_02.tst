## mbdyn_post_rigid_body_kinematics.tst:02
%!test
%! try
%! t = linspace(0, 10, 100);
%! x0 = 1;
%! v0 = 10;
%! a0 = 20;
%! a = repmat(a0, 1, length(t));
%! v = v0 + a0 * t;
%! x = x0 + v0 * t + 0.5 * a0 * t.^2;
%! n1 = [5; 7; 9];
%! n1 /= norm(n1);
%! o1 = [1; 2; 3];
%! res.t = t.';
%! res.trajectory{1} = [n1 * x;
%!                      zeros(3, length(t))].';
%! res.velocity{1} = [n1 * v;
%!                    zeros(3, length(t))].';
%! res.acceleration{1} = [n1 * a;
%!                        zeros(3, length(t))].';
%! res.node_id = 1234;
%! response(1).node_label = res.node_id;
%! response(1).offset = o1;
%! response(1).direction = n1;
%! log_dat.nodes(1).label = res.node_id;
%! log_dat.nodes(1).orientation_description = "euler123";
%! log_dat.nodes(1).X0 = zeros(3, 1);
%! log_dat.nodes(1).R0 = eye(3);
%! [a1, v1, s1] = mbdyn_post_rigid_body_kinematics(res, log_dat, response);
%! tol = eps;
%! assert_simple(a1{1}, a.', tol * norm(a));
%! assert_simple(v1{1}, v.', tol * norm(v));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
