## mbdyn_post_rigid_body_kinematics.tst:01
%!test
%! try
%! r1 = 15e-3;
%! omega0 = 100;
%! Phi0 = 30 * pi / 180;
%! t = linspace(0, abs(2 * pi / omega0), 100);
%! omegaP = 10;
%! Phi1 = Phi0 + omega0 * t + 0.5 * omegaP * t.^2;
%! omega1 = omega0 + omegaP * t;
%! a_rad = -r1 * omega1.^2;
%! a_tan = repmat(r1 * omegaP, 1, length(t));
%! v_tan = r1 * omega1;
%! n1 = [1; 0; 0];
%! o1 = [r1; 0; 0];
%! x1 = [r1 * cos(Phi1);
%!       r1 * sin(Phi1);
%!       zeros(1, length(t))];
%! res.t = t.';
%! res.trajectory{1} = [zeros(5, length(t));
%!                      Phi1].';
%! res.velocity{1} = [zeros(5, length(t));
%!                    omega1].';
%! res.acceleration{1} = [zeros(5, length(t));
%!                        repmat(omegaP, 1, length(t))].';
%! res.node_id = 1234;
%! response(1).node_label = res.node_id;
%! response(1).offset = o1;
%! response(1).direction = n1;
%! response(2).node_label = res.node_id;
%! response(2).offset = o1;
%! response(2).direction = [0; 1; 0];
%! response(3).node_label = res.node_id;
%! response(3).offset = o1;
%! response(3).direction = -n1;
%! response(4).node_label = res.node_id;
%! response(4).offset = o1;
%! response(4).direction = -[0; 1; 0];
%! log_dat.nodes(1).label = res.node_id;
%! log_dat.nodes(1).orientation_description = "euler123";
%! log_dat.nodes(1).X0 = zeros(3, 1);
%! log_dat.nodes(1).R0 = eye(3);
%! [a, v, s] = mbdyn_post_rigid_body_kinematics(res, log_dat, response);
%! tol = eps^0.8;
%! assert_simple(a{1}, a_rad.', tol * norm(a_rad));
%! assert_simple(v{1}, zeros(length(res.t), 1), tol * norm(v_tan));
%! assert_simple(a{2}, a_tan.', tol * norm(a_tan));
%! assert_simple(v{2}, v_tan.', tol * norm(v_tan));
%! assert_simple(a{3}, -a_rad.', tol * norm(a_rad));
%! assert_simple(v{3}, zeros(length(res.t), 1), tol * norm(v_tan));
%! assert_simple(a{4}, -a_tan.', tol * norm(a_tan));
%! assert_simple(v{4}, -v_tan.', tol * norm(v_tan));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
