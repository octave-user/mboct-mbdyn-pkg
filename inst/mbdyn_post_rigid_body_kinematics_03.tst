## mbdyn_post_rigid_body_kinematics.tst:03
%!test
%! try
%! r1 = 15e-3;
%! omega0 = 100;
%! Phi0 = 30 * pi / 180;
%! t = linspace(0, abs(2 * pi / omega0), 100);
%! omegaP = 500;
%! Phi1 = Phi0 + omega0 * t + 0.5 * omegaP * t.^2;
%! omega1 = omega0 + omegaP * t;
%! o1 = [r1; 0; 0];
%! x1 = r1 * [cos(Phi1);
%!            sin(Phi1);
%!            zeros(1, length(t))];
%! xP1 = r1 * [-sin(Phi1) .* omega1;
%!              cos(Phi1) .* omega1;
%!              zeros(1, length(t))];
%! xPP1 = r1 * [-cos(Phi1) .* omega1.^2 - sin(Phi1) * omegaP;
%!              -sin(Phi1) .* omega1.^2 + cos(Phi1) * omegaP;
%!              zeros(1, length(t))];
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
%! log_dat.nodes(1).label = res.node_id;
%! log_dat.nodes(1).orientation_description = "euler123";
%! log_dat.nodes(1).X0 = x1(:, 1)*2;
%! log_dat.nodes(1).R0 = euler123_to_rotation_matrix([0; 0; 2*Phi0]);
%! [a, v, s] = mbdyn_post_rigid_body_kinematics(res, log_dat, response);
%! tol = eps^0.8;
%! assert_simple(s{1}, (x1 - log_dat.nodes(1).X0 - log_dat.nodes(1).R0 * o1).', tol * norm(x1));
%! assert_simple(v{1}, xP1.', tol * norm(xP1));
%! assert_simple(a{1}, xPP1.', tol * norm(xPP1));
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
