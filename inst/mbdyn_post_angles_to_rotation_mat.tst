%!test
%! state = rand("state");
%! unwind_protect
%! rand("seed", 0);
%! N = 10;
%! orient = {"euler123", "euler321", "phi"};
%! func = {@rotation_matrix_to_euler123, ...
%!         @rotation_matrix_to_euler321, ...
%!         @rotation_matrix_to_rotation_vector};
%! res.node_id = 1:numel(orient);
%! for i=1:numel(orient)
%!   res.t = 1:N;
%!   res.trajectory{i} = [zeros(N, 3), (2 * rand(N, 3) - 1) * 0.5 * pi];
%!   log_dat.nodes(i).label = i;
%!   log_dat.nodes(i).orientation_description = orient{i};
%! endfor
%! R = mbdyn_post_angles_to_rotation_mat(res.node_id, res, log_dat);
%! for i=1:numel(R)
%!   Phi = feval(func{i}, R{i}).';
%!   assert_simple(Phi, res.trajectory{i}(:, 4:6), sqrt(eps) * pi);
%! endfor
%! unwind_protect_cleanup
%! rand("state", state);
%! end_unwind_protect

%!test
%! state = rand("state");
%! unwind_protect
%! rand("seed", 0);
%! N = 10;
%! orient = {"euler123", "euler321", "phi"};
%! func = {@rotation_matrix_to_euler123, ...
%!         @rotation_matrix_to_euler321, ...
%!         @rotation_matrix_to_rotation_vector};
%! res.node_id = 1:numel(orient);
%! for i=1:numel(orient)
%!   res.t = 1:N;
%!   res.trajectory{i} = [zeros(N, 3), (2 * rand(N, 3) - 1) * 0.5 * pi];
%!   log_dat.nodes(i).label = i;
%!   log_dat.nodes(i).orientation_description = orient{i};
%! endfor
%! R = mbdyn_post_angles_to_rotation_mat(res.node_id, res, log_dat, 1:2:N);
%! for i=1:numel(R)
%!   Phi = feval(func{i}, R{i}).';
%!   assert_simple(Phi, res.trajectory{i}(1:2:N, 4:6), sqrt(eps) * pi);
%! endfor
%! unwind_protect_cleanup
%! rand("state", state);
%! end_unwind_protect
