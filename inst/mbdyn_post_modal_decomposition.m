function [r, l, lambda, err, LAMBDA] = mbdyn_post_modal_decomposition(A, B, dCoef, sigma, k, opts)
  n = rows(A);

  if (n ~= columns(A) || n ~= rows(B) || n ~= columns(B))
    error("incompatible matrix size of A and B");
  endif

  pkg load mboct-fem-pkg;

  Bfact = fem_sol_factor(B, opts);

  oper_func = @(X) mbdyn_modal_operator_func(A, Bfact, n, X);

  opts.issym = false;
  opts.isreal = true;
  opts.p = min(2 * k, n);

  rndstate = rand("state");

  unwind_protect
    rand("seed", 0);
    [V, LAMBDA, flag] = eigs(oper_func, 2 * n, k, sigma, opts);
  unwind_protect_cleanup
    rand("state", rndstate);
  end_unwind_protect

  if (flag ~= 0)
    error("eigs failed with status %d", flag);
  endif

  LAMBDA = diag(LAMBDA)(1:end/2);

  r = V(1:n, 1:end/2);
  l = V((n + 1):2 * n, 1:end/2);

  s = zeros(1, columns(r));

  for i=1:columns(r)
    s(i) = l(:, i).' * B * r(:, i);
  endfor

  s = diag(1 ./ sqrt(s));

  r *= s;
  l *= s;

  err = 0;

  for i=1:columns(r)
    v1a = A * r(:, i);
    v2a = LAMBDA(i) * B * r(:, i);
    v1b = A.' * l(:, i);
    v2b = LAMBDA(i) * B.' * l(:, i);
    err = max([err, norm(v1a - v2a) / max(norm(v1a), norm(v2a)), norm(v1b - v2b) / max(norm(v1b),  norm(v2b))]);
  endfor

  lambda = (LAMBDA - 1) ./ (LAMBDA + 1) / dCoef;
endfunction

function Y = mbdyn_modal_operator_func(A, B, n, X)
  r = X(1:n, :);
  l = X((n + 1):2*n, :);

  Y = [B \  (A * r);
       ((l' * A) / B).'];
endfunction

%!test
%! function [R,lambda,f] = damped_eigenmode(M,D,S)
%!   C = [ -M \ D,           -M \ S;
%!         eye(length(M)),   zeros(length(M))];
%!   [R,lambda] = eig(C);
%!   lambda = diag(lambda);
%!   A = [abs(imag(lambda(:))), imag(lambda(:))];
%!   [~, idx] = sortrows(A, [-1,-2]);
%!   lambda = lambda(idx);
%!   R = R(:, idx);
%!   if ( nargout() >= 3 )
%!     f = imag(lambda(2*(1:length(lambda)/2)-1)) / ( 2 * pi );
%!   endif
%! endfunction
%! m1 = 10;
%! m2 = 100;
%! k1 = 10000;
%! k2 = 100000;
%! M = [1 / 3 * m1,        1 / 6 * m1,       0;
%!      1 / 6 * m1, 1 / 3 * (m1 + m2), 1 / 6 * m2;
%!               0,        1 / 6 * m2, 1 / 3 * m2];
%!
%! K = [ k1,      -k1,        0;
%!      -k1,  k1 + k2,      -k2;
%!        0,      -k2,  k1 + k2];
%!
%! beta = 0.2;
%!
%! D = beta * K;
%! f_constraint = true;
%! if (f_constraint)
%!   Phi = [-1, 0, 1];
%!   T = null(Phi);
%! else
%!   Phi = zeros(0, 3);
%!   T = eye(3);
%! endif
%! Mred = T.' * M * T;
%! Kred = T.' * K * T;
%! Dred = T.' * D * T;
%! [PHI_ref, lambda_ref, f_ref] = damped_eigenmode(Mred, Dred, Kred);

%! fmin = 0.1 * max(f_ref);
%! if (~fmin > 0)
%!   error("estimation for dCoef failed");
%! endif
%! dCoef = 5 / (2 * pi * fmin);
%!
%! df_dyP = [-M, zeros(size(M)), zeros(rows(M), rows(Phi));
%!            D, eye(columns(M)), zeros(rows(M), rows(Phi));
%!            zeros(rows(Phi), 2 * columns(M) + rows(Phi))];
%!
%! df_dy = [zeros(size(M)), eye(columns(M)), zeros(rows(M), rows(Phi));
%!          K, zeros(size(M)), -Phi.';
%!          Phi/dCoef, zeros(rows(Phi), columns(M) + rows(Phi))];
%!
%! Jac = @(dCoef) -df_dyP - dCoef * df_dy;
%! A = Jac(-dCoef);
%! B = Jac(dCoef);
%! A = sparse(A);
%! B = sparse(B);
%! sigma = "li";
%! k = 4;
%! opts.solver = "umfpack";
%! opts.refine_max_iter = 10;
%! opts.pre_scaling = true;
%! opts.maxit = 5000;
%! opts.tol = 0;
%! tol = sqrt(eps);
%! [VR, VL, lambda, err, LAMBDA] = mbdyn_post_modal_decomposition(A, B, dCoef, sigma, k, opts);
%! assert(isfinite(err));
%! assert(err < tol);
%! f = imag(lambda) / (2 * pi);
%! assert(VL.' * A * VR, diag(LAMBDA), tol*norm(LAMBDA));
%! assert(VL.' * B * VR, eye(2), tol);
%! assert(lambda, lambda_ref(1:2), tol*norm(lambda_ref(1:2)));
