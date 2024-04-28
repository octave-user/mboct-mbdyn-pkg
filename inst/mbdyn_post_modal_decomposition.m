function [r, l, lambda, err, LAMBDA] = mbdyn_post_modal_decomposition(A, B, dCoef, sigma, k, opts)
  n = rows(A);

  if (n ~= columns(A) || n ~= rows(B) || n ~= columns(B))
    error("incompatible matrix size of A and B");
  endif

  if (~(isreal(A) && isreal(B)))
    error("A and B must be real matrices");
  endif
  
  pkg load mboct-fem-pkg;

  Bfact = fem_sol_factor(B, opts);

  oper_func = @(X) mbdyn_modal_operator_func(A, Bfact, n, X);

  opts.issym = false;
  opts.isreal = true;

  if (~isfield(opts, "p"))
    opts.p = min(4 * k, n);
  endif

  rndstate = rand("state");

  unwind_protect
    rand("seed", 0);
    [V, LAMBDA, flag] = eigs(oper_func, 2 * n, 2 * k, sigma, opts);
  unwind_protect_cleanup
    rand("state", rndstate);
  end_unwind_protect

  if (flag ~= 0)
    error("eigs failed with status %d", flag);
  endif
  
  LAMBDA = diag(LAMBDA);
  
  r = V(1:n, :);
  l = V((n + 1):2 * n, :);

  s = zeros(1, columns(r));

  for i=1:columns(r)
    s(i) = l(:, i).' * B * r(:, i);
  endfor

  s = diag(1 ./ sqrt(s));

  r *= s;
  l *= s;

  ortho = true(1, columns(r));

  if (~isfield(opts, "tol_ortho"))
    opts.tol_ortho = eps^0.5;
  endif
  
  for i=1:columns(ortho)
    for j=i + 1:columns(ortho)
      if (ortho(j))
        ortho(j) = abs(l(:, i).' * B * r(:, j)) < opts.tol_ortho;
      endif
    endfor
  endfor

  idx = find(ortho);

  LAMBDA = LAMBDA(idx);
  
  r = r(:, idx);
  l = l(:, idx);
  
  err = 0;

  for i=1:columns(r)
    v1a = A * r(:, i);
    v2a = LAMBDA(i) * B * r(:, i);
    v1b = A.' * l(:, i);
    v2b = LAMBDA(i) * B.' * l(:, i);
    err = max([err, norm(v1a - v2a) / max(norm(v1a), norm(v2a)), norm(v1b - v2b) / max(norm(v1b),  norm(v2b))]);
  endfor

  lambda = (LAMBDA - 1) ./ (LAMBDA + 1) / dCoef;

  O = [abs(imag(lambda)), imag(lambda), real(lambda)];

  [~, idx] = sortrows(O, [1, -2, 3]);

  lambda = lambda(idx);
  LAMBDA = LAMBDA(idx);
  l = l(:, idx);
  r = r(:, idx);
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
%!   [~, idx] = sortrows(A, [1, -2]);
%!   lambda = lambda(idx);
%!   R = R(:, idx);
%!   if ( nargout() >= 3 )
%!     f = imag(lambda(2*(1:length(lambda)/2)-1)) / ( 2 * pi );
%!   endif
%! endfunction
%! n = 8;
%! elno = [1, 2;
%!         2, 3;
%!         3, 4,
%!         4, 5;
%!         5, 6;
%!         6, 7;
%!         7, 8;
%!         8, 0];
%! elm = [1000;
%!        1000;
%!        1000;
%!        1000;
%!        1000;
%!        1000;
%!        1000;
%!        1000];
%! elk = [10000;
%!        10000;
%!        10000;
%!        10000;
%!        10000;
%!        10000;
%!        10000;
%!        10000];
%! M = zeros(n, n);
%! K = zeros(n, n);
%! for i=1:rows(elno)
%!   Mi = [1 / 3 * elm(i), 1 / 6 * elm(i);
%!         1 / 6 * elm(i), 1 / 3 * elm(i)];
%!   Ki = [elk(i), -elk(i);
%!         -elk(i), elk(i)];
%!   for j=1:2
%!     for k=1:2
%!       if (elno(i, j) > 0 && elno(i, k) > 0)
%!       M(elno(i, j), elno(i, k)) += Mi(j, k);
%!       K(elno(i, j), elno(i, k)) += Ki(j, k);
%!       endif
%!     endfor
%!   endfor
%! endfor
%! assert(isdefinite(M));
%! assert(isdefinite(K));
%!
%! beta = 0.05;
%!
%! D = beta * K;
%! f_constraint = true;
%! if (f_constraint)
%!   Phi = zeros(1, n);
%!   Phi(1:2) = [1, -1];
%!   T = null(Phi);
%! else
%!   Phi = zeros(0, n);
%!   T = eye(3);
%! endif
%! Mred = T.' * M * T;
%! Kred = T.' * K * T;
%! Dred = T.' * D * T;
%! [PHI_ref, lambda_ref, f_ref] = damped_eigenmode(Mred, Dred, Kred);

%! fmin = 0.1 * min(f_ref);
%! if (~(fmin > 0))
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
%! k = 14;
%! opts.p = 2 * columns(A);
%! opts.solver = "umfpack";
%! opts.refine_max_iter = 20;
%! opts.pre_scaling = true;
%! opts.maxit = 5000;
%! opts.tol = 0;
%! tol = sqrt(eps);
%! [VR, VL, lambda, err, LAMBDA] = mbdyn_post_modal_decomposition(A, B, dCoef, sigma, k, opts);
%! assert(isfinite(err));
%! assert(err < tol);
%! f = imag(lambda) / (2 * pi);
%! assert(VL.' * A * VR, diag(LAMBDA), tol*norm(LAMBDA));
%! assert(VL.' * B * VR, eye(columns(VL)), tol);
%! assert(lambda, lambda_ref(1:rows(lambda)), tol*norm(lambda_ref(1:rows(lambda))));
