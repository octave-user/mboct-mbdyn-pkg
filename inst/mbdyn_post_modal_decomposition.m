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
  #lambda = (diag(LAMBDA) - 1) ./ (diag(LAMBDA) + 1) / dCoef
  ## for i=1:columns(V)
  ##   oper_func(V(:, i)) - LAMBDA(i,i) * V(:, i)
  ## endfor

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
%!   [lambda_imag,i_lambda] = sort(imag(lambda),'descend');
%!   lambda = lambda(i_lambda);
%!   R = R(:,i_lambda);
%!   [lambda_imag,i_lambda]     = sort(abs(imag(lambda)));
%!   lambda = lambda(i_lambda);
%!   R = R(:,i_lambda);
%!   if ( nargout() >= 3 )
%!     f = imag(lambda(2*(1:length(lambda)/2)-1)) / ( 2 * pi );
%!   endif
%! endfunction
%!
%! M = [10.5,     0, 0;
%!          0, 20.8, 0;
%!          0,    0,  120];
%!
%! K = [ 1000, 0,     0;
%!          0,  2000,   0;
%!          0,  0,  2000];
%!
%! beta = 0.1;
%!
%! D = beta * K;
%! f_constraint = false;
%! if (f_constraint)
%!   Phi = [-1, 1, 0];
%!   T = [Phi;
%!        zeros(1, 2), 1].';
%! else
%!   Phi = zeros(0, 3);
%!   T = eye(3);
%! endif
%! Mred = T.' * M * T;
%! Kred = T.' * K * T;
%! Dred = T.' * D * T;
%! [PHI_ref, lambda_ref, f_ref] = damped_eigenmode(Mred, Dred, Kred);
%! fmin = 0.1 * max(f_ref);
%! dCoef = 5 / (2 * pi * fmin);
%!
%! df_dyP = [-M, zeros(3, 3), zeros(3, rows(Phi));
%!            D, eye(3), zeros(3, rows(Phi));
%!            zeros(rows(Phi), 6 + rows(Phi))];
%!
%! df_dy = [zeros(3, 3), eye(3), zeros(3, rows(Phi));
%!          K, zeros(3, 3), -Phi.';
%!          Phi/dCoef, zeros(rows(Phi), 3 + rows(Phi))];
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
%! opts.pre_scaling = false;
%! opts.maxit = 5000;
%! opts.tol = 0;
%! tol = sqrt(eps);
%! [VR, VL, lambda, err, LAMBDA] = mbdyn_post_modal_decomposition(A, B, dCoef, sigma, k, opts);
%! assert(isfinite(err));
%! assert(err < tol);
%! f = imag(lambda) / (2 * pi);
%! assert(VL.' * A * VR, diag(LAMBDA), tol*norm(LAMBDA));
%! assert(VL.' * B * VR, eye(2), tol);
## A = [Jp, zeros(size(Jp));
##      zeros(size(Jp)), Jp.'];
## B = [Jm, zeros(size(Jm));
##      zeros(size(Jm)), Jm.'];
## [V2, LAMBDA2] = eig(B, A);
## LAMBDA2 = diag(LAMBDA2);
## lambda2 = (LAMBDA2 - 1) ./ (LAMBDA2 + 1) / dCoef;
## [~, idx] = sort(imag(lambda2));
## LAMBDA2 = LAMBDA2(idx);
## V2 = V2(:, idx);
## n = rows(Jp);
## VR2 = V2(1:n, 1:4);
## VL2 = V2((n+1):(2*n), 1:4);
## s = zeros(1, columns(VR2));
## for i=1:columns(VR2)
##  s(i) = VL2(:, i).' * Jm * VR2(:, i);
## endfor
## VR2 *= diag(s);
## VL2 *= diag(s);
## VL2.'*Jp*VR2
