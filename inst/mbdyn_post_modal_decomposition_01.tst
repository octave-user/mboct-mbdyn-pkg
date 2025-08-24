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
