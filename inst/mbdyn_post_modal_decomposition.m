## Copyright (C) 2011(-2023) Reinhard <octave-user@a1.net>
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*-
## @deftypefn {Function File} [@var{VR}, @var{VL}, @var{lambda}, @var{err}, @var{LAM}] = mbdyn_post_modal_decomposition(@var{A}, @var{B}, @var{dCoef}, @var{sigma}, @var{k}, @var{opts})
##
## Computes eigenvalues and corresponding left and right eigenvectors of the generalized eigenvalue problem
## @var{A} * @var{VR} = @var{LAM} * @var{B} * @var{VR}
## and
## @var{VL}.' * @var{A} = @var{LAM} * @var{VL}.' * @var{B}
## @end deftypefn

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
