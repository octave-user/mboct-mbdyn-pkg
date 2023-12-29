## Copyright (C) 2014(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{y} = mbdyn_post_frequency_response_helper(@var{i}, @var{df_dy}, @var{df_dy_dot}, @var{Tp}, @var{p}, @var{Tu}, @var{omega}, @var{options})
##
## Internal helper function for mbdyn_post_frequency_response.
##
## @seealso{mbdyn_post_frequency_response}
## @end deftypefn

function y = mbdyn_post_frequency_response_helper(i, df_dy, df_dy_dot, Tp, p, Tu, omega, options)
  if (nargin ~= 8)
    print_usage();
  endif

  if (columns(p) > 1)
    p = p(:, i); ## reduce memory usage: expand only one load vector per time
  endif

  b = Tp.' * p;

  A = df_dy + 1j * omega(i) * df_dy_dot;

  y = Tu * options.solver_func(A, b, options);
endfunction
