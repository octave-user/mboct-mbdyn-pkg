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
## @deftypefn {Function File} mbdyn_tests()
##
## Execute all tests of the Octave MBDyn package.
##
## @end deftypefn

function mbdyn_tests()
  pkg_info = pkg("describe", "mboct-mbdyn-pkg");
  num_passed = 0;
  num_total = 0;
  for k=1:numel(pkg_info)
    for j=1:numel(pkg_info{k}.provides)
      func = pkg_info{k}.provides{j}.functions;
      for i=1:numel(func)
        switch (exist(func{i}))
          case 2
            [N, NMAX] = test(func{i}, "quiet");
          case 3
            if (exist([func{i}, ".cpp"]))
              [N, NMAX] = test([func{i}, ".cpp"], "quiet");
            endif
          otherwise
            N = 0;
            NMAX = 0;
        endswitch
        num_passed += N;
        num_total += NMAX;
      endfor
    endfor
  endfor

  fprintf(stderr, "%d of %d tests passed\n", num_passed, num_total);
endfunction


