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
## @deftypefn {Function File} mbdyn_pre_beam_check_rotation_matrix(@var{R})
##
## Internal helper function for checking the rotation matrix of a curved beam.
##
## @end deftypefn

function mbdyn_pre_beam_check_rotation_matrix(R)
  if (nargin ~= 1)
    print_usage();
  endif
  
  if (~min(min(isfinite(R))))
    error("rotation matrix is not finite");
  endif
  
  f1 = max(max(abs(R.' * R - eye(3))));
  f2 = max(max(abs(R * R.' - eye(3))));
  
  if (max(f1, f2) > sqrt(eps))
    error("rotation matrix is not orthogonal");
  endif
endfunction
