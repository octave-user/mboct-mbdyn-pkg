## Copyright (C) 2015(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{beam}, @var{X}, @var{shape}] = mbdyn_pre_beam_helical_compute(@var{shape_options}, @var{number_of_elements_per_coil}, @var{interpolation_points})
## Build a beam model for a helical spring.
## @seealso{mbdyn_pre_beam_helical_shape}
## @end deftypefn

function [beam, X, shape] = mbdyn_pre_beam_helical_compute(shape_options, number_of_elements_per_coil, interpolation_points)
    [X, shape] = mbdyn_pre_beam_helical_shape(shape_options);
    number_of_coils = abs( shape.Phi(end) - shape.Phi(1) ) / ( 2 * pi );
    N = ceil(double(number_of_elements_per_coil) * number_of_coils);
    beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
endfunction
