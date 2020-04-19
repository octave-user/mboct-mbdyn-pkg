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
## @deftypefn {Function File} @var{Rg} = mbdyn_pre_beam_gauss_point_orientation(@var{beam}, @var{beam_id})
##
## Return the rotation matrix at Gauss points for beam <@var{beam_id}>.
##
## @end deftypefn

function Rg = mbdyn_pre_beam_gauss_point_orientation(beam, beam_id)
  if (nargin ~= 2 || nargout > 1)
    print_usage();
  endif
  
  Rg = beam.Rg(:, :, beam.beams(beam_id).gidx);
endfunction
