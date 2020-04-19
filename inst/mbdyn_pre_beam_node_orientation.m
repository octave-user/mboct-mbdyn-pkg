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
## @deftypefn {Function File} @var{Rn} = mbdyn_pre_beam_node_orientation(@var{beam}, @var{beam_id})
##
## Returns the orientation matrices of all nodes from beam <@var{beam_id}>.
##
## @end deftypefn

function Rn = mbdyn_pre_beam_node_orientation(beam, beam_id)
  if (nargin ~= 2 || nargout > 1)
    print_usage();
  endif
  
  Rn = beam.Rn(:, :, beam.beams(beam_id).nidx);
endfunction
