## Copyright (C) 2020(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{x}] = mbdyn_post_ehd_interp_grid(@var{bearing}, @var{node_type})
## Build a grid for interpolation of nodal data from MBDyn hydrodynamic bearings.
##
## @var{bearing} @dots{} Return value from mbdyn_post_load_log.
##
## @var{node_type} @dots{} One of "hydro", "thermal", "flux_x", "flux_z"
##
## @seealso{mbdyn_post_ehd_load_output, mbdyn_post_load_log}
## @end deftypefn

function x = mbdyn_post_ehd_interp_grid(bearing, node_type)
  if (nargin < 1 || nargin > 2 || nargout > 1)
    print_usage();
  endif

  if (nargin < 2)
    node_type = "hydro";
  endif
  
  inum_nodes_corner = int32(0);
  x = zeros(2, length(bearing.nodes));
  
  for k=1:columns(x)
    switch (bearing.nodes(k).type)
      case node_type
        x(:, ++inum_nodes_corner) = bearing.nodes(k).x;
    endswitch
  endfor

  x = x(:, 1:inum_nodes_corner);
endfunction
