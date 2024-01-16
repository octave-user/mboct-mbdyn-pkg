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
## @deftypefn {Function File} [@var{node_id}, @var{prm_value}] = mbdyn_post_load_output_prm(@var{mbdyn_filename}, @var{filter_node_id})
##
## Loads the MBDyn .prm file named "<@var{mbdyn_filename}.prm>".
##
## @var{filter_node_id} @dots{} Load only nodes with node numbers in array <@var{filter_node_id}>. If empty, load everything.
##
## @var{node_id} @dots{} Node number of the loaded nodes.
##
## @end deftypefn

function [node_id, prm_value] = mbdyn_post_load_output_prm(mbdyn_filename, filter_node_id)
  if (nargin < 1 || nargout > 2)
    print_usage();
  endif

  if (nargin < 2)
    filter_node_id = [];
  endif

  prm_filename = mbdyn_post_output_filename(mbdyn_filename, ".prm");

  column_count = 1;

  [node_id, prm_value] = mbdyn_post_load_output(prm_filename, column_count, filter_node_id);
endfunction
