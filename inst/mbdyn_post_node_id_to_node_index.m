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
## @deftypefn {Function File} @var{node_idx} = mbdyn_post_node_id_to_node_index(@var{node_id}, @var{vars}, @var{node_id_prefix}, @var{node_idx_prefix}, @var{node_idx_in})
##
## Converts node id's to node indices (e.g. @var{vars}.@var{node_id_}* => @var{vars}.@var{node_idx_}*).
##
## @end deftypefn

function node_idx = mbdyn_post_node_id_to_node_index(node_id, vars, node_id_prefix = "node_id_", node_idx_prefix = "node_idx_", node_idx_in)
  if (nargin >= 5)
    node_idx = node_idx_in;
  else
    node_idx = struct();
  endif

  names = fieldnames(vars);

  for i=1:length(names)
    if (strncmp(names{i}, node_id_prefix, length(node_id_prefix)))
      id_value = getfield(vars, names{i});
      idx_value = find(node_id == id_value);
      if (length(idx_value) > 0)
        idx_name = strcat(node_idx_prefix, names{i}((length(node_id_prefix) + 1):end));
        node_idx = setfield(node_idx, idx_name, idx_value);
      endif
    endif
  endfor
endfunction
