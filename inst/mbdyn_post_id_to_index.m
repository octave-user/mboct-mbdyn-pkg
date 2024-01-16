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
## @deftypefn {Function File} @var{vars} = mbdyn_post_id_to_index(@var{res},@var{vars_in})
##
## Convert node and element id numbers to index values (e.g. vars.node_id_* => vars.node_idx_*).
##
## @var{res} @dots{} Data structure with id numbers to be converted.
##
## @var{res}.node_id @dots{} Vector of node id's.
##
## @var{res}.force_id @dots{} Vector of force id's.
##
## @var{res}.torque_id @dots{} Vector of torque id's.
##
## @var{res}.joint_id @dots{} Vector of joint id's.
##
## @var{res}.elem_id @dots{} Vector of element id's.
##
## @var{vars_in} @dots{} Data structure returned from mbdyn_post_load_log_vars.
##
## @var{vars} @dots{} Data structure with the same information like <@var{vars_in}>
## plus the requested indices (e.g. node_id_xxx -> node_idx_xxx).
##
## @end deftypefn

function vars = mbdyn_post_id_to_index(res, vars_in)
  if (nargin ~= 2 || nargout > 1)
    print_usage();
  endif

  vars = vars_in;

  id_names = {"node",       "node";
              "abs_node",   "abs_node";
              "elec_node",  "elec_node";
              "prm_node",   "prm_node";
              "hyd_node",   "hyd_node";
              "therm_node", "therm_node";
              "force",      "force";
              "force",      "torque";
              "joint",      "joint";
              "elem",       "elem";
              "genel",      "genel";
              "drive",      "drive";
              "trace",      "trace"};

  for i=1:rows(id_names)
    if (isfield(res, [id_names{i, 1}, "_id"]))
      vars = mbdyn_post_node_id_to_node_index(getfield(res, [id_names{i, 1}, "_id"]), vars, [id_names{i, 2}, "_id_"], [id_names{i, 2}, "_idx_"], vars);
    endif
  endfor
endfunction
