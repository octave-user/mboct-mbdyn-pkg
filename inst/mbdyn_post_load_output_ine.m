## Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{node_id}, @var{beta}, @var{gamma}, @var{beta_dot}, @var{gamma_dot}] = mbdyn_post_load_output_ine(@var{mbdyn_filename}, @var{filter_node_id})
##
## Loads the MBDyn .ine file named "<@var{mbdyn_filename}.ine>".
##
## @var{node_id} @dots{} The node identifier
##
## @var{beta} @dots{} Three components of the momentum in the absolute reference frame
##
## @var{gamma} @dots{} The three components of the momenta moment in the absolute reference frame,
## with respect to the coordinates of the node, thus to a moving frame
##
## @var{beta_dot} @dots{} The three components of the derivative of the momentum
##
## @var{gamma_dot} @dots{} The three components of the derivative of the momentum moment
##
## @end deftypefn

function [node_id, beta, gamma, beta_dot, gamma_dot] = mbdyn_post_load_output_ine(mbdyn_filename, filter_node_id, append_rows)
  if (nargin < 1 || nargin > 3 || nargout > 5)
    print_usage();
  endif

  if (nargin < 2)
    filter_node_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  ine_filename = mbdyn_post_output_filename(mbdyn_filename, ".ine");

  column_count = (nargout - 1) * 3;

  [node_id, data] = mbdyn_post_load_output(ine_filename, column_count, filter_node_id, append_rows);

  for i=1:length(data)
    if (nargout >= 2)
      beta{i} = data{i}(:, 1:3);
    endif

    if (nargout >= 3)
      gamma{i} = data{i}(:, 4:6);
    endif

    if (nargout >= 4)
      beta_dot{i} = data{i}(:, 7:9);
    endif

    if (nargout >= 5)
      gamma_dot{i} = data{i}(:, 10:12);
    endif
  endfor
endfunction
