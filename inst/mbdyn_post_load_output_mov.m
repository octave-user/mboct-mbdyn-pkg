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
## @deftypefn {Function File} [@var{node_id}, @var{trajectory}, @var{velocity}, @var{acceleration}, @var{orientation_description}, @var{deformation}]=mbdyn_post_load_output_mov (@var{mbdyn_filename}, @var{filter_node_id}, @var{append_rows})
##
## Loads the MBDyn .mov file named "<@var{mbdyn_filename}.mov>".
##
## @var{filter_node_id} @dots{} Load only the data of the corresponding nodes if present.
## If empty, the data of all nodes is loaded.
##
## @var{append_rows} @dots{} This parameter is used in order to optimize the memory allocation and should be set to the expected number of rows per node.
##
## @var{append_columns} @dots{} This parameter is used in order to optimize the memory allocation and should be set to the expected number of columns.
##
## @var{node_id} @dots{} Node number of structural nodes.
##
## @var{trajectory} @dots{} Matrix of position and orientation of structural nodes.
##
## @var{velocity} @dots{} Matrix of velocity and angular velocity of structural nodes.
##
## @var{acceleration} @dots{} Matrix of acceleration and angular acceleration.
##
## @var{orientation_description} @dots{} Cell array of strings. One of ("euler123", "euler321", "euler313", "phi").
##
## @var{deformation} @dots{} The difference between actual position/orientation and initial position/orientation.
##
## @end deftypefn

function [node_id, trajectory, velocity, acceleration, orientation_description, deformation] = mbdyn_post_load_output_mov(mbdyn_filename, filter_node_id, append_rows, append_columns, auto_resize_rows)
  if (nargin < 1 || nargout > 6)
    print_usage();
  endif

  if (nargin < 2)
    filter_node_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  if (nargin < 4)
    append_columns = 6;
  endif

  if (nargin < 5)
    auto_resize_rows = true;
  endif

  if (~ischar(mbdyn_filename))
    error("mbdyn_filename must be a string!");
  endif

  if (~isscalar(append_rows))
    error("append_rows must be a scalar!");
  endif

  mov_filename = mbdyn_post_output_filename(mbdyn_filename, ".mov");
  log_filename = mbdyn_post_output_filename(mbdyn_filename, ".log");

  column_count = 0;

  if (nargout >= 2)
    column_count = 6;
  endif

  if (nargout >= 3)
    column_count = 12;
  endif

  nodes = mbdyn_post_load_log_node(log_filename);

  node_labels = [nodes.label];

  [node_id, data] = mbdyn_post_load_output(mov_filename, column_count, filter_node_id, append_rows, append_columns, 1, auto_resize_rows);

  node_label_idx = zeros(1,length(node_id));

  for i=1:length(node_id)
    node_label_idx_i = find(node_id(i) == node_labels);

    if (length(node_label_idx_i) == 1)
      node_label_idx(i) = node_label_idx_i;
    else
      error("node_id %d not found in file \"%s\"!", node_id(i), log_filename);
    endif
  endfor

  orientation_description = { nodes(node_label_idx).orientation_description };

  if (length(data) == 0)
    if (nargout >= 2)
      trajectory = {};
    endif

    if (nargout >=  3)
      velocity = {};
    endif

    if (nargout >= 4)
      acceleration = {};
    endif

    if (nargout >= 6)
      deformation = {};
    endif
  else
    for i=1:length(data)
      if (nargout >= 2)
        trajectory{i} = data{i}(:, 1:6);

        if (nargout >=  3)
          velocity{i} = data{i}(:, 7:12);
        endif

        if (nargout >= 4)
          if (columns(data{i}) >= 18)
            acceleration{i} = data{i}(:, 13:18);
          else
            acceleration{i} = zeros(rows(data{i}), 0);
          endif
        endif

        switch (orientation_description{i})
          case {"euler123", "euler321", "euler313"}
            trajectory{i}(:,4:6) *= pi / 180; % euler angles are output in degrees
          case "phi"
            ## rotation vectors are already in radians
        endswitch
      endif

      if (nargout >= 6)
        node_idx_i = find(node_id(i) == [nodes.label]);

        if (length(node_idx_i) ~= 1)
          error("node_id(%d)=%d must appear exactly one times in the .log file!",i,node_id(i));
        endif

        deformation{i} = trajectory{i} - repmat([nodes(node_idx_i).X0.', ...
                                                 nodes(node_idx_i).Phi0.'], ...
                                                rows(trajectory{i}),1);

      endif
    endfor
  endif
endfunction
