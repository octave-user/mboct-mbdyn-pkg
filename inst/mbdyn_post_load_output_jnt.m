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
## @deftypefn {Function File} [@var{joint_id}, @var{local_reaction}, @var{global_reaction}] = mbdyn_post_load_output_jnt(@var{mbdyn_filename}, @var{filter_joint_id}, @var{append_rows})
##
## Loads the MBDyn .jnt file named "<@var{mbdyn_filename}.jnt>".
##
## @var{joint_id} @dots{} The label of joint element.
##
## @var{local_reaction} @dots{} The three components of the reaction force and reaction couple in a local reference frame.
##
## @var{global_reaction} @dots{} The three components of the reaction force and reaction couple in the global frame.
##
## @end deftypefn

function [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(mbdyn_filename, filter_joint_id, append_rows, auto_resize_rows)
  if (nargin < 1 || nargin > 4 || nargout > 3)
    print_usage();
  endif

  if (nargin < 2)
    filter_joint_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  if (nargin < 4)
    auto_resize_rows = true;
  endif

  jnt_filename = mbdyn_post_output_filename(mbdyn_filename, ".jnt");

  column_count = 0;

  if (nargout >= 2)
    column_count = 6;
  endif

  if (nargout >= 3)
    column_count = 12;
  endif

  [joint_id, data] = mbdyn_post_load_output(jnt_filename, column_count, filter_joint_id, append_rows, 12, 1, auto_resize_rows);

  if (length(joint_id) == 0)
    local_reaction = {};
    global_reaction = {};
  else
    for i=1:length(data)
      if (nargout >= 2)
        local_reaction{i} = data{i}(:, 1:6);
      endif
      if (nargout >= 3)
        global_reaction{i} = data{i}(:, 7:12);
      endif
    endfor
  endif
endfunction
