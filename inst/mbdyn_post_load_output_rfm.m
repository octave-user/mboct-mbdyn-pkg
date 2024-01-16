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
## @deftypefn {Function File} [ @var{ref_id}, @var{position}, @var{orientation}, @var{velocity}, @var{angular_velocity} ] = mbdyn_post_load_output_rfm(@var{mbdyn_filename},@var{filter_ref_id})
##
## Loads the MBDyn .rfm file named "<@var{mbdyn_filename}.rfm>".
##
## @var{ref_id} @dots{} Vector of reference frame numbers used in the MBDyn input file.
##
## @var{position}(@var{i}) @dots{} Position of the reference frame <@var{ref_id}(@var{i})>
##
## @var{orientation}(@var{i}) @dots{} Orientation of the reference frame <@var{ref_id}(@var{i})>
##
## @var{velocity}(@var{i}) @dots{} Velocity of the reference frame <@var{ref_id}(@var{i})>
##
## @var{angular_velocity}(@var{i}) @dots{} Angular velocity of the reference frame <@var{ref_id}(@var{i})>
##
## @var{filter_ref_id} @dots{} Load only the data of the corresponding reference frames if present.
## If empty, the data of all reference frames is loaded.
##
## @end deftypefn

function [ref_id, position, orientation, velocity, angular_velocity] = mbdyn_post_load_output_rfm(mbdyn_filename, filter_ref_id)
  if (nargin < 1 || nargin > 2 || nargout > 5)
    print_usage();
  endif

  if (nargin < 2)
    filter_ref_id = [];
  endif

  if (~ischar(mbdyn_filename))
    error("mbdyn_filename must be a string!");
  endif

  rfm_filename = mbdyn_post_output_filename(mbdyn_filename, ".rfm");

  def_col_size = 0;

  if (nargout >= 2)
    def_col_size = 3;
  endif

  if (nargout >= 3)
    def_col_size = 6;
  endif

  if (nargout >= 4)
    def_col_size = 9;
  endif

  if (nargout >= 5)
    def_col_size = 12;
  endif

  [ref_id, data] = mbdyn_post_load_output(rfm_filename, def_col_size, filter_ref_id);

  for i=1:length(data)
    if (nargout >= 2)
      position{i} = data{i}(:,1:3);
    endif

    if (nargout >= 3)
      ## FIXME: we have to consider the orientation description of each reference frame!
      orientation{i} = data{i}(:, 4:6) * pi / 180;
    endif

    if (nargout >= 4)
      velocity{i} = data{i}(:,7:9);
    endif

    if (nargout >= 5)
      angular_velocity{i} = data{i}(:,10:12);
    endif
  endfor
endfunction
