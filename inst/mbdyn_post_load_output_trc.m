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
## @deftypefn {Function File} [@var{drive_id}, @var{drive_value}] = mbdyn_post_load_output_trc(@var{mbdyn_filename}, @var{filter_drive_id}, @var{append_rows})
##
## Loads the MBDyn .trc file "<@var{mbdyn_filename}.trc>".
##
## @var{filter_drive_id} @dots{} Load only drives with drive numbers in <@var{filter_drive_id}>.
##
## @var{append_rows} @dots{} Hint for memory reallocation.
##
## @var{drive_id} @dots{} Numbers of the loaded drives.
##
## @var{drive_value} @dots{} Value and derivative of the drive.
##
## @end deftypefn

function [drive_id, drive_value] = mbdyn_post_load_output_trc(mbdyn_filename, filter_drive_id, append_rows)
  if (nargin < 1 || nargout > 2)
    print_usage();
  endif

  if (nargin < 2)
    filter_drive_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  trc_filename = mbdyn_post_output_filename(mbdyn_filename, ".trc");

  column_count = 1;

  [drive_id, drive_value] = mbdyn_post_load_output(trc_filename, column_count, filter_drive_id, append_rows);
endfunction
