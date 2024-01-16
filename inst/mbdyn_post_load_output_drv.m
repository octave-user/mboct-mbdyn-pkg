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
## @deftypefn {Function File} [@var{drive_id}, @var{drive_value}] = mbdyn_post_load_output_drv(@var{mbdyn_filename}, @var{filter_drive_id}, @var{append_rows})
##
## Loads output from drive callers from MBDyn output file "<@var{mbdyn_filename}.drv>".
##
## @var{filter_drive_id} @dots{} Load only drives with drive numbers in <@var{filter_drive_id}>. If empty, everything is loaded.
##
## @var{drive_id} @dots{} Drive identifier of the loaded drives.
##
## @var{append_rows} @dots{} The expected number of rows for which memory should be allocated.
##
## @end deftypefn

function [drive_id, drive_value] = mbdyn_post_load_output_drv(mbdyn_filename, filter_drive_id, append_rows, auto_resize_rows)
  if (nargin < 1 || nargout > 2)
    print_usage();
  endif

  if (nargin < 2)
    filter_drive_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  if (nargin < 4)
    auto_resize_rows = true;
  endif

  drv_filename = mbdyn_post_output_filename(mbdyn_filename, ".drv");

  column_count = 1;

  [drive_id, drive_value] = mbdyn_post_load_output(drv_filename, column_count, filter_drive_id, append_rows, 2, 1, auto_resize_rows);
endfunction
