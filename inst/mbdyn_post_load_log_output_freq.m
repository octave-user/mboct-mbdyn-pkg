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
## @deftypefn {Function File} @var{f} = mbdyn_post_load_log_output_freq(@var{mbdyn_filename})
## Read the output frequency from log file "<@var{mbdyn_filename}.log>".
##
## @var{mbdyn_filename} @dots{} Name of .log file without extension.
##
## @end deftypefn

function f = mbdyn_post_load_log_output_freq(mbdyn_filename)
  if (nargin ~= 1 || nargout > 1)
    print_usage();
  endif

  command = sprintf("exec awk -F ' ' '/output frequency:/ {print $3}' '%s.log'", mbdyn_convert_path(mbdyn_filename));

  mbdyn_path_init();

  [status, output] = shell(command, true, "sync");

  if (status ~= 0)
    error("command \"%s\" failed with status %d", command, status);
  endif

  [f, count] = sscanf(output, "%d", "C");

  if (count == 0)
    error("could not find variable \"%s\" in file \"%s\"!", variable_name, mbdyn_filename);
  endif
endfunction
