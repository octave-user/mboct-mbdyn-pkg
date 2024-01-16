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
## @deftypefn {Function File} @var{var} = mbdyn_post_load_log_var(@var{mbdyn_filename}, @var{variable}, @var{format}, @var{data_type})
## Read the value of variable name "<@var{variable}>" from MBDyn .log file "<@var{mbdyn_filename}.log>".
##
## @var{mbdyn_filename} @dots{} Name of MBDyn output file without extension .log
##
## @var{variable} @dots{} Variable name
##
## @var{format} @dots{} Format to be passed to fscanf
##
## @var{data_type} @dots{} MBDyn data type (e.g. "integer", "real", "string")
##
## @end deftypefn

function var = mbdyn_post_load_log_var(mbdyn_filename, variable_name, format, data_type)
  if (nargin < 2 || nargout > 1)
    print_usage();
  endif

  if (nargin < 3)
    format = "%d";
  endif

  if (nargin < 4)
    data_type = "";
  endif

  command = sprintf("exec awk -F ' ' '/%s %s =/{ print $4 }' '%s'", ...
                    data_type, ...
                    variable_name, ...
                    mbdyn_convert_path(mbdyn_post_output_filename(mbdyn_filename, ".log")));

  mbdyn_path_init();

  [status, output] = shell(command, true, "sync");

  if (status ~= 0)
    error("command \"%s\" failed with status %d", command, status);
  endif

  [var, count] = sscanf(output, format, "C");

  if (count == 0)
    error("could not find variable \"%s\" in file \"%s\"!", variable_name, mbdyn_filename);
  endif
endfunction
