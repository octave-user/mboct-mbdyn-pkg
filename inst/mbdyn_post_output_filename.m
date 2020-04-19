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
## @deftypefn {Function File} @var{output_filename} = mbdyn_post_output_filename(@var{input_filename}, @var{ext})
## Append the extension <@var{ext}> to filename <@var{input_filename}>.
## If "<@var{input_filename}>" has an extension, it will be treated in the same way MBDyn does.
## @end deftypefn

function output_filename = mbdyn_post_output_filename(input_filename, ext)
  if (nargin < 1 || nargin > 2 || nargout > 1)
    print_usage();
  endif
  
  if (nargin < 2)
    ext = "";
  endif
  
  [dir, name, ext_out] = fileparts(input_filename);

  switch (ext_out)
    case {".mbd", ".mbdyn", ".mov", ext}
      ext_out = ext;
    otherwise
      ext_out = strcat(ext_out, ext);
  endswitch

  output_filename = fullfile(dir, cstrcat(name, ext_out));
endfunction
