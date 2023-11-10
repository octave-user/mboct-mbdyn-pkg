## Copyright (C) 2015(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} mbdyn_pre_input_file_copy(@var{input_files}, @var{working_directory})
## @deftypefnx {} mbdyn_pre_input_file_copy(@var{input_files}, @var{working_directory}, @var{input_directory})
## @deftypefnx {} mbdyn_pre_input_file_copy(@var{input_files}, @var{working_directory}, @var{input_directory}, @var{fcopy})
##
## Creates symbolic links or copies of files <@var{input_files}> inside <@var{working_directory}>.
## This is especially useful, if only a subset of input files are generated by a user program
## and additional hand written input files are needed. Those files must be located inside the
## same directory like the main input file, so MBDyn can find it.
##
## @var{input_files} @dots{} Cell array of character strings of filenames.
##
## @var{working_directory} @dots{} Path where the symbolic links or copies are created.
##
## @var{input_directory} @dots{} Path name where <@var{input_files}> are located.
##
## @var{fcopy} @dots{} If true, copy the files. Otherwise create symbolic links.
##
## @end deftypefn

function mbdyn_pre_input_file_copy(input_files, working_directory, input_directory, fcopy)
  if (nargin < 2)
    print_usage();
  endif

  if (nargin < 3)
    input_directory = pwd();
  endif

  if (nargin < 4)
    fcopy = false;
  endif
  
  for i=1:numel(input_files)
    old_file = make_absolute_filename(fullfile(input_directory, input_files{i}));
    [DIR, NAME, EXT] = fileparts(input_files{i});
    input_file_i = strcat(NAME, EXT);
    new_file = make_absolute_filename(fullfile(working_directory, input_file_i));

    if (~strcmp(old_file, new_file))
      [~] = unlink(new_file);

      if (fcopy || ~isunix())
        [err, msg] = copyfile(old_file, new_file);

        if (err ~= 1)
          error("copyfile(\"%s\",\"%s\") failed with code %d: %s", old_file, new_file, err, msg);
        endif
      else
        [err, msg] = symlink(old_file, new_file);

        if (err ~= 0)
          error("symlink(\"%s\",\"%s\") failed with code %d: %s", old_file, new_file, err, msg);
        endif
      endif
    endif
  endfor
endfunction
