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
## @deftypefn {Function File} mbdyn_pre_write_param_file(@var{output_file}, @var{param},  @var{options})
##
## Write all scalar variables in struct <@var{param}> to file "<@var{output_file}>".
##
## @var{output_file} @dots{} Output file which can be included in an MBDyn input file.
##
## @var{param} @dots{} Data structure with scalar variables.
##
## @var{options}.open_mode @dots{} Argument to be passed to fopen.
##
## @end deftypefn

function mbdyn_pre_write_param_file(output_file, param, options)
  if (nargin < 2 || nargin > 3)
    print_usage();
  endif

  if (~ischar(output_file) && ~isscalar(output_file))
    error("output_file must be a filename or a file descriptor!");
  endif

  if (~isstruct(param))
    error("param must be a struct!");
  endif

  if (nargin < 3)
    options = struct();
  endif

  if (~isstruct(options))
    error("options must be a struct!");
  endif

  if (~isfield(options, "open_mode"))
    options.open_mode = "wt";
  endif

  fout = -1;
  owns_fd = false;

  unwind_protect
    if (ischar(output_file))
      owns_fd = true;

      [fout, msg] = fopen(output_file, options.open_mode);

      if (fout == -1)
        error("could not open file \"%s\": %s",  output_file, msg);
      endif
    else
      fout = output_file;
    endif

    param_names = fieldnames(param);

    for i=1:length(param_names)
      param_val = getfield(param, param_names{i});
      if (isscalar(param_val))
        if (isinteger(param_val))
          fprintf(fout, "set: integer %s = %d;\n", param_names{i}, param_val);
        elseif (isbool(param_val))
          fprintf(fout, "set: bool %s = %s;\n",  param_names{i},  {"FALSE", "TRUE"}{param_val + 1});
        elseif (isreal(param_val))
          if (param_val > realmax)
            param_val = realmax;
          elseif (param_val < -realmax)
            param_val = -realmax;
          endif
          fprintf(fout, "set: real %s = %.16e;\n", param_names{i}, param_val);
        endif
      elseif (ischar(param_val))
        fprintf(fout, "set: string %s = \"%s\";\n", param_names{i}, param_val);
      endif
    endfor
  unwind_protect_cleanup
    if (owns_fd && fout ~= -1)
      fclose(fout);
    endif
  end_unwind_protect
endfunction
