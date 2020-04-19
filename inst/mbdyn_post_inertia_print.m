## Copyright (C) 2014(-2016) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} mbdyn_post_inertia_print(@var{inertia}, @var{mass_log_file})
##
## Prints inertia properties to the file <@var{mass_log_file}>.
##
## @var{inertia} @dots{} Struct array with information about the bodies in the log file.
##
## @var{mass_log_file} @dots{} String name or file descriptor of the output file.
##
## @end deftypefn

function mbdyn_post_inertia_print(inertia, mass_log_file)
  if (nargin < 1 || nargin > 2)
    print_usage();
  endif

  if (nargin < 2)
    mass_log_file = stdout;
  endif
  
  fid = -1;

  unwind_protect
    if (ischar(mass_log_file))
      [fid, msg] = fopen(mass_log_file, "wt");

      if (fid == -1)
        error("could not open file \"%s\":%s", mass_log_file, msg);
      endif
    elseif (isscalar(mass_log_file))
      fid = mass_log_file;
    else
      error("invalid_argument!");
    endif
    
    fprintf(fid,"name\tdm\tXgc1\tXgc2\tXgc3\tJ11\tJ22\tJ33\tJ12\tJ13\tJ23\n");

    for i=1:length(inertia)
      fprintf(fid, "%s\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n",
              inertia(i).name, 
              inertia(i).dm,
              inertia(i).Xgc(1),
              inertia(i).Xgc(2),
              inertia(i).Xgc(3),
              inertia(i).J(1,1),
              inertia(i).J(2,2),
              inertia(i).J(3,3),
              inertia(i).J(1,2),
              inertia(i).J(1,3),
              inertia(i).J(2,3));
    endfor
  unwind_protect_cleanup                                                                                                         
    if ((fid ~= -1) && ischar(mass_log_file))
      fclose(fid);
    endif
  end_unwind_protect  
endfunction
