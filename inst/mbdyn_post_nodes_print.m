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
## @deftypefn {Function File} mbdyn_post_nodes_print(@var{nodes}, @var{log_file})
## @deftypefnx {} mbdyn_post_nodes_print(@var{nodes})
##
## Print a list of nodes from struct array <@var{nodes}> to file "<@var{log_file}>" or to stdout.
##
## @end deftypefn

function mbdyn_post_nodes_print(nodes, log_file)
  if (nargin < 1 || nargin > 2)
    print_usage();
  endif

  if (nargin < 2)
    log_file = stdout;
  endif

  fid = -1;

  unwind_protect
    if (ischar(log_file))
      [fid, msg] = fopen(log_file, "wt");

      if (fid == -1)
        error("could not open file \"%s\":%s", log_file, msg);
      endif
    elseif (isscalar(log_file))
      fid = log_file;
    else
      error("invalid_argument!");
    endif

    fprintf(fid, "label\t");

    for i=1:3
      fprintf(fid, "X%d\t", i);
    endfor

    for i=1:3
      for j=1:3
        fprintf(fid, "R%d%d\t", i, j);
      endfor
    endfor

    fprintf(fid, "\n");

    for i=1:length(nodes)
      fprintf(fid, "%d\t", nodes(i).label);
      for j=1:3
        fprintf(fid, '%e\t', nodes(i).X0(j));
      endfor
      for j=1:3
        for k=1:3
          fprintf(fid, '%e\t', nodes(i).R0(j, k));
        endfor
      endfor
      fprintf(fid, "\n");
    endfor
  unwind_protect_cleanup
    if ((fid ~= -1) && ischar(log_file))
      fclose(fid);
    endif
  end_unwind_protect
endfunction
