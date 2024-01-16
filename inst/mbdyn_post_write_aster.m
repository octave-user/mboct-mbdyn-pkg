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
## @deftypefn {Function File} mbdyn_post_write_aster(@var{data}, @var{mesh_file})
##
## Converts a MBDyn model to Code_Aster .mail format.
##
## @var{data} @dots{} The result obtained from mbdyn_post_load_log.
##
## @var{mesh_file} @dots{} The name of the Code_Aster .mail file.
##
## @end deftypefn

function mbdyn_post_write_aster(data, mesh_file)
  if (nargin < 1 || nargin > 2)
    print_usage();
  endif

  if (nargin < 2)
    mesh_file = stdout;
  endif

  owns_fd = false;
  fd = -1;

  unwind_protect
    if (ischar(mesh_file))
      owns_fd = true;

      [fd, msg] = fopen(mesh_file, "wt");

      if (fd == -1)
        error("failed to open file \"%s\": %s", mesh_file, msg);
      endif
    else
      fd = mesh_file;

      if (~isscalar(fd))
        error("mesh_file must be a open file descriptor or a file name");
      endif
    endif


    fprintf(fd, "%% --------------------------------------------------------------------------------\n");
    fprintf(fd, " TITRE\n");
    fprintf(fd, "%s\n", ctime(time()));
    fprintf(fd, " FINSF\n");

    if (length(data.nodes) > 0)
      fprintf(fd, " %%\n");
      fprintf(fd, " COOR_3D\n");

      for i=1:length(data.nodes)
        fprintf(fd, " N%-8d ", data.nodes(i).label);

        for j=1:3
          fprintf(fd, "%.14e ", data.nodes(i).X0(j));
        endfor

        fprintf(fd, "\n");
      endfor

      fprintf(fd, " FINSF\n");
    endif

    if (length(data.beams2) > 0)
      fprintf(fd, " %%\n");

      fprintf(fd, " SEG2\n");

      for i=1:length(data.beams2)
        fprintf(fd, " M%-8d ", data.beams2(i).label);

        for j=1:2
          fprintf(fd, "N%-8d ", data.beams2(i).nodes(j).label);
        endfor

        fprintf(fd, "\n");
      endfor

      fprintf(fd, " FINSF\n");
    endif

    if (length(data.beams3) > 0)
      fprintf(fd, " %%\n");

      fprintf(fd, " SEG3\n");

      for i=1:length(data.beams3)
        fprintf(fd, " M%-8d ", data.beams3(i).label);

        for j=1:3
          fprintf(fd, "N%-8d ", data.beams3(i).nodes([1,3,2])(j).label);
        endfor

        fprintf(fd, "\n");
      endfor

      fprintf(fd, " FINSF\n");
    endif

    fprintf(fd, " FIN\n");
  unwind_protect_cleanup
    if (fd ~= -1 && owns_fd)
      fclose(fd);
    endif
  end_unwind_protect
endfunction
