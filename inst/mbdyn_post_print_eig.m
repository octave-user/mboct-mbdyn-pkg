## Copyright (C) 2016 Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} mbdyn_post_print_eig(@var{modal}, @var{log_dat}, @var{output_file})
##
## Print all eigenvalues from an eigenanalysis to a file.
##
## @var{modal} @dots{} Return value of mbdyn_post_load_output_eig.
##
## @var{log_dat} @dots{} Return value of mbdyn_post_load_log.
##
## @var{output_file} @dots{} File name for results.
##
## @end deftypefn

function mbdyn_post_print_eig(modal, log_dat, varargin)
  fd = -1;
  owns_fd = false;

  if (length(varargin) >= 1)
    if (ischar(varargin{1}))
      owns_fd = true;
      [fd, msg] = fopen(varargin{1}, "wt");
      if (fd == -1)
        error("failed to open file \"%s\": %s", varargin{1}, msg);
      endif
    else
      fd = varargin{1};
    endif
  else
    fd = stdout;
  endif

  unwind_protect
    if (isfield(log_dat.dof_info, "struct_node_labels"))
      for j=1:columns(modal.VR)
        fprintf(fd, "f=%.3fHz\t", modal.f(j));
        for i=1:length(log_dat.dof_info.struct_node_labels)
          fprintf(fd, "node %d\t", log_dat.dof_info.struct_node_labels(i));
        endfor
        fprintf(fd, "\n");

        norm_VR = max(abs(modal.VR(:, j)));

        for k=1:6
          switch (k)
            case {1, 2, 3}
              fprintf(fd, "X%d\t", k);
            case {4, 5, 6}
              fprintf(fd, "g%d\t", k - 3);
          endswitch
          for i=1:length(log_dat.dof_info.struct_node_dofs)
            fprintf(fd, "%.6f\t", abs(modal.VR(log_dat.dof_info.struct_node_dofs(i) + k, j)) / norm_VR);
          endfor
          fprintf(fd, "\n");
        endfor
      endfor
    endif
  unwind_protect_cleanup
    if (fd ~= -1 && owns_fd)
      fclose(fd);
    endif
  end_unwind_protect
endfunction
