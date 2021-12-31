## Copyright (C) 2014(-2021) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{M}, @var{D}, @var{K}] = mbdyn_post_load_fem_data(@var{filename})
## @var{M} @dots{} Reduced order mass matrix
##
## @var{D} @dots{} Reduced order damping matrix
##
## @var{K} @dots{} Reduced order stiffness matrix
##
## Load data from a FEM input file for MBDyn's model element.
## @end deftypefn

function [M, D, K] = mbdyn_post_load_fem_data(filename)
  if (nargin ~= 1 || ~ischar(filename))
    print_usage();
  endif
  
  fd = -1;
  NODES = int32(0);
  MODES = int32(0);
  M = [];
  D = [];
  K = [];
  
  unwind_protect
    [fd, msg] = fopen(filename, "rt");

    if (fd == -1)
      error("failed to open file \"%s\": %s", filename, msg);
    endif
    
    while (true)
      line = fgetl(fd);

      if (~ischar(line))
        break;
      endif

      switch (line)
        case "** RECORD GROUP 1, HEADER"
          line = fgetl(fd);
          
          if (~ischar(line))
            break;
          endif
          
          [REV, NODES, MODES, COUNT] = fscanf(fd, "%s %d %d %d %d %d", "C");        
        case "** RECORD GROUP 9, MODAL MASS MATRIX"
          M = fem_read_matrix(fd, MODES, MODES);
        case "** RECORD GROUP 10, MODAL STIFFNESS MATRIX"
          K = fem_read_matrix(fd, MODES, MODES);
        case "** RECORD GROUP 13, MODAL DAMPING MATRIX"
          D = fem_read_matrix(fd, MODES, MODES);
      endswitch
    endwhile
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction

function A = fem_read_matrix(fd, NR, NC)
  A = fscanf(fd, "%g", [NR, NC]);
endfunction
