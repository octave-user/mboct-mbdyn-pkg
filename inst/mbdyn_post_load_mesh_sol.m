## Copyright (C) 2025(-2025) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{mesh}] = mbdyn_post_load_mesh_sol(@var{output_file})
##
## Loads mesh data for solid and beam elements from an MBDyn output file "<@var{output_file}>".
##
## @var{output_file} @dots{} Name of the MBDyn output file without extension.
##
## @var{mesh} @dots{} Finite element mesh data structure.
##
## @seealso{fem_post_mesh_export}
## @end deftypefn

function [mesh, elem_types] = mbdyn_post_load_mesh_sol(output_file)
  if (nargin ~= 1 || ~ischar(output_file))
    print_usage();
  endif

  mesh = struct();

  log_file = [output_file, ".log"];

  [info, err, msg] = stat(log_file);

  if (err ~= 0)
    error("file \"%s\" not found: %s", log_file, msg);
  endif

  awk_cmd = sprintf("awk -F ' ' '/^structural node:/{printf(\"%%d %%g %%g %%g %%g %%g %%g %%d\\n\", $3, $4, $5, $6, $8, $9, $10, ($7 == \"euler123\") ? 123 : ($7 == \"phi\") ? 100 : 999);}' \"%s\"", log_file);

  fd = -1;

  unwind_protect
    fd = popen(awk_cmd, "r");

    if (fd == -1)
      error("failed to open file \"%s\"", log_file);
    endif

    [val, count, msg] = fscanf(fd, "%d %g %g %g %g %g %g %d\n", [8, inf]);
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect

  if (count == 0)
    error("no nodes found in file \"%s\"", log_file);
  endif

  mesh.node_id = val(1, :).';
  mesh.node_id_transl(mesh.node_id) = 1:rows(mesh.node_id);
  mesh.nodes = val(2:7, :).';
  mesh.node_type = val(8, :).';

  elem_types = mbdyn_post_elem_types();

  mesh.elements = struct();
  mesh.elem_id = struct();
  mesh.elem_node_offset = struct();

  for i=1:numel(elem_types)
    offset_col = sum(isspace(elem_types(i).elem_tag));
    arg_fmt = sprintf(", $%d", 2 + offset_col);

    for j=1:elem_types(i).num_cols
      arg_fmt = [arg_fmt, sprintf(", $%d", j + 2 + offset_col)];
    endfor

    awk_cmd = sprintf("awk -F ' ' '/^%s:/{printf(\"%%d%s\\n\"%s);}' \"%s\"", elem_types(i).elem_tag, repmat(" %g", 1, elem_types(i).num_cols), arg_fmt, log_file);

    fd = -1;

    unwind_protect
      fd = popen(awk_cmd, "r");

      if (fd == -1)
        error("failed to open file \"%s\"", log_file);
      endif

      [val, count, msg] = fscanf(fd, ["%d", repmat(" %g", 1, elem_types(i).num_cols), "\n"], [elem_types(i).num_cols + 1, inf]);

      if (count)
        elem_id = val(1, :).';
        elem_nodes = mesh.node_id_transl(val(1+elem_types(i).node_cols, :).');

        if (~isempty(elem_types(i).node_offset))
          elem_node_offset = reshape(val(1+elem_types(i).node_offset, :).', rows(elem_nodes), 3, columns(elem_nodes));
          mesh.elem_node_offset = setfield(mesh.elem_node_offset, elem_types(i).elem_type, elem_node_offset);
        endif

        if (~isfield(mesh.elem_id, elem_types(i).elem_type))
          mesh.elem_id = setfield(mesh.elem_id, elem_types(i).elem_type, []);
          mesh.elements = setfield(mesh.elements, elem_types(i).elem_type, []);
        endif

        mesh.elem_id = setfield(mesh.elem_id, elem_types(i).elem_type, [getfield(mesh.elem_id, elem_types(i).elem_type); elem_id]);
        mesh.elements = setfield(mesh.elements, elem_types(i).elem_type, [getfield(mesh.elements, elem_types(i).elem_type); elem_nodes]);
      endif
    unwind_protect_cleanup
      if (fd ~= -1)
        fclose(fd);
      endif
    end_unwind_protect
  endfor

  ## FIXME: beam3 is not supported by Gmsh
  if (isfield(mesh.elements, "beam3"))
    if (isfield(mesh.elements, "beam2"))
      elem_nodes = mesh.elements.beam2;
      elem_node_offset = mesh.elem_node_offset.beam2;
    else
      elem_nodes = zeros(0, 2);
    endif

    mesh.elements.beam2 = [elem_nodes;
                           [mesh.elements.beam3(:, 1:2);
                            mesh.elements.beam3(:, 2:3)]];

    mesh.elements = rmfield(mesh.elements, "beam3");
  endif

  mesh.elements.point1 = int32(1:rows(mesh.nodes))(:);
endfunction

function elem_types_out = mbdyn_post_elem_types()
  persistent elem_types = [];

  if (isempty(elem_types))
    empty_cell = cell(1, 28);

    elem_types = struct("elem_type", empty_cell, "elem_tag", empty_cell, "num_cols", empty_cell, "node_cols", empty_cell, "node_offset", empty_cell);

    elem_types(1).elem_type = "iso8";
    elem_types(1).node_cols = [1,2,3,4,5,6,7,8];
    elem_types(1).node_offset = [];
    elem_types(1).num_cols = 8;
    elem_types(1).elem_tag = "hexahedron8";

    elem_types(2).elem_type = "iso8upc";
    elem_types(2).node_cols = [1,2,3,4,5,6,7,8];
    elem_types(2).node_offset = [];
    elem_types(2).num_cols = 8;
    elem_types(2).elem_tag = "hexahedron8upc";

    elem_types(3).elem_type = "iso20";
    elem_types(3).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
    elem_types(3).node_offset = [];
    elem_types(3).num_cols = 20;
    elem_types(3).elem_tag = "hexahedron20";

    elem_types(4).elem_type = "iso20upc";
    elem_types(4).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
    elem_types(4).node_offset = [];
    elem_types(4).num_cols = 20;
    elem_types(4).elem_tag = "hexahedron20upc";

    elem_types(5).elem_type = "iso20upcr";
    elem_types(5).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
    elem_types(5).node_offset = [];
    elem_types(5).num_cols = 20;
    elem_types(5).elem_tag = "hexahedron20upcr";

    elem_types(6).elem_type = "iso27";
    elem_types(6).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27];
    elem_types(6).node_offset = [];
    elem_types(6).num_cols = 27;
    elem_types(6).elem_tag = "hexahedron27";

    elem_types(7).elem_type = "iso27upc";
    elem_types(7).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27];
    elem_types(7).node_offset = [];
    elem_types(7).num_cols = 27;
    elem_types(7).elem_tag = "hexahedron27upc";

    elem_types(8).elem_type = "iso20r";
    elem_types(8).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
    elem_types(8).node_offset = [];
    elem_types(8).num_cols = 20;
    elem_types(8).elem_tag = "hexahedron20r";

    elem_types(9).elem_type = "penta15";
    elem_types(9).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
    elem_types(9).node_offset = [];
    elem_types(9).num_cols = 15;
    elem_types(9).elem_tag = "pentahedron15";

    elem_types(10).elem_type = "penta15upc";
    elem_types(10).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
    elem_types(10).node_offset = [];
    elem_types(10).num_cols = 15;
    elem_types(10).elem_tag = "pentahedron15upc";

    elem_types(11).elem_type = "tet10h";
    elem_types(11).node_cols = [1,2,3,4,5,6,7,8,9,10];
    elem_types(11).node_offset = [];
    elem_types(11).num_cols = 10;
    elem_types(11).elem_tag = "tetrahedron10";

    elem_types(12).elem_type = "tet10upc";
    elem_types(12).node_cols = [1,2,3,4,5,6,7,8,9,10];
    elem_types(12).node_offset = [];
    elem_types(12).num_cols = 10;
    elem_types(12).elem_tag = "tetrahedron10upc";

    elem_types(13).elem_type = "beam2";
    elem_types(13).node_cols = [1,5];
    elem_types(13).node_offset = [2,3,4;
                                  6,7,8];
    elem_types(13).num_cols = 8;
    elem_types(13).elem_tag = "beam2";

    elem_types(14).elem_type = "beam3";
    elem_types(14).node_cols = [1,5,9];
    elem_types(14).node_offset = [2,3,4;
                                  6,7,8;
                                  10,11,12];
    elem_types(14).num_cols = 12;
    elem_types(14).elem_tag = "beam3";

    elem_types(15).elem_type = "iso8f";
    elem_types(15).node_cols = [1,2,3,4,5,6,7,8];
    elem_types(15).node_offset = [];
    elem_types(15).num_cols = 8;
    elem_types(15).elem_tag = "hexahedron8f";

    elem_types(16).elem_type = "iso20f";
    elem_types(16).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
    elem_types(16).node_offset = [];
    elem_types(16).num_cols = 20;
    elem_types(16).elem_tag = "hexahedron20f";

    elem_types(17).elem_type = "iso27f";
    elem_types(17).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27];
    elem_types(17).node_offset = [];
    elem_types(17).num_cols = 27;
    elem_types(17).elem_tag = "hexahedron27f";

    elem_types(18).elem_type = "iso20fr";
    elem_types(18).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
    elem_types(18).node_offset = [];
    elem_types(18).num_cols = 20;
    elem_types(18).elem_tag = "hexahedron20fr";

    elem_types(19).elem_type = "penta15f";
    elem_types(19).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
    elem_types(19).node_offset = [];
    elem_types(19).num_cols = 15;
    elem_types(19).elem_tag = "pentahedron15f";

    elem_types(20).elem_type = "tet10hf";
    elem_types(20).node_cols = [1,2,3,4,5,6,7,8,9,10];
    elem_types(20).node_offset = [];
    elem_types(20).num_cols = 10;
    elem_types(20).elem_tag = "tetrahedron10f";

    elem_types(21).elem_type = "tet20";
    elem_types(21).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
    elem_types(21).node_offset = [];
    elem_types(21).num_cols = 20;
    elem_types(21).elem_tag = "tetrahedron20";

    elem_types(22).elem_type = "tet20f";
    elem_types(22).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20];
    elem_types(22).node_offset = [];
    elem_types(22).num_cols = 20;
    elem_types(22).elem_tag = "tetrahedron20f";

    elem_types(23).elem_type = "line2";
    elem_types(23).node_cols = [1,23];
    elem_types(23).node_offset = [];
    elem_types(23).num_cols = 56;
    elem_types(23).elem_tag = "totaljoint";

    elem_types(24).elem_type = "line2";
    elem_types(24).node_cols = [1,14];
    elem_types(24).node_offset = [];
    elem_types(24).num_cols = 20;
    elem_types(24).elem_tag = "journal bearing";

    elem_types(25).elem_type = "line2";
    elem_types(25).node_cols = [1,14];
    elem_types(25).node_offset = [];
    elem_types(25).num_cols = 25;
    elem_types(25).elem_tag = "inline friction";

    elem_types(26).elem_type = "line2";
    elem_types(26).node_cols = [1,14];
    elem_types(26).node_offset = [];
    elem_types(26).num_cols = 25;
    elem_types(26).elem_tag = "deformablejoint";

    elem_types(27).elem_type = "penta18";
    elem_types(27).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
    elem_types(27).node_offset = [];
    elem_types(27).num_cols = 18;
    elem_types(27).elem_tag = "pentahedron18";

    elem_types(28).elem_type = "penta18f";
    elem_types(28).node_cols = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18];
    elem_types(28).node_offset = [];
    elem_types(28).num_cols = 18;
    elem_types(28).elem_tag = "pentahedron18f";
  endif

  elem_types_out = elem_types;
endfunction
