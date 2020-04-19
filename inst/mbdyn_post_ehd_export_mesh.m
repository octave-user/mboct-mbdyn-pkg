## Copyright (C) 2020(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} mbdyn_post_ehd_export_mesh(@var{mesh}, @var{filename})
## Export 3D mesh from hydrodynamic plain bearings Gmsh format
##
## @var{mesh} @dots{} Return value from mbdyn_post_ehd_create_mesh.
##
## @var{filename} @dots{} Character string name of output file in Gmsh format.
##
## @seealso{mbdyn_post_ehd_export_data, mbdyn_post_ehd_create_mesh}
## @end deftypefn

function mbdyn_post_ehd_export_mesh(mesh, filename)
  if (nargin ~= 2 || nargout > 0)
    print_usage();
  endif

  fd = -1;
  
  unwind_protect
    [fd] = fopen(filename, "wt");

    if (fd == -1)
      error("failed to open file \"%s\"", filename);
    endif

    fputs(fd, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

    inumelem = int32(0);
    elem_types = fieldnames(mesh.elements);

    for i=1:numel(elem_types)
      switch elem_types{i}
	case {"quad4"}
	  inumelem += rows(getfield(mesh.elements, elem_types{i}));
      endswitch
    endfor

    fputs(fd, "$PhysicalNames\n");
    fprintf(fd, "%d\n", numel(mesh.bearings));

    for i=1:numel(mesh.bearings)
      group_dim = 2;
      fprintf(fd, ...
	      "%d %d \"bearing[%d]\"\n", ...
	      group_dim, ...
	      mesh.bearings(i).label, ...
	      mesh.bearings(i).label);
    endfor

    fputs(fd, "$EndPhysicalNames\n");

    fprintf(fd, "$Nodes\n%d\n", rows(mesh.nodes));
    fprintf(fd, "%d %.16e %.16e %.16e\n", [1:rows(mesh.nodes); mesh.nodes.']);
    fputs(fd, "$EndNodes\n");

    fprintf(fd, "$Elements\n%d\n", inumelem);

    inumelem = int32(0);

    for i=1:numel(elem_types)
      switch elem_types{i}
	case "quad4"
	  elem_type_id = 3;
	  elem_node_order = 1:4;
	otherwise
	  continue
      endswitch

      elem_nodes = getfield(mesh.elements, elem_types{i});
      elem_groups = zeros(rows(elem_nodes), 1, "int32");
      elem_tags = [repmat([elem_type_id; 2], 1, rows(elem_nodes)); zeros(2, rows(elem_nodes))];

      for j=1:numel(mesh.bearings)
	elem_tags(3:4, mesh.bearings(j).elements) = mesh.bearings(j).label;
      endfor

      numcols = rows(elem_tags) + columns(elem_nodes) + 1;
      format = "%d";

      for j=1:numcols - 1
	format = [format, " %d"];
      endfor

      format = [format, "\n"];
      fprintf(fd, format, [inumelem + (1:rows(elem_nodes));
			   elem_tags;
			   elem_nodes(:, elem_node_order).']);
      inumelem += rows(elem_nodes);
    endfor

    fputs(fd, "$EndElements\n");
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif
  end_unwind_protect
endfunction
