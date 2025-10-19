## Copyright (C) 2023(-2025) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{options} = mbdyn_pre_solid_write_nodes(@var{mesh}, @var{nodes_file}, @var{options})
##
## Generate an MBDyn input file <@var{nodes_file}> containing all nodes from <@var{mesh}>.
##
## @var{mesh} @dots{} Finite Element mesh
##
## @var{nodes_file} @dots{} filename of output file
##
## @var{options}.struct_nodes.type @dots{} type of structural nodes to be created (e.g. MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP)
##
## @end deftypefn

function options = mbdyn_pre_solid_write_nodes(mesh, nodes_file, options)
  if (~(nargin == 3 && isstruct(mesh) && ischar(nodes_file) && isstruct(options)))
    print_usage();
  endif

  if (~isfield(options, "struct_nodes"))
    options.struct_nodes = struct();
  endif

  if (~isfield(options, "hydraulic_nodes"))
    options.hydraulic_nodes = struct();
  endif

  if (~isfield(options.struct_nodes, "reference_frame"))
    options.struct_nodes.reference_frame = "global";
  endif

  if (~isfield(options.struct_nodes, "number"))
    options.struct_nodes.number = 0;
  endif

  if (~isfield(options.hydraulic_nodes, "number"))
    options.hydraulic_nodes.number = 0;
  endif

  if (~isfield(options.hydraulic_nodes, "value"))
    options.hydraulic_nodes.value = 0;
  endif

  if (~isfield(options.struct_nodes, "type"))
    options.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
  endif

  if (~isfield(options, "orientation_description"))
    options.orientation_description = "euler123";
  endif

  fd = -1;

  unwind_protect
    [fd, msg] = fopen(nodes_file, "wt");

    if (fd == -1)
      error("failed to open file \"%s\": %s", nodes_file, msg);
    endif

    node_types = [MBDYN_NODE_TYPE_STATIC_STRUCT, ...
                  MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, ...
                  MBDYN_NODE_TYPE_DYNAMIC_STRUCT, ...
                  MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP];

    pressure_node_upc = false(rows(mesh.nodes));

    elem_types = fieldnames(mesh.elements);

    for j=1:numel(elem_types)
      switch (elem_types{j})
        case "iso8upc"
          elem_node_idx_upc = int32(1);
        case {"iso20upc", "iso20upcr"}
          elem_node_idx_upc = int32(1:8);
        case "iso27upc"
          elem_node_idx_upc = int32(1:27);
        case "penta15upc"
          elem_node_idx_upc = int32(1:6);
        case "tet10upc"
          elem_node_idx_upc = int32(1:4);
        otherwise
          elem_node_idx_upc = [];
      endswitch

      if (~isempty(elem_node_idx_upc))
        pressure_node_upc(getfield(mesh.elements, elem_types{j})(:, elem_node_idx_upc)) = true;
      endif
    endfor

    for j=node_types
      idx_node = find(options.struct_nodes.type == j);

      if (isempty(idx_node))
        continue;
      endif

      switch (j)
        case MBDYN_NODE_TYPE_STATIC_STRUCT
          struct_nodes_type = "static";
        case MBDYN_NODE_TYPE_STATIC_STRUCT_DISP
          struct_nodes_type = "static displacement";
        case MBDYN_NODE_TYPE_DYNAMIC_STRUCT
          struct_nodes_type = "dynamic";
        case MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP
          struct_nodes_type = "dynamic displacement";
        otherwise
          error("unexpected value for structural node type");
      endswitch

      switch (j)
        case {MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP}
          node_pos_idx = 1:3;
          nodes_format = sprintf("structural: %%d, %s, reference, %s, %%.16e,  %%.16e, %%.16e, reference, %s, null;\n", ...
                                 struct_nodes_type, ...
                                 options.struct_nodes.reference_frame,
                                 options.struct_nodes.reference_frame);
        case {MBDYN_NODE_TYPE_STATIC_STRUCT, MBDYN_NODE_TYPE_DYNAMIC_STRUCT}
          node_pos_idx = 1:6;
          nodes_format = sprintf("structural: %%d, %s, reference, %s, %%.16e,  %%.16e, %%.16e, reference, %s, %s, %%.16e, %%.16e, %%.16e, reference, %s, null, reference, %s, null;\n", ...
                                 struct_nodes_type, ...
                                 options.struct_nodes.reference_frame, ...
                                 options.struct_nodes.reference_frame, ...
                                 options.orientation_description, ...
                                 options.struct_nodes.reference_frame, ...
                                 options.struct_nodes.reference_frame);
        otherwise
          error("unexpected value for structural node type");
      endswitch

      node_id = options.struct_nodes.number + idx_node;

      node_data = [node_id.';
                   mesh.nodes(idx_node, node_pos_idx).'];

      fprintf(fd, nodes_format, node_data);

      idx_node_upc = find(pressure_node_upc(idx_node));

      if (~isempty(idx_node_upc))
        fprintf(fd, "hydraulic: %d, %g;\n", [node_id(idx_node_upc).';
                                             repmat(options.hydraulic_nodes.value, 1, numel(idx_node_upc))]);
      endif

      options.hydraulic_nodes.number += numel(idx_node_upc);
    endfor

    options.struct_nodes.number += rows(mesh.nodes);
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif

    fd = -1;
  end_unwind_protect
endfunction
