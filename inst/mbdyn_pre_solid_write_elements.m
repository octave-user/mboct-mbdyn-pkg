## Copyright (C) 2023(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} @var{options} = mbdyn_pre_solid_write_elements(@var{mesh}, @var{load_case_dof}, @var{load_case}, @var{elem_file}, @var{options})
##
## Generate an MBDyn input file <@var{elem_file}> containing all solid elements, constraints and loads from <@var{mesh}>.
##
## @var{mesh} @dots{} Finite Element mesh
##
## @var{elem_file} @dots{} filename of output file
##
## @var{options}.forces.time_function @dots{} expression for the time variation of all applied loads
##
## @end deftypefn

function options = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, options)
  if (~(nargin == 5 && isstruct(mesh) && isscalar(mesh) && isstruct(load_case_dof) && isscalar(load_case_dof) && isstruct(load_case) && ischar(elem_file) && isstruct(options)))
    print_usage();
  endif

  if (~isfield(options, "solids"))
    options.solids = struct();
  endif

  if (~isfield(options, "genels"))
    options.genels = struct();
  endif

  if (~isfield(options.solids, "number"))
    options.solids.number = 0;
  endif

  if (~isfield(options.genels, "number"))
    options.genels.number = 0;
  endif

  if (~isfield(options, "forces"))
    options.forces = struct();
  endif

  if (~isfield(options.forces, "number"))
    options.forces.number = 0;
  endif

  if (~isfield(options.forces, "time_function"))
    options.forces.time_function = "time";
  endif

  if (~iscell(options.forces.time_function))
    if (~ischar(options.forces.time_function))
      error("class of options.forces.time_function must be char or cell");
    endif

    cell_time_func = cell(1, numel(load_case));

    for i=1:numel(load_case)
      cell_time_func{i} = options.forces.time_function;
    endfor

    options.forces.time_function = cell_time_func;
    clear cell_time_func;
  endif

  if (~isfield(options, "joints"))
    options.joints = struct();
  endif

  if (~isfield(options.joints, "number"))
    options.joints.number = 0;
  endif

  if (~isfield(options, "surface_loads"))
    options.surface_loads = struct();
  endif

  if (~isfield(options.surface_loads, "number"))
    options.surface_loads.number = 0;
  endif

  if (~isfield(options.surface_loads, "time_function"))
    options.surface_loads.time_function = "time";
  endif

  if (~iscell(options.surface_loads.time_function))
    if (~ischar(options.surface_loads.time_function))
      error("class of options.surface_loads.time_function must be cell or char");
    endif

    cell_time_func = cell(1, numel(load_case));

    for i=1:numel(load_case)
      cell_time_func{i} = options.surface_loads.time_function;
    endfor

    options.surface_loads.time_function = cell_time_func;

    clear cell_time_func;
  endif

  if (~isfield(options, "drive_callers"))
    options.drive_callers = struct();
  endif

  if (~isfield(options.drive_callers, "number"))
    options.drive_callers.number = 0;
  endif

  fd = -1;

  unwind_protect
    [fd, msg] = fopen(elem_file, "wt");

    if (fd == -1)
      error("failed to open file \"%s\": %s", elem_file, msg);
    endif

    elem_type_solid = {"iso8", "iso8upc", "iso20", "iso20upc", "iso20upcr", "iso20r", "iso27", "iso27upc", "penta15", "penta15upc", "tet10h", "tet10upc"};

    for i=1:numel(elem_type_solid)
      if (~isfield(mesh.elements, elem_type_solid{i}))
        continue;
      endif

      elem_node_idx_upc = zeros(1, 0, "int32");

      switch (elem_type_solid{i})
        case "iso8"
          elem_name = "hexahedron8";
          num_colloc_points = 2^3;
        case "iso8upc"
          elem_name = "hexahedron8upc";
          num_colloc_points = 2^3;
          elem_node_idx_upc = int32(1);
        case "iso20"
          elem_name = "hexahedron20";
          num_colloc_points = 3^3;
        case "iso20upc"
          elem_name = "hexahedron20upc";
          num_colloc_points = 3^3;
          elem_node_idx_upc = int32(1:8);
        case "iso20upcr"
          elem_name = "hexahedron20upcr";
          num_colloc_points = 2^3;
          elem_node_idx_upc = int32(1:8);
        case "iso27"
          elem_name = "hexahedron27";
          num_colloc_points = 3^3;
        case "iso27upc"
          elem_name = "hexahedron27upc";
          num_colloc_points = 3^3;
          elem_node_idx_upc = int32(1:27);
        case "iso20r"
          elem_name = "hexahedron20r";
          num_colloc_points = 2^3;
        case "penta15"
          elem_name = "pentahedron15";
          num_colloc_points = 7 * 3;
        case "penta15upc"
          elem_name = "pentahedron15upc";
          num_colloc_points = 7 * 3;
          elem_node_idx_upc = int32(1:6);
        case "tet10h"
          elem_name = "tetrahedron10";
          num_colloc_points = 5;
        case "tet10upc"
          elem_name = "tetrahedron10upc";
          num_colloc_points = 5;
          elem_node_idx_upc = int32(1:4);
      endswitch

      elem_nodes = getfield(mesh.elements, elem_type_solid{i}) + (options.struct_nodes.number - rows(mesh.nodes));
      elem_mat = getfield(mesh.materials, elem_type_solid{i});

      rho = zeros(columns(elem_nodes), rows(elem_nodes));

      for j=1:columns(elem_nodes)
        rho(j, :) = [mesh.material_data(elem_mat).rho];
      endfor

      ## Assume that the same node labels are used for UPC pressure nodes
      elem_nodes_upc = [elem_nodes, elem_nodes(:, elem_node_idx_upc)];

      ## Need to convert everything to double; otherwise we are truncating rho!
      elem_data = [double((1:rows(elem_nodes)) + options.solids.number);
                   double(elem_nodes_upc).';
                   rho;
                   double(repmat(elem_mat.' + (options.const_laws.number - numel(mesh.material_data)), num_colloc_points, 1))];

      if (~isfield(options.solids, elem_type_solid{i}))
        options.solids = setfield(options.solids, elem_type_solid{i}, struct());
      endif

      if (~isfield(getfield(options.solids, elem_type_solid{i}), "flags"))
        options.solids = setfield(options.solids, elem_type_solid{i}, setfield(getfield(options.solids, elem_type_solid{i}), "flags", zeros(rows(elem_nodes), 1, "uint32")));
      endif

      elem_flags = getfield(options.solids, elem_type_solid{i}).flags;

      elem_flag_names = {"", ", static", ", lumped mass"};
      elem_flag_bits = [uint32(0), MBDYN_ELEMENT_FLAG_STATIC, MBDYN_ELEMENT_FLAG_LUMPED_MASS];

      for l=1:numel(elem_flag_names)
        idx_elem_flag = find(elem_flags == elem_flag_bits(l));

        if (isempty(idx_elem_flag))
          continue;
        endif

        elem_format = sprintf("%s: %%d%s%s%s%s;\n", ...
                              elem_name, ...
                              elem_flag_names{l},
                              repmat(", %d", 1, columns(elem_nodes_upc)), ...
                              repmat(", %.16e", 1, columns(elem_nodes)), ...
                              repmat(", reference, %d", 1, num_colloc_points));

        fprintf(fd, elem_format, elem_data(:, idx_elem_flag));
      endfor
      options.solids.number += rows(elem_nodes);
    endfor

    if (options.solids.number == 0)
      warning("file \"%s\": no solid elements were generated", elem_file);
    endif

    [nidx, dofidx] = find(load_case_dof.locked_dof(:, 1:3)); ## FIXME: cannot use genels for rotation parameters

    if (~isempty(nidx))
      genel_format = "genel: %d, clamp, %d, structural, %d, algebraic, from node;\n";

      fprintf(fd, genel_format, [(1:rows(nidx)) + options.genels.number; nidx.'; dofidx.';]);

      options.genels.number += rows(nidx);
    endif

    if (isfield(mesh.elements, "rbe3"))
      for i=1:numel(mesh.elements.rbe3)
        fprintf(fd, "joint: %d, rigid body displacement joint,\n", ++options.joints.number);
        fprintf(fd, "\t%d, %d", mesh.elements.rbe3(i).nodes(1), numel(mesh.elements.rbe3(i).nodes) - 1);

        for j=1:numel(mesh.elements.rbe3(i).nodes) - 1
          fprintf(fd, ", %d, weight, %.16e", mesh.elements.rbe3(i).nodes(j + 1), mesh.elements.rbe3(i).weight(j));
        endfor

        fprintf(fd, ";\n\n");
      endfor
    endif

    if (isfield(mesh.elements, "rbe2"))
      for i=1:numel(mesh.elements.rbe2)
        for j=2:numel(mesh.elements.rbe2(i).nodes)
          fprintf(fd, "joint: %d, offset displacement joint, %d, %d;\n", ++options.joints.number, mesh.elements.rbe2(i).nodes(1), mesh.elements.rbe2(i).nodes(j));
        endfor
      endfor
    endif

    if (isfield(load_case, "loads"))
      for j=1:numel(load_case)
        if (isempty(load_case(j).loads))
          continue;
        endif
        fprintf(fd, "drive caller: %d, name, \"forces for load_case(%d)\", %s;\n", ++options.drive_callers.number, j, options.forces.time_function{j});
        force_format = sprintf("force: %%d, abstract, %%d, structural, %%d, mult, const, %%.16e, reference, %d;\n", options.drive_callers.number);

        for i=1:6
          cond_load = load_case(j).loads(:, i);

          if (i > 3)
            loaded_node_type = options.struct_nodes.type(load_case(j).loaded_nodes);
            cond_load &= ((loaded_node_type == MBDYN_NODE_TYPE_STATIC_STRUCT) | (loaded_node_type == MBDYN_NODE_TYPE_DYNAMIC_STRUCT));
          endif

          ridx = find(cond_load);

          if (isempty(ridx))
            continue;
          endif

          force_data = [(1:rows(ridx)) + options.forces.number;
                        double(load_case(j).loaded_nodes(ridx).');
                        repmat(i, 1, rows(ridx));
                        load_case(j).loads(ridx, i).'];

          fprintf(fd, force_format, force_data);
          options.forces.number += rows(ridx);
        endfor
      endfor
    endif

    if (isfield(load_case, "pressure") || isfield(load_case, "traction"))
      elem_type_surf_load = {"iso4", "quad8", "quad9", "tria6h", "quad8r"};
      surf_load_type = {"pressure", "traction", "traction_abs"};
      surf_load_elem_prefix = {"pressure", "traction", "traction"};

      for j=1:numel(load_case)
        fprintf(fd, "drive caller: %d, name, \"time variation for surfaces loads in load_case(%d)\", %s;\n", ++options.drive_callers.number, j, options.surface_loads.time_function{j});

        for l=1:numel(surf_load_type)
          if (~isfield(load_case(j), surf_load_type{l}))
            continue;
          endif

          for i=1:numel(elem_type_surf_load)
            if (~isfield(getfield(load_case(j), surf_load_type{l}), elem_type_surf_load{i}))
              continue;
            endif

            switch (elem_type_surf_load{i})
              case "iso4"
                elem_name = "q4";
                num_nodes = 4;
              case "quad8"
                elem_name = "q8";
                num_nodes = 8;
              case "quad9"
                elem_name = "q9";
                num_nodes = 9;
              case "quad8r"
                elem_name = "q8r";
                num_nodes = 8;
              case "tria6h"
                elem_name = "t6";
                num_nodes = 6;
            endswitch

            elem_name = [surf_load_elem_prefix{l}, elem_name];

            elem_data = getfield(getfield(load_case(j), surf_load_type{l}), elem_type_surf_load{i});
            elem_nodes = elem_data.elements + (options.struct_nodes.number - rows(mesh.nodes));

            switch (surf_load_type{l})
              case "pressure"
                elem_press = elem_data.p;

                elem_format = sprintf("%s: %%d%s, from drives%s;\n", ...
                                      elem_name, ...
                                      repmat(", %d", 1, columns(elem_nodes)), ...
                                      repmat(sprintf(", mult, const, %%.16e, reference, %d", ...
                                                     options.drive_callers.number), 1, columns(elem_nodes)));

                ## Need to convert everything to double; otherwise we are truncating the pressure!
                elem_data = [double((1:rows(elem_nodes)) + options.surface_loads.number);
                             double(elem_nodes.');
                             elem_press.'];
              case {"traction", "traction_abs"}
                switch (surf_load_type{l})
                  case "traction_abs"
                    traction_type = "absolute";
                  otherwise
                    traction_type = "";
                endswitch

                if (~isempty(traction_type))
                  traction_type = [", ", traction_type];
                endif

                drv_caller_comp_fmt = sprintf(", mult, const, %%.16e, reference, %d", options.drive_callers.number);

                elem_format = sprintf("%s: %%d%s%s, from drives%s;\n", ...
                                      elem_name, ...
                                      traction_type, ...
                                      repmat(", %d", 1, columns(elem_nodes)), ...
                                      repmat([", component", repmat(drv_caller_comp_fmt, 1, 3)], 1, columns(elem_nodes)));

                f = zeros(3 * columns(elem_nodes), rows(elem_nodes));

                for k=1:3
                  f((1:3:3*columns(elem_nodes)) + k - 1, :) = elem_data.f(:, :, k).';
                endfor

                elem_data = [double((1:rows(elem_nodes)) + options.surface_loads.number);
                             double(elem_nodes.');
                             f];
              otherwise
                error("surface load type \"%s\" not implemented", surf_load_type{l});
            endswitch

            fprintf(fd, elem_format, elem_data);

            options.surface_loads.number += rows(elem_nodes);
          endfor
        endfor
      endfor
    endif
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif

    fd = -1;
  end_unwind_protect
endfunction
