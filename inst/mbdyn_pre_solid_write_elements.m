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

%!test
%! ## TEST 1
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.E = 55e6 / SI_unit_pascal;
%! param.delta = 0.0;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.p1 = 0.006e6 / SI_unit_pascal;
%! param.bound_cond = "edge";
%! param.analysis = "plain strain";
%! param.transfinite = true;
%! options.verbose = false;
%! elem_types = {"tet10h", "tet10upc", "iso8upc", "penta15", "penta15upc", "iso20", "iso20upc", "iso20upcr", "iso20r", "iso27"};
%! nu_val = [0.3, 0.499];
%! for idx_elem_type=1:numel(elem_types)
%! param.elem_type = elem_types{idx_elem_type};
%! for idx_nu = 1:numel(nu_val)
%! param.nu = nu_val(idx_nu);
%! switch (idx_nu)
%! case 2
%! switch (param.elem_type)
%! case {"tet10h", "iso8", "penta15", "iso20", "iso27"}
%!   continue;
%! endswitch
%! endswitch
%! switch (param.elem_type)
%! case {"iso8", "iso8upc"}
%! param.h = 10e-3 / 32 / SI_unit_meter;
%! otherwise
%! param.h = 10e-3 / 16 / SI_unit_meter;
%! endswitch
%! switch (idx_nu)
%! case 1
%!   param.N = 3;
%! otherwise
%! switch (param.elem_type)
%! case {"tet10upc", "iso8upc", "penta15upc", "iso20upc", "iso20upcr"};
%! param.N = 3;
%! otherwise
%! param.N = 20;
%! endswitch
%! endswitch
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fprintf(fd, "SetFactory(\"Built-in\");\n");
%!     fprintf(fd, "h = %g;\n", param.h);
%!     fprintf(fd, "Point(1) = {0,0,20};\n");
%!     fprintf(fd, "Point(2) = {5,0,20};\n");
%!     fprintf(fd, "Point(3) = {10,0,20};\n");
%!     fprintf(fd, "Point(4) = {10,0,10};\n");
%!     fprintf(fd, "Point(5) = {15 - 5 * Cos(Pi/4), 0, 10 - 5 * Sin(Pi/4)};\n");
%!     fprintf(fd, "Point(6) = {15, 0, 5};\n");
%!     fprintf(fd, "Point(7) = {65, 0, 5};\n");
%!     fprintf(fd, "Point(8) = {65, 0, 0};\n");
%!     fprintf(fd, "Point(9) = {65, 0, -5};\n");
%!     fprintf(fd, "Point(10) = {15, 0, -5};\n");
%!     fprintf(fd, "Point(11) = {15 - 5 * Cos(Pi/4), 0, -(10 - 5 * Sin(Pi/4))};\n");
%!     fprintf(fd, "Point(12) = {10, 0, -10};\n");
%!     fprintf(fd, "Point(13) = {10,0,-20};\n");
%!     fprintf(fd, "Point(14) = {5, 0, -20};\n");
%!     fprintf(fd, "Point(15) = {0, 0, -20};\n");
%!     fprintf(fd, "Point(16) = {0, 0, -10};\n");
%!     fprintf(fd, "Point(17) = {0, 0, 0};\n");
%!     fprintf(fd, "Point(18) = {0, 0, 10};\n");
%!     fprintf(fd, "Point(19) = {5, 0, 10};\n");
%!     fprintf(fd, "Point(20) = {5, 0, 0};\n");
%!     fprintf(fd, "Point(21) = {15, 0, 0};\n");
%!     fprintf(fd, "Point(22) = {5, 0, -10};\n");
%!     fprintf(fd, "Point(23) = {15,0,10};\n");
%!     fprintf(fd, "Point(25) = {15,0,-10};\n");
%!     fprintf(fd, "Line(1) = {1,2};\n");
%!     fprintf(fd, "Line(2) = {2,3};\n");
%!     fprintf(fd, "Line(3) = {3,4};\n");
%!     fprintf(fd, "Circle(4) = {4,23,5};\n");
%!     fprintf(fd, "Circle(5) = {5,23,6};\n");
%!     fprintf(fd, "Line(6) = {6,7};\n");
%!     fprintf(fd, "Line(7) = {7,8};\n");
%!     fprintf(fd, "Line(8) = {8,9};\n");
%!     fprintf(fd, "Line(9) = {9,10};\n");
%!     fprintf(fd, "Circle(10)={10,25,11};\n");
%!     fprintf(fd, "Circle(11)={11,25,12};\n");
%!     fprintf(fd, "Line(12) = {12,13};\n");
%!     fprintf(fd, "Line(13) = {13,14};\n");
%!     fprintf(fd, "Line(14) = {14,15};\n");
%!     fprintf(fd, "Line(15) = {15,16};\n");
%!     fprintf(fd, "Line(16) = {16,17};\n");
%!     fprintf(fd, "Line(17) = {17,18};\n");
%!     fprintf(fd, "Line(18) = {18,1};\n");
%!     fprintf(fd, "Line(19) = {2,19};\n");
%!     fprintf(fd, "Line(20) = {19,20};\n");
%!     fprintf(fd, "Line(21) = {20,22};\n");
%!     fprintf(fd, "Line(22) = {22,14};\n");
%!     fprintf(fd, "Line(23) = {18,19};\n");
%!     fprintf(fd, "Line(24) = {19, 4};\n");
%!     fprintf(fd, "Line(25) = {17,20};\n");
%!     fprintf(fd, "Line(26) = {16,22};\n");
%!     fprintf(fd, "Line(27) = {22,12};\n");
%!     fprintf(fd, "Line(28) = {20,21};\n");
%!     fprintf(fd, "Line(29) = {21,8};\n");
%!     fprintf(fd, "Line(30) = {6,21};\n");
%!     fprintf(fd, "Line(31) = {21,10};\n");
%!     fprintf(fd, "Line(32) = {20,5};\n");
%!     fprintf(fd, "Line(33) = {20,11};\n");
%!     if (param.transfinite)
%!       fprintf(fd, "Transfinite Curve(1) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(2) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(3) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(4) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(5) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(6) = Max(1, Round(50 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(7) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(8) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(9) = Max(1, Round(50 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(10) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(11) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(12) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(13) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(14) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(15) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(16) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(17) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(18) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(19) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(20) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(21) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(22) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(23) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(24) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(25) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(26) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(27) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(28) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(29) = Max(1, Round(50 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(30) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(31) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(32) = Max(1, Round(5 / h)) + 1;\n");
%!       fprintf(fd, "Transfinite Curve(33) = Max(1, Round(5 / h)) + 1;\n");
%!     endif
%!     fprintf(fd, "Line Loop(1) = {1,19,-23,18};\n");
%!     fprintf(fd, "Line Loop(2) = {2,3,-24,-19};\n");
%!     fprintf(fd, "Line Loop(3) = {23,20,-25,17};\n");
%!     fprintf(fd, "Line Loop(4) = {24,4,-32,-20};\n");
%!     fprintf(fd, "Line Loop(5) = {32,5,30,-28};\n");
%!     fprintf(fd, "Line Loop(6) = {6,7,-29,-30};\n");
%!     fprintf(fd, "Line Loop(7) = {25,21,-26,16};\n");
%!     fprintf(fd, "Line Loop(8) = {33,11,-27,-21};\n");
%!     fprintf(fd, "Line Loop(9) = {28,31,10,-33};\n");
%!     fprintf(fd, "Line Loop(10) = {29,8,9,-31};\n");
%!     fprintf(fd, "Line Loop(11) = {26, 22, 14, 15};\n");
%!     fprintf(fd, "Line Loop(12) = {27, 12, 13, -22};\n");
%!     fprintf(fd, "Plane Surface(1) = {1};\n");
%!     fprintf(fd, "Plane Surface(2) = {2};\n");
%!     fprintf(fd, "Plane Surface(3) = {3};\n");
%!     fprintf(fd, "Plane Surface(4) = {4};\n");
%!     fprintf(fd, "Plane Surface(5) = {5};\n");
%!     fprintf(fd, "Plane Surface(6) = {6};\n");
%!     fprintf(fd, "Plane Surface(7) = {7};\n");
%!     fprintf(fd, "Plane Surface(8) = {8};\n");
%!     fprintf(fd, "Plane Surface(9) = {9};\n");
%!     fprintf(fd, "Plane Surface(10) = {10};\n");
%!     fprintf(fd, "Plane Surface(11) = {11};\n");
%!     fprintf(fd, "Plane Surface(12) = {12};\n");
%!     if (param.transfinite)
%!       fprintf(fd, "Transfinite Surface(1) = {PointsOf{Surface{1};}};\n");
%!       fprintf(fd, "Transfinite Surface(2) = {PointsOf{Surface{2};}};\n");
%!       fprintf(fd, "Transfinite Surface(3) = {PointsOf{Surface{3};}};\n");
%!       fprintf(fd, "Transfinite Surface(4) = {PointsOf{Surface{4};}};\n");
%!       fprintf(fd, "Transfinite Surface(5) = {PointsOf{Surface{5};}};\n");
%!       fprintf(fd, "Transfinite Surface(6) = {PointsOf{Surface{6};}};\n");
%!       fprintf(fd, "Transfinite Surface(7) = {PointsOf{Surface{7};}};\n");
%!       fprintf(fd, "Transfinite Surface(8) = {PointsOf{Surface{8};}};\n");
%!       fprintf(fd, "Transfinite Surface(9) = {PointsOf{Surface{9};}};\n");
%!       fprintf(fd, "Transfinite Surface(10) = {PointsOf{Surface{10};}};\n");
%!       fprintf(fd, "Transfinite Surface(11) = {PointsOf{Surface{11};}};\n");
%!       fprintf(fd, "Transfinite Surface(12) = {PointsOf{Surface{12};}};\n");
%!     endif
%!     switch (param.elem_type)
%!     case {"tet10h", "tet10upc"}
%!       extrude_opt = "";
%!     otherwise
%!       extrude_opt = "Layers{1}; Recombine;";
%!     endswitch
%!     fprintf(fd, "v1[] = Extrude {0, h, 0}{ Surface{1}; %s };\n", extrude_opt);
%!     fprintf(fd, "v2[] = Extrude {0, h, 0}{ Surface{2}; %s };\n", extrude_opt);
%!     fprintf(fd, "v3[] = Extrude {0, h, 0}{ Surface{3}; %s };\n", extrude_opt);
%!     fprintf(fd, "v4[] = Extrude {0, h, 0}{ Surface{4}; %s };\n", extrude_opt);
%!     fprintf(fd, "v5[] = Extrude {0, h, 0}{ Surface{5}; %s };\n", extrude_opt);
%!     fprintf(fd, "v6[] = Extrude {0, h, 0}{ Surface{6}; %s };\n", extrude_opt);
%!     fprintf(fd, "v7[] = Extrude {0, h, 0}{ Surface{7}; %s };\n", extrude_opt);
%!     fprintf(fd, "v8[] = Extrude {0, h, 0}{ Surface{8}; %s };\n", extrude_opt);
%!     fprintf(fd, "v9[] = Extrude {0, h, 0}{ Surface{9}; %s };\n", extrude_opt);
%!     fprintf(fd, "v10[] = Extrude {0, h, 0}{ Surface{10}; %s };\n", extrude_opt);
%!     fprintf(fd, "v11[] = Extrude {0, h, 0}{ Surface{11}; %s };\n", extrude_opt);
%!     fprintf(fd, "v12[] = Extrude {0, h, 0}{ Surface{12}; %s };\n", extrude_opt);
%!     switch (param.elem_type)
%!     case {"iso8", "iso8upc", "iso20", "iso20upc", "iso20upcr", "iso20r", "iso27"}
%!       fprintf(fd, "Recombine Surface{1, v1[0]};\n");
%!       fprintf(fd, "Recombine Surface{2, v2[0]};\n");
%!       fprintf(fd, "Recombine Surface{3, v3[0]};\n");
%!       fprintf(fd, "Recombine Surface{4, v4[0]};\n");
%!       fprintf(fd, "Recombine Surface{5, v5[0]};\n");
%!       fprintf(fd, "Recombine Surface{6, v6[0]};\n");
%!       fprintf(fd, "Recombine Surface{7, v7[0]};\n");
%!       fprintf(fd, "Recombine Surface{8, v8[0]};\n");
%!       fprintf(fd, "Recombine Surface{9, v9[0]};\n");
%!       fprintf(fd, "Recombine Surface{10, v10[0]};\n");
%!       fprintf(fd, "Recombine Surface{11, v11[0]};\n");
%!       fprintf(fd, "Recombine Surface{12, v12[0]};\n");
%!     endswitch
%!     fprintf(fd, "Physical Volume(\"volume\", 1) = {1, 2, 3, 4, 7, 8, 9, 5, 10, 6, 11, 12};\n");
%!     fprintf(fd, "Physical Surface(\"clamp\", 2) = {54, 98, 186, 274};\n");
%!     fprintf(fd, "Physical Surface(\"pressure\", 3) = {68, 112, 134, 152};\n");
%!     fprintf(fd, "Physical Surface(\"displacement\", 4) = {156, 244};\n");
%!     fprintf(fd, "Physical Surface(\"stress\", 5) = {112, 134};\n");
%!     fprintf(fd, "Physical Surface(\"top\", 6) = {1,2,3,4,5,6,7,8,9,10,11,12};\n");
%!     fprintf(fd, "Physical Surface(\"bottom\", 7) = {v1[0],v2[0],v3[0],v4[0],v5[0],v6[0],v7[0],v8[0],v9[0],v10[0],v11[0],v12[0]};\n");
%!     fprintf(fd, "Physical Curve(\"clampedge\", 8) = {40, 269};\n");
%!     switch (param.elem_type)
%!     case {"tet10h", "tet10upc"}
%!       f_use_mesh_size = true;
%!     otherwise
%!       f_use_mesh_size = ~param.transfinite;
%!     endswitch
%!     if (f_use_mesh_size)
%!       fputs(fd, "MeshSize{PointsOf{Volume{1,2,3,4,5,6,7,8,9,10,11,12};}} = h;\n");
%!     endif
%!     switch (param.elem_type)
%!     case {"penta15", "penta15upc", "iso20", "iso20upc", "iso20upcr", "iso20r"}
%!       fprintf(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     otherwise
%!       fprintf(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!     endswitch
%!     if (~param.transfinite)
%!     switch (param.elem_type)
%!     case {"tet10h", "tet10upc", "penta15", "penta15upc"}
%!       fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!     endswitch
%!     endif
%!     switch (param.elem_type)
%!     case {"iso8", "iso8upc"}
%!       fprintf(fd, "Mesh.ElementOrder = 1;\n");
%!     otherwise
%!       fprintf(fd, "Mesh.ElementOrder = 2;\n");
%!     endswitch
%!     fprintf(fd, "Coherence Mesh;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   [~] = unlink([filename, ".msh"]);
%!   #spawn_wait(spawn("gmsh", {[filename, ".geo"]}));
%!   pid = spawn("gmsh", {"-format", "msh2", "-3", [filename, ".geo"]});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   [~] = unlink([filename, ".geo"]);
%!   switch (param.elem_type)
%!   case {"tet10h", "tet10upc"}
%!     param.elem_type_surf = {"tria6h"};
%!   case {"iso8", "iso8upc"}
%!     param.elem_type_surf = {"iso4"};
%!   case {"iso20", "iso20upc"}
%!     param.elem_type_surf = {"quad8"};
%!   case {"iso20upcr", "iso20r"}
%!     param.elem_type_surf = {"quad8r"};
%!   case {"iso27"}
%!     param.elem_type_surf = {"quad9"};
%!   case {"penta15", "penta15upc"}
%!     param.elem_type_surf = {"quad8", "tria6h"};
%!   endswitch
%!   switch (param.elem_type)
%!   case {"iso8", "iso8upc"}
%!     param.elem_type_line = "line2";
%!   otherwise
%!     param.elem_type_line = "line3";
%!   endswitch
%!   param.material = "hookean linear elastic isotropic";
%!   opt_msh.elem_type = {param.elem_type, param.elem_type_surf{:}, "line2", "line3", "point1"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   [~] = unlink([filename, ".msh"]);
%!   load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!   load_case.pressure = struct();
%!   mesh.materials = struct();
%!   for i=1:numel(param.elem_type_surf)
%!   switch (param.bound_cond)
%!   case "surface"
%!     grp_id_clamp = find([[getfield(mesh.groups, param.elem_type_surf{i})].id] == 2);
%!     if (~isempty(grp_id_clamp))
%!       load_case_dof.locked_dof(getfield(mesh.groups, param.elem_type_surf{i})(grp_id_clamp).nodes, 1:3) = true;
%!     endif
%!   case "edge"
%!     grp_id_clamp = find([[getfield(mesh.groups, param.elem_type_line)].id] == 8);
%!     if (~isempty(grp_id_clamp))
%!       load_case_dof.locked_dof(getfield(mesh.groups, param.elem_type_line)(grp_id_clamp).nodes, 1:3) = true;
%!     endif
%!   endswitch
%!   switch (param.analysis)
%!   case "plain strain"
%!   if (isfield(mesh.groups, param.elem_type_surf{i}))
%!   grp_id_top = find([[getfield(mesh.groups, param.elem_type_surf{i})].id] == 6);
%!   grp_id_bottom = find([[getfield(mesh.groups, param.elem_type_surf{i})].id] == 7);
%!   if (~isempty(grp_id_top))
%!     load_case_dof.locked_dof(getfield(mesh.groups, param.elem_type_surf{i})(grp_id_top).nodes, 2) = true;
%!   endif
%!   if (~isempty(grp_id_bottom))
%!     load_case_dof.locked_dof(getfield(mesh.groups, param.elem_type_surf{i})(grp_id_bottom).nodes, 2) = true;
%!   endif
%!   endif
%!   endswitch
%!   grp_id_p1 = find([[getfield(mesh.groups, param.elem_type_surf{i})].id] == 3);
%!   if (~isempty(grp_id_p1))
%!   elem_id_p1 = getfield(mesh.groups, param.elem_type_surf{i})(grp_id_p1).elements;
%!   elno_p1 = getfield(mesh.elements, param.elem_type_surf{i})(elem_id_p1, :);
%!   press_load.elements = elno_p1;
%!   press_load.p = [repmat(param.p1, rows(elno_p1), columns(elno_p1))];
%!   load_case.pressure = setfield(load_case.pressure, param.elem_type_surf{i}, press_load);
%!   endif
%!   endfor
%!   mesh.materials = setfield(mesh.materials, param.elem_type, ones(rows(getfield(mesh.elements, param.elem_type)), 1, "int32"));
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).delta = param.delta;
%!   mesh.material_data(1).type = param.material;
%!   mesh.material_data(1).rho = param.rho;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!   if (mesh.material_data(1).nu == 0.5)
%!     ++opt_mbd_mesh.genels.number;
%!   endif
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fprintf(fd, " set: real N = %d;\n", param.N);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-6, test, norm, 1e-6, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 20, linear solver max iterations, 40, minimum step, 1e-12, recovery step, 1e-12, verbose, 3;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, t1, forever, t1 / 10;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    hydraulic nodes: %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!     fprintf(fd, "    surface loads: %d;\n", opt_mbd_mesh.surface_loads.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "    genels: %d;\n", opt_mbd_mesh.genels.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     if (mesh.material_data(1).nu == 0.5)
%!       ## Because the pressure will have no impact on strain if the material is fully incompressible
%!       fprintf(fd, "genel: %d, clamp, %d, hydraulic, null;\n", opt_mbd_mesh.genels.number, getfield(mesh.elements, param.elem_type)(1,1));
%!     endif
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   opt_mbdyn.mbdyn_command = "mbdyn -C";
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol_stat] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   delta = nan(1, numel(param.elem_type_surf));
%!   sigma1_max = nan(1, numel(param.elem_type_surf));
%!   for k=1:numel(param.elem_type_surf)
%!   grp_id_displacement = find([[getfield(mesh.groups, param.elem_type_surf{k})].id] == 4);
%!   if (~isempty(grp_id_displacement))
%!   if (isfield(mesh.groups, param.elem_type_surf{k}))
%!     elem_id_displacement = getfield(mesh.groups, param.elem_type_surf{k})(grp_id_displacement).elements;
%!     elno_id_displacement = getfield(mesh.elements, param.elem_type_surf{k})(elem_id_displacement, :);
%!     delta(k) = mean(sol_stat.def(elno_id_displacement, 3, end));
%!   endif
%!   endif
%!   grp_id_stress = find([[getfield(mesh.groups, param.elem_type_surf{k})].id] == 5);
%!   if(~isempty(grp_id_stress))
%!   elem_id_stress = getfield(mesh.groups, param.elem_type_surf{k})(grp_id_stress).elements;
%!   elno_id_stress = getfield(mesh.elements, param.elem_type_surf{k})(elem_id_stress, :);
%!   taum = zeros(6, numel(elno_id_stress));
%!   taum_n = zeros(1, numel(elno_id_stress));
%!   for i=1:numel(elno_id_stress)
%!     [ridx, cidx] = find(getfield(mesh.elements, param.elem_type) == elno_id_stress(i));
%!     for j=1:numel(ridx)
%!       taum(:, i) += reshape(getfield(sol_stat.stress.taum, param.elem_type)(ridx(j), cidx(j), :, end), 6, 1);
%!       ++taum_n(i);
%!     endfor
%!   endfor
%!   taum *= diag(1 ./ taum_n);
%!   for i=1:columns(taum)
%!     TAU = [taum(1, i), taum(4, i), taum(6, i);
%!            taum(4, i), taum(2, i), taum(5, i);
%!            taum(6, i), taum(5, i), taum(3, i)];
%!     sigma1_max(k) = max(sigma1_max(k), max(eig(TAU)));
%!   endfor
%!   endif
%!   endfor
%!   delta = mean(delta(isfinite(delta)));
%!   sigma1_max = mean(sigma1_max(isfinite(sigma1_max)));
%!   fprintf(stderr, "mesh size=%.1f\n", param.h);
%!   fprintf(stderr, "max(sigma1)=%.3f [MPa]\n", sigma1_max);
%!   fprintf(stderr, "delta=%.3f [mm]\n", delta);
%!   switch (param.nu)
%!   case 0.3
%!     ## K.J.Bathe page 329 chaper 4.4
%!     sigma1_max_ref = 0.6056e6 / SI_unit_pascal;
%!     delta_ref = -1.669e-3 / SI_unit_meter;
%!   case 0.499
%!     ## K.J.Bathe page 333 chaper 4.4
%!     sigma1_max_ref = 0.5998e6 / SI_unit_pascal;
%!     delta_ref = -1.393e-3 / SI_unit_meter;
%!   endswitch
%!   tol_sigma = 2.5e-2;
%!   tol_delta = 1.25e-2;
%!   fprintf(stderr, "difference(sigam1_max)=%.2f%%\n", (sigma1_max / sigma1_max_ref - 1) * 100);
%!   fprintf(stderr, "difference(delta)=%.2f%%\n", (delta / delta_ref - 1) * 100);
%!   assert_simple(sigma1_max, sigma1_max_ref, tol_sigma * abs(sigma1_max_ref));
%!   assert_simple(delta, delta_ref, tol_delta * abs(delta_ref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
%! endfor
%! endfor

%!test
%! ## TEST 2
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.113)
%! pkg load mboct-fem-pkg;
%! close all;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.C10 = param.G / 2;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 2.5e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"quad8", "iso20", "penta15", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"quad8", "tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"quad8", "tria6h"});
%!   if (isfield(mesh.elements, "iso20"))
%!     mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!     mesh.materials.iso20([mesh.groups.iso20(find([[mesh.groups.iso20.id] == 1])).elements]) = 1;
%!   endif
%!   if (isfield(mesh.elements, "penta15"))
%!     mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!     mesh.materials.penta15([mesh.groups.penta15(find([[mesh.groups.penta15.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * param.C10 * D; ## (4.113)
%!   Nref = -1/2 * pi * Ra^4 * param.C10 * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-3;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 3
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.113)
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.C10 = param.G / 2;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 0.8333e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"iso4", "iso8", "tria3", "penta6"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"iso4"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"iso4"});
%!   if (isfield(mesh.elements, "iso8"))
%!     mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!     mesh.materials.iso8([mesh.groups.iso8(find([[mesh.groups.iso8.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * param.C10 * D; ## (4.113)
%!   Nref = -1/2 * pi * Ra^4 * param.C10 * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 1.5e-2;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 4
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.113)
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.C10 = param.G / 2;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 1.25e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fprintf(fd, "h = %e;\n", param.h);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "Point(2) = {l1, 0, 0};\n");
%!     fputs(fd, "Point(3) = {l1, 0.5 * d1, 0};\n");
%!     fputs(fd, "Point(4) = {0, 0.5 * d1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Curve Loop(1) = {1, 2, 3, 4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "v[] = Extrude {{1, 0, 0}, {0, 0, 0}, 2*Pi}{ Surface{1}; };\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{v[1]}; } } = h;\n");
%!     fputs(fd, "Physical Volume(1) = {v[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {4};\n");
%!     fputs(fd, "Physical Surface(2) = {2};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"tria6h", "tet10h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"tria6h"});
%!   if (isfield(mesh.elements, "tet10h"))
%!     mesh.materials.tet10h = zeros(rows(mesh.elements.tet10h), 1, "int32");
%!     mesh.materials.tet10h([mesh.groups.tet10h(find([[mesh.groups.tet10h.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 50;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"0.5 * pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * param.C10 * D; ## (4.113)
%!   Nref = -1/2 * pi * Ra^4 * param.C10 * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-3;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 5
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.113)
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.C10 = param.G / 2;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 2.5e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"quad8", "iso20r", "penta15", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"quad8", "tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"quad8", "tria6h"});
%!   if (isfield(mesh.elements, "iso20r"))
%!     mesh.materials.iso20r = zeros(rows(mesh.elements.iso20r), 1, "int32");
%!     mesh.materials.iso20r([mesh.groups.iso20r(find([[mesh.groups.iso20r.id] == 1])).elements]) = 1;
%!   endif
%!   if (isfield(mesh.elements, "penta15"))
%!     mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!     mesh.materials.penta15([mesh.groups.penta15(find([[mesh.groups.penta15.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * param.C10 * D; ## (4.113)
%!   Nref = -1/2 * pi * Ra^4 * param.C10 * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-3;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 6
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! SI_unit_joule = SI_unit_newton * SI_unit_meter;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 10e6 / SI_unit_pascal;
%! param.nu = 0.3;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.beta = 1e-3 / SI_unit_second;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.h = 0.625e-3/2 / SI_unit_meter;
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.d2 = param.d1 - 2 * param.h;
%! param.l1 = param.h;
%! param.z0 = 5e-3 / SI_unit_meter;
%! param.vz0 = 0 / (SI_unit_meter / SI_unit_second);
%! param.g = 9.81 / (SI_unit_meter / SI_unit_second^2);
%! param.Zeta = -80 * pi / 180;
%! param.t1 = 0.1 / SI_unit_second;
%! param.N = 100;
%! param.mus = 0.5;
%! param.muc = 0.5;
%! param.vs = 1 / (SI_unit_meter / SI_unit_second);
%! param.kv = 0 / (SI_unit_pascal / (SI_unit_meter / SI_unit_second));
%! param.i = 1;
%! param.sigma0x = 1e4 / SI_unit_meter^-1;
%! param.sigma0y = 1e4 / SI_unit_meter^-1;
%! param.sigma1x = 0 / (SI_unit_second * SI_unit_meter^-1);
%! param.sigma1y = 0 / (SI_unit_second * SI_unit_meter^-1);
%! param.m = (param.d1^2 - param.d2^2) * pi / 4 * param.l1 * param.rho;
%! param.gamma = 1e-16;
%! param.sn = param.m * param.g / (param.gamma * param.d1);
%! param.z0 -= param.m * param.g * sin(abs(param.Zeta)) / (4 * param.sn);
%! options.verbose = false;
%! options.do_plot = false;
%! options.var_step_size = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0,  0, z0 + 0.5 * d2};\n");
%!     fputs(fd, "Point(2) = {0, l1, z0 + 0.5 * d2};\n");
%!     fputs(fd, "Point(3) = {0, l1, z0 + 0.5 * d1};\n");
%!     fputs(fd, "Point(4) = {0,  0, z0 + 0.5 * d1};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Curve Loop(1) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(l1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round((d1 - d2) / 2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(l1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round((d1 - d2) / 2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Surface(1) = {1,2,3,4};\n");
%!     fputs(fd, "v1[] = Extrude {{0, 1, 0}, {0, 0, z0}, Pi}{ Surface{1}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{0, 1, 0}, {0, 0, z0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Recombine Surface{v1[0],v2[0]};\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {4,9};\n");
%!     fputs(fd, "Mesh.ElementOrder = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"iso8", "iso4", "tet4", "penta6", "tria3"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_id_ground = rows(mesh.nodes) + 1;
%!   mesh.nodes(node_id_ground, 1:3) = zeros(1, 3);
%!   mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!   mesh.materials.iso8([mesh.groups.iso8(find([[mesh.groups.iso8.id] == 1])).elements]) = 1;
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean viscoelastic";
%!   mesh.material_data(1).beta = param.beta;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_id_ground) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_sphere";
%!   opt_mbd_mesh.joints.number = 1;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   node_id_cont = mesh.groups.iso4(1).nodes;
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fputs(fd, "set: integer ref_id_sphere = 1000;\n");
%!     fputs(fd, "set: integer drive_id_E = 2000;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fx = 2001;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fy = 2002;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fz = 2003;\n");
%!     fprintf(fd, " set: integer node_id_ground = %d;\n", node_id_ground);
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     if (options.var_step_size)
%!       fputs(fd, "    min time step: t1 / N / 100;\n");
%!       fputs(fd, "    max time step: t1 / N;\n");
%!       fputs(fd, " strategy: factor, 0.5, 1, 1.25, 3, 3, 5;\n");
%!     endif
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-6, test, norm, 1e-6, test, norm;\n");
%!     fputs(fd, "    linear solver: pardiso, grad, max iterations, 100;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: mcp newton min fb;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    #output meter: closest next, 0., forever, t1 / 100;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "    loadable elements: %d;\n", numel(node_id_cont));
%!     fputs(fd, "      gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, "reference: ref_id_sphere,\n");
%!     fputs(fd, "  position, reference, global, null,\n");
%!     fputs(fd, "  orientation, reference, global, eye,\n");
%!     fputs(fd, "  velocity, reference, global, 0., 0., vz0,\n");
%!     fputs(fd, "  angular velocity, reference, global, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     for i=1:numel(node_id_cont)
%!       fprintf(fd, "user defined: %d, unilateral disp in plane,\n", i);
%!       fprintf(fd, "node1, %d,\n", node_id_cont(i));
%!       fputs(fd, "    offset, 1,\n");
%!       fputs(fd, "    reference, node, null,\n");
%!       fputs(fd, "    stiffness, sn,\n");
%!       fputs(fd, "    enable mcp, yes,\n");
%!       fputs(fd, " node2, node_id_ground,\n");
%!       fputs(fd, " offset, reference, node, null,\n");
%!       fputs(fd, " hinge, 3, 0, 0, 1,\n");
%!       fputs(fd, "        1, 1, 0, 0,\n");
%!       fputs(fd, " coulomb friction coefficient, muc,\n");
%!       fputs(fd, " static friction coefficient, mus,\n");
%!       fputs(fd, " sliding velocity coefficient, vs,\n");
%!       fputs(fd, " sliding velocity exponent, i,\n");
%!       fputs(fd, " micro slip stiffness x, sigma0x,\n");
%!       fputs(fd, " micro slip stiffness y, sigma0y,\n");
%!       fputs(fd, " micro slip damping x, sigma1x,\n");
%!       fputs(fd, " micro slip damping y, sigma1y,\n");
%!       fputs(fd, " viscous friction coefficient, kv;\n");
%!     endfor
%!     fputs(fd, "  joint: 1, clamp, node_id_ground, node, node;\n");
%!     fputs(fd, "  gravity: uniform, component, g * cos(Zeta), 0., g * sin(Zeta);\n");
%!     fprintf(fd, "drive caller: drive_id_E, array, %d", rows(mesh.elements.iso8));
%!     for i=1:rows(mesh.elements.iso8)
%!       fprintf(fd, ",\n    element, %d, solid, string, \"E\", direct", i);
%!     endfor
%!     fputs(fd, ",\n  output, yes;\n\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fx, element, 1, joint, string, \"Fx\", direct, output, yes;\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fy, element, 1, joint, string, \"Fy\", direct, output, yes;\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fz, element, 1, joint, string, \"Fz\", direct, output, yes;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [drive_id, drive_data] = mbdyn_post_load_output_drv(opt_mbdyn.output_file);
%!   Wkin_res = drive_data{1};
%!   Fx_res = drive_data{2};
%!   Fy_res = drive_data{3};
%!   Fz_res = drive_data{4};
%!   ra = 0.5 * param.d1;
%!   ri = 0.5 * param.d2;
%!   h = param.l1;
%!   m = param.rho * pi * (ra^2 - ri^2) * h;
%!   J = m * (ra^2 + ri^2) / 2;
%!   mred = m + J / ra^2;
%!   Phi = param.Zeta + pi / 2;
%!   qddot = m * param.g * sin(Phi) / mred;
%!   qdot = qddot * sol.t;
%!   q = 0.5 * qddot * sol.t.^2;
%!   Wkin_ref = 0.5 * mred * qdot.^2;
%!   Fz_ref = repmat(-m * param.g * cos(Phi), size(sol.t));
%!   Fx_ref = repmat(m * param.g * sin(Phi) * (1 - m / mred), size(sol.t));
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Wkin_ref * SI_unit_joule, "-x;Wkin_r_e_f;0");
%!     plot(sol.t * SI_unit_second, Wkin_res * SI_unit_joule, "-o;Wkin_r_e_s;1");
%!     xlabel("t [s]");
%!     ylabel("W [J]");
%!     grid on;
%!     grid minor on;
%!     title("kinetic energy versus time");
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Fz_ref * SI_unit_newton, "-x;Fz_r_e_f;0");
%!     plot(sol.t * SI_unit_second, Fz_res * SI_unit_newton, "-o;Fz_r_e_s;1");
%!     xlabel("t [s]");
%!     ylabel("F [N]");
%!     grid on;
%!     grid minor on;
%!     title("normal force");
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Fx_ref * SI_unit_newton, "-x;Fx_r_e_f;0");
%!     plot(sol.t * SI_unit_second, Fx_res * SI_unit_newton, "-o;Fx_r_e_s;1");
%!     xlabel("t [s]");
%!     ylabel("F [N]");
%!     grid on;
%!     grid minor on;
%!     title("tangent force");
%!   endif
%!   tol = 2e-2;
%!   assert_simple(Wkin_res, Wkin_ref, tol * max(abs(Wkin_res)));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 7
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! SI_unit_joule = SI_unit_newton * SI_unit_meter;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 10e6 / SI_unit_pascal;
%! param.nu = 0.3;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.beta = 1e-3 / SI_unit_second;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.h = 0.625e-3 / SI_unit_meter;
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.d2 = param.d1 - 2 * param.h;
%! param.l1 = param.h;
%! param.z0 = 5e-3 / SI_unit_meter;
%! param.vz0 = 0 / (SI_unit_meter / SI_unit_second);
%! param.g = 9.81 / (SI_unit_meter / SI_unit_second^2);
%! param.Zeta = -80 * pi / 180;
%! param.t1 = 0.1 / SI_unit_second;
%! param.N = 100;
%! param.mus = 0.5;
%! param.muc = 0.5;
%! param.vs = 1 / (SI_unit_meter / SI_unit_second);
%! param.kv = 0 / (SI_unit_pascal / (SI_unit_meter / SI_unit_second));
%! param.i = 1;
%! param.sigma0x = 1e4 / SI_unit_meter^-1;
%! param.sigma0y = 1e4 / SI_unit_meter^-1;
%! param.sigma1x = 0 / (SI_unit_second * SI_unit_meter^-1);
%! param.sigma1y = 0 / (SI_unit_second * SI_unit_meter^-1);
%! param.m = (param.d1^2 - param.d2^2) * pi / 4 * param.l1 * param.rho;
%! param.gamma = 1e-16;
%! param.sn = param.m * param.g / (param.gamma * param.d1);
%! param.z0 -= param.m * param.g * sin(abs(param.Zeta)) / (8 * param.sn);
%! options.verbose = false;
%! options.do_plot = false;
%! options.var_step_size = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0,  0, z0 + 0.5 * d2};\n");
%!     fputs(fd, "Point(2) = {0, l1, z0 + 0.5 * d2};\n");
%!     fputs(fd, "Point(3) = {0, l1, z0 + 0.5 * d1};\n");
%!     fputs(fd, "Point(4) = {0,  0, z0 + 0.5 * d1};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Curve Loop(1) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(l1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round((d1 - d2) / 2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(l1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round((d1 - d2) / 2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Surface(1) = {1,2,3,4};\n");
%!     fputs(fd, "v1[] = Extrude {{0, 1, 0}, {0, 0, z0}, Pi}{ Surface{1}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{0, 1, 0}, {0, 0, z0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Recombine Surface{v1[0],v2[0]};\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {4,9};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"iso20", "quad8"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_id_ground = rows(mesh.nodes) + 1;
%!   mesh.nodes(node_id_ground, 1:3) = zeros(1, 3);
%!   mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!   mesh.materials.iso20([mesh.groups.iso20(find([[mesh.groups.iso20.id] == 1])).elements]) = 1;
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean viscoelastic";
%!   mesh.material_data(1).beta = param.beta;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_id_ground) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_sphere";
%!   opt_mbd_mesh.joints.number = 1;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   node_id_cont = mesh.groups.quad8(1).nodes;
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fputs(fd, "set: integer ref_id_sphere = 1000;\n");
%!     fputs(fd, "set: integer drive_id_E = 2000;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fx = 2001;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fy = 2002;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fz = 2003;\n");
%!     fprintf(fd, " set: integer node_id_ground = %d;\n", node_id_ground);
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     if (options.var_step_size)
%!       fputs(fd, "    min time step: t1 / N / 100;\n");
%!       fputs(fd, "    max time step: t1 / N;\n");
%!       fputs(fd, " strategy: factor, 0.5, 1, 1.25, 3, 3, 5;\n");
%!     endif
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-6, test, norm, 1e-6, test, norm;\n");
%!     fputs(fd, "    linear solver: pardiso, grad, max iterations, 100;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: mcp newton min fb;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    #output meter: closest next, 0., forever, t1 / 100;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "    loadable elements: %d;\n", numel(node_id_cont));
%!     fputs(fd, "      gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, "reference: ref_id_sphere,\n");
%!     fputs(fd, "  position, reference, global, null,\n");
%!     fputs(fd, "  orientation, reference, global, eye,\n");
%!     fputs(fd, "  velocity, reference, global, 0., 0., vz0,\n");
%!     fputs(fd, "  angular velocity, reference, global, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     for i=1:numel(node_id_cont)
%!       fprintf(fd, "user defined: %d, unilateral disp in plane,\n", i);
%!       fprintf(fd, "node1, %d,\n", node_id_cont(i));
%!       fputs(fd, "    offset, 1,\n");
%!       fputs(fd, "    reference, node, null,\n");
%!       fputs(fd, "    stiffness, sn,\n");
%!       fputs(fd, "    enable mcp, yes,\n");
%!       fputs(fd, " node2, node_id_ground,\n");
%!       fputs(fd, " offset, reference, node, null,\n");
%!       fputs(fd, " hinge, 3, 0, 0, 1,\n");
%!       fputs(fd, "        1, 1, 0, 0,\n");
%!       fputs(fd, " coulomb friction coefficient, muc,\n");
%!       fputs(fd, " static friction coefficient, mus,\n");
%!       fputs(fd, " sliding velocity coefficient, vs,\n");
%!       fputs(fd, " sliding velocity exponent, i,\n");
%!       fputs(fd, " micro slip stiffness x, sigma0x,\n");
%!       fputs(fd, " micro slip stiffness y, sigma0y,\n");
%!       fputs(fd, " micro slip damping x, sigma1x,\n");
%!       fputs(fd, " micro slip damping y, sigma1y,\n");
%!       fputs(fd, " viscous friction coefficient, kv;\n");
%!     endfor
%!     fputs(fd, "  joint: 1, clamp, node_id_ground, node, node;\n");
%!     fputs(fd, "  gravity: uniform, component, g * cos(Zeta), 0., g * sin(Zeta);\n");
%!     fprintf(fd, "drive caller: drive_id_E, array, %d", rows(mesh.elements.iso20));
%!     for i=1:rows(mesh.elements.iso20)
%!       fprintf(fd, ",\n    element, %d, solid, string, \"E\", direct", i);
%!     endfor
%!     fputs(fd, ",\n  output, yes;\n\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fx, element, 1, joint, string, \"Fx\", direct, output, yes;\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fy, element, 1, joint, string, \"Fy\", direct, output, yes;\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fz, element, 1, joint, string, \"Fz\", direct, output, yes;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [drive_id, drive_data] = mbdyn_post_load_output_drv(opt_mbdyn.output_file);
%!   Wkin_res = drive_data{1};
%!   Fx_res = drive_data{2};
%!   Fy_res = drive_data{3};
%!   Fz_res = drive_data{4};
%!   ra = 0.5 * param.d1;
%!   ri = 0.5 * param.d2;
%!   h = param.l1;
%!   m = param.rho * pi * (ra^2 - ri^2) * h;
%!   J = m * (ra^2 + ri^2) / 2;
%!   mred = m + J / ra^2;
%!   Phi = param.Zeta + pi / 2;
%!   qddot = m * param.g * sin(Phi) / mred;
%!   qdot = qddot * sol.t;
%!   q = 0.5 * qddot * sol.t.^2;
%!   Wkin_ref = 0.5 * mred * qdot.^2;
%!   Fz_ref = repmat(-m * param.g * cos(Phi), size(sol.t));
%!   Fx_ref = repmat(m * param.g * sin(Phi) * (1 - m / mred), size(sol.t));
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Wkin_ref * SI_unit_joule, "-x;Wkin_r_e_f;0");
%!     plot(sol.t * SI_unit_second, Wkin_res * SI_unit_joule, "-o;Wkin_r_e_s;1");
%!     xlabel("t [s]");
%!     ylabel("W [J]");
%!     grid on;
%!     grid minor on;
%!     title("kinetic energy versus time");
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Fz_ref * SI_unit_newton, "-x;Fz_r_e_f;0");
%!     plot(sol.t * SI_unit_second, Fz_res * SI_unit_newton, "-o;Fz_r_e_s;1");
%!     xlabel("t [s]");
%!     ylabel("F [N]");
%!     grid on;
%!     grid minor on;
%!     title("normal force");
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Fx_ref * SI_unit_newton, "-x;Fx_r_e_f;0");
%!     plot(sol.t * SI_unit_second, Fx_res * SI_unit_newton, "-o;Fx_r_e_s;1");
%!     xlabel("t [s]");
%!     ylabel("F [N]");
%!     grid on;
%!     grid minor on;
%!     title("tangent force");
%!   endif
%!   tol = 2e-2;
%!   assert_simple(Wkin_res, Wkin_ref, tol * max(abs(Wkin_res)));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 8
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.114)
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! M_alpha = 0.25;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.kappa = param.E / (3 * (1 - 2 * param.nu));
%! param.C1 = param.G / (2 * (1 + M_alpha));
%! param.C2 = M_alpha * param.C1;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 2.5e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"quad8", "iso20", "penta15", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"quad8", "tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"quad8", "tria6h"});
%!   if (isfield(mesh.elements, "iso20"))
%!     mesh.materials.iso20 = zeros(rows(mesh.elements.iso20), 1, "int32");
%!     mesh.materials.iso20([mesh.groups.iso20(find([[mesh.groups.iso20.id] == 1])).elements]) = 1;
%!   endif
%!   if (isfield(mesh.elements, "penta15"))
%!     mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!     mesh.materials.penta15([mesh.groups.penta15(find([[mesh.groups.penta15.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).C1 = param.C1;
%!   mesh.material_data(1).C2 = param.C2;
%!   mesh.material_data(1).kappa = param.kappa;
%!   mesh.material_data(1).type = "mooney rivlin elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * (param.C1 + param.C2) * D; ## (4.114)
%!   Nref = -1/2 * pi * Ra^4 * (param.C1 + 2 * param.C2) * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-3;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 9
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.114)
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! M_alpha = 0.25;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.kappa = param.E / (3 * (1 - 2 * param.nu));
%! param.C1 = param.G / (2 * (1 + M_alpha));
%! param.C2 = M_alpha * param.C1;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 0.8333e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"iso4", "iso8", "penta6", "iso4"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"iso4", "iso4"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"iso4", "iso4"});
%!   if (isfield(mesh.elements, "iso8"))
%!     mesh.materials.iso8 = zeros(rows(mesh.elements.iso8), 1, "int32");
%!     mesh.materials.iso8([mesh.groups.iso8(find([[mesh.groups.iso8.id] == 1])).elements]) = 1;
%!   endif
%!   if (isfield(mesh.elements, "penta6"))
%!     mesh.materials.penta6 = zeros(rows(mesh.elements.penta6), 1, "int32");
%!     mesh.materials.penta6([mesh.groups.penta6(find([[mesh.groups.penta6.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).C1 = param.C1;
%!   mesh.material_data(1).C2 = param.C2;
%!   mesh.material_data(1).kappa = param.kappa;
%!   mesh.material_data(1).type = "mooney rivlin elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * (param.C1 + param.C2) * D; ## (4.114)
%!   Nref = -1/2 * pi * Ra^4 * (param.C1 + 2 * param.C2) * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-3;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 10
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.114)
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! M_alpha = 0.25;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.kappa = param.E / (3 * (1 - 2 * param.nu));
%! param.C1 = param.G / (2 * (1 + M_alpha));
%! param.C2 = M_alpha * param.C1;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 2.5e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"quad8", "iso20r", "penta15", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"quad8", "tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"quad8", "tria6h"});
%!   if (isfield(mesh.elements, "iso20r"))
%!     mesh.materials.iso20r = zeros(rows(mesh.elements.iso20r), 1, "int32");
%!     mesh.materials.iso20r([mesh.groups.iso20r(find([[mesh.groups.iso20r.id] == 1])).elements]) = 1;
%!   endif
%!   if (isfield(mesh.elements, "penta15"))
%!     mesh.materials.penta15 = zeros(rows(mesh.elements.penta15), 1, "int32");
%!     mesh.materials.penta15([mesh.groups.penta15(find([[mesh.groups.penta15.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).C1 = param.C1;
%!   mesh.material_data(1).C2 = param.C2;
%!   mesh.material_data(1).kappa = param.kappa;
%!   mesh.material_data(1).type = "mooney rivlin elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * (param.C1 + param.C2) * D; ## (4.114)
%!   Nref = -1/2 * pi * Ra^4 * (param.C1 + 2 * param.C2) * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-3;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 11
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.114)
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.499;
%! M_alpha = 0.25;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.kappa = param.E / (3 * (1 - 2 * param.nu));
%! param.C1 = param.G / (2 * (1 + M_alpha));
%! param.C2 = M_alpha * param.C1;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 1.25e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fprintf(fd, "h = %e;\n", param.h);
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "Point(2) = {l1, 0, 0};\n");
%!     fputs(fd, "Point(3) = {l1, 0.5 * d1, 0};\n");
%!     fputs(fd, "Point(4) = {0, 0.5 * d1, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Curve Loop(1) = {1, 2, 3, 4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "v[] = Extrude {{1, 0, 0}, {0, 0, 0}, 2*Pi}{ Surface{1}; };\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{v[1]}; } } = h;\n");
%!     fputs(fd, "Physical Volume(1) = {v[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {4};\n");
%!     fputs(fd, "Physical Surface(2) = {2};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!     fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"tria6h", "tet10h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"tria6h"});
%!   if (isfield(mesh.elements, "tet10h"))
%!     mesh.materials.tet10h = zeros(rows(mesh.elements.tet10h), 1, "int32");
%!     mesh.materials.tet10h([mesh.groups.tet10h(find([[mesh.groups.tet10h.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).C1 = param.C1;
%!   mesh.material_data(1).C2 = param.C2;
%!   mesh.material_data(1).kappa = param.kappa;
%!   mesh.material_data(1).type = "mooney rivlin elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 50;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-3, test, norm, 1e-3, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 30;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"0.5 * pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * (param.C1 + param.C2) * D; ## (4.114)
%!   Nref = -1/2 * pi * Ra^4 * (param.C1 + 2 * param.C2) * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-3;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 12
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.113)
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.5;
%! param.delta = 0.25;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.C10 = param.G / (2 * (1 + param.delta));
%! param.C01 = param.C10 * param.delta;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 2.5e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"quad8", "iso20upcr", "penta15upc", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"quad8", "tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"quad8", "tria6h"});
%!   if (isfield(mesh.elements, "iso20upcr"))
%!     mesh.materials.iso20upcr = zeros(rows(mesh.elements.iso20upcr), 1, "int32");
%!     mesh.materials.iso20upcr([mesh.groups.iso20upcr(find([[mesh.groups.iso20upcr.id] == 1])).elements]) = 1;
%!   endif
%!   if (isfield(mesh.elements, "penta15upc"))
%!     mesh.materials.penta15upc = zeros(rows(mesh.elements.penta15upc), 1, "int32");
%!     mesh.materials.penta15upc([mesh.groups.penta15upc(find([[mesh.groups.penta15upc.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).delta = param.delta;
%!   mesh.material_data(1).type = "mooney rivlin elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   if (mesh.material_data(1).nu == 0.5)
%!     ++opt_mbd_mesh.genels.number;
%!   endif
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-6, test, norm, 1e-6, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 10;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    hydraulic nodes: %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "    genels: %d;\n", opt_mbd_mesh.genels.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     if (mesh.material_data(1).nu == 0.5)
%!       ## Because the pressure will have no impact on strain if the material is fully incompressible
%!       fprintf(fd, "genel: %d, clamp, %d, hydraulic, null;\n", opt_mbd_mesh.genels.number, mesh.elements.iso20upcr(1,1));
%!     endif
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * (param.C10 + param.C01) * D; ## (4.113)
%!   Nref = -1/2 * pi * Ra^4 * (param.C10 + 2 * param.C01) * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-4;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 13
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.113)
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.49999; ## FIXME: nu=0.5 is possible but the condition number during the first iteration may become very bad!
%! param.delta = 0.25;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.C10 = param.G / (2 * (1 + param.delta));
%! param.C01 = param.C10 * param.delta;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 1.25e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"iso4", "iso8upc", "penta6upc", "tria3"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"iso4", "tria3"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"iso4", "tria3"});
%!   if (isfield(mesh.elements, "iso8upc"))
%!     mesh.materials.iso8upc = zeros(rows(mesh.elements.iso8upc), 1, "int32");
%!     mesh.materials.iso8upc([mesh.groups.iso8upc(find([[mesh.groups.iso8upc.id] == 1])).elements]) = 1;
%!   endif
%!   if (isfield(mesh.elements, "penta6upc"))
%!     mesh.materials.penta6upc = zeros(rows(mesh.elements.penta6upc), 1, "int32");
%!     mesh.materials.penta6upc([mesh.groups.penta6upc(find([[mesh.groups.penta6upc.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).delta = param.delta;
%!   mesh.material_data(1).type = "mooney rivlin elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   if (mesh.material_data(1).nu == 0.5)
%!     ++opt_mbd_mesh.genels.number;
%!   endif
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-6, test, norm, 1e-6, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4, 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 10;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    hydraulic nodes: %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "    genels: %d;\n", opt_mbd_mesh.genels.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     if (mesh.material_data(1).nu == 0.5)
%!       ## Because the pressure will have no impact on strain if the material is fully incompressible
%!       fprintf(fd, "genel: %d, clamp, %d, hydraulic, null;\n", opt_mbd_mesh.genels.number, mesh.elements.iso8upc(1,1));
%!     endif
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * (param.C10 + param.C01) * D; ## (4.113)
%!   Nref = -1/2 * pi * Ra^4 * (param.C10 + 2 * param.C01) * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 1e-2;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 14
%! ## torsion-tension coupling effect of an incompressible cylinder
%! ## Viskoelastisches Materialverhalten von Elastomerwerkstoffen
%! ## Experimentelle Untersuchung und Modellbildung
%! ## Konstantin Sedlan 2000
%! ## page 102
%! ## equation (4.113)
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 5e6 / SI_unit_pascal;
%! param.nu = 0.5;
%! param.delta = 0.25;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.C10 = param.G / (2 * (1 + param.delta));
%! param.C01 = param.C10 * param.delta;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.l1 = 20e-3 / SI_unit_meter;
%! param.h = 2.5e-3 / SI_unit_meter;
%! options.verbose = false;
%! options.do_plot = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0, 0, 0};\n");
%!     fputs(fd, "L[] = Extrude {0, 0.5 * d1, 0}{ Point{1}; Layers{Ceil(0.5 * d1 / h)}; };\n");
%!     fputs(fd, "s1[] = Extrude {l1, 0, 0}{ Line{L[1]}; Layers{Ceil(l1 / h)}; Recombine; };\n");
%!     fputs(fd, "v1[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{s1[1]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{1, 0, 0}, {0, 0, 0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {3,7};\n");
%!     fputs(fd, "Physical Surface(2) = {4,8};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 1;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"quad8", "iso20upc", "penta15upc", "tria6h"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_idx_bearing1 = rows(mesh.nodes) + 1;
%!   node_idx_bearing2 = rows(mesh.nodes) + 2;
%!   X1 = [0; 0; 0];
%!   X2 = [param.l1; 0; 0];
%!   mesh.nodes(node_idx_bearing2, 1:3) = X2.';
%!   mesh.nodes(node_idx_bearing1, 1:3) = X1.';
%!   mesh.elements.rbe2(2) = fem_pre_mesh_rbe3_from_surf(mesh, 2, node_idx_bearing2, {"quad8", "tria6h"});
%!   mesh.elements.rbe2(1) = fem_pre_mesh_rbe3_from_surf(mesh, 1, node_idx_bearing1, {"quad8", "tria6h"});
%!   if (isfield(mesh.elements, "iso20upc"))
%!     mesh.materials.iso20upc = zeros(rows(mesh.elements.iso20upc), 1, "int32");
%!     mesh.materials.iso20upc([mesh.groups.iso20upc(find([[mesh.groups.iso20upc.id] == 1])).elements]) = 1;
%!   endif
%!   if (isfield(mesh.elements, "penta15upc"))
%!     mesh.materials.penta15upc = zeros(rows(mesh.elements.penta15upc), 1, "int32");
%!     mesh.materials.penta15upc([mesh.groups.penta15upc(find([[mesh.groups.penta15upc.id] == 1])).elements]) = 1;
%!   endif
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).delta = param.delta;
%!   mesh.material_data(1).type = "mooney rivlin elastic";
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.joints.number = 2;
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing1) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.type(node_idx_bearing2) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   if (mesh.material_data(1).nu == 0.5)
%!     ++opt_mbd_mesh.genels.number;
%!   endif
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fprintf(fd, " set: real t1 = %g;\n", 1 / SI_unit_second);
%!     fputs(fd, " set: real N = 40;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     fputs(fd, "    min time step: t1 / N / 100;\n");
%!     fputs(fd, "    max time step: t1 / N;\n");
%!     fputs(fd, " strategy: factor, 0.8, 1, 1.25, 3, 5, 10;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-6, test, norm, 1e-6, test, norm;\n");
%!     fputs(fd, "    linear solver: umfpack, grad, scale, iterative, always, max iterations, 0;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: nox, modified, 100, keep jacobian matrix, jacobian operator, newton krylov, forcing term, type 2, forcing term min tolerance, 1e-10, forcing term max tolerance, 1e-6, inner iterations before assembly, 10;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    model: static;\n");
%!     fputs(fd, "    output meter: closest next, 0., forever, t1 / 20;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    hydraulic nodes: %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "    genels: %d;\n", opt_mbd_mesh.genels.number);
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "    joint: 1, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing1);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing1, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fputs(fd, "    joint: 2, total pin joint,\n");
%!     fprintf(fd, "                    %d,\n", node_idx_bearing2);
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fprintf(fd, "                            position,          reference, global, %s\n", sprintf("%.16e, ", mesh.nodes(node_idx_bearing2, 1:3)));
%!     fputs(fd, "                            position orientation,   reference, global, eye,\n");
%!     fputs(fd, "                            rotation orientation,   reference, global, eye,\n");
%!     fputs(fd, "                    position constraint, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active, \n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                    orientation constraint,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                                    active,\n");
%!     fputs(fd, "                            component, string, \"pi * Time / t1\",\n");
%!     fputs(fd, "                                       const, 0.,\n");
%!     fputs(fd, "                                       const, 0.;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     if (mesh.material_data(1).nu == 0.5)
%!       ## Because the pressure will have no impact on strain if the material is fully incompressible
%!       fprintf(fd, "genel: %d, clamp, %d, hydraulic, null;\n", opt_mbd_mesh.genels.number, mesh.elements.iso20upc(1,1));
%!     endif
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [joint_id, local_reaction, global_reaction] = mbdyn_post_load_output_jnt(opt_mbdyn.output_file);
%!   Phi1 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing1), 4, :)(:), 2 * pi);
%!   Phi2 = mod(2 * pi + sol.def(find(mesh_sol.node_id == node_idx_bearing2), 4, :)(:), 2 * pi);
%!   x1 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing1), 1, :)(:);
%!   x2 = mesh_sol.nodes(find(mesh_sol.node_id == node_idx_bearing2), 1, :)(:);
%!   D = (Phi2 - Phi1) ./ (x2 - x1);
%!   Ra = param.d1 / 2;
%!   Mref = pi * Ra^4 * (param.C10 + param.C01) * D; ## (4.113)
%!   Nref = -1/2 * pi * Ra^4 * (param.C10 + 2 * param.C01) * D.^2;
%!   N = -global_reaction{2}(:, 1); ## negative sign because of MBDyn's conventions
%!   M = -global_reaction{2}(:, 4);
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, N * SI_unit_newton, "-;N(D);r");
%!     plot(D * SI_unit_meter^-1, Nref * SI_unit_newton, "-;Nref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("axial force [N]");
%!     title("axial force versus twist");
%!     grid minor on;
%!     figure("visible", "off");
%!     hold on;
%!     plot(D * SI_unit_meter^-1, M * SI_unit_newton * SI_unit_meter, "-;M(D);r");
%!     plot(D * SI_unit_meter^-1, Mref * SI_unit_newton * SI_unit_meter, "-;Mref(D);b");
%!     xlabel("twist [rad/m]");
%!     ylabel("torque [Nm]");
%!     title("torque versus twist");
%!     grid minor on;
%!   endif
%!   tol = 5e-4;
%!   assert_simple(N, Nref, tol * norm(Nref));
%!   assert_simple(M, Mref, tol * norm(Mref));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! ## TEST 15
%! close all;
%! pkg load mboct-fem-pkg;
%! ## Define the unit system
%! SI_unit_meter = 1e-3;
%! SI_unit_second = 1;
%! SI_unit_kilogram = 1e3;
%! SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%! SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%! SI_unit_joule = SI_unit_newton * SI_unit_meter;
%! param.alpha = 0 / (1 / SI_unit_second);
%! param.beta = 0 / (SI_unit_second);
%! param.E = 10e6 / SI_unit_pascal;
%! param.nu = 0.3;
%! param.G = param.E / (2 * (1 + param.nu));
%! param.beta = 1e-3 / SI_unit_second;
%! param.rho = 1000 / (SI_unit_kilogram / SI_unit_meter^3);
%! param.h = 0.625e-3 / SI_unit_meter;
%! param.d1 = 10e-3 / SI_unit_meter;
%! param.d2 = param.d1 - 2 * param.h;
%! param.l1 = param.h;
%! param.z0 = 5e-3 / SI_unit_meter;
%! param.vz0 = 0 / (SI_unit_meter / SI_unit_second);
%! param.g = 9.81 / (SI_unit_meter / SI_unit_second^2);
%! param.Zeta = -80 * pi / 180;
%! param.t1 = 0.1 / SI_unit_second;
%! param.N = 100;
%! param.mus = 0.5;
%! param.muc = 0.5;
%! param.vs = 1 / (SI_unit_meter / SI_unit_second);
%! param.kv = 0 / (SI_unit_pascal / (SI_unit_meter / SI_unit_second));
%! param.i = 1;
%! param.sigma0x = 1e4 / SI_unit_meter^-1;
%! param.sigma0y = 1e4 / SI_unit_meter^-1;
%! param.sigma1x = 0 / (SI_unit_second * SI_unit_meter^-1);
%! param.sigma1y = 0 / (SI_unit_second * SI_unit_meter^-1);
%! param.m = (param.d1^2 - param.d2^2) * pi / 4 * param.l1 * param.rho;
%! param.gamma = 1e-16;
%! param.sn = param.m * param.g / (param.gamma * param.d1);
%! param.z0 -= param.m * param.g * sin(abs(param.Zeta)) / (8 * param.sn);
%! options.verbose = false;
%! options.do_plot = false;
%! options.var_step_size = false;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   fd = -1;
%!   mbdyn_file = [filename, ".mbdyn"];
%!   elem_file = [filename, ".elm"];
%!   nodes_file = [filename, ".nod"];
%!   csl_file = [filename, ".csl"];
%!   geometry_file = [filename, ".geo"];
%!   unwind_protect
%!     [fd, msg] = fopen(geometry_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", geometry_file);
%!     endif
%!     fn = fieldnames(param);
%!     for i=1:length(fn)
%!       fprintf(fd, "%s = %.3e;\n", fn{i}, getfield(param, fn{i}));
%!     endfor
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Point(1) = {0,  0, z0 + 0.5 * d2};\n");
%!     fputs(fd, "Point(2) = {0, l1, z0 + 0.5 * d2};\n");
%!     fputs(fd, "Point(3) = {0, l1, z0 + 0.5 * d1};\n");
%!     fputs(fd, "Point(4) = {0,  0, z0 + 0.5 * d1};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 1};\n");
%!     fputs(fd, "Curve Loop(1) = {1,2,3,4};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "Transfinite Curve(1) = Max(1, Round(l1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(2) = Max(1, Round((d1 - d2) / 2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(3) = Max(1, Round(l1 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Curve(4) = Max(1, Round((d1 - d2) / 2 / h)) + 1;\n");
%!     fputs(fd, "Transfinite Surface(1) = {1,2,3,4};\n");
%!     fputs(fd, "v1[] = Extrude {{0, 1, 0}, {0, 0, z0}, Pi}{ Surface{1}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "v2[] = Extrude {{0, 1, 0}, {0, 0, z0}, Pi}{ Surface{v1[0]}; Layers{Ceil(Pi * d1 / 2 / h)}; Recombine; };\n");
%!     fputs(fd, "Recombine Surface{1,v1[0]};\n");
%!     fputs(fd, "Recombine Surface{v1[0],v2[0]};\n");
%!     fputs(fd, "Physical Volume(1) = {v1[1],v2[1]};\n");
%!     fputs(fd, "Physical Surface(1) = {4,9};\n");
%!     fputs(fd, "Mesh.ElementOrder = 2;\n");
%!     fputs(fd, "Mesh.SecondOrderIncomplete = 0;\n");
%!     fputs(fd, "Mesh 3;\n");
%!     fputs(fd, "Coherence Mesh;\n");
%!     fputs(fd, "Mesh.Format = 1;\n");
%!     fprintf(fd, "Save \"%s.msh\";\n", filename);
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   fprintf(stderr, "meshing ...\n");
%!   pid = spawn("gmsh", {"-0", "-format", "msh2", geometry_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     warning("gmsh failed with status %d", status);
%!   endif
%!   fprintf(stderr, "loading mesh ...\n");
%!   opt_msh.elem_type = {"iso27", "quad9"};
%!   mesh = fem_pre_mesh_reorder(fem_pre_mesh_import([filename, ".msh"], "gmsh", opt_msh));
%!   node_id_ground = rows(mesh.nodes) + 1;
%!   mesh.nodes(node_id_ground, 1:3) = zeros(1, 3);
%!   mesh.materials.iso27 = zeros(rows(mesh.elements.iso27), 1, "int32");
%!   mesh.materials.iso27([mesh.groups.iso27(find([[mesh.groups.iso27.id] == 1])).elements]) = 1;
%!   mesh.material_data(1).rho = param.rho;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).type = "neo hookean viscoelastic";
%!   mesh.material_data(1).beta = param.beta;
%!   opt_mbd_mesh = struct();
%!   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!   opt_mbd_mesh.struct_nodes.type(node_id_ground) = MBDYN_NODE_TYPE_STATIC_STRUCT;
%!   opt_mbd_mesh.struct_nodes.reference_frame = "ref_id_sphere";
%!   opt_mbd_mesh.joints.number = 1;
%!   opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!   opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!   load_case_dof.locked_dof = false(size(mesh.nodes));
%!   load_case_empty = struct();
%!   opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case_empty, elem_file, opt_mbd_mesh);
%!   node_id_cont = mesh.groups.quad9(1).nodes;
%!   unwind_protect
%!     [fd, msg] = fopen(mbdyn_file, "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", mbdyn_file);
%!     endif
%!     fputs(fd, "set: integer ref_id_sphere = 1000;\n");
%!     fputs(fd, "set: integer drive_id_E = 2000;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fx = 2001;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fy = 2002;\n");
%!     fputs(fd, "set: integer drive_id_ground_Fz = 2003;\n");
%!     fprintf(fd, " set: integer node_id_ground = %d;\n", node_id_ground);
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "    problem: initial value; # the default\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "    initial time: 0;\n");
%!     fputs(fd, "    final time: t1;\n");
%!     fputs(fd, "    time step: t1 / N;\n");
%!     if (options.var_step_size)
%!       fputs(fd, "    min time step: t1 / N / 100;\n");
%!       fputs(fd, "    max time step: t1 / N;\n");
%!       fputs(fd, " strategy: factor, 0.5, 1, 1.25, 3, 3, 5;\n");
%!     endif
%!     fputs(fd, "    max iterations: 10000;\n");
%!     fputs(fd, "    tolerance: 1.e-6, test, norm, 1e-6, test, norm;\n");
%!     fputs(fd, "    linear solver: pardiso, grad, max iterations, 100;\n");
%!     fputs(fd, "    method: implicit euler;\n");
%!     fputs(fd, "         derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: 1e-9, auto;\n");
%!     fputs(fd, "         output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "    nonlinear solver: mcp newton min fb;\n");
%!     fprintf(fd, "    threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "    threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "    #output meter: closest next, 0., forever, t1 / 100;\n");
%!     fputs(fd, "        use automatic differentiation;\n");
%!     fprintf(fd, "    structural nodes: %d;\n", opt_mbd_mesh.struct_nodes.number);
%!     fprintf(fd, "    joints: %d;\n", opt_mbd_mesh.joints.number);
%!     fprintf(fd, "    solids: %d;\n", opt_mbd_mesh.solids.number);
%!     fprintf(fd, "    loadable elements: %d;\n", numel(node_id_cont));
%!     fputs(fd, "      gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, "reference: ref_id_sphere,\n");
%!     fputs(fd, "  position, reference, global, null,\n");
%!     fputs(fd, "  orientation, reference, global, eye,\n");
%!     fputs(fd, "  velocity, reference, global, 0., 0., vz0,\n");
%!     fputs(fd, "  angular velocity, reference, global, null;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", nodes_file);
%!     fputs(fd, " end: nodes;\n");
%!     fprintf(fd, "include: \"%s\";\n", csl_file);
%!     fputs(fd, " begin: elements;\n");
%!     fprintf(fd, "include: \"%s\";\n", elem_file);
%!     for i=1:numel(node_id_cont)
%!       fprintf(fd, "user defined: %d, unilateral disp in plane,\n", i);
%!       fprintf(fd, "node1, %d,\n", node_id_cont(i));
%!       fputs(fd, "    offset, 1,\n");
%!       fputs(fd, "    reference, node, null,\n");
%!       fputs(fd, "    stiffness, sn,\n");
%!       fputs(fd, "    enable mcp, yes,\n");
%!       fputs(fd, " node2, node_id_ground,\n");
%!       fputs(fd, " offset, reference, node, null,\n");
%!       fputs(fd, " hinge, 3, 0, 0, 1,\n");
%!       fputs(fd, "        1, 1, 0, 0,\n");
%!       fputs(fd, " coulomb friction coefficient, muc,\n");
%!       fputs(fd, " static friction coefficient, mus,\n");
%!       fputs(fd, " sliding velocity coefficient, vs,\n");
%!       fputs(fd, " sliding velocity exponent, i,\n");
%!       fputs(fd, " micro slip stiffness x, sigma0x,\n");
%!       fputs(fd, " micro slip stiffness y, sigma0y,\n");
%!       fputs(fd, " micro slip damping x, sigma1x,\n");
%!       fputs(fd, " micro slip damping y, sigma1y,\n");
%!       fputs(fd, " viscous friction coefficient, kv;\n");
%!     endfor
%!     fputs(fd, "  joint: 1, clamp, node_id_ground, node, node;\n");
%!     fputs(fd, "  gravity: uniform, component, g * cos(Zeta), 0., g * sin(Zeta);\n");
%!     fprintf(fd, "drive caller: drive_id_E, array, %d", rows(mesh.elements.iso27));
%!     for i=1:rows(mesh.elements.iso27)
%!       fprintf(fd, ",\n    element, %d, solid, string, \"E\", direct", i);
%!     endfor
%!     fputs(fd, ",\n  output, yes;\n\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fx, element, 1, joint, string, \"Fx\", direct, output, yes;\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fy, element, 1, joint, string, \"Fy\", direct, output, yes;\n");
%!     fputs(fd, "drive caller: drive_id_ground_Fz, element, 1, joint, string, \"Fz\", direct, output, yes;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!     fd = -1;
%!   end_unwind_protect
%!   if (~options.verbose)
%!     opt_mbdyn.logfile = [filename, ".stdout"];
%!   endif
%!   opt_mbdyn.output_file = [filename, "_mbdyn"];
%!   info_mbdyn = mbdyn_solver_run(mbdyn_file, opt_mbdyn);
%!   [mesh_sol, sol] = mbdyn_post_load_output_sol(opt_mbdyn.output_file);
%!   [drive_id, drive_data] = mbdyn_post_load_output_drv(opt_mbdyn.output_file);
%!   Wkin_res = drive_data{1};
%!   Fx_res = drive_data{2};
%!   Fy_res = drive_data{3};
%!   Fz_res = drive_data{4};
%!   ra = 0.5 * param.d1;
%!   ri = 0.5 * param.d2;
%!   h = param.l1;
%!   m = param.rho * pi * (ra^2 - ri^2) * h;
%!   J = m * (ra^2 + ri^2) / 2;
%!   mred = m + J / ra^2;
%!   Phi = param.Zeta + pi / 2;
%!   qddot = m * param.g * sin(Phi) / mred;
%!   qdot = qddot * sol.t;
%!   q = 0.5 * qddot * sol.t.^2;
%!   Wkin_ref = 0.5 * mred * qdot.^2;
%!   Fz_ref = repmat(-m * param.g * cos(Phi), size(sol.t));
%!   Fx_ref = repmat(m * param.g * sin(Phi) * (1 - m / mred), size(sol.t));
%!   if (options.do_plot)
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Wkin_ref * SI_unit_joule, "-x;Wkin_r_e_f;0");
%!     plot(sol.t * SI_unit_second, Wkin_res * SI_unit_joule, "-o;Wkin_r_e_s;1");
%!     xlabel("t [s]");
%!     ylabel("W [J]");
%!     grid on;
%!     grid minor on;
%!     title("kinetic energy versus time");
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Fz_ref * SI_unit_newton, "-x;Fz_r_e_f;0");
%!     plot(sol.t * SI_unit_second, Fz_res * SI_unit_newton, "-o;Fz_r_e_s;1");
%!     xlabel("t [s]");
%!     ylabel("F [N]");
%!     grid on;
%!     grid minor on;
%!     title("normal force");
%!     figure("visible", "off");
%!     hold on;
%!     plot(sol.t * SI_unit_second, Fx_ref * SI_unit_newton, "-x;Fx_r_e_f;0");
%!     plot(sol.t * SI_unit_second, Fx_res * SI_unit_newton, "-o;Fx_r_e_s;1");
%!     xlabel("t [s]");
%!     ylabel("F [N]");
%!     grid on;
%!     grid minor on;
%!     title("tangent force");
%!   endif
%!   tol = 2e-2;
%!   assert_simple(Wkin_res, Wkin_ref, tol * max(abs(Wkin_res)));
%! unwind_protect_cleanup
%!   if (numel(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
