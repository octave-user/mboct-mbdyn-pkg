## Copyright (C) 2014(-2023) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{log_dat}] = mbdyn_post_load_log(@var{mbdyn_filename})
##
## Parses the MBDyn log-file "<@var{mbdyn_filename}>.log" and returns information about variables, nodes, elements, and structural node degrees of freedom.
##
## @var{mbdyn_filename} @dots{} MBDyn input file without extension .mbdyn.
##
## @var{log_dat}.nodes @dots{} Struct array with information about all nodes including dummy nodes in ascending order.
##
## @var{log_dat}.nodes(@var{i}).X0 @dots{} Initial position of the structural node @var{i}.
##
## @var{log_dat}.nodes(@var{i}).Phi0 @dots{} Angles of rotation of the structural node @var{i}.
##
## @var{log_dat}.nodes(@var{i}).R0 @dots{} Initial rotation matrix of the structural node @var{i}.
##
## @var{log_dat}.nodes(@var{i}).orientation_description @dots{} Orientation_description of the node: one of  ("euler123", "euler321", "euler313","phi","mat").
##
## @var{log_dat}.dof_info.struct_node_dofs @dots{} Index of the first degree of freedom minus one of the node in the global vector of degrees of freedom (e.g X_i = X(struct_node_dofs(i) + 1:3); g_i = X(struct_node_dofs(i) + 4:6)).
##
## @var{log_dat}.dof_info.struct_node_force_index @dots{} Index of the first equation of the force and moment equilibrium in the global residual vector (e.g F_i = f(struct_node_force_index(i) + 1:3); M_i = f(struct_node_force_index(i) + 4:6)).
##
## @var{log_dat}.dof_info.struct_node_labels @dots{} Labels of the corresponding structural nodes.
##
## @var{X}(struct_node_dofs(i) + 1:6) @dots{} Is the vector of degrees of freedom of a static or dynamic structural node with label @var{struct_node_labels}(i).
##
##
## @end deftypefn

function [log_dat] = mbdyn_post_load_log(mbdyn_filename, options)
  if (nargin < 1 || nargout > 1)
    print_usage();
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (~isfield(options, "dof_info"))
    options.dof_info = true;
  endif

  if (~isfield(options, "nodes"))
    options.nodes = true;
  endif

  if (~isfield(options, "hydrodynamic_bearings"))
    options.hydrodynamic_bearings = true;
  endif

  if (~isfield(options, "vars"))
    options.vars = true;
  endif

  if (~isfield(options, "beams2"))
    options.beams2 = true;
  endif

  if (~isfield(options, "beams3"))
    options.beams3 = true;
  endif

  if (~isfield(options, "modal"))
    options.modal = true;
  endif

  if (~isfield(options, "joints"))
    options.joints = true;
  endif

  if (~isfield(options, "node_id_prefix"))
    options.node_id_prefix = "node_id_";
  endif

  if (~isfield(options, "node_idx_prefix"))
    options.node_idx_prefix = "node_idx_";
  endif

  if (~isfield(options, "force_id_prefix"))
    options.force_id_prefix = "force_id_";
  endif

  if (~isfield(options, "force_idx_prefix"))
    options.force_idx_prefix = "force_idx_";
  endif

  log_filename = mbdyn_post_output_filename(mbdyn_filename, ".log");

  fid = -1;

  unwind_protect
    [fid, msg] = fopen(log_filename, "rt");

    if (fid == -1)
      error("mbdyn:post", "could not open file \"%s\": %s", log_filename,msg);
    endif

    iNode = 0;
    iBeam2 = 0;
    iBeam3 = 0;
    iBearing = 0;
    iModal = 0;
    iJoint = 0;

    empty_cell = cell(0, 0);
    log_dat.nodes = struct("label", empty_cell, "X0", empty_cell, "R0", empty_cell, ...
                           "Phi0", empty_cell, "orientation_description", empty_cell);
    log_dat.vars = struct();
    log_dat.dof_info = struct();
    log_dat.beams2 = struct("label", empty_cell, "nodes", empty_cell);
    log_dat.beams3 = struct("label", empty_cell, "nodes", empty_cell);
    log_dat.bearings = struct("label", empty_cell, "type", empty_cell, "B", empty_cell, ...
                              "d", empty_cell, "eta", empty_cell);
    log_dat.joints = struct("label", empty_cell, "type", empty_cell, "nodes", empty_cell);

    iLine = 0;

    try

      while (1)
        line = fgets(fid);

        if (~ischar(line) && line == -1)
          break;
        endif

        ++iLine;

        tag_end = find(line == ':');

        if (length(tag_end) >= 1)
          tag_end = tag_end(1);
          tag = line(1:tag_end);
          data = line(length(tag)+1:end);

          switch (tag)
            case {"structural node:", "relative frame struct node:", "relative frame structural node:"}
              if (options.nodes)
                [node_label, x, y, z, orientation_description, a11, a12, a13, a21, a22, a23, a31, a32, a33, count] = sscanf(data,"%d %g %g %g %s %g %g %g %g %g %g %g %g %g","C");

                X0 = [ x; y; z ];

                switch (orientation_description)
                  case {"euler123", "euler321", "euler313"}
                    if (count ~= 8)
                      error("while reading logfile \"%s\": \"%s\"",log_filename,line);
                    endif

                    Phix = a11 * pi / 180;
                    Phiy = a12 * pi / 180;
                    Phiz = a13 * pi / 180;

                    Phi0 = [ Phix; Phiy; Phiz ];

                    switch (orientation_description)
                      case "euler123"
                        R0 = euler123_to_rotation_matrix(Phi0);
                      case "euler313"
                        R0  = euler313_to_rotation_matrix(Phi0);
                      case "euler321"
                        R0  = euler321_to_rotation_matrix(Phi0);
                    endswitch
                  case "phi"
                    if (count ~= 8)
                      error("while reading logfile \"%s\": \"%s\"",log_filename,line);
                    endif
                    Phix = a11;
                    Phiy = a12;
                    Phiz = a13;

                    Phi0 = [ Phix; Phiy; Phiz ];

                    R0  = rotation_vector_to_rotation_matrix(Phi0);
                  case "mat"
                    if (count ~= 14)
                      error("while reading logfile \"%s\": \"%s\"",log_filename,line);
                    endif

                    R0 = [ a11, a12, a13;
                           a21, a22, a23;
                           a31, a32, a33 ];

                    Phi0 = rotation_matrix_to_rotation_vector(R0);
                  case ""
                    if (count ~= 4)
                      error("while reading logfile \"%s\": \"%s\"", log_filename, line);
                    endif
                    orientation_description = "none";
                    R0 = eye(3);
                    Phi0 = zeros(3, 1);
                  otherwise
                    error("orientation %s is not implemented!",orientation_description);
                endswitch

                log_dat.nodes(++iNode).label = int32(node_label);
                log_dat.nodes(iNode).X0 = X0;
                log_dat.nodes(iNode).Phi0 = Phi0;
                log_dat.nodes(iNode).R0 = R0;
                log_dat.nodes(iNode).orientation_description = orientation_description;
              endif
            case "struct node dofs:"
              if (options.dof_info)
                log_dat.dof_info.struct_node_dofs = int32(sscanf(data, "%d ")).';
              endif

            case "struct node eqs:"
              if (options.dof_info)
                log_dat.dof_info.struct_node_force_index = int32(sscanf(data, "%d ")).';
              endif

            case "struct node labels:"
              if (options.dof_info)
                log_dat.dof_info.struct_node_labels = int32(sscanf(data, "%d ")).';
              endif

            case "beam2:"
              if (options.beams2)
                [beamLabel, data] = next_token(data, "%d");
                log_dat.beams2(++iBeam2).label = int32(beamLabel);

                for i=1:2
                  [nodeLabel, data] = next_token(data, "%d");
                  log_dat.beams2(iBeam2).nodes(i).label = int32(nodeLabel);

                  for j=1:3
                    [log_dat.beams2(iBeam2).nodes(i).offset(j), data] = next_token(data, "%g");
                  endfor
                endfor
              endif

            case "beam3:"
              if (options.beams3)
                [beamLabel, data] = next_token(data, "%d");
                log_dat.beams3(++iBeam3).label = int32(beamLabel);

                for i=1:3
                  [nodeLabel, data] = next_token(data, "%d");
                  log_dat.beams3(iBeam3).nodes(i).label = int32(nodeLabel);

                  for j=1:3
                    [log_dat.beams3(iBeam3).nodes(i).offset(j), data] = next_token(data, "%g");
                  endfor
                endfor
              endif
            case "modal:"
              if (options.modal)
                [log_dat.modal(++iModal).label, data] = next_token(data, "%d");
                [log_dat.modal(iModal).modal_node, data] = next_token(data, "%d");
              endif
            case "hydrodynamic_plain_bearing_with_offset:"
              if (options.hydrodynamic_bearings)
                data = sscanf(data, "%g");
                idx = 0;
                log_dat.bearings(++iBearing).type = "butenschoen";
                log_dat.bearings(iBearing).label = int32(data(++idx));

                for i=1:2
                  log_dat.bearings(iBearing).nodes(i).label = int32(data(++idx));
                  log_dat.bearings(iBearing).nodes(i).o = data(idx + (1:3));
                  idx += 3;
                endfor

                log_dat.bearings(iBearing).B = data(++idx);
                log_dat.bearings(iBearing).d = data(++idx);
                log_dat.bearings(iBearing).Psi = data(++idx);
                log_dat.bearings(iBearing).eta = data(++idx);
              endif
            case "hydrodynamic plain bearing:"
              if (options.hydrodynamic_bearings)
                [log_dat.bearings(++iBearing).label, data]= next_token(data, "%d");
                log_dat.bearings(iBearing).type = "simple cylindrical";
                [log_dat.bearings(iBearing).d, data] = next_token(data, "%g");
                [log_dat.bearings(iBearing).D, data] = next_token(data, "%g");
                [log_dat.bearings(iBearing).B, data] = next_token(data, "%g");
                [log_dat.bearings(iBearing).eta, data] = next_token(data, "%g");
                [log_dat.bearings(iBearing).M, data] = next_token(data, "%d");

                for i=1:log_dat.bearings(iBearing).M
                  [log_dat.bearings(iBearing).z(i), data] = next_token(data, "%g");
                endfor

                [log_dat.bearings(iBearing).N, data] = next_token(data, "%d");

                for i=1:log_dat.bearings(iBearing).N
                  [log_dat.bearings(iBearing).Phi(i), data] = next_token(data, "%g");
                endfor

                for i=1:2
                  [log_dat.bearings(iBearing).nodes(i).label, data] = next_token(data, "%d");
                  for j=1:3
                    [log_dat.bearings(iBearing).nodes(i).o(j, 1), data] = next_token(data, "%g");
                  endfor

                  for j=1:3
                    for k=1:3
                      [log_dat.bearings(iBearing).nodes(i).Rb(j, k), data] = next_token(data, "%g");
                    endfor
                  endfor
                endfor

                [log_dat.bearings(iBearing).output_force, data] = next_token(data, "%d");
                [log_dat.bearings(iBearing).output_pressure, data] = next_token(data, "%d");
                [log_dat.bearings(iBearing).output_clearance, data] = next_token(data, "%d");
                [log_dat.bearings(iBearing).output_clearance_derivative, data] = next_token(data, "%d");
                [log_dat.bearings(iBearing).output_velocity, data] = next_token(data, "%d");
                [log_dat.bearings(iBearing).output_stress, data] = next_token(data, "%d");
              endif
            case "hydrodynamic plain bearing2:"
              if (options.hydrodynamic_bearings)
                data = sscanf(data, "%g");
                idx = 0;
                log_dat.bearings(++iBearing).label = int32(data(++idx));
                log_dat.bearings(iBearing).type = "generic";
                log_dat.bearings(iBearing).eta = data(++idx);
                num_flags = int32(data(++idx));

                if (num_flags < 13)
                  error("incompatible file version: \"%s\"", log_filename);
                endif

                switch (num_flags)
                  case {13, 15}
                    node_types = {"hydro", "thermal"};
                  case {19, 21}
                    node_types = {"hydro", "thermal", "flux_x", "flux_z"};
                  otherwise
                    error("incompatible file version: \"%s\"", log_filename);
                endswitch

                for i=1:num_flags
                  bflags(i) = logical(data(++idx));
                endfor

                output_columns = {struct("name", "p", "index", 1, "type", "hydro", "size", 1), ...
                                  struct("name", "pc", "index", 2, "type", "hydro", "size", 1), ...
                                  struct("name", "rho", "index", 3, "type", "hydro", "size", 1), ...
                                  struct("name", "h", "index", 4, "type", "hydro", "size", 1), ...
                                  struct("name", "dh_dt", "index", 5, "type", "hydro", "size", 1), ...
                                  struct("name", "U1x", "index", 6, "type", "hydro", "size", 1), ...
                                  struct("name", "U1z", "index", 6, "type", "hydro", "size", 1), ...
                                  struct("name", "U2x", "index", 6, "type", "hydro", "size", 1), ...
                                  struct("name", "U2z", "index", 6, "type", "hydro", "size", 1), ...
                                  struct("name", "tau_xy_0", "index", 7, "type", "hydro", "size", 1), ...
                                  struct("name", "tau_yz_0", "index", 7, "type", "hydro", "size", 1), ...
                                  struct("name", "tau_xy_h", "index", 7, "type", "hydro", "size", 1), ...
                                  struct("name", "tau_yz_h", "index", 7, "type", "hydro", "size", 1), ...
                                  struct("name", "tauc_xy_0", "index", 8, "type", "hydro", "size", 1), ...
                                  struct("name", "tauc_yz_0", "index", 8, "type", "hydro", "size", 1), ...
                                  struct("name", "wtot", "index", 9, "type", "hydro", "size", 1), ...
                                  struct("name", "dwtot_dt", "index", 9, "type", "hydro", "size", 1), ...
                                  struct("name", "T", "index", 10, "type", "thermal", "size", 1), ...
                                  struct("name", "dT_dt", "index", 10, "type", "thermal", "size", 1), ...
                                  struct("name", "F1", "index", 11, "type", "once", "size", 3), ...
                                  struct("name", "M1", "index", 11, "type", "once", "size", 3), ...
                                  struct("name", "F2", "index", 11, "type", "once", "size", 3), ...
                                  struct("name", "M2", "index", 11, "type", "once", "size", 3), ...
                                  struct("name", "Pff", "index", 12, "type", "once", "size", 1), ...
                                  struct("name", "Pfc", "index", 12, "type", "once", "size", 1), ...
                                  struct("name", "qx", "index", 13, "type", "flux_x", "size", 1), ...
                                  struct("name", "mdotx", "index", 14, "type", "flux_x", "size", 1), ...
                                  struct("name", "qz", "index", 15, "type", "flux_z", "size", 1), ...
                                  struct("name", "mdotz", "index", 16, "type", "flux_z", "size", 1), ...
                                  struct("name", "Qx", "index", 17, "type", "flux_x", "size", 1), ...
                                  struct("name", "Qz", "index", 18, "type", "flux_z", "size", 1), ...
                                  struct("name", "w1", "index", 19, "type", "hydro", "size", 1), ...
                                  struct("name", "w2", "index", 20, "type", "hydro", "size", 1)
                                 };

                bcol_present = false(1, length(output_columns));

                for i=1:length(output_columns)
                  switch (output_columns{i}.type)
                    case {node_types{:}, "once"}
                      bcol_present(i) = true;
                  endswitch
                endfor

                output_columns = {output_columns{find(bcol_present)}};

                log_dat.bearings(iBearing).column_output = struct()([]);
                log_dat.bearings(iBearing).nodes = struct()([]);

                num_nodes = zeros(length(node_types), 1, "int32");
                inode = int32(0);

                for l=1:length(num_nodes)
                  num_nodes(l) = int32(data(++idx));

                  for i=1:num_nodes(l)
                    log_dat.bearings(iBearing).nodes(++inode).number = int32(data(++idx));

                    for j=1:2
                      log_dat.bearings(iBearing).nodes(inode).x(j, 1) = data(++idx);
                    endfor

                    log_dat.bearings(iBearing).nodes(inode).index = int32(data(++idx));
                    log_dat.bearings(iBearing).nodes(inode).type = node_types{l};
                  endfor
                endfor

                j = int32(0);
                icol_node = zeros(length(node_types), 1, "int32");
                icol_other = int32(0);
                inum_per_node_cols = zeros(length(node_types), 1, "int32");
                inum_other_cols = int32(0);

                for i=1:length(bflags)
                  if (bflags(i))
                    for k=1:length(output_columns)
                      if (i == output_columns{k}.index)
                        isize = output_columns{k}.size;

                        switch (output_columns{k}.type)
                          case node_types
                            for l=1:length(node_types)
                              if (strcmp(output_columns{k}.type, node_types{l}))
                                break;
                              endif
                            endfor
                            ++inum_per_node_cols(l);
                            isize *= num_nodes(l);
                            icol = ++icol_node(l);
                          case "once"
                            ++inum_other_cols;
                            icol = ++icol_other;
                        endswitch

                        log_dat.bearings(iBearing).column_output(++j).name = output_columns{k}.name;
                        log_dat.bearings(iBearing).column_output(j).icol = icol;
                        log_dat.bearings(iBearing).column_output(j).size = isize;
                        log_dat.bearings(iBearing).column_output(j).type = output_columns{k}.type;
                      endif
                    endfor
                  endif
                endfor

                ilast_col = int32(0);
                ioffset = int32(0);

                for l=1:length(node_types)
                  for i=1:length(log_dat.bearings(iBearing).column_output)
                    switch (log_dat.bearings(iBearing).column_output(i).type)
                      case node_types{l}
                        log_dat.bearings(iBearing).column_output(i).column_start = log_dat.bearings(iBearing).column_output(i).icol + ioffset;
                        log_dat.bearings(iBearing).column_output(i).column_step = inum_per_node_cols(l);
                        ilast_col = log_dat.bearings(iBearing).column_output(i).icol + ioffset + (num_nodes(l) - 1) * inum_per_node_cols(l);
                        log_dat.bearings(iBearing).column_output(i).column_end = ilast_col;
                    endswitch
                  endfor
                  ioffset = ilast_col;
                endfor

                for i=1:length(log_dat.bearings(iBearing).column_output)
                  switch (log_dat.bearings(iBearing).column_output(i).type)
                    case "once"
                      log_dat.bearings(iBearing).column_output(i).column_start = ++ilast_col;
                      log_dat.bearings(iBearing).column_output(i).column_step = 1;
                      ilast_col += log_dat.bearings(iBearing).column_output(i).size - 1;
                      log_dat.bearings(iBearing).column_output(i).column_end = ilast_col;
                  endswitch
                endfor

                num_elem = int32(data(++idx));

                for i=1:num_elem
                  log_dat.bearings(iBearing).elements(i).number = int32(data(++idx));
                  num_nodes_elem = int32(data(++idx));

                  for j=1:num_nodes_elem
                    log_dat.bearings(iBearing).elements(i).nodes(j) = int32(data(++idx));
                  endfor
                endfor

                mesh_id = uint32(data(++idx));

                switch (bitand(mesh_id, uint32(0x1)))
                  case 1
                    ## cylindrical rigid body bearing
                    log_dat.bearings(iBearing).cylindrical.d = data(++idx);
                    log_dat.bearings(iBearing).cylindrical.D = data(++idx);
                    log_dat.bearings(iBearing).cylindrical.B = data(++idx);

                    switch (bitand(mesh_id, uint32(0xF0)))
                      case 0x10
                        log_dat.bearings(iBearing).cylindrical.dm = log_dat.bearings(iBearing).cylindrical.d;
                        log_dat.bearings(iBearing).cylindrical.mesh_pos = "journal";
                      case 0x20
                        log_dat.bearings(iBearing).cylindrical.dm = log_dat.bearings(iBearing).cylindrical.D;
                        log_dat.bearings(iBearing).cylindrical.mesh_pos = "shell";
                      otherwise
                        ## old version of MBDyn
                    endswitch

                    for i=1:2
                      log_dat.bearings(iBearing).cylindrical.nodes(i).label = int32(data(++idx));

                      for j=1:3
                        log_dat.bearings(iBearing).cylindrical.nodes(i).o(j, 1) = data(++idx);
                      endfor

                      for j=1:9
                        R(j) = data(++idx);
                      endfor

                      log_dat.bearings(iBearing).cylindrical.nodes(i).Rb = reshape(R, 3, 3).';
                    endfor
                  otherwise
                    ## unkown bearing geometry
                endswitch
              endif
            case "totaljoint:"
              if (options.joints)
                [log_dat.joints(++iJoint).label, data] = next_token(data, "%d");
                log_dat.joints(iJoint).type = "totaljoint";
                log_dat.joints(iJoint).nodes = struct("label", cell(1, 2), "offset", cell(1, 2), "orientation", cell(1, 2));

                for k=1:2
                  [log_dat.joints(iJoint).nodes(k).label, data] = next_token(data, "%d");

                  log_dat.joints(iJoint).nodes(k).offset = zeros(1, 3);

                  for i=1:3
                    [log_dat.joints(iJoint).nodes(k).offset(i), data] = next_token(data, "%g");
                  endfor

                  log_dat.joints(iJoint).nodes(k).orientation = zeros(3, 3, 2);

                  for l=1:2
                    for i=1:3
                      for j=1:3
                        [log_dat.joints(iJoint).nodes(k).orientation(i, j, l), data] = next_token(data, "%g");
                      endfor
                    endfor
                  endfor
                endfor
              endif
            case "journal bearing:"
              if (options.joints)
                [log_dat.joints(++iJoint).label, data] = next_token(data, "%d");
                log_dat.joints(iJoint).type = "journal bearing";
                log_dat.joints(iJoint).nodes = struct("label", cell(1, 2), "offset", cell(1, 2), "orientation", cell(1, 2));

                for k=1:2
                  [log_dat.joints(iJoint).nodes(k).label, data] = next_token(data, "%d");
                  log_dat.joints(iJoint).nodes(k).offset = zeros(1, 3);

                  for i=1:3
                    [log_dat.joints(iJoint).nodes(k).offset(i), data] = next_token(data, "%g");
                  endfor

                  if (k == 1)
                    log_dat.joints(iJoint).nodes(k).orientation = zeros(3, 3);

                    for i=1:3
                      for j=1:3
                        [log_dat.joints(iJoint).nodes(k).orientation(i, j), data] = next_token(data, "%g");
                      endfor
                    endfor
                  endif
                endfor
              endif
            case "inline friction:"
              if (options.joints)
                [log_dat.joints(++iJoint).label, data] = next_token(data, "%d");
                log_dat.joints(iJoint).type = "inline friction";
                log_dat.joints(iJoint).nodes = struct("label", cell(1, 2), "offset", cell(1, 2), "orientation", cell(1, 2));

                for k=1:2
                  [log_dat.joints(iJoint).nodes(k).label, data] = next_token(data, "%d");
                  log_dat.joints(iJoint).nodes(k).offset = zeros(1, 3);

                  for i=1:3
                    [log_dat.joints(iJoint).nodes(k).offset(i), data] = next_token(data, "%g");
                  endfor

                  if (k == 1)
                    log_dat.joints(iJoint).nodes(k).orientation = zeros(3, 3);

                    for i=1:3
                      for j=1:3
                        [log_dat.joints(iJoint).nodes(k).orientation(i, j), data] = next_token(data, "%g");
                      endfor
                    endfor
                  endif
                endfor
              endif
          endswitch
        else
          [type, name, value, count] = sscanf(line, "%s %s = %s", "C");

          if (count == 3)
            switch (type)
              case "real"
                value = sscanf(value, "%g", "C");
              case "integer"
                value = int32(sscanf(value, "%d", "C"));
              case "bool"
                value = logical(sscanf(value, "%d", "C"));
              case "string"
                idx1 = strfind(line, " = ");
                if (numel(idx1) >= 1 && numel(line) >= idx1(1) + 4 && numel(line) > 2)
                  value = line((idx1(1) + 4):numel(line) - 2);
                else
                  error("invalid syntax in log file detected");
                endif
            endswitch

            log_dat.vars = setfield(log_dat.vars, name, value);
          endif
        endif
      endwhile
    catch
      err = lasterror();
      err.message = sprintf("error at line %d while parsing file \"%s\": %s", iLine, log_filename, err.message);
      rethrow(err);
    end_try_catch
  unwind_protect_cleanup
    if (fid ~= -1)
      fclose(fid);
    endif
  end_unwind_protect

  if (options.nodes)
    [node_labels, label_idx] = sort([log_dat.nodes.label]);
    log_dat.nodes = log_dat.nodes(label_idx);
  endif

  if (options.vars)
    node_idx = mbdyn_post_node_id_to_node_index([log_dat.nodes.label], log_dat.vars, options.node_id_prefix, options.node_idx_prefix);
    names = fieldnames(node_idx);

    for i=1:length(names)
      log_dat.vars = setfield(log_dat.vars, names{i}, getfield(node_idx, names{i}));
    endfor
  endif

  if (isfield(options, "struct_force_id") && options.vars)
    force_idx = mbdyn_post_node_id_to_node_index(options.struct_force_id, log_dat.vars, options.force_id_prefix, options.force_idx_prefix);

    names = fieldnames(force_idx);

    for i=1:length(names)
      log_dat.vars = setfield(log_dat.vars, names{i}, getfield(force_idx, names{i}));
    endfor
  endif
endfunction

function [val, rema] = next_token(text, format)
  [tok, rema] = strtok(text, " ");

  if (length(tok) == 0)
    error("missing token in string \"%s\"", text);
  endif

  [val, count] = sscanf(tok, format, "C");

  if (count ~= 1)
    error("token \"%s\" cannot be converted", tok);
  endif
endfunction
