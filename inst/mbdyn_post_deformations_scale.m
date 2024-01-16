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
## @deftypefn {Function File} mbdyn_post_deformations_scale(@var{log_file}, @var{input_mov_file}, @var{output_mov_file}, @var{node_groups}, @var{idx_t_ref }, @var{filter_node_id}, @var{output_frequency}, @var{append_rows})
##
## Scales up the deformations of selected nodes.
##
## @var{log_file} @dots{} MBDyn .log file which contains information about the orientation description of structural nodes.
##
## @var{input_mov_file} @dots{} Input MBDyn .mov file which should be scaled.
##
## @var{output_mov_file} @dots{} Output MBDyn .mov file where the scaled deformations are written to.
##
## @var{node_groups} @dots{} Struct or struct array that contains information about the nodes who's deformations should be scaled.
##
## @var{node_groups}(i).R0ref_node_id(1) @dots{} MBDyn node id of the origin point of the reference frame which is used in order to compute the rigid body motion.
##
## @var{node_groups}(i).R0ref_node_id(2) @dots{} MBDyn node id of a point at the 1-axis of the reference frame which is used in order to compute the rigid body motion.
##
## @var{node_groups}(i).R0ref_node_id(3) @dots{} MBDyn node id of a point in the 12-plane of the reference frame which is used in order to compute the rigid body motion.
##
## @var{node_groups}(i).X0_node_id @dots{} MBDyn node id of the reference point used in order to compute the rigid body motion.
##
## @var{node_groups}(i).scale @dots{} Scale factor for the elastic deformations.
##
## @var{idx_t_ref} @dots{} Index of the time step that represents the reference configuration.
##
## @var{filter_node_id} @dots{} Vector of MBDyn node id's which should be processed. If not present, all nodes are processed.
##
## @var{output_frequency} @dots{} Output interval for the file <@var{output_mov_file}>.
##
## @var{append_rows} @dots{} Used to optimize memory allocation. This parameter should be set equal to the number of time steps in the file <@var{input_mov_file}>.
##
## @end deftypefn

function mbdyn_post_deformations_scale(log_file, input_mov_file, output_mov_file, node_groups, idx_t_ref = -1, filter_node_id = [], output_frequency = 1, append_rows = 1024, extra_nodes = struct(), log_dat_nodes = [])
  if (~exist("log_file", "var") ||
      ~exist("input_mov_file", "var") ||
      ~exist("output_mov_file", "var") ||
      ~exist("node_groups", "var"))
    print_usage();
  endif

  if (~isstruct(node_groups))
    error("node_groups is not a valid struct array");
  endif

  if (idx_t_ref == -1)
    if (length(log_dat_nodes) == 0)
      [dir_mbdyn, name_mbdyn, ext_mbdyn] = fileparts(input_mov_file);
      log_filename = mbdyn_post_output_filename(fullfile(dir_mbdyn, name_mbdyn), ".log");
      log_dat = mbdyn_post_load_log(log_filename, struct("nodes", true));
    else
      log_dat.nodes = log_dat_nodes;
    endif
  else
    log_dat = struct();
  endif

  [node_id, X, R, velocity, acceleration, orientation_description] = load_deformation_input_file(input_mov_file, filter_node_id, append_rows, idx_t_ref);

  N = length(X);

  if (isfield(extra_nodes, "node_id"))
    node_id(end + (1:length(extra_nodes.node_id))) = extra_nodes.node_id;

    for i=1:length(extra_nodes.node_id)
      X{N + i} = extra_nodes.X{i};
      R{N + i} = extra_nodes.R{i};
    endfor
  endif

  if (isfield(node_groups, "stages"))
    for i=1:length(node_groups.stages)
      [X, R] = scale_deformation(node_id, X, R, node_groups.stages(i).node_groups, idx_t_ref, log_dat);
    endfor
  else
    [X, R] = scale_deformation(node_id, X, R, node_groups,idx_t_ref, log_dat);
  endif

  write_output_file(node_id, X, R, velocity, acceleration, orientation_description, output_mov_file, output_frequency, N);
endfunction

function [node_id, X, R, velocity, acceleration, orientation_description] = load_deformation_input_file(input_mov_file, filter_node_id, append_rows, idx_t_ref)
  [node_id, trajectory, velocity, acceleration, orientation_description]= mbdyn_post_load_output_mov(input_mov_file, filter_node_id, append_rows);

  R = cell(1, numel(trajectory));

  for i=1:length(trajectory)
    if (idx_t_ref > rows(trajectory{i}))
      error("idx_t_ref = %d > rows(trajectory{%d}) = %d", idx_t_ref, i, rows(trajectory{i}));
    endif

    X{i} = trajectory{i}(:,1:3).';

    phi_ij = trajectory{i}(:, 4:6).';

    switch (orientation_description{i})
      case "euler123"
        R_ij = euler123_to_rotation_matrix(phi_ij);
      case "euler313"
        R_ij  = euler313_to_rotation_matrix(phi_ij);
      case "euler321"
        R_ij  = euler321_to_rotation_matrix(phi_ij);
      case "phi"
        R_ij  = rotation_vector_to_rotation_matrix(phi_ij);
      otherwise
        error("orientation description \"%s\" is not implemented!", orientation_description{i});
    endswitch

    R{i} = R_ij;
  endfor
endfunction

function write_output_file(node_id, X, R, velocity, acceleration, orientation_description, output_mov_file, output_frequency, N)
  fid = -1;

  unwind_protect
    [fid, msg] = fopen(output_mov_file, "wt");

    if (fid == -1)
      error("could not open file \"%s\": %s", output_mov_file, msg);
    endif

    if (length(X) > 0)
      Phi = cell(1, N);

      for j=1:N
        switch(orientation_description{j})
          case "euler123"
            Phi{j} = rotation_matrix_to_euler123(R{j}) * (180 / pi);
          case "euler313"
            Phi{j} = rotation_matrix_to_euler313(R{j}) * (180 / pi);
          case "euler321"
            Phi{j} = rotation_matrix_to_euler321(R{j}) * (180 / pi);
          case "phi"
            Phi{j} = rotation_matrix_to_rotation_vector(R{j});
          otherwise
            error("orientation description \"%s\" not supported!", orientation_description{j});
        endswitch
      endfor

      for i=1:output_frequency:columns(X{1})
        for j=1:N
          fprintf(fid, "%d ", node_id(j));
          fprintf(fid, "%e ", X{j}(:, i));
          fprintf(fid, "%e ", Phi{j}(:, i));
          fprintf(fid, "%e ", velocity{j}(i, :));

          if (columns(acceleration{j}) > 0)
            fprintf(fid, "%e ", acceleration{j}(i, :));
          endif
          fprintf(fid, "\n");
        endfor
      endfor
    endif
  unwind_protect_cleanup
    if (fid ~= -1)
      fclose(fid);
    endif
  end_unwind_protect
endfunction

function [X1, R1] = scale_deformation(node_id, X, R, node_groups, idx_t_ref, log_dat)
  X1 = X;
  R1 = R;

  if (~isfield(node_groups, "R0ref_node_id") ||
      ~isfield(node_groups, "X0_node_id") ||
      ~isfield(node_groups, "node_id_scale") ||
      ~isfield(node_groups, "scale"))
    error("node_groups is not a valid struct array");
  endif

  if (idx_t_ref == -1)
    node_id_log = [log_dat.nodes(:).label];
  endif

  for i=1:length(node_groups)
    X0_node_idx = find(node_groups(i).X0_node_id == node_id);

    if (length(X0_node_idx) == 0)
      error("node_groups(%d).X0_node_id=%d not found!", i, node_groups(i).X0_node_id);
    endif

    switch (length(node_groups(i).R0ref_node_id))
      case 1
        R0ref_node_idx = find(node_groups(i).R0ref_node_id == node_id);

        if (idx_t_ref == -1)
          R0ref = log_dat.nodes(find(node_groups(i).R0ref_node_id == node_id_log)).R0;
        else
          R0ref = R{R0ref_node_idx}(:, :, idx_t_ref);
        endif
      case 3
        for j=1:3
          R0ref_node_idx_j = find(node_groups(i).R0ref_node_id(j) == node_id);

          if (length(R0ref_node_idx_j) == 0)
            error("node_groups(%d).R0ref_node_id(%d) not found!", i, j);
          endif

          R0ref_node_idx(j) = R0ref_node_idx_j;

          if (idx_t_ref == -1)
            XR0ref(:, j) = log_dat.nodes(find(node_groups(i).R0ref_node_id(j) == node_id_log)).X0;
          else
            XR0ref(:, j) = X{R0ref_node_idx(j)}(:, idx_t_ref);
          endif
        endfor

        R0ref = rotation_matrix_from_three_points(XR0ref);
      otherwise
        error("length(node_groups(%d).R0ref_node_id) is invalid!", i);
    endswitch

    if (idx_t_ref == -1)
      X0ref = log_dat.nodes(find(node_groups(i).X0_node_id == node_id_log)).X0;
    else
      X0ref = X{X0_node_idx}(:, idx_t_ref);
    endif

    scale = node_groups(i).scale;

    for j=1:length(node_groups(i).node_id_scale)
      node_idx_j = find(node_groups(i).node_id_scale(j) == node_id);

      if (length(node_idx_j) == 0)
        error("node_groups(%d).node_id_scale(%d)=%d not found!", i, j, node_groups(i).node_id_scale(j));
      endif

      if (idx_t_ref == -1)
        node_idx = find(node_groups(i).node_id_scale(j) == node_id_log);

        if (length(node_idx) == 0)
          error("node_id %d not found in log_dat.nodes", node_groups(i).node_id_scale(j));
        endif

        Rref = log_dat.nodes(node_idx).R0;
        Xref = log_dat.nodes(node_idx).X0;
      else
        Rref = R{node_idx_j}(:, :, idx_t_ref);
        Xref = X{node_idx_j}(:, idx_t_ref);
      endif

      for k=1:columns(X{node_idx_j})
        X0 = X{X0_node_idx}(:, k);

        switch (length(R0ref_node_idx))
          case 1
            R0 = R{R0ref_node_idx}(:, :, k);
          case 3
            for l=1:3
              XR0(:,l) = X{R0ref_node_idx(l)}(:, k);
            endfor
            R0 = rotation_matrix_from_three_points(XR0);
          otherwise
            error("length(R0ref_node_idx) is invalid");
        endswitch

        X1{node_idx_j}(:, k) = scale * X{node_idx_j}(:, k) + (1 - scale) * (X0 + R0 * R0ref.' * (Xref - X0ref));
        DeltaR = Rref.' * R0ref * R0.' * R{node_idx_j}(:,:,k);
        DeltaPhi = scale * rotation_matrix_to_rotation_vector(DeltaR);
        DeltaR = rotation_vector_to_rotation_matrix(DeltaPhi);
        R1{node_idx_j}(:, :, k) = R0 * R0ref.' * Rref * DeltaR;
      endfor
    endfor
  endfor
endfunction

function R = rotation_matrix_from_three_points(X)
  persistent A = 1;
  persistent B = 2;
  persistent C = 3;

  a = X(:, B) - X(:, A);
  b = X(:, C) - X(:, A);

  R = hinge_to_rotation_matrix(a, b);
endfunction

function R = hinge_to_rotation_matrix(a, b)
  e1 = a;
  e3 = cross(a, b);
  e2 = cross(e3, e1);
  R = [e1 / norm(e1), e2 / norm(e2), e3 / norm(e3)];
endfunction
