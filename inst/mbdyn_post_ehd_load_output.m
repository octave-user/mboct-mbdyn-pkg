## Copyright (C) 2013(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{bearings}] = mbdyn_post_ehd_load_output(@var{mbdyn_filename}, @var{log_dat}, @var{options})
##
## Loads output data from MBDyn's hydrodynamic plain bearings.
##
## @var{mbdyn_filename} @dots{} Location of MBDyn's output file excluding extension (e.g. the -o option).
##
## @var{log_dat} @dots{} Return value from mbdyn_post_load_log.
##
## @var{options} @dots{} Scalar struct with several options.
##
## @var{options}.num_steps @dots{} Total number of time steps in the output file.
##
## @var{options}.output_index @dots{} Optionally load only a subset of time steps indicated by output_index.
##
## @var{options}.interpolate_mesh @dots{} If true, return each field as a three dimensional matrix which can be passed directly to contourf.
##
## @var{options}.loaded_fields @dots{} Optionally load only a subset of available fields indicated by this Cell array of character strings.
##
## @var{options}.verbose @dots{} Enable verbose output.
##
## @end deftypefn

function bearings = mbdyn_post_ehd_load_output(mbdyn_filename, log_dat, options)
  if (nargin < 2 || nargin > 3 || nargout > 1)
    print_usage();
  endif

  if (nargin < 3)
    options = struct();
  endif

  fnames = fieldnames(options);

  for i=1:length(fnames)
    switch (fnames{i})
      case {"num_steps", ...
            "output_index", ...
            "interpolate_mesh", ...
            "assemble_stress_matrix", ...
            "loaded_fields", ...
            "verbose", ...
            "output_step"}
      otherwise
        error("invalid option \"%s\" provided in struct options", fnames{i});
    endswitch
  endfor

  if (~isfield(options, "num_steps"))
    options.num_steps = 1024;
  endif

  if (~isfield(options, "output_index"))
    options.output_index = [];
  endif

  if (~isfield(options, "output_step"))
    options.output_step = 1;
  endif

  if (~isfield(options, "interpolate_mesh"))
    options.interpolate_mesh = true;
  endif

  if (~isfield(options, "assemble_stress_matrix"))
    options.assemble_stress_matrix = false;
  endif

  if (~isfield(options, "loaded_fields"))
    options.loaded_fields = {};
  endif

  if (~isfield(options, "verbose"))
    options.verbose = false;
  endif

  if (~isstruct(log_dat) && isfield(log_dat, "bearings"))
    error("log_dat must be a struct which contains a struct array named bearings");
  endif

  [dir, name, ext] = fileparts(mbdyn_filename);
  mbdyn_filename = fullfile(dir, name);
  usr_filename = mbdyn_post_output_filename(mbdyn_filename, ".usr");
  input_mov_file = mbdyn_post_output_filename(mbdyn_filename, ".mov");

  empty_cell = cell(1, length(log_dat.bearings));

  bearings = struct("label", empty_cell, "type", empty_cell, "cylindrical", empty_cell);

  for i=1:length(log_dat.bearings)
    bearings(i).label = log_dat.bearings(i).label;
    bearings(i).type = log_dat.bearings(i).type;

    switch (log_dat.bearings(i).type)
      case "butenschoen"
        [node_id, X, R, XP, omega] = load_trajectory(input_mov_file, ...
                                                     [log_dat.bearings(i).nodes.label], ...
                                                     options.num_steps, ...
                                                     options.output_index);
        col_count = 16;
        [elem_id, usr_data] = mbdyn_post_load_output(usr_filename, col_count, log_dat.bearings(i).label, options.num_steps, 16, options.output_step);

        if (length(usr_data) > 0)
          bearings(i).Pf = zeros(rows(usr_data{1}), 1);

          for j=1:floor((columns(usr_data{1}))/col_count)
            offset = (j - 1) * col_count;
            e_R2 = usr_data{1}(:, (1:2) + offset);
            e_dot_R2 = usr_data{1}(:, (3:4) + offset);
            e_R1 = zeros(rows(e_R2), columns(e_R2));
            e_dot_R1 = zeros(rows(e_dot_R2), columns(e_dot_R2));
            gamma_R2 = zeros(rows(e_R2), columns(e_R2));

            for k=1:rows(e_R1)
              e_R1(k, :) = (R{1}(:, 1:2, k).' * (R{2}(:, 1:2, k) * e_R2(k, :).')).';
              e_dot_R1(k, :) = (R{1}(:, 1:2, k).' * (R{2}(:, 1:2, k) * e_dot_R2(k, :).')).';
              gamma_R2(k, :) = R{2}(:, 1:2, k).' * R{1}(:, 3, k);
            endfor

            bearings(i).cylindrical.e_R1(:, (1:2) + (j - 1) * 2) = e_R1;
            bearings(i).cylindrical.e_dot_R1(:, (1:2) + (j - 1) * 2) = e_dot_R1;
            bearings(i).cylindrical.e_R2(:, (1:2) + (j - 1) * 2) = e_R2;
            bearings(i).cylindrical.e_dot_R2(:, (1:2) + (j - 1) * 2) = e_dot_R2;
            bearings(i).cylindrical.gamma_R2(:, (1:2) + (j - 1) * 2) = gamma_R2;
            bearings(i).cylindrical.omega1z(:, j) = usr_data{1}(:, 5 + offset);
            bearings(i).cylindrical.omega2z(:, j) = usr_data{1}(:, 6 + offset);
            bearings(i).cylindrical.epsilon(:, j) = usr_data{1}(:, 7 + offset);
            bearings(i).cylindrical.epsilon_dot(:, j) = usr_data{1}(:, 8 + offset);
            bearings(i).cylindrical.delta(:, j) = usr_data{1}(:, 9 + offset);
            bearings(i).F2_R2(:, (1:2) + (j - 1) * 2) = usr_data{1}(:, (10:11) + offset);
            bearings(i).M2z_R2(:, j) = usr_data{1}(:, 12 + offset);
            bearings(i).cylindrical.SoD(:, j) = usr_data{1}(:, 13 + offset);
            bearings(i).cylindrical.SoV(:, j) = usr_data{1}(:, 14 + offset);
            bearings(i).cylindrical.mu(:, j) = usr_data{1}(:, 15 + offset);
            bearings(i).cylindrical.beta(:, j) = usr_data{1}(:, 16 + offset);
            bearings(i).Pf += bearings(i).M2z_R2(:, j) .* (bearings(i).cylindrical.omega2z(:, j) ...
                                                           - bearings(i).cylindrical.omega1z(:, j));
          endfor
        endif
      case "generic"
        node_types = {"hydro", "thermal", "flux_x", "flux_z"};
        node_idx = false(length(log_dat.bearings(i).nodes), length(node_types));

        for j=1:length(log_dat.bearings(i).nodes)
          for k=1:length(node_types)
            switch (log_dat.bearings(i).nodes(j).type)
              case node_types{k}
                node_idx(j, k) = true;
            endswitch
          endfor
        endfor

        for j=1:length(node_types)
          x{j} = [log_dat.bearings(i).nodes(find(node_idx(:, j))).x];
        endfor

        [bearings(i).xi, bearings(i).zi] = meshgrid(unique(x{1}(1, :)), unique(x{1}(2, :))); ## Use corner based grid for output

        inum_columns = int32(0);

        for j=1:length(log_dat.bearings(i).column_output)
          inum_columns = max(inum_columns, log_dat.bearings(i).column_output(j).column_end);
        endfor

        ++inum_columns;

        clear bearing_data;

        if (options.verbose)
          fprintf(stderr, "loading bearing %d (rows = %d, cols = %d)...\n", log_dat.bearings(i).label, options.num_steps, inum_columns);
        endif

        [bearing_label, bearing_data] = mbdyn_post_load_output(usr_filename, ...
                                                               inum_columns, ...
                                                               log_dat.bearings(i).label, ...
                                                               options.num_steps, ...
                                                               inum_columns, ...
                                                               options.output_step);

        if (numel(bearing_data) > 0)
          if (~all(bearing_label == log_dat.bearings(i).label))
            error("labels do not match");
          endif

          if (length(options.output_index) == 0)
            options.output_index = 1:rows(bearing_data{1});
          endif

          cols = struct();

          for j=1:length(log_dat.bearings(i).column_output)
            if (length(options.loaded_fields) == 0)
              f_load_field = true;
            else
              f_load_field = false;

              for k=1:length(options.loaded_fields)
                if (strcmp(options.loaded_fields{k}, log_dat.bearings(i).column_output(j).name))
                  f_load_field = true;
                  break;
                endif
              endfor
            endif

            if (~f_load_field)
              continue;
            endif

            if (options.verbose)
              fprintf(stderr, "\textracting field %s ...\n", log_dat.bearings(i).column_output(j).name);
            endif

            icol_start = log_dat.bearings(i).column_output(j).column_start;
            icol_step = log_dat.bearings(i).column_output(j).column_step;
            icol_end = log_dat.bearings(i).column_output(j).column_end;
            z = bearing_data{1}(options.output_index, [icol_start:icol_step:icol_end]);

            col_name = log_dat.bearings(i).column_output(j).name;

            switch (log_dat.bearings(i).column_output(j).type)
              case {"hydro", "thermal", "flux_x", "flux_z"}
                if (options.interpolate_mesh)
                  l = -1;
                  for m=1:length(node_types)
                    if (strcmp(log_dat.bearings(i).column_output(j).type, node_types{m}))
                      l = m;
                      break;
                    endif
                  endfor

                  if (columns(x{l}) == 0)
                    continue;
                  endif

                  zi = zeros(rows(bearings(i).xi), columns(bearings(i).xi), rows(z));

                  tri_grid = delaunay(x{l}(1, :), x{l}(2, :));

                  for k=1:rows(z)
                    zi(:, :, k) = griddata_prepared(x{l}(1, :), x{l}(2, :), z(k, :), bearings(i).xi, bearings(i).zi, tri_grid);

                    if (any(any(isnan(zi(:, :, k)))))
                      idx_row = ones(1, size(zi, 1), "int32");
                      idx_col = ones(1, size(zi, 2), "int32");

                      for n=1:length(idx_row)
                        if (all(isnan(zi(n, :, k))))
                          idx_row(n) = 0;
                        endif
                      endfor

                      for n=1:length(idx_col)
                        if (all(isnan(zi(:, n, k))))
                          idx_col(n) = 0;
                        endif
                      endfor

                      idx_row = find(idx_row);
                      idx_col = find(idx_col);

                      zi(:, :, k) = interp2(bearings(i).xi(idx_row, idx_col), ...
                                            bearings(i).zi(idx_row, idx_col), ...
                                            zi(idx_row, idx_col, k), ...
                                            bearings(i).xi, ...
                                            bearings(i).zi, ...
                                            "linear");
                    endif
                  endfor

                  cols = setfield(cols, col_name, zi);
                endif

                cols = setfield(cols, cstrcat(col_name, "_n"), z);
              case "once"
                cols = setfield(cols, col_name, z);
            endswitch
          endfor

          if (isfield(cols, "Pff") && isfield(cols, "Pfc"))
            bearings(i).Pf = cols.Pff + cols.Pfc;
          endif

          bearings(i).columns = cols;
          clear cols;

          if (options.assemble_stress_matrix)
            cidx = data = zeros(4 * length(log_dat.bearings(i).elements), 1);
            dataidx = int32(0);

            for j=1:length(log_dat.bearings(i).elements)
              switch (length(log_dat.bearings(i).elements(j).nodes))
                case 4
                  N = length(log_dat.bearings(i).elements(j).nodes);
                  cidx(dataidx + int32(1:N), 1) = log_dat.bearings(i).elements(j).nodes.';
                  xe = [log_dat.bearings(i).nodes(log_dat.bearings(i).elements(j).nodes).x](1, :);
                  ze = [log_dat.bearings(i).nodes(log_dat.bearings(i).elements(j).nodes).x](2, :);
                  Ae =((xe(3)-xe(1))*ze(4)+(ze(1)-ze(3))*xe(4)+xe(2)*ze(3)-ze(2)*xe(3)+xe(1)*ze(2)-ze(1)*xe(2))/2;
                  data(dataidx + int32(1:N)) = 0.25 * Ae;
                  dataidx += N;
              endswitch
            endfor

            cidx = cidx(1:dataidx);
            data = data(1:dataidx);
            ridx = ones(length(cidx), 1);

            bearings(i).matrices.A_tau = full(sparse(ridx, cidx, data));
          endif

          if (isfield(log_dat.bearings(i), "cylindrical") && isstruct(log_dat.bearings(i).cylindrical))
            bearings(i).Vn = 0.5 * (log_dat.bearings(i).cylindrical.D - log_dat.bearings(i).cylindrical.d) ...
                             * log_dat.bearings(i).cylindrical.d * pi * log_dat.bearings(i).cylindrical.B;

            if (isfield(bearings(i).columns, "h"))
              bearings(i).V = zeros(size(bearings(i).columns.h, 3), 1);
              for k=1:size(bearings(i).columns.h, 3)
                bearings(i).V(k) = cylindrical_bearing_integrate(bearings(i).xi, bearings(i).zi, bearings(i).columns.h(:, :, k));
              endfor
            endif

            if (isfield(bearings(i).columns, "dh_dt"))
              bearings(i).dV_dt = zeros(size(bearings(i).columns.dh_dt, 3), 1);
              for k=1:size(bearings(i).columns.dh_dt, 3)
                bearings(i).dV_dt(k) = cylindrical_bearing_integrate(bearings(i).xi, bearings(i).zi, bearings(i).columns.dh_dt(:, :, k));
              endfor
            endif

            node_id_bearing = [log_dat.bearings(i).cylindrical.nodes.label];

            [node_id, X, R, XP, omega] = load_trajectory(input_mov_file, ...
                                                         node_id_bearing, ...
                                                         options.num_steps, ...
                                                         options.output_index);

            idx_shaft = find(node_id == node_id_bearing(1));
            idx_bearing = find(node_id == node_id_bearing(2));

            forces = {"F1", "M1", "F2", "M2"};

            have_forces = true;

            for j=1:numel(forces)
              if (~isfield(bearings(i).columns, forces{j}))
                have_forces = false;
                break;
              endif
            endfor

            if (have_forces)
              [bearings(i).cylindrical.e_R1, ...
               bearings(i).cylindrical.e_dot_R1, ...
               bearings(i).cylindrical.e_R2, ...
               bearings(i).cylindrical.e_dot_R2, ...
               bearings(i).cylindrical.epsilon, ...
               bearings(i).cylindrical.epsilon_dot, ...
               bearings(i).cylindrical.delta, ...
               bearings(i).cylindrical.delta_dot, ...
               bearings(i).cylindrical.gamma_R2, ...
               bearings(i).cylindrical.theta_R2, ...
               bearings(i).cylindrical.omega1z, ...
               bearings(i).cylindrical.omega2z, ...
               bearings(i).cylindrical.omega_res, ...
               bearings(i).P] = ...
              cylindrical_bearing_kinematics(log_dat.bearings(i).cylindrical.nodes(1).Rb, ...
                                             log_dat.bearings(i).cylindrical.nodes(1).o, ...
                                             X{idx_shaft}, ...
                                             R{idx_shaft}, ...
                                             XP{idx_shaft}, ...
                                             omega{idx_shaft}, ...
                                             log_dat.bearings(i).cylindrical.nodes(2).Rb, ...
                                             log_dat.bearings(i).cylindrical.nodes(2).o, ...
                                             X{idx_bearing}, ...
                                             R{idx_bearing}, ...
                                             XP{idx_bearing}, ...
                                             omega{idx_bearing}, ...
                                             bearings(i).columns.F1.', ...
                                             bearings(i).columns.M1.', ...
                                             bearings(i).columns.F2.', ...
                                             bearings(i).columns.M2.', ...
                                             log_dat.bearings(i).cylindrical.D,
                                             log_dat.bearings(i).cylindrical.d);
            endif
          endif
        else
          warning("no valid records found for bearing %d in file \"%s\"", log_dat.bearings(i).label, usr_filename);
        endif
    endswitch
  endfor
endfunction

function [field, col_out] = extract_field(data, col, bearing)
  field_size = bearing.M * bearing.N;
  field = nan(bearing.M, bearing.N, rows(data));
  for row=1:rows(data)
    field(:, :, row) = reshape(data(row, col + (1:field_size)), bearing.N, bearing.M).';
  endfor
  col_out = col + field_size;
endfunction

function [node_id, X, R, XP, omega] = load_trajectory(input_mov_file, filter_node_id, append_rows, output_index)
  [node_id, trajectory, velocity, acceleration, orientation_description]= mbdyn_post_load_output_mov(input_mov_file, filter_node_id, append_rows);

  if (length(output_index) == 0)
    output_index = 1:rows(trajectory{1});
  endif

  for i=1:length(trajectory)
    X{i} = trajectory{i}(output_index, 1:3).';
    R{i} = zeros(3,3,length(output_index));
    XP{i} = velocity{i}(output_index, 1:3).';
    omega{i} = velocity{i}(output_index, 4:6).';

    phi_ij = trajectory{i}(output_index, 4:6).';

    switch( orientation_description{i} )
      case "euler123"
        R_ij = euler123_to_rotation_matrix(phi_ij);
      case "euler313"
        R_ij  = euler313_to_rotation_matrix(phi_ij);
      case "euler321"
        R_ij  = euler321_to_rotation_matrix(phi_ij);
      case "phi"
        R_ij  = rotation_vector_to_rotation_matrix(phi_ij);
      otherwise
        error("orientation description \"%s\" is not implemented!",orientation_description{i});
    endswitch

    R{i} = R_ij;
  endfor
endfunction

function [e_R1, e_dot_R1, ...
          e_R2, e_dot_R2, ...
          epsilon, epsilon_dot, ...
          delta, delta_dot, ...
          gamma_R2, ...
          theta_R2, ...
          omega1z, omega2z, ...
          omega_res, Pf] = cylindrical_bearing_kinematics(Rb1, o1, X1, R1, XP1, omega1, ...
                                                          Rb2, o2, X2, R2, XP2, omega2, ...
                                                          F1, M1, F2, M2, D, d)
  N = columns(X1);
  e = zeros(N, 2);
  e_dot = zeros(N, 2);
  delta = zeros(N, 1);
  delta_dot = zeros(N, 1);
  epsilon = zeros(N, 1);
  epsilon_dot = zeros(N, 1);
  gamma_R2 = zeros(N, 2);
  theta_R2 = zeros(N, 3);
  omega1z = zeros(N, 1);
  omega2z = zeros(N, 1);
  Pf = zeros(N, 1);
  s = D - d;

  for i=1:N
    e = X1(:, i) + R1(:, :, i) * o1 - X2(:, i) - R2(:, :, i) * o2;

    V1 = XP1(:, i) + cross(omega1(:, i), R1(:, :, i) * o1);
    V2 = XP2(:, i) + cross(omega2(:, i), R2(:, :, i) * o2);

    e_dot = V1 - V2;

    e_R1(i, :) = (Rb1(:, 1:2).' * (R1(:, :, i).' * e));
    e_dot_R1(i, :) = (Rb1(:, 1:2).' * (R1(:, :, i).' * e_dot));

    e_R2(i, :) = (Rb2(:, 1:2).' * (R2(:, :, i).' * e));
    e_dot_R2(i, :) = (Rb2(:, 1:2).' * (R2(:, :, i).' * e_dot));

    gamma_R2(i, :) = (Rb2(:, 1:2).' * (R2(:, :, i).' * (R1(:, :, i) * Rb1(:, 3)))).';

    theta_R2(i, :) = rotation_matrix_to_rotation_vector(Rb2.' * (R2(:, :, i).' * (R1(:, :, i) * Rb1))).';

    %% angle of the position with minimum clearance between shaft and bearing
    delta(i) = atan2(e_R2(i, 2), e_R2(i, 1));

    %% angle of the velocity vector of the shaft relative to the bearing
    phi = atan2(e_dot_R2(i, 2), e_dot_R2(i, 1));

    %% angle between velocity vector and minimum clearance
    kappa = phi - delta(i);

    %% absolute value of the eccentricity of the shaft inside the bearing
    abs_e = norm(e_R2(i, :));

    %% absolute value of the velocity of the shaft relative to the bearing
    abs_e_dot = norm(e_dot_R2(i, :));

    %% time derivative of the relative eccentricity of the shaft
    epsilon_dot(i) = 2 * cos(kappa) * abs_e_dot / s;

    %% relative eccentricity of the shaft
    epsilon(i) = 2 * abs_e / s;

    norm_e_R2 = ( e_R2(i, 2)^2 + e_R2(i, 1)^2 );

    if (norm_e_R2 == 0)
      delta_dot(i) = 0;
    else
      delta_dot(i) = ( e_R2(i, 1) * e_dot_R2(i, 2) - e_R2(i, 2) * e_dot_R2(i, 1) ) / norm_e_R2;
    endif

    omega1z(i) = (Rb2(:, 3).' * R2(:, :, i).' * omega1(:, i));
    omega2z(i) = (Rb2(:, 3).' * R2(:, :, i).' * omega2(:, i));
    Pf(i) = XP1(:, i).' * F1(:, i) + omega1(:, i).' * M1(:, i) ...
            + XP2(:, i).' * F2(:, i) + omega2(:, i).' * M2(:, i);
  endfor

  omega_res = omega1z + omega2z - 2 * delta_dot;
endfunction

function v = cylindrical_bearing_integrate(x, z, h)
  v = sum(sum(0.25 * (h(1:end - 1, 2:end - 1) + h(2:end, 2:end - 1) + h(2:end, 3:end) + h(1:end - 1, 3:end)) ...
              * 0.5 .* (x(1:end - 1, 3:end) + x(2:end, 3:end) - x(1:end - 1, 2:end - 1) - x(2:end, 2:end - 1)) ...
              * 0.5 .* (z(2:end, 2:end - 1) + z(2:end, 3:end) - z(1:end - 1, 2:end - 1) - z(1:end - 1, 3:end))));
endfunction
