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

  if (~isfield(options, "num_rows"))
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
          assert(bearing_label, log_dat.bearings(i).label);

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

%!test
%! ## References:
%! ## Hans Juergen Butenschoen
%! ## Das hydrodynamische, zylindrische Gleitlager endlicher Breite unter instationaerer Belastung, Karlsruhe 1976
%! ## Test case for rotation
%! B_d_r = [1, 1/2, 1/2.5, 1/3, 1/4, 1/5, 1/6, 1/8];
%!
%! epsilon_r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99];
%!
%! So_r =  [0.1196 0.0368  0.0243  0.0171  0.0098  0.0063  0.0044  0.0025
%!          0.2518 0.0783  0.0518  0.0366  0.0210  0.0136  0.0095  0.0054
%!          0.4091 0.1304  0.0867  0.0615  0.0354  0.0229  0.0160  0.0091
%!          0.6108 0.2026  0.1357  0.0968  0.0560  0.0364  0.0254  0.0144
%!          0.8903 0.3124  0.2117  0.1522  0.0888  0.0579  0.0406  0.0231
%!          1.3146 0.4982  0.3435  0.2496  0.1476  0.0969  0.0682  0.0390
%!          2.0432 0.8595  0.6075  0.4492  0.2708  0.1797  0.1274  0.0732
%!          3.5663 1.7339  1.2756  0.9687  0.6043  0.4085  0.2930  0.1706
%!          8.4352 5.0881  4.0187  3.2201  2.1595  1.5261  1.1263  0.6776
%!          18.7895 13.3083 11.2225 9.4993  6.9248  5.1833  3.9831  2.5202
%!          24.1172 17.7934 15.2823 13.1525 9.8517  7.5257  5.8712  3.7907
%!          33.1297 25.5920 22.4503 19.7034 15.2697 11.9804 9.5425  6.3397
%!          51.4774 41.9230 37.7412 33.9523 27.5040 22.3556 18.3874 12.7695
%!          107.7868 93.7881 87.2906 81.1597 70.0359 60.3874 52.1425 39.2568];
%!
%! beta_r = [79.410 81.767  82.119  82.299  82.481  82.561  82.608  82.653
%!           73.858 75.142  75.283  75.344  75.287  75.406  75.414  75.422
%!           68.264 68.493  68.423  68.368  68.204  68.261  68.238  68.211
%!           62.571 61.778  61.544  61.382  61.208  61.115  61.062  61.007
%!           56.705 54.993  54.603  54.348  54.069  53.922  53.863  53.784
%!           50.536 48.049  47.521  47.191  46.825  46.647  46.554  46.449
%!           43.859 40.803  40.156  39.756  39.326  39.108  38.983  38.865
%!           36.235 32.938  32.216  31.761  31.249  30.988  30.840  30.692
%!           26.482 23.566  22.849  22.368  21.790  21.476  21.289  21.089
%!           19.450 17.265  16.648  16.207  15.632  15.291  15.075  14.832
%!           17.609 15.660  15.089  14.669  14.109  13.768  13.547  13.262
%!           15.484 13.818  13.313  12.929  12.299  12.062  11.838  11.570
%!           12.903 11.598  11.178  10.850  10.375  10.057  9.835   9.558
%!            9.416  8.587   8.301   8.066   7.703   7.440   7.242   6.975] * pi / 180;
%!
%! Q_r = [0.1603  0.0939  0.0768  0.0648  0.0492  0.0396  0.0331  0.0249
%!        0.3196  0.1878  0.1536  0.1296  0.0984  0.0792  0.0662  0.0498
%!        0.4765  0.2816  0.2304  0.1944  0.1476  0.1188  0.0993  0.0747
%!        0.6318  0.3755  0.3072  0.2592  0.1968  0.1583  0.1324  0.0996
%!        0.7852  0.4694  0.384   0.324   0.246   0.198   0.1655  0.1245
%!        0.9374  0.5634  0.461   0.3889  0.2953  0.2376  0.1986  0.1494
%!        1.0888  0.6578  0.5383  0.454   0.3446  0.2772  0.2317  0.1743
%!        1.2392  0.7529  0.6159  0.5193  0.394   0.3169  0.2648  0.1992
%!        1.3898  0.8453  0.6944  0.5852  0.4436  0.3567  0.298   0.2241
%!        1.4655  0.8984  0.7343  0.6186  0.4687  0.3767  0.3147  0.2366
%!        1.4808  0.9084  0.7425  0.6254  0.4737  0.3807  0.3181  0.2391
%!        1.4968  0.9185  0.7507  0.6322  0.4788  0.3848  0.3214  0.2417
%!        1.5124  0.9287  0.7585  0.6391  0.484   0.3889  0.3248  0.2442
%!        1.5277  0.9351  0.7674  0.6461  0.4892  0.393   0.3282  0.2467];
%!
%! So = beta = mu = Q = dQ = mu_r = nan(numel(epsilon_r), numel(B_d_r));
%! test_freq = 4;
%! verbose = false;
%! for j=1:test_freq:numel(epsilon_r)
%! for k=1:test_freq:numel(B_d_r)
%! param.d1 = 10e-3;
%! param.Psi = 1e-3;
%! param.M = int32(20);
%! param.N = int32(85);
%! param.output_bearing_data = true;
%! param.B = param.d1 * B_d_r(k);
%! param.epsilon = epsilon_r(j);
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer hyd_node_id_inlet = 2003;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: integer genel_id_inlet = 4005;\n");
%!     fputs(fd, "set: real omega1z = 2 * pi * 50;\n");
%!     fputs(fd, "set: real omega2z = 2 * pi * 15;\n");
%!     fputs(fd, "set: real v1z = 0;\n");
%!     fputs(fd, "set: real v2z = 0;\n");
%!     fputs(fd, "set: real n = 1e-3;\n");
%!     fputs(fd, "set: integer K = 36000;\n");
%!     fputs(fd, "set: real t1 = abs(2 * pi * n / omega1z);\n");
%!     fputs(fd, "set: real dt = t1 / (n * K);\n");
%!     fputs(fd, "set: real D2 = d1 / (1. - Psi);\n");
%!     fputs(fd, "set: real eta = 1e-3;\n");
%!     fputs(fd, "set: real rho = 3600. * 1000.; # [l/h]\n\n");
%!     fputs(fd, "set: real epsilon_t0 = epsilon;\n");
%!     fputs(fd, "set: real epsilon_t1 = epsilon;\n");
%!     fputs(fd, "set: real delta_t0 = pi;\n");
%!     fputs(fd, "set: real delta_t1 = pi + omega2z * t1;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real p_pB2 = 0;\n");
%!     fputs(fd, "set: real p_mB2 = 0;\n");
%!     fputs(fd, "set: real p_in = 0;\n");
%!     fputs(fd, "set: real pmax = 500000.;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 30;\n");
%!     fputs(fd, "        tolerance: 1e-10, test, minmax;\n");
%!     fputs(fd, "        linear solver: umfpack, map, colamd, scale, row max column max, always, max iterations, 100;\n");
%!     fputs(fd, "        nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e5;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: bdf;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 10;\n");
%!     fputs(fd, "        derivatives coefficient: 110;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "        structural nodes: 2;\n");
%!     fputs(fd, "        joints: 2;\n");
%!     fputs(fd, "        loadable elements: 1;\n");
%!     fputs(fd, "        hydraulic nodes: 3;\n");
%!     fputs(fd, "        genels: 3;\n");
%!     fputs(fd, "        print: dof stats, to file;\n");
%!     fputs(fd, "        print: equation description, to file;\n");
%!     fputs(fd, "        print: dof description, to file;\n");
%!     fputs(fd, "        print: element connection, to file;\n");
%!     fputs(fd, "        print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 8;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        v1z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing, null,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        v2z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_inlet, p_in;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            mult, const, v1z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, v2z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!     fputs(fd, "    genel: genel_id_inlet, clamp,\n");
%!     fputs(fd, "        hyd_node_id_inlet, hydraulic, p_in;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!     fputs(fd, "                density, rho,\n");
%!     fputs(fd, "                1.,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                viscosity, eta,\n");
%!     fputs(fd, "                temperature, 0,\n");
%!     fputs(fd, "            viscosity vapor, eta,\n");
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at bearing,\n");
%!     fputs(fd, "                bearing width, B,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure coupling conditions radial, 1,\n");
%!     fputs(fd, "                position, 0.5 * D2 * (delta_t0 + pi), 0.,\n");
%!     fputs(fd, "                rectangle, width, D2 * pi / (N - 1), height, B,\n");
%!     fputs(fd, "                pressure node, hyd_node_id_inlet,\n");
%!     fputs(fd, "            pressure dof scale, pmax,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, no,\n");
%!     fputs(fd, "            output velocity, no,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   if (~verbose)
%!     opt_sol.logfile = [output_file, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%!   d = res.log_dat.bearings.cylindrical.d;
%!   D = res.log_dat.bearings.cylindrical.D;
%!   B = res.log_dat.bearings.cylindrical.B;
%!   s = (D - d);
%!   Psi = (D - d) / D;
%!   eta = res.log_dat.bearings.eta;
%!   node_idx_2 = find(res.node_id == res.log_dat.bearings.cylindrical.nodes(2).label);
%!   Phi2 = res.trajectory{node_idx_2}(:, 4:6).';
%!   R2 = euler123_to_rotation_matrix(Phi2);
%!   Rb1 = res.log_dat.bearings.cylindrical.nodes(1).Rb;
%!   Rb2 = res.log_dat.bearings.cylindrical.nodes(2).Rb;
%!   F1 = -res.bearings.columns.F1.';
%!   M1 = -res.bearings.columns.M1.';
%!   omega1z = res.bearings.cylindrical.omega1z;
%!   omega2z = res.bearings.cylindrical.omega2z;
%!   omega_res = res.bearings.cylindrical.omega_res;
%!   epsilon = res.bearings.cylindrical.epsilon;
%!   delta = res.bearings.cylindrical.delta;
%!   Q_out1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end) / res.log_dat.vars.rho;
%!   Q_out2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end) / res.log_dat.vars.rho;
%!   Q_in = -res.genel_data{res.log_dat.vars.genel_idx_inlet}(end) / res.log_dat.vars.rho;
%!   P = norm(Rb2(:, 1:2).' * R2.' * F1);
%!   alpha = atan2(Rb2(:, 2).' * R2.' * F1, Rb2(:, 1).' * R2.' * F1);
%!   beta(j, k) = sign(omega_res) * (delta - alpha);
%!   beta(j, k) = mod(beta(j, k) + pi, 2 * pi) - pi;
%!   mu(j, k) = 2 * Rb2(:, 3).' * R2.' * M1 / (P * d) * sign(omega1z - omega2z);
%!   So(j, k) = P * Psi^2 / (B * D * eta * abs(omega_res));
%!   Q(j, k) = (Q_out1 + Q_out2) / ((0.5 * D)^3 * Psi * abs(omega_res));
%!   dQ(j, k) = (Q_out1 + Q_out2 + Q_in) / Q_in;
%!   mu_r(j, k) = Psi * (abs((omega1z - omega2z) / omega_res) * pi / (sqrt(1 - epsilon^2) * So_r(j, k)) + sin(beta_r(j, k)) * abs(epsilon) / 2);
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
%! endfor
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   semilogy(epsilon_r(1:test_freq:end), So_r(1:test_freq:end, i), "-;So_r(epsilon);0");
%!   semilogy(epsilon_r(1:test_freq:end), So(1:test_freq:end, i), "-;So(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("So [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Sommerfeld number for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta_r(1:test_freq:end, i), "-;beta_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta(1:test_freq:end, i), "-;beta(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("beta [deg]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Angular offset of reaction force for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), mu_r(1:test_freq:end, i), "-;mu_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), mu(1:test_freq:end, i), "-;mu(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("mu [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Coefficient of friction for pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), Q_r(1:test_freq:end, i), "-;Q_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), Q(1:test_freq:end, i), "-;Q(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("Q [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Non-dimensional oil flow pure rotation B/d=%.2f", B_d_r(i)));
%! endfor
%! assert(mean(mean(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.03);
%! assert(mean(mean(abs(beta(1:test_freq:end, 1:test_freq:end) ./ beta_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.02);
%! assert(mean(mean(abs(mu(1:test_freq:end, 1:test_freq:end) ./ mu_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.02);
%! assert(mean(mean(abs(Q(1:test_freq:end, 1:test_freq:end) ./ Q_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.07);
%! assert(max(max(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.05);
%! assert(max(max(abs(beta(1:test_freq:end, 1:test_freq:end) ./ beta_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.04);
%! assert(max(max(abs(mu(1:test_freq:end, 1:test_freq:end) ./ mu_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.05);
%! assert(max(max(abs(Q(1:test_freq:end, 1:test_freq:end) ./ Q_r(1:test_freq:end, 1:test_freq:end) - 1))) < 0.08);

%!test
%! ## References:
%! ## Hans Juergen Butenschoen
%! ## Das hydrodynamische, zylindrische Gleitlager endlicher Breite unter instationaerer Belastung, Karlsruhe 1976
%! ## Test case for squeeze flow
%! 
%! B_d_r = [1, 1/2, 1/2.5, 1/3, 1/4, 1/5, 1/6, 1/8];
%! epsilon_r = [-0.99, -0.98, -0.97, -0.96, -0.95, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0, ...
%!                0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99];

%! So_r = [0.4925  0.1303  0.0842  0.0588  0.0333  0.0214  0.0149  0.0084
%!         0.4986  0.1319  0.0853  0.0596  0.0338  0.0217  0.0151  0.0085
%!         0.5049  0.1337  0.0864  0.0604  0.0342  0.022  0.0153  0.0086
%!         0.5113  0.1354  0.0875  0.0612  0.0347  0.0223  0.0155  0.0087
%!         0.5177  0.1372  0.0887  0.062  0.0351  0.0226  0.0157  0.0085
%!         0.5518  0.1466  0.0948  0.0663  0.0376  0.0241  0.0168  0.0095
%!         0.6293  0.1685  0.1085  0.0761  0.0432  0.0277  0.0193  0.0109
%!         0.7242  0.1944  0.126  0.0881  0.05  0.0321  0.0223  0.0126
%!         0.8342  0.2265  0.1469  0.1028  0.0583  0.0375  0.0261  0.0147
%!         0.971  0.2664  0.173  0.1211  0.0688  0.0442  0.0308  0.0174
%!         1.1391  0.3164  0.2058  0.1443  0.082  0.0527  0.0367  0.0207
%!         1.3494  0.3803  0.2479  0.1739  0.0989  0.0637  0.0444  0.025
%!         1.6134  0.4632  0.3027  0.2127  0.1212  0.078  0.0544  0.0307
%!         1.9579  0.5732  0.3757  0.2645  0.1509  0.0973  0.0678  0.0383
%!         2.4076  0.7224  0.4754  0.3355  0.1919  0.1238  0.0863  0.0488
%!         3.0122  0.9814  0.6159  0.4357  0.2499  0.1614  0.1127  0.0637
%!         3.8485  1.2236  0.8205  0.5826  0.3354  0.2171  0.1517  0.0859
%!         5.0479  1.6904  1.1329  0.8081  0.4675  0.3033  0.2122  0.1203
%!         6.8503  2.4202  1.6392  1.1756  0.6845  0.4455  0.3123  0.1774
%!         9.7319  3.675  2.5206  1.8236  1.0715  0.7005  0.4923  0.2803
%!         14.769  6.0649  4.2362  3.1002  1.8455  1.2147  0.857  0.4899
%!         24.842  11.373  8.156  6.0733  3.6891  2.4543  1.7424  1.0027
%!         50.349  26.676  15.94  15.282  9.6161  6.5227  4.6863  2.7312
%!         160.25  104.75  84.39  68.484  46.559  33.102  24.5  14.77
%!         490.18  371.98  320.59  275.76  205.15  155.23  119.99  76.315
%!         699.59  549.66  492.49  422.2  323.24  250.15  196.59  127.78
%!         1098.6  900.95  806.73  720.35  571.88  455.42  366.1  245.48
%!         2073.5  1774.3  1629.4  1490.8  1239.4  1027.2  853.42  600.66
%!         6057.9  5486.9  5171.8  4884.9  4335.3  3813.7  2331  2560.1];
%!
%! Q_r = [3.7975	1.9601	1.5761	1.3172	0.9917	0.7954	0.6638	0.4988
%!        3.7949	1.9601	1.5761	1.3172	0.9917	0.7954	0.6638	0.4988
%!        3.7927	1.9600	1.5760	1.3172	0.9917	0.7954	0.6638	0.4988
%!        3.7901	1.9600	1.5760	1.3172	0.9916	0.7954	0.6638	0.4988
%!        3.7878	1.9599	1.5759	1.3171	0.9916	0.7953	0.6638	0.4988
%!        3.7749	1.9594	1.5755	1.3171	0.9915	0.7953	0.6638	0.4988
%!        3.7498	1.9577	1.5743	1.3170	0.9914	0.7952	0.6637	0.4988
%!        3.7216	1.9549	1.5723	1.3164	0.9912	0.7951	0.6637	0.4988
%!        3.6937	1.9512	1.5706	1.3156	0.9909	0.7950	0.6636	0.4988
%!        3.6648	1.9470	1.5688	1.3144	0.9906	0.7949	0.6635	0.4988
%!        3.6338	1.9422	1.5664	1.3133	0.9903	0.7948	0.6634	0.4987
%!        3.5994	1.9369	1.5639	1.3119	0.9899	0.7947	0.6634	0.4987
%!        3.5585	1.9308	1.5611	1.3104	0.9895	0.7945	0.6633	0.4987
%!        3.5124	1.9231	1.5581	1.3090	0.9890	0.7943	0.6632	0.4987
%!        3.4640	1.9152	1.5549	1.3072	0.9885	0.7941	0.6632	0.4986
%!        3.4002	1.9063	1.5512	1.3051	0.9879	0.7939	0.6631	0.4986
%!        3.3313	1.8961	1.5476	1.3031	0.9873	0.7937	0.6630	0.4986
%!        3.2523	1.8859	1.5431	1.3010	0.9866	0.7935	0.6629	0.4985
%!        3.1603	1.8738	1.5381	1.2987	0.9860	0.7932	0.6627	0.4985
%!        3.0839	1.8603	1.5328	1.2959	0.9852	0.7929	0.6626	0.4985
%!        2.9292	1.8455	1.5266	1.2927	0.9843	0.7925	0.6625	0.4984
%!        2.8114	1.8257	1.5188	1.2889	0.9832	0.7921	0.6623	0.4984
%!        2.6815	1.8013	1.5086	1.2841	0.9819	0.7917	0.6622	0.4984
%!        2.5363	1.7723	1.4934	1.2770	0.9803	0.7911	0.6620	0.4983
%!        2.4560	1.7543	1.4843	1.2735	0.9793	0.7908	0.6619	0.4983
%!        2.4407	1.7503	1.4822	1.2725	0.9791	0.7907	0.6618	0.4983
%!        2.4244	1.7463	1.4799	1.2715	0.9788	0.7906	0.6618	0.4982
%!        2.4089	1.7416	1.4773	1.2702	0.9785	0.7905	0.6618	0.4982
%!        2.3923	1.7362	1.4746	1.2690	0.9782	0.7904	0.6617	0.4982];
%! So = beta = mu = Q = dQ = nan(numel(epsilon_r), numel(B_d_r));
%! cavitation = "mass conserving";
%! test_freq = 4;
%! verbose = false;
%! for j=1:test_freq:numel(epsilon_r)
%! for k=1:test_freq:numel(B_d_r)
%! param.d1 = 10e-3;
%! param.D2 = 10.01e-3;
%! param.M = int32(20);
%! param.N = int32(85);
%! param.output_bearing_data = true;
%! param.B = param.d1 * B_d_r(k);
%! param.epsilon = epsilon_r(j);
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer hyd_node_id_inlet = 2003;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: integer genel_id_inlet = 4005;\n");
%!     fputs(fd, "set: real omega1z = 2 * pi * 0;\n");
%!     fputs(fd, "set: real omega2z = 2 * pi * 0;\n");
%!     fputs(fd, "set: real v1z = 0;\n");
%!     fputs(fd, "set: real v2z = 0;\n");
%!     fputs(fd, "set: real n = 1e-3;\n");
%!     fputs(fd, "set: integer K = 36;\n");
%!     fputs(fd, "set: real t1 = abs(2 * pi * n / max(1.,abs(omega1z)));\n");
%!     fputs(fd, "set: real dt = t1 / K;\n");
%!     fputs(fd, "set: real eta = 1e-3;\n");
%!     fputs(fd, "set: real rho = 3600. * 1000.; # [l/h]\n\n");
%!     fputs(fd, "set: real epsilon_t0 = (abs(epsilon) * (1 - sign(epsilon) * n)) * (abs(epsilon) > 0) + n * (abs(epsilon) == 0);\n");
%!     fputs(fd, "set: real epsilon_t1 = abs(epsilon);\n");
%!     fputs(fd, "set: real delta_t0 = pi;\n");
%!     fputs(fd, "set: real delta_t1 = pi + (omega1z + omega2z) / 2 * t1;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real p_pB2 = 0;\n");
%!     fputs(fd, "set: real p_mB2 = 0;\n");
%!     fputs(fd, "set: real p_in = 0;\n");
%!     fputs(fd, "set: real pmax = 500000.;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 30;\n");
%!     fputs(fd, "        tolerance: 1e-10, test, minmax;\n");
%!     fputs(fd, "        linear solver: umfpack, map, colamd, scale, row max column max, always, max iterations, 100;\n");
%!     fputs(fd, "        nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e7;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: bdf;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 10;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-3;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "    structural nodes: 2;\n");
%!     fputs(fd, "    joints: 2;\n");
%!     fputs(fd, "    loadable elements: 1;\n");
%!     fputs(fd, "    hydraulic nodes: 3;\n");
%!     fputs(fd, "    genels: 3;\n");
%!     fputs(fd, "    print: dof stats, to file;\n");
%!     fputs(fd, "    print: equation description, to file;\n");
%!     fputs(fd, "    print: dof description, to file;\n");
%!     fputs(fd, "    print: element connection, to file;\n");
%!     fputs(fd, "    print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 16;\n");
%!     fputs(fd, "    default output: reference frames;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.5 * (D2 - d1) * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        0.5 * (D2 - d1) * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        v1z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing, null,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        v2z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_inlet, p_in;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            mult, const, v1z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, v2z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!     fputs(fd, "    genel: genel_id_inlet, clamp,\n");
%!     fputs(fd, "        hyd_node_id_inlet, hydraulic, p_in;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     switch (cavitation)
%!     case "non mass conserving"
%!       fputs(fd, "            hydraulic fluid, incompressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                pressure, 0.,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!     case "mass conserving"
%!       fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                1.,\n");
%!       fputs(fd, "                0.,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!       fputs(fd, "            viscosity vapor, eta,\n");
%!     otherwise
%!       error("unknown cavitation model: %s", cavitation);
%!     endswitch
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at bearing,\n");
%!     fputs(fd, "                bearing width, B,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure coupling conditions radial, 1,\n");
%!     fputs(fd, "                position, 0.5 * D2 * (delta_t0 + pi * (epsilon_dot > 0)), 0.,\n");
%!     fputs(fd, "                rectangle, width, 1.5 * D2 * pi / (N - 1), height, B,\n");
%!     fputs(fd, "                pressure node, hyd_node_id_inlet,\n");
%!     fputs(fd, "            pressure dof scale, pmax,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, yes,\n");
%!     fputs(fd, "            output velocity, yes,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   if (~verbose)
%!     opt_sol.logfile = [output_file, ".stdout"];
%!   endif
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%!   d = res.log_dat.bearings.cylindrical.d;
%!   D = res.log_dat.bearings.cylindrical.D;
%!   B = res.log_dat.bearings.cylindrical.B;
%!   s = (D - d);
%!   Psi = (D - d) / D;
%!   eta = res.log_dat.bearings.eta;
%!   node_idx_2 = find(res.node_id == res.log_dat.bearings.cylindrical.nodes(2).label);
%!   Phi2 = res.trajectory{node_idx_2}(end, 4:6).';
%!   R2 = euler123_to_rotation_matrix(Phi2);
%!   Rb1 = res.log_dat.bearings.cylindrical.nodes(1).Rb;
%!   Rb2 = res.log_dat.bearings.cylindrical.nodes(2).Rb;
%!   F1 = -res.bearings.columns.F1(end, :).';
%!   M1 = -res.bearings.columns.M1(end, :).';
%!   omega1z = res.bearings.cylindrical.omega1z(end, :);
%!   omega2z = res.bearings.cylindrical.omega2z(end, :);
%!   omega_res = res.bearings.cylindrical.omega_res(end);
%!   epsilon = res.bearings.cylindrical.epsilon(end);
%!   delta = res.bearings.cylindrical.delta(end);
%!   e1_dot = res.bearings.cylindrical.e_dot_R2(end, :);
%!   epsilon_dot = res.bearings.cylindrical.epsilon_dot(end);
%!   Q_out1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end) / res.log_dat.vars.rho;
%!   Q_out2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end) / res.log_dat.vars.rho;
%!   Q_in = -res.genel_data{res.log_dat.vars.genel_idx_inlet}(end) / res.log_dat.vars.rho;
%!   P = norm(Rb2(:, 1:2).' * R2.' * F1);
%!   alpha = atan2(Rb2(:, 2).' * R2.' * F1, Rb2(:, 1).' * R2.' * F1);
%!   beta(j, k) = sign(omega_res) * (delta - alpha);
%!   beta(j, k) = mod(beta(j, k) + pi, 2 * pi) - pi;
%!   mu(j, k) = 2 * Rb2(:, 3).' * R2.' * M1 / (P * d) * sign(omega1z - omega2z);
%!   So(j, k) = P * Psi^2 / (B * D * eta * abs(epsilon_dot));
%!   Q(j, k) = (Q_out1 + Q_out2) / ((0.5 * D)^3 * Psi * abs(epsilon_dot));
%!   dQ(j, k) = (Q_out1 + Q_out2 + Q_in) / Q_in;
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
%! endfor
%! endfor
%! for i=1:columns(beta)
%!   beta(beta(:, i) < -0.9 * pi, i) += 2 * pi;
%! endfor
%! mu_r = zeros(numel(epsilon_r), numel(B_d_r));
%! beta_r = zeros(numel(epsilon_r), numel(B_d_r));
%! beta_r(find(epsilon_r < 0), :) = pi;
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   semilogy(epsilon_r(1:test_freq:end), So_r(1:test_freq:end, i), "-;So_r(epsilon);0");
%!   semilogy(epsilon_r(1:test_freq:end), So(1:test_freq:end, i), "-;So(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("So [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Sommerfeld number for pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta_r(1:test_freq:end, i), "-;beta_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), 180 / pi * beta(1:test_freq:end, i), "-;beta(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("beta [deg]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Angular offset of reaction force for pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), mu_r(1:test_freq:end, i), "-;mu_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), mu(1:test_freq:end, i), "-;mu(epsilon);1");
%!   xlabel("epsilon [1]");
%!   ylabel("mu [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Coefficient of friction for pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! for i=1:test_freq:numel(B_d_r)
%!   figure("visible", "off");
%!   hold("on");
%!   plot(epsilon_r(1:test_freq:end), Q_r(1:test_freq:end, i), "-;Q_r(epsilon);0");
%!   plot(epsilon_r(1:test_freq:end), Q(1:test_freq:end, i), "-;Q(epsilon);1");
%!   ylim([0, ylim()(2)]);
%!   xlabel("epsilon [1]");
%!   ylabel("Q [1]");
%!   grid on;
%!   grid minor on;
%!   title(sprintf("Non-dimensional oil flow pure displacement B/d=%.2f", B_d_r(i)));
%! endfor
%! assert(mean(mean(abs(So(1:test_freq:end, 1:test_freq:end) ./ So_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.02);
%! assert(mean(mean(abs(beta(1:test_freq:end,1:test_freq:end) - beta_r(1:test_freq:end, 1:test_freq:end)))) < 1e-6 * pi / 180);
%! assert(mean(mean(abs(mu(1:test_freq:end,1:test_freq:end)))) < 1e-8);
%! assert(mean(mean(abs(Q(1:test_freq:end,1:test_freq:end) ./ Q_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.06);
%! assert(max(max(abs(So(1:test_freq:end,1:test_freq:end) ./ So_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.11);
%! assert(max(max(abs(beta(1:test_freq:end,1:test_freq:end) - beta_r(1:test_freq:end,1:test_freq:end)))) < 1e-6 * pi / 180);
%! assert(max(max(abs(mu(1:test_freq:end,1:test_freq:end)))) < 1e-8);
%! assert(max(max(abs(Q(1:test_freq:end,1:test_freq:end) ./ Q_r(1:test_freq:end,1:test_freq:end) - 1))) < 0.07);

%!test
%! function [p, mdotz] = reference_solution(h0, B, dh_dz, p0, p1, eta, U1z, U2z, rho, dx, dh_dt, z)
%!    p = "(dh_dz*h0^2*p1*z^2+2*B*dh_dz^2*h0*p1*z^2+B^2*dh_dz^3*p1*z^2-dh_dz*h0^2*p0*z^2+6*B*U2z*dh_dz*eta*z^2+6*B*U1z*dh_dz*eta*z^2+12*B*dh_dt*eta*z^2+2*h0^3*p1*z+4*B*dh_dz*h0^2*p1*z+2*B^2*dh_dz^2*h0*p1*z-2*h0^3*p0*z-6*B^2*U2z*dh_dz*eta*z-6*B^2*U1z*dh_dz*eta*z-12*B^2*dh_dt*eta*z+2*B*h0^3*p0+B^2*dh_dz*h0^2*p0)/(B*(2*h0+B*dh_dz)*(dh_dz*z+h0)^2)";
%!    mdotz = "-((12*B*dh_dt*eta*h0+6*B^2*dh_dt*dh_dz*eta)*rho*z+((h0^4+2*B*dh_dz*h0^3+B^2*dh_dz^2*h0^2)*p1+(-h0^4-2*B*dh_dz*h0^3-B^2*dh_dz^2*h0^2)*p0+(-6*B*U2z-6*B*U1z)*eta*h0^2+((-6*B^2*U2z-6*B^2*U1z)*dh_dz-6*B^2*dh_dt)*eta*h0)*rho)/(12*B*eta*h0+6*B^2*dh_dz*eta)";
%!    p = eval(vectorize(p));
%!    mdotz = eval(mdotz) * dx;
%! endfunction
%! param.d1 = 10e-3;
%! param.h = 1e-5;
%! param.M = int32(750);
%! param.N = int32(4);
%! param.output_bearing_data = true;
%! param.B = 25e-3;
%! cavitation = "mass conserving";
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: real omega1z = 0;\n");
%!     fputs(fd, "set: real omega2z = 0;\n");
%!     fputs(fd, "set: real v1z = 15;\n");
%!     fputs(fd, "set: real v2z = -7;\n");
%!     fputs(fd, "set: integer K = 36;\n");
%!     fputs(fd, "set: real t1 = 1e-2;\n");
%!     fputs(fd, "set: real dt = t1 / K;\n");
%!     fputs(fd, "set: real D2 = d1 + 2 * h;\n");
%!     fputs(fd, "set: real dy = 1.5 * h;\n");
%!     fputs(fd, "set: real eta = 1e-3;\n");
%!     fputs(fd, "set: real rho = 3600. * 1000.; # [l/h]\n\n");
%!     fputs(fd, "set: real epsilon_t0 = 0;\n");
%!     fputs(fd, "set: real epsilon_t1 = 0;\n");
%!     fputs(fd, "set: real delta_t0 = 0;\n");
%!     fputs(fd, "set: real delta_t1 = 0;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real p_pB2 = 10e5;\n");
%!     fputs(fd, "set: real p_mB2 = 0.8e5;\n");
%!     fputs(fd, "set: real pmax = 500000.;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 30;\n");
%!     fputs(fd, "        tolerance: 1e-8, test, minmax;\n");
%!     fputs(fd, "        linear solver: umfpack, map, colamd, scale, row max column max, always, max iterations, 100;\n");
%!     fputs(fd, "        nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e5;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: bdf;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 10;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-6;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "        structural nodes: 2;\n");
%!     fputs(fd, "        joints: 2;\n");
%!     fputs(fd, "        loadable elements: 1;\n");
%!     fputs(fd, "        hydraulic nodes: 2;\n");
%!     fputs(fd, "        genels: 2;\n");
%!     fputs(fd, "        print: dof stats, to file;\n");
%!     fputs(fd, "        print: equation description, to file;\n");
%!     fputs(fd, "        print: dof description, to file;\n");
%!     fputs(fd, "        print: element connection, to file;\n");
%!     fputs(fd, "        print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 8;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        -v1z * t1,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        v1z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        -v2z * t1,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        v2z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            array, 2, -v1z * t1,\n");
%!     fputs(fd, "                      mult, const, v1z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            array, 2, -v2z * t1, mult, const, v2z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     switch (cavitation)
%!     case "non mass conserving"
%!       fputs(fd, "            hydraulic fluid, incompressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                pressure, 0.,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!     case "mass conserving"
%!       fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!       fputs(fd, "                density, rho,\n");
%!       fputs(fd, "                1.,\n");
%!       fputs(fd, "                0.,\n");
%!       fputs(fd, "                viscosity, eta,\n");
%!       fputs(fd, "                temperature, 0,\n");
%!       fputs(fd, "            viscosity vapor, eta,\n");
%!     otherwise
%!       error("unknown cavitation model: %s", cavitation);
%!     endswitch
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at bearing,\n");
%!     fputs(fd, "                bearing width, B,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                pockets,\n");
%!     fputs(fd, "                  shaft, 1,\n");
%!     fputs(fd, "                    position, d1 * pi * 0.5, 0.,\n");
%!     fputs(fd, "                    complete surface,\n");
%!     fputs(fd, "                    pocket height,\n");
%!     fputs(fd, "                    linear,\n");
%!     fputs(fd, "                      x, 0., d1 * pi,\n");
%!     fputs(fd, "                      z, -0.5 * B, 0.5 * B,\n");
%!     fputs(fd, "                      delta y, 0., -dy,\n");
%!     fputs(fd, "                               0., -dy,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure dof scale, pmax,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, no,\n");
%!     fputs(fd, "            output velocity, no,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   opt_sol.logfile = [output_file, ".stdout"];
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%!   mdot1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end);
%!   mdot2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end);
%!   d = res.log_dat.bearings.cylindrical.d;
%!   D = res.log_dat.bearings.cylindrical.D;
%!   B = res.log_dat.bearings.cylindrical.B;
%!   h0 = 0.5 * (D - d);
%!   dh_dz = res.log_dat.vars.dy / res.log_dat.vars.B;
%!   p0 = res.log_dat.vars.p_mB2;
%!   p1 = res.log_dat.vars.p_pB2;
%!   eta = res.log_dat.vars.eta;
%!   U1z = res.log_dat.vars.v1z - res.log_dat.vars.v2z;
%!   U2z = 0;
%!   rho = res.log_dat.vars.rho;
%!   dx = res.log_dat.bearings.cylindrical.dm * pi;
%!   dh_dt = -dh_dz * (U1z - U2z);
%!   z = res.bearings.zi(:, 1) + 0.5 * B;
%!   [p_ref, mdotz_ref] = reference_solution(h0, B, dh_dz, p0, p1, eta, U1z, U2z, rho, dx, dh_dt, z);
%!   figure("visible", "off");
%!   hold on;
%!   plot(res.bearings.zi(:, 1) + 0.5 * B, res.bearings.columns.p(:, 1), "-;p(z);1");
%!   plot(z, p_ref, "-;p_r;0");
%!   xlabel("z [m]");
%!   ylabel("p [Pa]");
%!   grid on;
%!   grid minor on;
%!   title("axial pressure distribution");
%!   figure("visible", "off");
%!   hold on;
%!   plot(res.bearings.zi(:, 1) + 0.5 * B, res.bearings.columns.h(:, 1), "-;h(z);1");
%!   plot(z, h0 + dh_dz * z, "-;h_r(z);0");
%!   xlabel("z [m]");
%!   ylabel("h [m]");
%!   grid on;
%!   grid minor on;
%!   title("radial clearance versus time");
%!   assert(res.bearings.columns.p(:, 1), p_ref, 1e-4 * max(abs(p_ref)));
%!   assert(-mdot1, mdotz_ref(1), 0.5e-2 * abs(mdotz_ref(1)));
%!   assert(mdot2, mdotz_ref(end), 0.5e-2 * abs(mdotz_ref(end)));
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! close all;
%! ## References:
%! ## Hans Juergen Butenschoen
%! ## Das hydrodynamische, zylindrische Gleitlager endlicher Breite unter instationaerer Belastung, Karlsruhe 1976
%!
%! B_d_r = [1, 1/2, 1/2.5, 1/3, 1/4, 1/5, 1/6, 1/8];
%!
%! epsilon_r = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.96, 0.97, 0.98, 0.99, 0.995, 0.999];
%!
%! So_r =  [0.1196 0.0368  0.0243  0.0171  0.0098  0.0063  0.0044  0.0025
%!          0.2518 0.0783  0.0518  0.0366  0.0210  0.0136  0.0095  0.0054
%!          0.4091 0.1304  0.0867  0.0615  0.0354  0.0229  0.0160  0.0091
%!          0.6108 0.2026  0.1357  0.0968  0.0560  0.0364  0.0254  0.0144
%!          0.8903 0.3124  0.2117  0.1522  0.0888  0.0579  0.0406  0.0231
%!          1.3146 0.4982  0.3435  0.2496  0.1476  0.0969  0.0682  0.0390
%!          2.0432 0.8595  0.6075  0.4492  0.2708  0.1797  0.1274  0.0732
%!          3.5663 1.7339  1.2756  0.9687  0.6043  0.4085  0.2930  0.1706
%!          8.4352 5.0881  4.0187  3.2201  2.1595  1.5261  1.1263  0.6776
%!          18.7895 13.3083 11.2225 9.4993  6.9248  5.1833  3.9831  2.5202
%!          24.1172 17.7934 15.2823 13.1525 9.8517  7.5257  5.8712  3.7907
%!          33.1297 25.5920 22.4503 19.7034 15.2697 11.9804 9.5425  6.3397
%!          51.4774 41.9230 37.7412 33.9523 27.5040 22.3556 18.3874 12.7695
%!          107.7868 93.7881 87.2906 81.1597 70.0359 60.3874 52.1425 39.2568
%!          223.8850 203.2450 193.3490 183.8040 166.1540 149.8690 134.8910 109.6090
%!          1174.5400 1124.6200 1102.0700 1078.3100 1032.8700 989.1500 945.6700 864.7400];
%!
%! beta_r = [79.410 81.767  82.119  82.299  82.481  82.561  82.608  82.653
%!           73.858 75.142  75.283  75.344  75.287  75.406  75.414  75.422
%!           68.264 68.493  68.423  68.368  68.204  68.261  68.238  68.211
%!           62.571 61.778  61.544  61.382  61.208  61.115  61.062  61.007
%!           56.705 54.993  54.603  54.348  54.069  53.922  53.863  53.784
%!           50.536 48.049  47.521  47.191  46.825  46.647  46.554  46.449
%!           43.859 40.803  40.156  39.756  39.326  39.108  38.983  38.865
%!           36.235 32.938  32.216  31.761  31.249  30.988  30.840  30.692
%!           26.482 23.566  22.849  22.368  21.790  21.476  21.289  21.089
%!           19.450 17.265  16.648  16.207  15.632  15.291  15.075  14.832
%!           17.609 15.660  15.089  14.669  14.109  13.768  13.547  13.262
%!           15.484 13.818  13.313  12.929  12.299  12.062  11.838  11.570
%!           12.903 11.598  11.178  10.850  10.375  10.057  9.835   9.558
%!            9.416  8.587   8.301   8.066   7.703   7.440   7.242   6.975
%!            6.829  6.325   6.143   5.987   5.733   5.526   5.380   5.151
%!            3.196  3.048   2.989   2.940   2.848   2.769   2.708   2.599] * pi / 180;
%!
%! Q_r = [0.1603  0.0939  0.0768  0.0648  0.0492  0.0396  0.0331  0.0249
%!        0.3196  0.1878  0.1536  0.1296  0.0984  0.0792  0.0662  0.0498
%!        0.4765  0.2816  0.2304  0.1944  0.1476  0.1188  0.0993  0.0747
%!        0.6318  0.3755  0.3072  0.2592  0.1968  0.1583  0.1324  0.0996
%!        0.7852  0.4694  0.384   0.324   0.246   0.198   0.1655  0.1245
%!        0.9374  0.5634  0.461   0.3889  0.2953  0.2376  0.1986  0.1494
%!        1.0888  0.6578  0.5383  0.454   0.3446  0.2772  0.2317  0.1743
%!        1.2392  0.7529  0.6159  0.5193  0.394   0.3169  0.2648  0.1992
%!        1.3898  0.8453  0.6944  0.5852  0.4436  0.3567  0.298   0.2241
%!        1.4655  0.8984  0.7343  0.6186  0.4687  0.3767  0.3147  0.2366
%!        1.4808  0.9084  0.7425  0.6254  0.4737  0.3807  0.3181  0.2391
%!        1.4968  0.9185  0.7507  0.6322  0.4788  0.3848  0.3214  0.2417
%!        1.5124  0.9287  0.7585  0.6391  0.484   0.3889  0.3248  0.2442
%!        1.5277  0.9351  0.7674  0.6461  0.4892  0.393   0.3282  0.2467
%!        1.5354  0.9441  0.7716  0.6496  0.4918  0.3951  0.33    0.248
%!        1.5409  0.9482  0.7745  0.6525  0.494   0.3968  0.3314  0.2491];
%!
%! param.d1 = 10e-3;
%! param.Psi = 1e-3;
%! param.M = int32(20);
%! param.N = int32(75);
%! param.output_bearing_data = true;
%! param.B = 7e-3;
%! param.epsilon = 0.6;
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "mbdyn_post_ehd_load_output_XXXXXX"));
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, "set: integer node_id_rotor = 1001;\n");
%!     fputs(fd, "set: integer node_id_stator = 1002;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet1 = 2001;\n");
%!     fputs(fd, "set: integer hyd_node_id_outlet2 = 2002;\n");
%!     fputs(fd, "set: integer hyd_node_id_inlet = 2003;\n");
%!     fputs(fd, "set: integer joint_id_rotor = 3001;\n");
%!     fputs(fd, "set: integer joint_id_stator = 3002;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 4001;\n");
%!     fputs(fd, "set: integer ref_id_bearing = 5000;\n");
%!     fputs(fd, "set: integer ref_id_rotor = 5001;\n");
%!     fputs(fd, "set: integer ref_id_stator = 5002;\n");
%!     fputs(fd, "set: integer genel_id_outlet1 = 4003;\n");
%!     fputs(fd, "set: integer genel_id_outlet2 = 4004;\n");
%!     fputs(fd, "set: integer genel_id_inlet = 4005;\n");
%!     fputs(fd, "set: real omega1z = 2 * pi * 50;\n");
%!     fputs(fd, "set: real omega2z = 2 * pi * 15;\n");
%!     fputs(fd, "set: real v1z = 0;\n");
%!     fputs(fd, "set: real v2z = 0;\n");
%!     fputs(fd, "set: real n = 3;\n");
%!     fputs(fd, "set: integer K = 36;\n");
%!     fputs(fd, "set: real t1 = abs(2 * pi * n / omega1z);\n");
%!     fputs(fd, "set: real dt = t1 / (n * K);\n");
%!     fputs(fd, "set: real D2 = d1 / (1. - Psi);\n");
%!     fputs(fd, "set: real eta = 1e-3;\n");
%!     fputs(fd, "set: real rho = 3600. * 1000.; # [l/h]\n\n");
%!     fputs(fd, "set: real epsilon_t0 = epsilon;\n");
%!     fputs(fd, "set: real epsilon_t1 = epsilon;\n");
%!     fputs(fd, "set: real delta_t0 = pi;\n");
%!     fputs(fd, "set: real delta_t1 = pi + omega2z * t1;\n");
%!     fputs(fd, "set: real gamma_t0 = 0;\n");
%!     fputs(fd, "set: real gamma_t1 = 0;\n");
%!     fputs(fd, "set: real p_pB2 = 0;\n");
%!     fputs(fd, "set: real p_mB2 = 0;\n");
%!     fputs(fd, "set: real p_in = 0;\n");
%!     fputs(fd, "set: real pmax = 500000.;\n");
%!     fputs(fd, "set: real Phi1x_t0 = atan(gamma_t0 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real Phi1x_t1 = atan(gamma_t1 * (D2 - d1) / B);\n");
%!     fputs(fd, "set: real omega1x = (Phi1x_t1 - Phi1x_t0) / t1;\n");
%!     fputs(fd, "set: real epsilon_dot = (epsilon_t1 - epsilon_t0) / t1;\n");
%!     fputs(fd, "set: real delta_dot = (delta_t1 - delta_t0) / t1;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value; # the default\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        max iterations: 30;\n");
%!     fputs(fd, "        tolerance: 1e-10, test, minmax;\n");
%!     fputs(fd, "        linear solver: umfpack, map, colamd, scale, row max column max, always, max iterations, 100;\n");
%!     fputs(fd, "        nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e5;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: bdf;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 10;\n");
%!     fputs(fd, "        derivatives coefficient: 110;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    skip initial joint assembly;\n");
%!     fputs(fd, "    output meter: closest next, t1 - 0.5 * dt, forever, const, dt;\n");
%!     fputs(fd, "    use: loadable elements, in assembly;\n");
%!     fputs(fd, "    default orientation: euler123;\n");
%!     fputs(fd, "        structural nodes: 2;\n");
%!     fputs(fd, "        joints: 2;\n");
%!     fputs(fd, "        loadable elements: 1;\n");
%!     fputs(fd, "        hydraulic nodes: 3;\n");
%!     fputs(fd, "        genels: 3;\n");
%!     fputs(fd, "        print: dof stats, to file;\n");
%!     fputs(fd, "        print: equation description, to file;\n");
%!     fputs(fd, "        print: dof description, to file;\n");
%!     fputs(fd, "        print: element connection, to file;\n");
%!     fputs(fd, "        print: node connection, to file;\n");
%!     fputs(fd, "    output precision: 8;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_bearing,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, eye,\n");
%!     fputs(fd, "        reference, global, null,\n");
%!     fputs(fd, "        reference, global, null;");
%!     fputs(fd, "reference: ref_id_rotor,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * cos(delta_t0),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * epsilon_t0 * sin(delta_t0),\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing, euler123,\n");
%!     fputs(fd, "        Phi1x_t0,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * cos(delta_t0) - epsilon_t0 * sin(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        (D2 - d1) / 2 * (epsilon_dot * sin(delta_t0) + epsilon_t0 * cos(delta_t0) * delta_dot),\n");
%!     fputs(fd, "        v1z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        omega1x,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega1z;\n");
%!     fputs(fd, "reference: ref_id_stator,\n");
%!     fputs(fd, "        reference, ref_id_bearing, null,\n");
%!     fputs(fd, "        reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        v2z,\n");
%!     fputs(fd, "        reference, ref_id_bearing,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        0.,\n");
%!     fputs(fd, "        omega2z;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_rotor, static,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                reference, ref_id_rotor, null;\n");
%!     fputs(fd, "        structural: node_id_stator, static,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                reference, ref_id_stator, null,\n");
%!     fputs(fd, "                reference, ref_id_stator, null;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet1, p_mB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_outlet2, p_pB2;\n");
%!     fputs(fd, "    hydraulic: hyd_node_id_inlet, p_in;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_rotor, total pin joint,\n");
%!     fputs(fd, "                node_id_rotor,\n");
%!     fputs(fd, "                        position, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * cos(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            string, \"(D2 - d1) / 2 * (epsilon_t0 + epsilon_dot * Time) * sin(delta_t0 + delta_dot * Time)\",\n");
%!     fputs(fd, "            mult, const, v1z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            array, 2,\n");
%!     fputs(fd, "                const, Phi1x_t0,\n");
%!     fputs(fd, "                mult, const, omega1x, time,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega1z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "        joint: joint_id_stator, total pin joint,\n");
%!     fputs(fd, "                node_id_stator,\n");
%!     fputs(fd, "                        position, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                        position orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                        rotation orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_bearing, eye,\n");
%!     fputs(fd, "        position constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, v2z,\n");
%!     fputs(fd, "                  time,\n");
%!     fputs(fd, "        orientation constraint,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "            active,\n");
%!     fputs(fd, "        component,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            const, 0.,\n");
%!     fputs(fd, "            mult, const, omega2z,\n");
%!     fputs(fd, "                  time;\n");
%!     fputs(fd, "    genel: genel_id_outlet1, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet1, hydraulic, p_mB2;\n");
%!     fputs(fd, "    genel: genel_id_outlet2, clamp,\n");
%!     fputs(fd, "        hyd_node_id_outlet2, hydraulic, p_pB2;\n");
%!     fputs(fd, "    genel: genel_id_inlet, clamp,\n");
%!     fputs(fd, "        hyd_node_id_inlet, hydraulic, p_in;\n");
%!     fputs(fd, "    user defined: elem_id_bearing,\n");
%!     fputs(fd, "        hydrodynamic plain bearing2,\n");
%!     fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!     fputs(fd, "                density, rho,\n");
%!     fputs(fd, "                1.,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                viscosity, eta,\n");
%!     fputs(fd, "                temperature, 0,\n");
%!     fputs(fd, "            viscosity vapor, eta,\n");
%!     fputs(fd, "            mesh, linear finite difference,\n");
%!     fputs(fd, "            geometry, cylindrical,\n");
%!     fputs(fd, "                mesh position, at bearing,\n");
%!     fputs(fd, "                bearing width, B,\n");
%!     fputs(fd, "                shaft diameter, d1,\n");
%!     fputs(fd, "                bearing diameter, D2,\n");
%!     fputs(fd, "                shaft node,\n");
%!     fputs(fd, "                    node_id_rotor,\n");
%!     fputs(fd, "                    offset, reference, ref_id_rotor, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_rotor, eye,\n");
%!     fputs(fd, "                bearing node,\n");
%!     fputs(fd, "                    node_id_stator,\n");
%!     fputs(fd, "                    offset, reference, ref_id_stator, null,\n");
%!     fputs(fd, "                    orientation, reference, ref_id_stator, eye,\n");
%!     fputs(fd, "            number of nodes z, M,\n");
%!     fputs(fd, "            number of nodes Phi, N,\n");
%!     fputs(fd, "            pressure coupling conditions axial,\n");
%!     fputs(fd, "                hyd_node_id_outlet1,\n");
%!     fputs(fd, "                hyd_node_id_outlet2,\n");
%!     fputs(fd, "            pressure coupling conditions radial, 1,\n");
%!     fputs(fd, "                position, 0.5 * D2 * (delta_t0 + pi), 0.,\n");
%!     fputs(fd, "                rectangle, width, D2 * pi / (N - 1), height, B,\n");
%!     fputs(fd, "                pressure node, hyd_node_id_inlet,\n");
%!     fputs(fd, "            pressure dof scale, pmax,\n");
%!     fputs(fd, "            output pressure, yes,\n");
%!     fputs(fd, "            output density, yes,\n");
%!     fputs(fd, "            output clearance, yes,\n");
%!     fputs(fd, "            output clearance derivative, no,\n");
%!     fputs(fd, "            output velocity, no,\n");
%!     fputs(fd, "            output stress, no,\n");
%!     fputs(fd, "            output reaction force, yes,\n");
%!     fputs(fd, "            output friction loss, yes,\n");
%!     fputs(fd, "            output mesh, yes,\n");
%!     fputs(fd, "            output, output_bearing_data;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!       fd = -1;
%!     endif
%!   end_unwind_protect
%!   opt_sol.output_file = output_file;
%!   opt_sol.verbose = false;
%!   opt_sol.logfile = [output_file, ".stdout"];
%!   info = mbdyn_solver_run(output_file, opt_sol);
%!   res.log_dat = mbdyn_post_load_log(opt_sol.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id] = mbdyn_post_load_output_struct(opt_sol.output_file);
%!   [res.genel_id, res.genel_data] = mbdyn_post_load_output([opt_sol.output_file, ".gen"], 1);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   opt_load.verbose = false;
%!   opt_load.num_steps = numel(res.t);
%!   opt_load.output_index = 1:numel(res.t);
%!   opt_load.loaded_fields = {};
%!   opt_load.interpolate_mesh = true;
%!   res.bearings = mbdyn_post_ehd_load_output(opt_sol.output_file, res.log_dat, opt_load);
%!   d = res.log_dat.bearings.cylindrical.d;
%!   D = res.log_dat.bearings.cylindrical.D;
%!   B = res.log_dat.bearings.cylindrical.B;
%!   s = (D - d);
%!   Psi = (D - d) / D;
%!   eta = res.log_dat.bearings.eta;
%!   node_idx_2 = find(res.node_id == res.log_dat.bearings.cylindrical.nodes(2).label);
%!   Phi2 = res.trajectory{node_idx_2}(:, 4:6).';
%!   R2 = euler123_to_rotation_matrix(Phi2);
%!   Rb1 = res.log_dat.bearings.cylindrical.nodes(1).Rb;
%!   Rb2 = res.log_dat.bearings.cylindrical.nodes(2).Rb;
%!   F1 = -res.bearings.columns.F1.';
%!   M1 = -res.bearings.columns.M1.';
%!   omega1z = res.bearings.cylindrical.omega1z;
%!   omega2z = res.bearings.cylindrical.omega2z;
%!   omega_res = res.bearings.cylindrical.omega_res;
%!   epsilon = res.bearings.cylindrical.epsilon;
%!   delta = res.bearings.cylindrical.delta;
%!   Q_out1 = -res.genel_data{res.log_dat.vars.genel_idx_outlet1}(end) / res.log_dat.vars.rho;
%!   Q_out2 = -res.genel_data{res.log_dat.vars.genel_idx_outlet2}(end) / res.log_dat.vars.rho;
%!   Q_in = -res.genel_data{res.log_dat.vars.genel_idx_inlet}(end) / res.log_dat.vars.rho;
%!   P = norm(Rb2(:, 1:2).' * R2.' * F1);
%!   alpha = atan2(Rb2(:, 2).' * R2.' * F1, Rb2(:, 1).' * R2.' * F1);
%!   beta = sign(omega_res) * (delta - alpha);
%!   beta = mod(beta + pi, 2 * pi) - pi;
%!   mu = 2 * Rb2(:, 3).' * R2.' * M1 / (P * d) * sign(omega1z - omega2z);
%!   So = P * Psi^2 / (B * D * eta * abs(omega_res));
%!   Q = (Q_out1 + Q_out2) / ((0.5 * D)^3 * Psi * abs(omega_res));
%!   dQ = (Q_out1 + Q_out2 + Q_in) / Q_in;
%!   So_ref = interp2(B_d_r, epsilon_r, So_r, B / d, epsilon, "linear");
%!   beta_ref = interp2(B_d_r, epsilon_r, beta_r, B / d, epsilon, "linear");
%!   mu_ref = Psi * (abs((omega1z - omega2z) / omega_res) * pi / (sqrt(1 - epsilon^2) * So_ref) + sin(beta_ref) * abs(epsilon) / 2);
%!   Q_ref = interp2(B_d_r, epsilon_r, Q_r, B / d, epsilon, "linear");
%!   figure("visible", "off");
%!   set(plot(res.bearings.xi(floor(end/2), :) / (0.5 * D) * 180 / pi, res.bearings.columns.p(floor(end/2), :)), "linewidth", 3);
%!   grid on;
%!   grid minor on;
%!   xticks([0:30:360]);
%!   xlabel("Phi [deg]");
%!   ylabel("p [Pa]");
%!   title("midplane pressure distribution");
%!   figure("visible", "off");
%!   set(plot(res.bearings.xi(floor(end/2), :) / (0.5 * D) * 180 / pi, res.bearings.columns.h(floor(end/2), :) / (0.5 * (D - d))), "linewidth", 3);
%!   grid on;
%!   grid minor on;
%!   xlabel("Phi [deg]");
%!   ylabel("h/h0 [1]");
%!   title("relative radial clearance");
%!   figure("visible", "off");
%!   set(plot(res.bearings.xi(floor(end/2), :) / (0.5 * D) * 180 / pi, res.bearings.columns.rho(floor(end/2), :) / res.log_dat.vars.rho), "linewidth", 3);
%!   grid on;
%!   grid minor on;
%!   xlabel("Phi [deg]");
%!   ylabel("rho/rho0 [1]");
%!   title("volumetric filling ratio");
%!   figure_list();
%!   fprintf(stderr, "So / So_ref - 1 = %.2f\n", So / So_ref - 1);
%!   fprintf(stderr, "beta / beta_ref - 1 = %.2f\n", beta / beta_ref - 1);
%!   fprintf(stderr, "mu / mu_ref - 1 = %.2f\n", mu / mu_ref - 1);
%!   fprintf(stderr, "Q / Q_ref - 1 = %.2f\n", Q / Q_ref - 1);
%!   assert(So, So_ref, 0.03 * So_ref);
%!   assert(beta, beta_ref, 0.02 * beta_ref);
%!   assert(mu, mu_ref, 0.03 * mu_ref);
%!   assert(Q, Q_ref, 0.07 * Q_ref);
%!   assert(dQ, 0, 1e-3);
%! unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       fn_i = fullfile(fn(i).folder, fn(i).name);
%!       status = unlink(fn_i);
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn_i);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
