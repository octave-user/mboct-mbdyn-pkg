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

%!test
%! f_print_input_file = false;
%! F = 10;
%! d = 1e-3;
%! A = d^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = d^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! D = 10e-3;
%! L = 1e-3;
%! n = 1;
%! N = n * 7;
%! Phi = linspace(0, 2 * pi * n, n * 36);
%! X = [0.5 * D * cos(Phi);
%!      0.5 * D * sin(Phi);
%!      linspace(0, L, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 10);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_deformation_scale_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "joint: 2, total pin joint,\n");
%!     fprintf(fd, " %d,\n", columns(beam.Xn));
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fputs(fd, "     position constraint, active, active, inactive, null,\n");
%!     fputs(fd, "     orientation constraint, active, active, active, null;\n");
%!     fprintf(fd, "force: 1, absolute, %d, position, reference, %d, null, 0., 0., -1., mult, time, F;\n", columns(beam.Xn), columns(beam.Xn));
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   if (f_print_input_file)
%!     spawn_wait(spawn("cat", {fname}));
%!   endif
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   for k=1:2
%!     node_groups.X0_node_id = log_dat.nodes(1).label;
%!     switch (k)
%!     case 1
%!       node_groups.R0ref_node_id = log_dat.nodes(1).label;
%!     case 2
%!       node_groups.R0ref_node_id = [log_dat.nodes([1, ceil(end / 2), end]).label];
%!     endswitch
%!     node_groups.scale = 100;
%!     node_groups.node_id_scale = [log_dat.nodes.label];
%!     mbdyn_post_deformations_scale([fname, ".log"], [fname, ".mov"], [fname, "_scaled.mov"], node_groups);
%!     copyfile([fname, ".log"], [fname, "_scaled.log"]);
%!     copyfile([fname, ".out"], [fname, "_scaled.out"]);
%!     [scaled{k}.t, scaled{k}.trajectory, scaled{k}.deformation] = mbdyn_post_load_output_struct([fname, "_scaled"]);
%!   endfor
%!   R = G * d^4 / (8 * n * D^3);
%!   wref = -F / R;
%!   tol = 1e-2;
%!   assert(deformation{end}(end, 3), wref, tol * abs(wref));
%!   for k=1:numel(scaled)
%!     assert(scaled{k}.deformation{end}(end, 3), node_groups.scale * wref, tol * abs(node_groups.scale * wref));
%!   endfor
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! f_print_input_file = false;
%! F = 10;
%! d = 1e-3;
%! A = d^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = d^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! D = 10e-3;
%! L = 1e-3;
%! n = 1;
%! N = n * 7;
%! Phi = linspace(0, 2 * pi * n, n * 36);
%! X = [0.5 * D * cos(Phi);
%!      0.5 * D * sin(Phi);
%!      linspace(0, L, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 10);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_deformation_scale_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 2;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, "joint: 2, total pin joint,\n");
%!     fprintf(fd, " %d,\n", columns(beam.Xn));
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fprintf(fd, "   position, reference, %d, null,\n", columns(beam.Xn));
%!     fputs(fd, "     position orientation, reference, global, eye,\n");
%!     fputs(fd, "     rotation orientation, reference, global, eye,\n");
%!     fputs(fd, "     position constraint, active, active, inactive, null,\n");
%!     fputs(fd, "     orientation constraint, active, active, active, null;\n");
%!     fprintf(fd, "force: 1, absolute, %d, position, reference, %d, null, 0., 0., -1., mult, time, F;\n", columns(beam.Xn), columns(beam.Xn));
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   if (f_print_input_file)
%!     spawn_wait(spawn("cat", {fname}));
%!   endif
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   for k=1:2
%!     node_groups.X0_node_id = log_dat.nodes(1).label;
%!     switch (k)
%!     case 1
%!       node_groups.R0ref_node_id = log_dat.nodes(1).label;
%!     case 2
%!       node_groups.R0ref_node_id = [log_dat.nodes([1, ceil(end / 2), end]).label];
%!     endswitch
%!     node_groups.scale = 100;
%!     node_groups.node_id_scale = [log_dat.nodes.label];
%!     mbdyn_post_deformations_scale([fname, ".log"], [fname, ".mov"], [fname, "_scaled.mov"], node_groups);
%!     copyfile([fname, ".log"], [fname, "_scaled.log"]);
%!     copyfile([fname, ".out"], [fname, "_scaled.out"]);
%!     [scaled{k}.t, scaled{k}.trajectory, scaled{k}.deformation] = mbdyn_post_load_output_struct([fname, "_scaled"]);
%!   endfor
%!   R = G * d^4 / (8 * n * D^3);
%!   wref = -F / R;
%!   tol = 1e-2;
%!   assert(deformation{end}(end, 3), wref, tol * abs(wref));
%!   for k=1:numel(scaled)
%!     assert(scaled{k}.deformation{end}(end, 3), node_groups.scale * wref, tol * abs(node_groups.scale * wref));
%!   endfor
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
