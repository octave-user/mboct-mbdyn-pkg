## Copyright (C) 2016 Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} mbdyn_post_eig_to_mov_file(@var{input_file}, @var{output_filename}, @var{options}, @var{modal})
##
## Converts the results from an direct eigenanalysis generated by MBDyn to a MBDyn compatible .mov file.
##
## @var{input_file} @dots{} MBDyn input file without extension .mbdyn
##
## @var{output_filename} @dots{} Name of the file where the output is written to
##
## @var{options}.@var{frames} @dots{} The number of time steps in the output file (default 100)
##
## @var{options}.@var{scale} @dots{} The scale factor of the eigenvector @var{modal}.@var{VR} (default 1)
##
## @var{options}.@var{mode_index} @dots{} The index of the eigenmode to be converted (default 1; limit 1...columns(@var{modal}.@var{VR}))
##
## @var{modal} @dots{} A data structure containing the output of an MBDyn eigenanalysis.
##
## @end deftypefn

function mode_index = mbdyn_post_eig_to_mov_file(input_file, output_filename_template, options, modal)
  if (nargin < 2)
    print_usage();
  endif

  output_filename_template = mbdyn_convert_path(output_filename_template);

  if (nargin < 3)
    options = struct();
  endif

  if (~isfield(options, "frames"))
    options.frames = 100;
  endif

  if (~isfield(options, "scale"))
    options.scale = 1;
  endif

  if (~isfield(options, "periods"))
    options.periods = 1;
  endif

  if (~isfield(options, "mode_index"))
    options.mode_index = 1:columns(modal.VR);
  endif

  if (~isfield(options, "tolerance_lambda"))
    options.tolerance_lambda = sqrt(eps);
  endif

  if (~isfield(options, "verbose"))
    options.verbose = false;
  endif

  log_file = canonicalize_file_name(mbdyn_post_output_filename(input_file, ".log"));

  [nodes, dof_info ] = mbdyn_post_load_log_node(input_file);

  node_labels = [nodes.label];

  if (nargin < 4)
    modal = mbdyn_post_load_output_eig(input_file);
  endif

  node_idx_modal = zeros(numel(nodes), 1, "int32");

  for i=1:numel(nodes)
    node_idx_modal_i = find(nodes(i).label == modal.labels);

    if (isempty(node_idx_modal_i))
      node_idx_modal(i) = -1;
    else
      node_idx_modal(i) = node_idx_modal_i;
    endif
  endfor

  dof_index = repmat(modal.idx, 1, 6);

  for i=1:6
    dof_index((i - 1) * length(modal.idx) + (1:length(modal.idx))) += i;
  endfor

  mode_index = zeros(columns(modal.VR), 1, "int32");
  mode_num = int32(0);

  for k=1:length(options.mode_index)
    if (options.mode_index(k) < 1 || options.mode_index(k) > columns(modal.VR))
      error("mode_index=%d exceeds limits [%d-%d]", options.mode_index(k), 1, columns(modal.VR));
    endif

    norm_VR = max(abs(modal.VR(:, options.mode_index(k))));
    norm_VR_struct = max(abs(modal.VR(dof_index, options.mode_index(k))));

    if (norm_VR_struct >= options.tolerance_lambda * norm_VR)
      mode_index(++mode_num) = options.mode_index(k);
    endif
  endfor

  if (mode_num == 0)
    return;
  endif

  mode_index = mode_index(1:mode_num);

  for k=1:length(mode_index)
    output_filename = sprintf(output_filename_template, mode_index(k));

    [out_dir, out_name, out_ext] = fileparts(output_filename);

    out_filename = fullfile(out_dir, [out_name, ".out"]);

    fd_mov = -1;
    fd_out = -1;

    unwind_protect
      [fd_mov, msg] = fopen(output_filename, "wt");

      if (fd_mov == -1)
        error("could not open file \"%s\": %s", output_filename, msg);
      endif

      [fd_out, msg] = fopen(out_filename, "wt");

      if (fd_out == -1)
        error("failed to open file \"%s\": %s", out_filename, msg);
      endif

      fprintf(fd_out, "# Derivatives solution step at time 0 performed in 0 iterations with 0 error\n");
      fprintf(fd_out, "# Key for lines starting with \"Step\":\n");
      fprintf(fd_out, "# Step Time TStep NIter ResErr SolErr SolConv Out\n");

      if (options.verbose)
        fprintf(stderr, "mode_index=%d file=\"%s\"\n", mode_index(k), output_filename);
      endif

      norm_VR = max(abs(modal.VR(dof_index(1:3 * length(modal.idx)), mode_index(k))));

      lambda = modal.lambda(mode_index(k));

      phase = 2 * pi * ((1:options.frames) - 1) * options.periods / (options.frames - 1);

      for i=1:numel(phase)
        omega = 2 * pi * modal.f(mode_index(k));

        fprintf(fd_out, "Step %d %e %e %d %e %e %d %d\n", i - 1, phase(i) / omega, (phase(2) - phase(1)) / omega, 0, 0, 0, 0, 1);

        for j=1:numel(nodes)
          fprintf(fd_mov, "%d ", nodes(j).label);

          if (node_idx_modal(j) > 0)
            if (~all(modal.labels(node_idx_modal(j)) == nodes(j).label))
              error("labels do not match");
            endif

            idx = modal.idx(node_idx_modal(j));

            VR = modal.VR(idx + (1:6), mode_index(k)) / norm_VR;

            dX = options.scale * real(VR * exp(1j * phase(i)));
            XP = options.scale * real(lambda * VR * exp(1j * phase(i)));

            X = modal.X0(node_idx_modal(j), 1:3).' + dX(1:3);
            Phi = modal.X0(node_idx_modal(j), 4:6).' + dX(4:6);
            R = rotation_vector_to_rotation_matrix(Phi);
          else
            X = nodes(j).X0;
            Phi = nodes(j).Phi0;
            R = nodes(j).R0;
            XP = zeros(6, 1);
          endif

          fprintf(fd_mov, "%e ", X);

          od = nodes(j).orientation_description;

          switch (od)
            case "euler123"
              Phi = rotation_matrix_to_euler123(R);
              fprintf(fd_mov, "%e ", Phi * 180 / pi);
            case "euler313"
              Phi = rotation_matrix_to_euler313(R);
              fprintf(fd_mov, "%e ", Phi * 180 / pi);
            case "euler321"
              Phi = rotation_matrix_to_euler321(R);
              fprintf(fd_mov, "%e ", Phi * 180 / pi);
            case "phi"
              fprintf(fd_mov, "%e ", Phi);
            case "mat"
              fprintf(fd_mov, "%e ", R.');
            case "none"
              ## must be ignored
            otherwise
              error("orientation description \"%s\" not supported!",od);
          endswitch

          fprintf(fd_mov, "%e ", XP);
          fprintf(fd_mov, "\n");
        endfor
      endfor
    unwind_protect_cleanup
      if (fd_mov ~= -1)
        fclose(fd_mov);
      endif

      fd_mov = -1;

      if (fd_out ~= -1)
        fclose(fd_out);
      endif

      fd_out = -1;
    end_unwind_protect

    if (isfield(options, "f_run_mbdyn2easyanim") && options.f_run_mbdyn2easyanim)
      input_file_k = mbdyn_post_output_filename(output_filename);
      log_file_k = mbdyn_post_output_filename(input_file_k, ".log");

      if (2 == exist(log_file_k, "file"))
        if (options.verbose)
          fprintf(stderr, "log file \"%s\" is removed ...\n", log_file_k);
        endif

        [err, msg] = unlink(log_file_k);

        if (err ~= 0)
          error("unlink(\"%s\") returned with status %d: %s", log_file_k, err, msg);
        endif
      endif

      if (options.verbose)
        fprintf(stderr, "creating file \"%s\" ...\n", log_file_k);
      endif

      [status, msg] = copyfile(log_file, log_file_k);

      if (status ~= 1)
        error("copyfile(\"%s\", \"%s\") failed with status %d: %s", log_file, log_file_k, status, msg);
      endif

      options_run.f_run_mbdyn2easyanim = options.f_run_mbdyn2easyanim;

      if (isfield(options, "every"))
        options_run.every = options.every;
      endif

      if (isfield(options, "showAll"))
        options_run.showAll = options.showAll;
      endif

      if (isfield(options, "f_run_EasyAnim"))
        options_run.f_run_EasyAnim = options.f_run_EasyAnim;
      endif

      options_run.f_run_mbdyn = false;

      if (options.verbose)
        fprintf(stderr, "eigenmode %d: %gHz file=\"%s\"\n", mode_index(k), imag(modal.lambda(mode_index(k))) / (2*pi), input_file_k);
      endif

      options_run.verbose = options.verbose;

      mbdyn_solver_run(input_file_k, options_run);
    endif
  endfor
endfunction
