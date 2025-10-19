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
## @deftypefn {Function File} @var{modal} = mbdyn_post_load_output_eig(@var{mbdyn_output_file})
## @deftypefnx {} @dots{} = mbdyn_post_load_output_eig(@var{mbdyn_output_file}, @var{options})
## @deftypefnx {} @dots{} = mbdyn_post_load_output_eig(@var{mbdyn_output_file}, @var{options}, @var{index})
##
## Loads data form an eigenanalysis from MBDyn output file "<@var{mbdyn_output_file}_%02d.m>".
##
## @var{mbdyn_output_file} @dots{} Name of the MBDyn output file without extension.
##
## @var{options}.positive_frequencies @dots{} Return only eigenvalues with imaginary values greater than zero.
##
## @var{options}.solve_qz @dots{} Solve the eigenvalue problem within Octave.
##
## @var{index} @dots{} If there are more than one eigenanalysis results (e.g. "mod_01.m", "mod_02.m", @dots{}) then load the file corresponding to the index <@var{index}>.
##
## @var{modal} @dots{} Structure that contains mode shapes, eigenvalues and node positions.
##
## @end deftypefn

function modal = mbdyn_post_load_output_eig(mbdyn_output_file, options, index)
  if (nargin < 1 || nargin > 3 || nargout > 1)
    print_usage();
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (nargin < 3)
    index = 0;
  endif

  if (~isfield(options, "positive_frequencies"))
    options.positive_frequencies = true;
  endif

  [inp_dir, inp_name, inp_ext] = fileparts(mbdyn_output_file);

  empty_cell = cell(1, 3);

  output_files = struct("name", empty_cell, "type", empty_cell);

  output_files(1).name = fullfile(inp_dir, [inp_name, ".nc"]);
  output_files(1).type = "netcdf";
  output_files(2).name = fullfile(inp_dir, sprintf("%s_%02d.m", inp_name, index));
  output_files(2).type = "m";
  output_files(3).name = fullfile(inp_dir, sprintf("%s_%d.m", inp_name, index));
  output_files(3).type = "m";  
  output_files(4).name = fullfile(inp_dir, [inp_name, ".m"]);
  output_files(4).type = "m";

  idx = false(size(output_files));

  for i=1:numel(output_files)
    [info, err, msg] = stat(output_files(i).name);
    idx(i) = (err == 0);
  endfor

  output_files = output_files(idx);

  clear idx;

  if (isempty(output_files))
    error("file not found \"%s\"", mbdyn_output_file);
  endif

  if (isfield(options, "use_netcdf"))
    for i=1:numel(output_files)
      switch (output_files(i).type)
        case "netcdf"
          if (options.use_netcdf)
            output_files = output_files(i);
            break;
          endif
        otherwise
          if (~options.use_netcdf)
            output_files = output_files(i);
            break;
          endif
      endswitch
    endfor
  endif

  fprintf(stderr, "loading file \"%s\":%d\n", output_files(1).name, index);

  switch (output_files(1).type)
    case "netcdf"
      modal = mbdyn_post_load_modal_data_nc(output_files(1).name, index);
    otherwise
      modal = mbdyn_post_load_modal_data_m(output_files(1).name);
  endswitch

  sparse_threshold = 0.1;

  if (isfield(modal, "Aplus"))
    if (nnz(modal.Aplus) < sparse_threshold * numel(modal.Aplus))
      modal.Aplus = sparse(modal.Aplus);
    endif
  endif

  if (isfield(modal, "Aminus"))
    if (nnz(modal.Aminus) < sparse_threshold * numel(modal.Aminus))
      modal.Aminus = sparse(modal.Aminus);
    endif
  endif

  if (isfield(modal, "alpha") && isfield(modal, "dCoef"))
    LAMBDA = (modal.alpha(:,1) + 1j * modal.alpha(:, 2)) ./ modal.alpha(:, 3);
    modal.lambda = 1 / modal.dCoef * (LAMBDA - 1) ./ (LAMBDA + 1);
    modal.f = imag(modal.lambda) / (2 * pi);
  else
    modal.lambda = [];
    modal.f = [];
  endif

  if (~isempty(modal.f))
    [modal.f, idx_f] = sort(modal.f);

    if (options.positive_frequencies)
      idx_gtz = find(modal.f > 0);
      modal.f = modal.f(idx_gtz);
      idx_f = idx_f(idx_gtz);
    endif

    modal.lambda = modal.lambda(idx_f);

    if (isfield(modal, "VR"))
      modal.VR = modal.VR(:, idx_f);
    endif

    if (isfield(modal, "VL"))
      modal.VL = modal.VL(:, idx_f);
    endif

    modal.alpha = modal.alpha(idx_f, :);
  endif

  if (isfield(modal, "lambda"))
    delta = -real(modal.lambda);
    omegad = imag(modal.lambda);
    a0 = (delta ./ omegad).^2;
    modal.D = sqrt(a0 ./ (1 + a0));
    modal.omega0 = omegad ./ sqrt(1 - modal.D.^2);
  endif

  if (~isfield(modal, "VR"))
    modal.VR = [];
  endif

  if (~isfield(modal, "VL"))
    modal.VL = [];
  endif
endfunction

function modal = mbdyn_post_load_modal_data_m(mbdyn_output_file)
  source(mbdyn_output_file);

  var_names = {"dTime", "dCoef", "lStep", "Aplus", "Aminus", "VR", "VL", "alpha", "X0", "idx", "labels"};

  modal = struct();

  for i=1:numel(var_names)
    if (1 == exist(var_names{i}, "var"))
      modal = setfield(modal, var_names{i}, eval(var_names{i}));
    endif
  endfor

  if (isfield(modal, "X0"))
    modal.X0 = reshape(modal.X0, 6, numel(modal.idx)).';
  endif
endfunction

function modal = mbdyn_post_load_modal_data_nc(mbdyn_output_file, index)
  pkg load netcdf;

  modal = struct();

  prefix = sprintf("eig.%d.", index);

  var_names = {[prefix, "time"],   "dTime";
               [prefix, "dCoef"],  "dCoef";
               [prefix, "step"],   "lStep";
               [prefix, "Aplus"],  "Aplus";
               [prefix, "Aminus"], "Aminus";
               [prefix, "VR"],     "VR";
               [prefix, "VL"],     "VL";
               [prefix, "alpha"],  "alpha";
               [prefix, "X0"],     "X0";
               "eig.idx",          "idx";
               "eig.labels",       "labels";
               "eig.joint.idx",    "joint_idx";
               "eig.joint.labels", "joint_labels";
               "eig.genel.idx",    "genel_idx";
               "eig.genel.labels", "genel_labels"};

  for i=1:rows(var_names)
    try
      value = ncread(mbdyn_output_file, var_names{i, 1});
    catch
      switch (lasterror.message)
        case "NetCDF: Variable not found"
          continue
        otherwise
          rethrow(lasterror());
      endswitch
    end_try_catch
    modal = setfield(modal, var_names{i, 2}, value);
  endfor

  if (isfield(modal, "time"))
    modal.dTime = modal.time;
    modal = rmfield(modal, "time");
  endif

  if (isfield(modal, "alpha"))
    modal.alpha = modal.alpha.';
  endif

  if (isfield(modal, "VR"))
    modal.VR = complex(modal.VR(:, :, 1), modal.VR(:, :, 2));
  endif

  if (isfield(modal, "VL"))
    modal.VL = complex(modal.VL(:, :, 1), modal.VL(:, :, 2));
  endif

  if (isfield(modal, "Aplus"))
    Aplus_attr = ncreadatt(mbdyn_output_file, [prefix, "Aplus"], "matrix type");
    Aminus_attr = ncreadatt(mbdyn_output_file, [prefix, "Aminus"], "matrix type");

    switch (Aplus_attr)
      case "sparse"
        modal.Aplus = spconvert(modal.Aplus.');
      case "dense"
        if (rows(modal.Aplus) ~= columns(modal.Aplus))
          error("invalid matrix format");
        endif
      otherwise
        error("unknown matrix type: \"%s\"", Aplus_attr);
    endswitch

    switch (Aminus_attr)
      case "sparse"
        modal.Aminus = spconvert(modal.Aminus.');
      case "dense"
        if (rows(modal.Aminus) ~= columns(modal.Aminus))
          error("invalid matrix format");
        endif
      otherwise
        error("unknown matrix type: \"%s\"", Aminus_attr);
    endswitch
  endif
endfunction
