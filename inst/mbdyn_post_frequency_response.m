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
## @deftypefn {Function File} [@var{F}] = mbdyn_post_frequency_response(@var{modal}, @var{dof_info}, @var{excitation}, @var{response}, @var{omega}, @var{p}, @var{options})
## Compute the frequency response function of the linearized equations of motion of a MBDyn model.
##
## @var{F} @dots{} Complex frequency response function
##
## @var{modal} @dots{} Data structure returned from mbdyn_post_load_output_eig()
##
## @var{dof_info} @dots{} Data structure returned form mbdyn_post_load_log_node()
##
## @var{excitation}.node_label @dots{} Label of the node that receives the force @var{p}
##
## @var{excitation}.offset @dots{} Offset between force and node
##
## @var{excitation}.direction @dots{} Direction of the force @var{p}
##
## @var{response}(@var{i}).node_label @dots{} Label of the node where the response @var{i} is measured
##
## @var{response}(@var{i}).offset @dots{} Offset between the point where the response @var{i} is measured and the node
##
## @var{response}(@var{i}).direction @dots{} Measurement direction of the response @var{i}
##
## @var{omega} @dots{} Angular velocity of the excitation force @var{p}
##
## @var{p} @dots{} Complex excitation force
##
## @var{options}.singular @dots{} If true, it is assumed that df/dy + j * omega * df/dy_dot is singular
##
## @var{options}.matrix_type @dots{} "wre" means F(omega_idx, response_idx, excitation_idx)
##
##                                   "rew" means F(response_idx, excitation_idx, omega_idx)
##
## @end deftypefn

function [F, Tp, Tu] = mbdyn_post_frequency_response(modal, dof_info, excitation, response, omega, p, options)
  if (nargin < 6 || nargin > 7 || nargout > 1)
    print_usage();
  endif

  if (nargin < 7)
    options = struct();
  endif

  if (~isfield(options,"singular"))
    options.singular = false;
  endif

  if (~isfield(options, "number_of_processors"))
    options.number_of_processors = 1;
  endif

  if (~isfield(options, "solver"))
    options.solver = "mldivide";
  endif

  if (~isfield(options, "symmetric"))
    options.symmetric = false; ## The default for MBDyn since A will be symmetric only in rare situations.
  endif

  if (~isfield(options, "matrix_format"))
    options.matrix_format = "wre";
  endif

  for i=1:length(excitation)
    if (~isfield(excitation, "offset"))
      excitation(i).offset = [];
    endif

    if (~isfield(excitation, "direction"))
      excitation(i).direction = [];
    endif

    if (~isfield(excitation, "component"))
      excitation(i).component = [];
    endif
  endfor

  for i=1:length(response)
    if (~isfield(response, "offset"))
      response(i).offset = [];
    endif

    if (~isfield(response, "direction"))
      response(i).direction = [];
    endif

    if (~isfield(response, "component"))
      response(i).component = [];
    endif
  endfor

  df_dy_dot = (modal.Aplus + modal.Aminus) / 2;
  df_dy = (modal.Aplus - modal.Aminus) / (2 * modal.dCoef);

  Tp = mbdyn_build_excitation_matrix(dof_info, columns(df_dy), excitation, "force");
  Tu = mbdyn_build_excitation_matrix(dof_info, columns(df_dy), response, "displacement");

  N_response = rows(Tu);

  if (columns(p) == 1)
    N_excitation = rows(Tp);
  else
    N_excitation = 1;
  endif

  switch(options.matrix_format)
    case "wre"
      F = zeros(length(omega), N_response, N_excitation);
    case "rew"
      F = zeros(N_response, N_excitation, length(omega));
  endswitch

  switch (options.solver)
    case "mldivide"
      options.solver_func = @mbdyn_post_linear_solver_mldivide;
    otherwise
      if (~exist("fem_sol_factor", "file"))
        pkg load mboct-fem-pkg;
      endif
      options.solver_func = @mbdyn_post_linear_solver_factor;
  endswitch

  if (options.number_of_processors == 1)
    for i=1:length(omega)
      y = mbdyn_post_frequency_response_helper(i, df_dy, df_dy_dot, Tp, p, Tu, omega, options);

      switch(options.matrix_format)
        case "wre"
          F(i, :, :) = y;
        case "rew"
          F(:, :, i) = y;
      endswitch
    endfor
  else
    options_par.number_of_processors = options.number_of_processors;
    options_par.number_of_parameters = length(omega);
    y = run_parallel(options_par, @mbdyn_post_frequency_response_helper, df_dy, df_dy_dot, Tp, p, Tu, omega, options);

    for i=1:length(omega)
      switch(options.matrix_format)
        case "wre"
          F(i, :, :) = y{i};
        case "rew"
          F(:, :, i) = y{i};
      endswitch
    endfor
  endif
endfunction

function Tp = mbdyn_build_excitation_matrix(dof_info, number_dofs, excitation, type)
  Tp_i = cell(numel(excitation), 1);
  number_of_rows = int32(0);

  for i=1:numel(excitation)
    Tp_i{i} = mbdyn_post_trans_mat_struct_node(dof_info, number_dofs, excitation(i).node_label, type, excitation(i).offset, excitation(i).direction, excitation(i).component);
    number_of_rows += rows(Tp_i{i});
  endfor

  Tp = sparse([], [], [], number_of_rows, number_dofs);

  number_of_rows = int32(0);

  for i=1:numel(excitation)
    Tp((number_of_rows + 1):(number_of_rows + rows(Tp_i{i})), :) = Tp_i{i};
    number_of_rows += rows(Tp_i{i});
  endfor
endfunction
