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
## @deftypefn {Function File} @var{options} = mbdyn_pre_solid_write_const_laws(@var{mesh}, @var{csl_file}, @var{options})
##
## Generate an MBDyn input file <@var{csl_file}> containing all constitutive laws from <@var{mesh}>.
##
## @var{mesh} @dots{} Finite Element mesh
##
## @var{csl_file} @dots{} filename of output file
##
## @var{options}.const_laws.number @dots{} starting index of the next constitutive law minus one
##
## @end deftypefn

function options = mbdyn_pre_solid_write_const_laws(mesh, csl_file, options)
  if (~(nargin == 3 && isstruct(mesh) && ischar(csl_file) && isstruct(options)))
    print_usage();
  endif

  if (~isfield(options, "const_laws"))
    options.const_laws = struct();
  endif

  if (~isfield(options.const_laws, "number"))
    options.const_laws.number = 0;
  endif

  fd = -1;

  unwind_protect
    [fd, msg] = fopen(csl_file, "wt");

    if (fd == -1)
      error("failed to open file \"%s\": %s", csl_file, msg);
    endif

    for i=1:numel(mesh.material_data)
      f_linear_elastic_generic = false;

      if (isfield(mesh.material_data, "type"))
        mat_type = mesh.material_data(i).type;
      else
        f_linear_elastic_generic = isfield(mesh.material_data, "C") && ~isempty(mesh.material_data(i).C);

        if (isfield(mesh.material_data, "beta") && ~isempty(mesh.material_data(i).beta) && mesh.material_data(i).beta ~= 0)
          if (f_linear_elastic_generic)
            mat_type = "linear viscoelastic generic";
          else
            mat_type = "hookean linear viscoelastic isotropic";
          endif
        else
          if (f_linear_elastic_generic)
            mat_type = "linear elastic generic";
          else
            mat_type = "hookean linear elastic isotropic";
          endif
        endif
      endif

      if (isfield(mesh.material_data, "extra_data") && ~isempty(mesh.material_data(i).extra_data))
        extra_data = mesh.material_data(i).extra_data;
      else
        extra_data = "";
      endif

      if (~ischar(extra_data))
        error("mesh.material_data(%d).extra_data must be a string", i);
      endif

      switch (mat_type)
        case {"neo hookean elastic", "neo hookean viscoelastic", "mooney rivlin elastic", "hookean linear elastic isotropic", "hookean linear viscoelastic isotropic", "linear elastic generic", "linear viscoelastic generic", "bilinear isotropic hardening", "linear viscoelastic maxwell1", "linear viscoelastic maxwelln", "mfront small strain", "mfront finite strain"}
          switch (mat_type)
            case {"linear elastic generic", "linear viscoelastic generic", "linear viscoelastic maxwell1", "linear viscoelastic maxwelln"}
              if (f_linear_elastic_generic)
                C = mesh.material_data(i).C;
              else
                C = fem_pre_mat_isotropic(mesh.material_data(i).E, mesh.material_data(i).nu);
              endif
          endswitch

          switch (mat_type)
            case "neo hookean elastic"
              extra_format = sprintf("E, %.16e, nu, %.16e", mesh.material_data(i).E, mesh.material_data(i).nu);
            case "neo hookean viscoelastic"
              extra_format = sprintf("E, %.16e, nu, %.16e, beta, %.16e", mesh.material_data(i).E, mesh.material_data(i).nu, mesh.material_data(i).beta);
            case {"mooney rivlin elastic"}
              if (isfield(mesh.material_data(i), "C1") && isfield(mesh.material_data(i), "C2") && isfield(mesh.material_data(i), "kappa"))
                extra_format = sprintf("C1, %.16e, C2, %.16e, kappa, %.16e", mesh.material_data(i).C1, mesh.material_data(i).C2, mesh.material_data(i).kappa);
              else
                if (isfield(mesh.material_data(i), "delta"))
                  delta = mesh.material_data(i).delta;
                else
                  delta = 0;
                endif
                extra_format = sprintf("E, %.16e, nu, %.16e, delta, %.16e", mesh.material_data(i).E, mesh.material_data(i).nu, delta);
              endif
            case "linear elastic generic"
              extra_format = sprintf("matr %s", sprintf(", %.16e", C));
            case "linear viscoelastic generic"
              extra_format = sprintf("matr %s, proportional, %.16e", sprintf(", %.16e", C), mesh.material_data(i).beta);
            case {"hookean linear elastic isotropic"}
              extra_format = sprintf("E, %.16e, nu, %.16e", ...
                                     mesh.material_data(i).E, ...
                                     mesh.material_data(i).nu);
            case "hookean linear viscoelastic isotropic"
              extra_format = sprintf("E, %.16e, nu, %.16e, beta, %.16e", ...
                                     mesh.material_data(i).E, ...
                                     mesh.material_data(i).nu, ...
                                     mesh.material_data(i).beta);
            case {"bilinear isotropic hardening"}
              extra_format = sprintf("E, %.16e, nu, %.16e, ET, %.16e, sigmayv, %.16e", ...
                                     mesh.material_data(i).E, ...
                                     mesh.material_data(i).nu, ...
                                     mesh.material_data(i).ET, ...
                                     mesh.material_data(i).sigmayv);
            case {"linear viscoelastic maxwell1", "linear viscoelastic maxwelln"}
              E0 = mesh.material_data(i).E;

              if (numel(mesh.material_data(i).theta) ~= numel(mesh.material_data(i).tau))
                error("numel(mesh.material_data(%d).theta) does not match numel(mesh.material_data(%d).tau)", i, i);
              endif

              E1 = E0 * mesh.material_data(i).theta;
              eta1 = E1 .* mesh.material_data(i).tau;

              switch (mat_type)
                case "linear viscoelastic maxwelln"
                  N = sprintf("%d, ", numel(E1));
                otherwise
                  N = "";

                  if (~isscalar(mesh.material_data(i).theta))
                    error("mesh.material_data(%d).theta must be a scalar", i);
                  endif
              endswitch

              extra_format = sprintf("E0, %.16e, %s%sC%s", ...
                                     E0, ...
                                     N, ...
                                     sprintf("E1, %.16e, eta1, %.16e, ", ...
                                             [reshape(E1, 1, numel(E1));
                                              reshape(eta1, 1, numel(eta1))]), ...
                                     sprintf(", %.16e", C.' / E0));
            case {"mfront small strain", "mfront finite strain"}
              extra_format = sprintf("library path, \"%s\", name, \"%s\"", mesh.material_data(i).library_path, mesh.material_data(i).name);

              var_names = {"parameters", "properties"};

              for j=1:numel(var_names)
                if (isfield(mesh.material_data(i), var_names{j}))
                  vars = getfield(mesh.material_data(i), var_names{j});
                  extra_format = [extra_format, sprintf(", %s, %d", var_names{j}, numel(vars))];
                  for j=1:numel(vars)
                    extra_format = [extra_format, ...
                                    sprintf(", \"%s\", %.16e", ...
                                            vars(j).name, ...
                                            vars(j).value)];
                  endfor
                endif
              endfor
          endswitch

          enable_csl_dim = false(MBDYN_CSL_COUNT, 1);

          f_enable_compressible_mat = true;
          f_enable_incompressible_mat = false;
          f_enable_green_lagrange = true;
          f_enable_deformation_gradient = false;

          switch (mat_type)
            case {"hookean linear elastic isotropic", "mooney rivlin elastic", "bilinear isotropic hardening"}
              if (isfield(mesh.material_data, "nu") && ~isempty(mesh.material_data(i).nu))
                f_enable_compressible_mat = mesh.material_data(i).nu < 0.5;
              elseif (isfield(mesh.material_data, "kappa") && ~isempty(mesh.material_data.kappa))
                f_enable_compressible_mat = mesh.material_data(i).kappa < inf;
              endif
              f_enable_incompressible_mat = true;
          endswitch

          switch (mat_type)
            case {"mooney rivlin elastic", "mfront finite strain"}
              f_enable_deformation_gradient = true;
          endswitch

          switch (mat_type)
            case "mfront finite strain"
              f_enable_green_lagrange = false;
          endswitch

          enable_csl_dim(MBDYN_CSL_IDX_COMPRESSIBLE) = f_enable_compressible_mat && f_enable_green_lagrange;
          enable_csl_dim(MBDYN_CSL_IDX_INCOMPRESSIBLE) = f_enable_incompressible_mat && f_enable_green_lagrange;
          enable_csl_dim(MBDYN_CSL_IDX_DEF_GRAD) = f_enable_compressible_mat && f_enable_deformation_gradient;

          csl_format = "constitutive law: %d, name, \"solid%d\", %d, %s, %s;\n";

          csl_dim = int32([MBDYN_CSL_DIM_COMPRESSIBLE;
                           MBDYN_CSL_DIM_INCOMPRESSIBLE;
                           MBDYN_CSL_DIM_DEF_GRAD]);

          csl_dim = csl_dim(enable_csl_dim);

          extra_format = [extra_format, extra_data];

          ++options.const_laws.number;

          for j=1:numel(csl_dim)
            ## MBDyn allows us to use the same constitutive law number for different dimensions
            fprintf(fd, csl_format, options.const_laws.number, i, csl_dim(j), mat_type, extra_format);
          endfor
        otherwise
          error("unknown material type: \"%s\"", mat_type);
      endswitch
    endfor
  unwind_protect_cleanup
    if (fd ~= -1)
      fclose(fd);
    endif

    fd = -1;
  end_unwind_protect
endfunction

function idx = MBDYN_CSL_IDX_COMPRESSIBLE()
  idx = int32(1);
endfunction

function idx = MBDYN_CSL_IDX_INCOMPRESSIBLE()
  idx = int32(2);
endfunction

function idx = MBDYN_CSL_IDX_DEF_GRAD()
  idx = int32(3);
endfunction

function idx = MBDYN_CSL_COUNT()
  idx = int32(3);
endfunction

function idx = MBDYN_CSL_DIM_COMPRESSIBLE()
  idx = int32(6);
endfunction

function idx = MBDYN_CSL_DIM_INCOMPRESSIBLE()
  idx = int32(7);
endfunction

function idx = MBDYN_CSL_DIM_DEF_GRAD()
  idx = int32(9);
endfunction
