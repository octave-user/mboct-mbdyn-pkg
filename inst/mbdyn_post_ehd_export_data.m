## Copyright (C) 2020(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{filenames}] = mbdyn_post_ehd_export_data(@var{mesh}, @var{res}, @var{outputfile}, @var{idx_t})
## Export hydrodynamic plain bearing data to Gmsh format.
##
## @var{mesh} @dots{} Return value form mbdyn_post_ehd_create_mesh.
##
## @var{res}.log_dat @dots{} Return value from mbdyn_post_load_log.
##
## @var{res}.trajectory @dots{} Return value from mbdyn_post_load_output_struct.
##
## @var{res}.bearings @dots{} Return value from mbdyn_post_ehd_load_output.
##
## @var{output_file} @dots{} Prefix for output filename.
##
## @var{idx_t} @dots{} If only a subset of available data should be exported, @var{idx_t} should contain the time step indices.
##
## @seealso{mbdyn_post_ehd_create_mesh, mbdyn_post_load_log, mbdyn_post_load_output_struct, mbdyn_post_ehd_load_output}
## @end deftypefn

function filenames = mbdyn_post_ehd_export_data(mesh, res, outputfile, idx_t, options)
  if (nargin < 3 || nargin > 5 || nargout > 1)
    print_usage();
  endif

  if (nargin < 4)
    idx_t = 1:numel(res.t);
  endif

  if (nargin < 5)
    options = struct();
  endif

  if (~isfield(options, "deformation_scale"))
    options.deformation_scale = 1000;
  endif

  filenames = cell(numel(idx_t), 1);

  [outputdir, outputname, outputext] = fileparts(outputfile);

  for k=1:numel(idx_t)
    filenames{k} = fullfile(outputdir, sprintf("%s_%03d.msh", outputname, k));

    fd = -1;

    unwind_protect
      fd = fopen(filenames{k}, "wt");

      if (fd == -1)
        error("failed to open file \"%s\"", filenames{k});
      endif

      disp_fields = {"nodal position", "deformation journal", "deformation shell"};

      U = zeros(3, rows(mesh.nodes), numel(disp_fields));

      for i=1:numel(mesh.bearings)
        switch (mesh.bearings(i).cylindrical.mesh_pos)
          case "journal"
            imnode = 1;
          case "shell"
            imnode = 2;
          otherwise
            error("invalid mesh position found in bearing data");
        endswitch

        idxnode = find([res.log_dat.nodes.label] == mesh.bearings(i).cylindrical.nodes(imnode).label);

        X_rb = res.trajectory{idxnode}(idx_t(k), 1:3).';
        Phi_rb = res.trajectory{idxnode}(idx_t(k), 4:6).';

        persistent rotfuncs = {"euler123", @euler123_to_rotation_matrix;
                               "euler313", @euler313_to_rotation_matrix;
                               "euler321", @euler321_to_rotation_matrix;
                               "phi", @rotation_vector_to_rotation_matrix};

        rotfunc = [];

        for j=1:numel(rotfuncs)
          if (strcmp(rotfuncs{j, 1}, res.log_dat.nodes(idxnode).orientation_description))
            rotfunc = rotfuncs{j, 2};
            break;
          endif
        endfor

        if (~numel(rotfunc))
          error("invalid orientation description");
        endif

        dm = mesh.bearings(i).cylindrical.dm;
        z = mesh.bearings(i).z;
        Phi = mesh.bearings(i).Phi;

        if (isfield(res.bearings(i).columns, "w1_n"))
          w1 = options.deformation_scale * res.bearings(i).columns.w1_n(idx_t(k), :);
          w1 -= (mesh.bearings(i).A * (mesh.bearings(i).A \ w1(:))).';
        else
          w1 = zeros(1, numel(mesh.bearings(i).nodeidx));
        endif

        if (isfield(res.bearings(i).columns, "w2_n"))
          w2 = options.deformation_scale * res.bearings(i).columns.w2_n(idx_t(k), :);
          w2 -= (mesh.bearings(i).A * (mesh.bearings(i).A \ w2(:))).';
        else
          w2 = zeros(1, numel(mesh.bearings(i).nodeidx));
        endif

        R_rb = rotfunc(Phi_rb);
        Rb = mesh.bearings(i).cylindrical.nodes(imnode).Rb;
        ob = mesh.bearings(i).cylindrical.nodes(imnode).o;
        Xw1 = [(0.5 * dm - w1) .* cos(Phi);
               (0.5 * dm - w1) .* sin(Phi);
               z];
        Xw2 = [(0.5 * dm + w2) .* cos(Phi);
               (0.5 * dm + w2) .* sin(Phi);
               z];
        Xn = mesh.bearings(i).Xn;
        X0 = mesh.nodes(mesh.noffset(i) + (1:columns(Xn)), :).';
        Urb = R_rb * (Rb * Xn + ob) + X_rb - X0;
        Uw1 = R_rb * (Rb * Xw1 + ob) + X_rb - X0;
        Uw2 = R_rb * (Rb * Xw2 + ob) + X_rb - X0;
        U(:, mesh.noffset(i) + (1:columns(Urb)), 1) = Urb;
        U(:, mesh.noffset(i) + (1:columns(Uw1)), 2) = Uw1;
        U(:, mesh.noffset(i) + (1:columns(Uw2)), 3) = Uw2;
      endfor

      fputs(fd, "$MeshFormat\n2.2 0 8\n$EndMeshFormat\n");

      for i=1:size(U, 3)
        fprintf(fd, "$NodeData\n1\n\"%s\"\n1\n%e\n3\n%d\n3\n%d\n", disp_fields{i}, res.t(idx_t(k)), k - 1, rows(mesh.nodes));
        fprintf(fd, "%d %e %e %e\n", [1:rows(mesh.nodes); U(:, :, i)]);
        fputs(fd, "$EndNodeData\n");
      endfor

      for i=1:numel(res.bearings)
        if (isfield(res.bearings(i).columns, "p_n") && isfield(res.bearings(i).columns, "pc_n"))
          res.bearings(i).columns.ptot_n = res.bearings(i).columns.p_n + res.bearings(i).columns.pc_n;
        endif
      endfor

      field_data = {"p_n",         "fluid pressure";
                    "pc_n",        "contact pressure";
                    "ptot_n",      "total pressure";
                    "rho_n",       "fluid density";
                    "h_n",         "radial clearance";
                    "dh_dt_n",     "radial clearance - time derivative";
                    "U1x_n",       "circumferential sliding velocity journal";
                    "U1z_n",       "axial sliding velocity journal";
                    "U2x_n",       "circumferential sliding velocity shell";
                    "U2z_n",       "axial sliding velocity shell";
                    "tau_xy_0_n",  "circumferential fluid shear stress shell";
                    "tau_yz_0_n",  "axial fluid shear stress shell";
                    "tau_xy_h_n",  "circumferential fluid shear stress journal";
                    "tau_yz_h_n",  "axial fluid shear stress journal";
                    "tauc_xy_0_n", "circumferential contact shear stress shell";
                    "tauc_yz_0_n", "axial contact shear stress shell";
                    "wtot_n",      "total radial deformation";
                    "dwtot_dt_n",  "total radial deformation - time derivative";
                    "w1_n",        "radial deformation journal";
                    "w2_n",        "radial deformation shell"};

      field_name = {field_data{:, 1}};
      field_type = {field_data{:, 2}};

      if (isfield(mesh.elements, "quad4"))
        for l=1:numel(field_name)
          inumelem = int32(0);
          valid_data = false(numel(mesh.bearings), 1);

          for i=1:numel(mesh.bearings)
            valid_data(i) = isfield(res.bearings(i), "columns") ...
                            && isfield(res.bearings(i).columns, field_name{l});
            if (valid_data(i))
              inumelem += numel(mesh.bearings(i).elements);
            endif
          endfor

          fprintf(fd, "$ElementNodeData\n1\n\"%s\"\n1\n%e\n3\n%d\n%d\n%d\n", field_type{l}, res.t(idx_t(k)), k - 1, 1, inumelem);

          for i=1:numel(mesh.bearings)
            if (valid_data(i))
              idxnode = int32([1:4]);

              val = getfield(res.bearings(i).columns, field_name{l});

              valout = zeros(2 + numel(idxnode), numel(mesh.bearings(i).elements));
              valout(1, :) = mesh.bearings(i).elements;
              valout(2, :) = numel(idxnode);

              for m=1:numel(idxnode)
                valout(2 + m, :) = val(idx_t(k), mesh.elements.quad4(mesh.bearings(i).elements, idxnode(m)) - mesh.noffset(i));
              endfor

              format = ["%d %d", repmat(" %e", 1, numel(idxnode)), "\n"];

              fprintf(fd, format, valout);
            endif
          endfor

          fputs(fd, "$EndElementNodeData\n");
        endfor
      endif
    unwind_protect_cleanup
      if (fd ~= -1)
        fclose(fd);
        fd = -1;
      endif
    end_unwind_protect
  endfor
endfunction
