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
        else
          w1 = zeros(1, numel(mesh.bearings(i).nodeidx));
        endif

        if (isfield(res.bearings(i).columns, "w2_n"))
          w2 = options.deformation_scale * res.bearings(i).columns.w2_n(idx_t(k), :);
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

%!demo
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
%!            3.196  3.048   2.989   2.940   2.848   2.769   2.708   2.599];
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
%! close all;
%! param.d1 = 12e-3;
%! param.Psi = 1.2e-3;
%! param.M = int32(30);
%! param.N = int32(300);
%! param.output_bearing_data = true;
%! param.B = 12e-3;
%! param.epsilon = 0.9;
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
%!     fputs(fd, "set: real n = 1;\n");
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
%!   opt_load.loaded_fields = {"p", "h", "rho", "F1", "M1", "F2", "M2"};
%!   opt_load.interpolate_mesh = false;
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
%!   beta_ref = pi / 180 * interp2(B_d_r, epsilon_r, beta_r, B / d, epsilon, "linear");
%!   mu_ref = Psi * (abs((omega1z - omega2z) / omega_res) * pi / (sqrt(1 - epsilon^2) * So_ref) + sin(beta_ref) * abs(epsilon) / 2);
%!   Q_ref = interp2(B_d_r, epsilon_r, Q_r, B / d, epsilon, "linear");
%!   post_pro_mesh_file = [output_file, "_hd_mesh.msh"];
%!   mesh = mbdyn_post_ehd_create_mesh(res.log_dat);
%!   mbdyn_post_ehd_export_mesh(mesh, post_pro_mesh_file);
%!   post_pro_data_files = mbdyn_post_ehd_export_data(mesh, res, [output_file, "_hd_data"]);
%!   post_pro_script_file = [output_file, "_hd_post_pro.geo"];
%!   fd = -1;
%!   unwind_protect
%!     [fd] = fopen(post_pro_script_file, "wt");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", post_pro_script_file);
%!     endif
%!     fprintf(fd, "Merge \"%s\";\n", post_pro_mesh_file);
%!     fputs(fd, "Mesh.SurfaceEdges = 0;\n");
%!     fputs(fd, "Mesh.SurfaceFaces = 0;\n");
%!     fputs(fd, "Mesh.SurfaceNumbers = 0;\n");
%!     fputs(fd, "Mesh.VolumeEdges = 0;\n");
%!     fputs(fd, "Mesh.VolumeFaces = 0;\n");
%!     for k=1:numel(post_pro_data_files)
%!       fprintf(fd, "Merge \"%s\";\n", post_pro_data_files{k});
%!     endfor
%!     fputs(fd, "View[0].Type = 1;\n");
%!     fputs(fd, "View[0].VectorType = 5;\n");
%!     fputs(fd, "View[0].Visible = 1;\n");
%!     fprintf(fd, "View[0].DisplacementFactor = %g;\n", 1);
%!     fputs(fd, "View[0].ShowTime = 1;\n");
%!     fputs(fd, "View[0].ShowElement = 1;\n");
%!     fputs(fd, "View[0].IntervalsType = 3;\n");
%!     fputs(fd, "View[0].NbIso = 20;\n");
%!     fputs(fd, "View[0].ExternalView = 3;\n");
%!     fputs(fd, "View[0].SaturateValues = 1;\n");
%!     fputs(fd, "iView = 0;\n");
%!     fputs(fd, "For (1:PostProcessing.NbViews - 1)\n");
%!     fputs(fd, "  iView++;\n");
%!     fputs(fd, "  View[iView].Visible = 0;\n");
%!     fputs(fd, "EndFor\n");
%!     fputs(fd, "General.Trackball = 0;\n");
%!     fputs(fd, "General.RotationX = 290;\n");
%!     fputs(fd, "General.RotationY = 0;\n");
%!     fputs(fd, "General.RotationZ = 130;\n");
%!     scale = 0.6;
%!     fprintf(fd, "General.ScaleX = %g;\n", scale);
%!     fprintf(fd, "General.ScaleY = %g;\n", scale);
%!     fprintf(fd, "General.ScaleZ = %g;\n", scale);
%!     fputs(fd, "General.Axes = 0;\n");
%!     fputs(fd, "General.Orthographic = 1;\n");
%!     fputs(fd, "General.RotationCenterGravity = 1;\n");
%!     fputs(fd, "General.DisplayBorderFactor = 0;\n");
%!     for k=1:numel(post_pro_data_files)
%!       fprintf(fd, "View[0].TimeStep = %d;\n", k - 1);
%!       fprintf(fd, "Print \"%s_hd_post_pro_%03d.jpg\";\n", output_file, k);
%!     endfor
%!     fputs(fd, "Exit;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {post_pro_script_file});
%!   status = spawn_wait(pid);
%!   if (status ~= 0)
%!     error("gmsh returned with status %d", status);
%!   endif
%!   fn = dir([output_file, "_hd_post_pro_*.jpg"]);
%!   for k=1:numel(fn)
%!     figure("visible", "off");
%!     [img, map, alpha] = imread(fullfile(fn(k).folder, fn(k).name));
%!     imshow(img, map);
%!     title(printable_title(fn(k).name));
%!   endfor
%!   figure_list();
%!   fprintf(stderr, "So / So_ref - 1 = %.2f\n", So / So_ref - 1);
%!   fprintf(stderr, "beta / beta_ref - 1 = %.2f\n", beta / beta_ref - 1);
%!   fprintf(stderr, "mu / mu_ref - 1 = %.2f\n", mu / mu_ref - 1);
%!   fprintf(stderr, "Q / Q_ref - 1 = %.2f\n", Q / Q_ref - 1);
%!   assert(So, So_ref, 0.04 * So_ref);
%!   assert(beta, beta_ref, 0.02 * beta_ref);
%!   assert(mu, mu_ref, 0.03 * mu_ref);
%!   assert(Q, Q_ref, 0.07 * Q_ref);
%!   assert(dQ, 0, 0.02);
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
