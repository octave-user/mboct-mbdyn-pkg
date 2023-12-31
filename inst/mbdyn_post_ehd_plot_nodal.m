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
## @deftypefn {Function File} mbdyn_post_ehd_plot_nodal(@var{t}, @var{bearings}, @var{options}, @var{log_dat})
## Create a surface or contour plot for output data from MBDyn's hydrodynamic bearings.
##
## @var{t} @dots{} Simulation time returned from mbdyn_post_load_output_out.
##
## @var{bearings} @dots{} Return value from mbdyn_post_ehd_load_output.
##
## @var{options} @dots{} Scalar struct with options.
##
## @var{options}.plot_func @dots{} One of "contour", "contourf", "surf", "surfl", "surfc".
##
## @var{options}.plot_title @dots{} Optional string to be prepended to the title.
##
## @var{options}.plot_range @dots{} If only a subset of available time steps should be plotted, this index vector indicates selected time steps.
##
## @var{options}.plot_bearings @dots{} If only a subset of available bearings should be plotted, this index vector indicates selected bearings.
##
## @var{options}.bearing_title @dots{} Optional cell array of character string names for each bearing to be used in the title.
##
## @var{options}.print_to_file @dots{} Character string name of the output file in Pdf format.
##
## @var{options}.f_show_all @dots{} Make all figures visible.
##
## @var{options}.width @dots{} Width of the figure used for printing.
##
## @var{options}.height @dots{} Height of the figure used for printing.
##
## @var{options}.f_plot_* @dots{} Flag to activate plotting a individual fields (e.g. f_plot_p for pressure, f_plot_h for clearance, ...).
##
## @var{log_dat} @dots{} Return value from mbdyn_post_load_log.
##
## @seealso{mbdyn_post_ehd_load_output, mbdyn_post_load_output_out}
## @end deftypefn

function mbdyn_post_ehd_plot_nodal(t, bearings, options, log_dat)
  if (nargin ~= 4 || nargout > 0)
    print_usage();
  endif

  fields = struct()([]);

  fields(end + 1).column = "p";
  fields(end).description = "fluid pressure";
  fields(end).unit = "Pa";

  fields(end + 1).column = "pc";
  fields(end).description = "contact pressure";
  fields(end).unit = "Pa";

  fields(end + 1).column = "rho";
  fields(end).description = "fluid density";
  fields(end).unit = "kg/m^3";

  fields(end + 1).column = "T";
  fields(end).description = "fluid temperature";
  fields(end).unit = "K";

  fields(end + 1).column = "h";
  fields(end).description = "clearance";
  fields(end).unit = "m";

  fields(end + 1).column = "dh_dt";
  fields(end).description = "clearance derivative";
  fields(end).unit = "m/s";

  fields(end + 1).column = "wtot";
  fields(end).description = "total deformation";
  fields(end).unit = "m";

  fields(end + 1).column = "dwtot_dt";
  fields(end).description = "total deformation derivative";
  fields(end).unit = "m/s";

  fields(end + 1).column = "w1";
  fields(end).description = "deformation at shaft";
  fields(end).unit = "m";

  fields(end + 1).column = "w2";
  fields(end).description = "deformation at bearing";
  fields(end).unit = "m";

  fields(end + 1).column = "U1x";
  fields(end).description = "surface velocity";
  fields(end).unit = "m/s";

  fields(end + 1).column = "U1z";
  fields(end).description = "surface velocity";
  fields(end).unit = "m/s";

  fields(end + 1).column = "U2x";
  fields(end).description = "surface velocity";
  fields(end).unit = "m/s";

  fields(end + 1).column = "U2z";
  fields(end).description = "surface velocity";
  fields(end).unit = "m/s";

  fields(end + 1).column = "tau_xy_0";
  fields(end).description = "fluid stress";
  fields(end).unit = "Pa";

  fields(end + 1).column = "tau_yz_0";
  fields(end).description = "fluid stress";
  fields(end).unit = "Pa";

  fields(end + 1).column = "tau_xy_h";
  fields(end).description = "fluid stress";
  fields(end).unit = "Pa";

  fields(end + 1).column = "tau_yz_h";
  fields(end).description = "fluid stress";
  fields(end).unit = "Pa";

  fields(end + 1).column = "tauc_xy_0";
  fields(end).description = "contact stress";
  fields(end).unit = "Pa";

  fields(end + 1).column = "tauc_yz_0";
  fields(end).description = "contact stress";
  fields(end).unit = "Pa";

  fields(end + 1).column = "qx";
  fields(end).description = "volume flux";
  fields(end).unit = "m^2/s";

  fields(end + 1).column = "qz";
  fields(end).description = "volume flux";
  fields(end).unit = "m^2/s";

  fields(end + 1).column = "mdotx";
  fields(end).description = "mass flux";
  fields(end).unit = "kg/(m s)";

  fields(end + 1).column = "mdotz";
  fields(end).description = "mass flux";
  fields(end).unit = "kg/(m s)";

  fields(end + 1).column = "Qx";
  fields(end).description = "heat flux";
  fields(end).unit = "W/m^2";

  fields(end + 1).column = "Qz";
  fields(end).description = "heat flux";
  fields(end).unit = "W/m^2";

  if (~isfield(options, "plot_title"))
    options.plot_title = "";
  endif

  if (~isfield(options, "plot_range"))
    options.plot_range = 1:length(t);
  endif

  if (~isfield(options, "plot_bearings"))
    options.plot_bearings = 1:length(bearings);
  endif

  if (~isfield(options, "bearing_title"))
    for i=1:length(bearings)
      options.bearing_title{i} = sprintf("%d", log_dat.bearings(i).label);
    endfor
  endif

  if (~isfield(options, "plot_func"))
    options.plot_func = "contourf";
  endif

  if (~isfield(options, "print_to_file"))
    options.print_to_file = "";
  endif

  if (~isfield(options, "f_show_all"))
    options.f_show_all = false;
  endif

  if (~isfield(options, "width"))
    options.width = 1024 / 2;
  endif

  if (~isfield(options, "height"))
    options.height = 768 / 2;
  endif

  ifig = 0;
  file_names_pdf = {};

  for i=1:length(options.plot_bearings)
    xh = mbdyn_post_ehd_interp_grid(log_dat.bearings(options.plot_bearings(i)), "hydro");
    xfx = mbdyn_post_ehd_interp_grid(log_dat.bearings(options.plot_bearings(i)), "flux_x");
    xfz = mbdyn_post_ehd_interp_grid(log_dat.bearings(options.plot_bearings(i)), "flux_z");

    dm = log_dat.bearings(options.plot_bearings(i)).cylindrical.dm;

    tri_grid_h = [];
    tri_grid_fx = [];
    tri_grid_fz = [];

    for k=1:length(fields)
      plot_flag = ["f_plot_", fields(k).column];

      if (isfield(options, plot_flag) && getfield(options, plot_flag))
        field_nodes = [fields(k).column, "_n"];

        if (isfield(bearings(options.plot_bearings(i)).columns, field_nodes))
          zt_n = getfield(bearings(options.plot_bearings(i)).columns, field_nodes);

          if (~numel(zt_n))
            continue;
          endif

          zt_mag = max(max(abs(zt_n)));

          for j=1:length(options.plot_range)
            if (~numel(options.print_to_file))
              nfig = figure("visible", "off");

              if (options.f_show_all)
                set(nfig, "visible", "on");
              endif
            else
              if (ifig == 0)
                nfig = figure("visible", "off");
              endif
              clf();
            endif

            switch (fields(k).column)
              case {"qx", "mdotx", "Qx"}
                x = xfx;

                if ~numel(tri_grid_fx)
                  tri_grid_fx = delaunay(xfx(1, :), xfx(2, :));
                endif

                tri_grid = tri_grid_fx;
              case {"qz", "mdotz", "Qz"}
                x = xfz;

                if ~numel(tri_grid_fz)
                  tri_grid_fz = delaunay(xfz(1, :), xfz(2, :));
                endif

                tri_grid = tri_grid_fz;
              otherwise
                x = xh;

                if ~numel(tri_grid_h)
                  tri_grid_h = delaunay(xh(1, :), xh(2, :));
                endif

                tri_grid = tri_grid_h;
            endswitch

            z = griddata_prepared(x(1, :), ...
                                  x(2, :), ...
                                  zt_n(options.plot_range(j), :), ...
                                  bearings(options.plot_bearings(i)).xi, ...
                                  bearings(options.plot_bearings(i)).zi, ...
                                  tri_grid);

            ++ifig;
            colormap("jet");

            if (isfield(options, [fields(k).column, "_plot_range"]))
              switch (options.plot_func)
                case {"contourf", "contour"}
                  z_range = {getfield(options, [fields(k).column, "_plot_range"])};
                otherwise
                  z_range = {};
              endswitch
            else
              z_range = {};
            endif

            xi = bearings(options.plot_bearings(i)).xi;
            zi = bearings(options.plot_bearings(i)).zi;

            switch (options.plot_func)
              case "deformed shape"
                phii = 2 * xi / dm;
                scale = 0.25 * dm / zt_mag;
                ri = 0.5 * dm + scale * z;
                surf(ri .* cos(phii), ri .* sin(phii), zi);
              otherwise
                feval(options.plot_func, xi, zi, z, z_range{:});
            endswitch

            if (isfield(options, "z_lim"))
              xlim(options.z_lim);
            endif

            if (isfield(options, "x_lim"))
              ylim(options.x_lim);
            endif

            switch (options.plot_func)
              case {"surf", "surfl", "surfc"}
                daspect([1, 1, daspect()(3)]);
              otherwise
                daspect([1, 1, 1]);
            endswitch

            xlabel('x [m]');
            ylabel('z [m]');

            switch (options.plot_func)
              case {"contour", "contourf"}
                colorbar();
              otherwise
                zlabel(sprintf('%s [%s]', fields(k).column, fields(k).unit));
            endswitch

            if (isfield(options, "plot_parameter") && isfield(options, "plot_parameter_format"))
              param = sprintf(options.plot_parameter_format, options.plot_parameter(options.plot_range(j)));
            else
              param = "";
            endif

            grid on;
            grid minor on;
            title(sprintf('%s %s %s [%s] of bearing %d (%s) at t=%g %s', ...
                          options.plot_title, ...
                          fields(k).description, ...
                          printable_title(fields(k).column), ...
                          fields(k).unit, ...
                          options.plot_bearings(i), ...
                          options.bearing_title{options.plot_bearings(i)}, ...
                          t(options.plot_range(j)), ...
                          param));

            if (length(options.print_to_file) > 0)
              file_name_prefix = sprintf("%s_%04d", options.print_to_file, ifig);
              file_name_pdf = [file_name_prefix, ".pdf"];

              unlink(file_name_pdf);

              print("-dpdf", "-color", "-portrait", sprintf("-S%d,%d", options.width, options.height), file_name_pdf);

              [info, err, msg] = stat(file_name_pdf);

              if (err ~= 0 || info.size == 0)
                error("failed to create file \"%s\"", file_name_pdf);
              endif

              file_names_pdf{end + 1} = file_name_pdf;
            endif
          endfor
        endif
      endif
    endfor
  endfor

  if (length(file_names_pdf) > 0)
    [output_dir, output_name, output_ext] = fileparts(options.print_to_file);
    pdf_merge(file_names_pdf, fullfile(output_dir, [output_name, ".pdf"]));
  endif
endfunction

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
%! param.d1 = 15e-3;
%! param.Psi = 0.8e-3;
%! param.M = int32(20);
%! param.N = int32(75);
%! param.output_bearing_data = true;
%! param.B = 8e-3;
%! param.epsilon = 0.8;
%! fd = -1;
%! output_file = "";
%! unwind_protect
%!   unwind_protect
%!     output_dir = tempdir();
%!     [fd, output_file] = mkstemp(fullfile(output_dir, "oct-mbdyn_post_ehd_load_output_XXXXXX"));
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
%!     fputs(fd, "set: real n = 2.75;\n");
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
%!     fputs(fd, "        enforce constraint equations: constraint violations, scale factor, 1e5;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "        method: bdf;\n");
%!     fputs(fd, "        output: messages;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-4;\n");
%!     fputs(fd, "        derivatives max iterations: 10;\n");
%!     fputs(fd, "        derivatives coefficient: 110;\n");
%!     fputs(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "        nonlinear solver: nox,\n");
%!     fputs(fd, "             jacobian operator, newton,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 100,\n");
%!     fputs(fd, "             krylov subspace size, 100;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "    use automatic differentiation;\n");
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
%!   beta_ref = pi / 180 * interp2(B_d_r, epsilon_r, beta_r, B / d, epsilon, "linear");
%!   mu_ref = Psi * (abs((omega1z - omega2z) / omega_res) * pi / (sqrt(1 - epsilon^2) * So_ref) + sin(beta_ref) * abs(epsilon) / 2);
%!   Q_ref = interp2(B_d_r, epsilon_r, Q_r, B / d, epsilon, "linear");
%!   opt_plot.plot_title = "hydrodynamic journal plain bearing";
%!   opt_plot.plot_func = "contourf";
%!   opt_plot.plot_range = numel(res.t);
%!   opt_plot.plot_bearings = 1;
%!   opt_plot.f_plot_p = true;
%!   opt_plot.f_plot_rho = true;
%!   opt_plot.f_plot_h = true;
%!   mbdyn_post_ehd_plot_nodal(res.t, res.bearings, opt_plot, res.log_dat);
%!   figure_list();
%!   fprintf(stderr, "So / So_ref - 1 = %.2f\n", So / So_ref - 1);
%!   fprintf(stderr, "beta / beta_ref - 1 = %.2f\n", beta / beta_ref - 1);
%!   fprintf(stderr, "mu / mu_ref - 1 = %.2f\n", mu / mu_ref - 1);
%!   fprintf(stderr, "Q / Q_ref - 1 = %.2f\n", Q / Q_ref - 1);
%!   assert_simple(So, So_ref, 0.03 * So_ref);
%!   assert_simple(beta, beta_ref, 0.02 * beta_ref);
%!   assert_simple(mu, mu_ref, 0.03 * mu_ref);
%!   assert_simple(Q, Q_ref, 0.07 * Q_ref);
%!   assert_simple(dQ, 0, 1e-3);
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
