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
