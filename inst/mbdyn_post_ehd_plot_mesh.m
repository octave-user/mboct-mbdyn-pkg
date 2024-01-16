## Copyright (C) 2018(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} mbdyn_post_ehd_plot_mesh(@var{bearing})
## @deftypefnx {} mbdyn_post_ehd_plot_mesh(@var{bearing}, @var{value_n})
## @deftypefnx {} mbdyn_post_ehd_plot_mesh(@var{bearing}, @var{value_n}, @var{idx_t})
## Plot hydrodynamic bearing data returned from mbdyn_post_ehd_load_output.
##
## @var{bearing} @dots{} Return value from mbdyn_post_ehd_load_output.
##
## @var{value_n} @dots{} Nodal field to be plotted (e.g. bearings(i).columns.p_n, bearings(i).columns.h_n, ...)
##
## @var{idx_t} @dots{} Time step index to be plotted (e.g. bearings(i).columns.p_n(idx_t, :)).
##
## @seealso{mbdyn_post_ehd_load_output, mbdyn_post_ehd_plot_nodal}
## @end deftypefn

function mbdyn_post_ehd_plot_mesh(bearing, value_n, idx_t)
  if (nargin < 1 || nargin > 3 || nargout > 0)
    print_usage();
  endif

  if (nargin < 2)
    value_n = zeros(1, columns(bearing.nodes));
  endif

  if (nargin < 3)
    idx_t = 1;
  endif

  switch (bearing.type)
    case "generic"
      hfig = figure("visible", "off");
      hold on;
      colormap jet;
      colors = get(hfig, "colormap");

      num_faces = int32(0);
      num_vertices = int32(0);

      for j=1:numel(bearing.elements)
        switch (numel(bearing.elements(j).nodes))
          case 4
            curr_faces = 2;
          case 9
            curr_faces = 8;
          otherwise
            continue
        endswitch

        num_faces += curr_faces;
        num_vertices += numel(bearing.elements(j).nodes);
      endfor

      face_data = zeros(num_faces, 3, "int32");
      vertex_data = zeros(num_vertices, 3);
      color_data = zeros(num_faces, 3);

      num_faces = int32(0);
      num_vertices = int32(0);

      for j=1:numel(bearing.elements)
        x = [bearing.nodes(bearing.elements(j).nodes).x];

        r = 0.5 * bearing.cylindrical.dm;

        Phi = x(1, :) / r;

        if (nargin < 2)
          dr = 0;
        else
          dr = (value_n(idx_t, bearing.elements(j).nodes) - min(min(value_n))) / (max(max(value_n)) - min(min(value_n))) * r;
        endif

        X = [(r + dr) .* cos(Phi);
             (r + dr) .* sin(Phi);
             x(2, :)];

        if (nargin < 2)
          color_val = ones(1, numel(bearing.elements(j).nodes));
        else
          color_val = (value_n(idx_t, bearing.elements(j).nodes) - min(min(value_n))) / (max(max(value_n)) - min(min(value_n)));
        endif

        switch (columns(x))
          case 4
            faces = int32([1, 2, 3;
                           1, 3, 4]);
          case 9
            faces = int32([1, 5, 9;
                           1, 9, 8;
                           5, 2, 6;
                           5, 6, 9;
                           9, 6, 3;
                           9, 3, 7;
                           8, 9, 7;
                           8, 7, 4]);
          otherwise
            continue;
        endswitch

        color_faces = colors(floor(mean(color_val(faces), 2) * (rows(colors) - 1)) + 1, :);
        face_data(num_faces + (1:rows(faces)), :) = faces + num_vertices;
        vertex_data(num_vertices + (1:columns(X)), :) = X.';
        color_data(num_faces + (1:rows(faces)), :) = color_faces;
        num_faces += rows(faces);
        num_vertices += columns(X);
      endfor

      patch('Faces', face_data, ...
            'Vertices', vertex_data, ...
            'FaceVertexCData', color_data, ...
            'FaceColor', 'interp');

      xlabel("x");
      ylabel("y");
      zlabel("z");
      title(sprintf("bearing number %d", bearing.label));
      grid on;
      grid minor on;
      daspect([1, 1, 1]);
  endswitch
endfunction
