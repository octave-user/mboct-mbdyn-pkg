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
## @deftypefn {Function File} mbdyn_pre_beam_plot(@var{beam}, @var{options})
##
## Generate a plot of a curved beam.
##
## @var{beam} @dots{} Return value from mbdyn_pre_beam_compute.
##
## @var{options}.Rn @dots{} If true, plot node orientations.
##
## @var{options}.Rg @dots{} If true, plot Gauss point orientations.
##
## @var{options}.X @dots{} If true, plot node positions.
##
## @var{options}.s @dots{} Size factor for rotation matrices.
##
## @var{options}.figure @dots{} Figure handle for plotting.
##
## @end deftypefn

function mbdyn_pre_beam_plot(beam, options)
  if (nargin < 1 || nargin > 2)
    print_usage();
  endif

  if (nargin < 2)
    options = struct();
  endif

  if (~isstruct(beam))
    error("beam must be a struct");
  endif

  if (~isfield(beam, "beams"))
    error("field beam.beams does not exist!");
  endif

  if (~isfield(options, "s"))
    options.s = 0.3;
  endif

  if (~isfield(options, "Rn"))
    options.Rn = false;
  endif

  if (~isfield(options, "Rg"))
    options.Rg = false;
  endif

  if (~isfield(options, "X"))
    options.X = true;
  endif

  if (isfield(options, "figure"))
    figure(options.figure, "visible", get(options.figure, "visible"));
  else
    figure("visible", "off");
  endif

  hold on;
  color_map_Xn = rainbow(length(beam.beams));
  color_map_Rn = eye(3);
  color_map_Rg = 0.8 * color_map_Rn;
  for i=1:length(beam.beams)
    Xn = mbdyn_pre_beam_node_position(beam, i);
    Xg = mbdyn_pre_beam_gauss_point_position(beam, i);
    Rn = mbdyn_pre_beam_node_orientation(beam, i);
    Rg = mbdyn_pre_beam_gauss_point_orientation(beam, i);

    if (options.X)
      X = [ Xn(:, 1),  Xg(:, 1),  Xn(:, 2),  Xg(:, 2),  Xn(:, 3) ];
      set(plot3(X(1, :), X(2, :), X(3, :), sprintf('x-;b%d;', i)), 'color', color_map_Xn(i, :));
    endif

    if (options.Rn || options.Rg)
      ds = options.s * norm(Xn(:, 1)-Xn(:, 3));

      if (options.Rn)
        for j=1:3
          text(Xn(1, j), Xn(2, j), Xn(3, j), sprintf('n%d', j));
          for k=1:3
            e_k = [ Xn(:, j),  Xn(:, j) + ds * Rn(:, k, j) ];
            set(plot3(e_k(1, :),  e_k(2, :),  e_k(3, :), '-'), 'color', color_map_Rn(k, :));
          endfor
        endfor
      endif

      if (options.Rg)
        for j=1:2
          text(Xg(1, j), Xg(2, j), Xg(3, j), sprintf('g%d', j));
          for k=1:3
            e_k = [ Xg(:, j),  Xg(:, j) + ds * Rg(:, k, j) ];
            set(plot3(e_k(1, :),  e_k(2, :),  e_k(3, :), '-'), 'color', color_map_Rg(k, :));
          endfor
        endfor
      endif
    endif
  endfor
  set(gca(), 'dataaspectratio', [1, 1, 1]);
  xlabel('x [m]');
  ylabel('y [m]');
  zlabel('z [m]');
  grid on;
  grid minor on;
  title('curved beam');
endfunction
