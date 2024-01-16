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
## @deftypefn {Function File} mbdyn_post_plot_model(@var{output_file}, @var{res}, @var{idx_t})
##
## Plot the deformed shape of an MBDyn model (currently only beam elements are plotted).
##
## @end deftypefn

function mbdyn_post_plot_model(output_file, res, idx_t)
  if (nargin < 2 || nargin > 3 || nargout > 0)
    print_usage();
  endif

  if (nargin < 3)
    idx_t = 1:numel(res.t);
  endif

  R = mbdyn_post_angles_to_rotation_mat([res.log_dat.nodes.label], res, res.log_dat);

  X = [res.log_dat.nodes.X0];

  lim = zeros(3, 2);

  for i=1:3
    lim(i, 1) = min(X(i, :));
    lim(i, 2) = max(X(i, :));
  endfor

  node_ids = [res.log_dat.nodes.label];

  for k=1:numel(idx_t)
    for i=1:numel(res.log_dat.beams3)
      Xi = zeros(3, numel(res.log_dat.beams3(i).nodes));

      for j=1:numel(res.log_dat.beams3(i).nodes)
        node_idx = find(node_ids == res.log_dat.beams3(i).nodes(j).label);
        Xi(:, j) = res.trajectory{node_idx}(idx_t(k), 1:3).' + R{node_idx}(:, :, idx_t(k)) * res.log_dat.beams3(i).nodes(j).offset(:);
      endfor

      for j=1:3
        lim(j, 1) = min(lim(j, 1), min(Xi(j, :)));
        lim(j, 2) = max(lim(j, 2), max(Xi(j, :)));
      endfor
    endfor
  endfor

  hfig = 0;

  unwind_protect
    hfig = figure("visible", "off");

    for k=1:numel(idx_t)
      clf;
      hold on;

      for i=1:numel(res.log_dat.beams3)
        Xi = zeros(3, numel(res.log_dat.beams3(i).nodes));

        for j=1:numel(res.log_dat.beams3(i).nodes)
          node_idx = find(node_ids == res.log_dat.beams3(i).nodes(j).label);
          Xi(:, j) = res.trajectory{node_idx}(idx_t(k), 1:3).' + R{node_idx}(:, :, idx_t(k)) * res.log_dat.beams3(i).nodes(j).offset(:);
        endfor

        line("xdata", Xi(1, :), "ydata", Xi(2, :), "zdata", Xi(3, :));
      endfor

      daspect(ones(1,3));
      view(30, 30);

      if (lim(1, 2) > lim(1, 1))
        xlim(lim(1, :));
      endif

      if (lim(2, 2) > lim(2, 1))
        ylim(lim(2, :));
      endif

      if (lim(3, 2) > lim(3, 1))
        zlim(lim(3, :));
      endif

      xlabel("x");
      ylabel("y");
      zlabel("z");
      print(hfig, "-dpng", sprintf("%s_%03d.png", output_file, k));
    endfor

    pid = spawn("ffmpeg", {"-i", [output_file, "_%03d.png"], "-y", [output_file, ".mpg"]});

    status = spawn_wait(pid);

    if (status ~= 0)
      warning("ffmpeg failed with status %d", status);
    endif
  unwind_protect_cleanup
    if (hfig ~= 0)
      close(hfig);
    endif
  end_unwind_protect
endfunction
