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
## @deftypefn {Function File} [@var{inertia}] = mbdyn_post_inertia_compute(@var{body_groups}, @var{bodies}, @var{nodes})
## @deftypefnx {} [@var{inertia}] = mbdyn_post_inertia_compute(@var{body_groups}, @var{mbdyn_filename})
##
## Computes the mass, center of gravity and momentum of inertia of groups of bodies.
##
## @var{body_groups}.name @dots{} Name of the group of bodies
##
## @var{body_groups}.labels @dots{} Array of labels of all bodies which belong to this group
##
## @var{bodies} @dots{} Data structure returned from mbdyn_post_load_log_body
##
## @var{nodes} @dots{} Data structure returned from mbdyn_post_load_log
##
## @var{mbdyn_filename} @dots{} Name of mbdyn output files
##
## @end deftypefn

function [inertia] = mbdyn_post_inertia_compute(body_groups = "all", varargin)
  if (nargin < 2 || nargin > 3 || nargout > 1)
    print_usage();
  endif

  inertia = struct();

  if (~isstruct(body_groups) && ~ischar(body_groups))
    error("body_groups must be a cell array or a string!");
  endif

  if (length(varargin) == 1 && ischar(varargin{1}))
    mbdyn_filename = varargin{1};
    bodies = mbdyn_post_load_log_body(mbdyn_filename);
    nodes = mbdyn_post_load_log_node(mbdyn_filename);
  elseif (length(varargin) == 2 && isstruct(varargin{1}) && isstruct(varargin{2}))
    bodies = varargin{1};
    nodes = varargin{2};
  else
    print_usage();
    return;
  endif

  N = length(bodies);

  if (ischar(body_groups))
    switch(body_groups)
      case 'all'
        body_groups = struct();
        body_groups(1).labels =  [ bodies.label ];
        body_groups(1).name = 'all';
      case 'every'
        body_groups = struct();
        for i=1:length(bodies)
          body_groups(i).labels = [bodies(i).label];
          body_groups(i).name = sprintf('body: %d',bodies(i).label);
        endfor
      otherwise
        error("invalid argument: body_groups=\"%s\"",body_groups);
    endswitch
  endif

  for i=1:length(body_groups)

    inertia(i).name = body_groups(i).name;
    inertia(i).dm = 0;
    inertia(i).Xgc = zeros(1,3);
    inertia(i).J = zeros(3,3);
    Mgc = zeros(3,1);

    for j=1:length(body_groups(i).labels)
      idx_body = find([bodies.label] == body_groups(i).labels(j));

      if (length(idx_body) ~= 1)
        error("bodies.label is not unique");
      endif

      if (length(idx_body) == 0)
        error("body %d not found!",body_groups(i).labels(j));
      endif

      idx_node = find([nodes.label] == bodies(idx_body).node);

      if (length(idx_body) ~= 1)
        error("body label not found or not unique!");
      endif

      if (length(idx_node)  ~= 1)
        error("node label not found or not unique!");
      endif

      m_j = bodies(idx_body).dm;
      J_j = bodies(idx_body).J;
      R_j = nodes(idx_node).R0;
      dXgc_j = R_j * bodies(idx_body).Xgc;
      Xgc_j = dXgc_j + nodes(idx_node).X0;

      Mgc +=  Xgc_j * m_j;
      inertia(i).dm += m_j;
      inertia(i).J += R_j * J_j * R_j.' + (skew(dXgc_j) * skew(dXgc_j) - skew(Xgc_j) * skew(Xgc_j)) * m_j;
    endfor
    inertia(i).Xgc = Mgc / inertia(i).dm;
    inertia(i).J += skew(inertia(i).Xgc) * skew(inertia(i).Xgc) * inertia(i).dm;
    inertia(i).bodies = sort(body_groups(i).labels);
  endfor
endfunction
