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
## @deftypefn {Function File} [@var{mesh}] = mbdyn_post_ehd_create_mesh(@var{log_dat})
## Build a 3D mesh from hydrodynamic plain bearing data for export to Gmsh.
##
## @var{log_dat} @dots{} Return value from mbdyn_post_load_log.
##
## @seealso{mbdyn_post_ehd_export_mesh, mbdyn_post_ehd_export_data, mbdyn_post_load_log}
## @end deftypefn

function mesh = mbdyn_post_ehd_create_mesh(log_dat)
  if (nargin ~= 1 || nargout > 1)
    print_usage();
  endif
  
  mesh.noffset = zeros(numel(log_dat.bearings), 1, "int32");
  mesh.eoffset = zeros(numel(log_dat.bearings), 1, "int32");
  
  noffset = int32(0);
  eoffset = int32(0);
  
  for i=1:numel(log_dat.bearings)
    mesh.noffset(i) = noffset;
    mesh.eoffset(i) = eoffset;
    
    inode = int32(0);
    mesh.bearings(i).cylindrical = log_dat.bearings(i).cylindrical;
    mesh.bearings(i).nodeidx = zeros(numel(log_dat.bearings(i).nodes), 1, "int32");
    mesh.bearings(i).label = log_dat.bearings(i).label;
    
    for j=1:numel(log_dat.bearings(i).nodes)
      switch (log_dat.bearings(i).nodes(j).type)
	case "hydro"
	  mesh.bearings(i).nodeidx(++inode) = j;
	  ++noffset;
      endswitch
    endfor

    mesh.bearings(i).nodeidx = mesh.bearings(i).nodeidx(1:inode);

    mesh.bearings(i).elemidx = zeros(numel(log_dat.bearings(i).elements), 1, "int32");
    mesh.bearings(i).elements = zeros(numel(log_dat.bearings(i).elements), 1, "int32");
    ielem = int32(0);
    
    for j=1:numel(log_dat.bearings(i).elements)
      switch (numel(log_dat.bearings(i).elements(j).nodes))
	case 4
	  mesh.bearings(i).elemidx(++ielem) = j;
	  mesh.bearings(i).elements(ielem) = ++eoffset;
      endswitch
    endfor

    mesh.bearings(i).elemidx = mesh.bearings(i).elemidx(1:ielem);
    mesh.bearings(i).elements = mesh.bearings(i).elements(1:ielem);
  endfor

  mesh.nodes = zeros(noffset, 3);

  for i=1:numel(log_dat.bearings)
    dm = log_dat.bearings(i).cylindrical.dm;
    x = [log_dat.bearings(i).nodes(mesh.bearings(i).nodeidx).x];
    
    mesh.bearings(i).Phi = x(1, :) / (0.5 * log_dat.bearings(i).cylindrical.dm);
    mesh.bearings(i).z = x(2, :);
    
    switch (log_dat.bearings(i).cylindrical.mesh_pos)
      case "journal"
	imnode = 1;
      case "shell"
	imnode = 2;
      otherwise
	error("invalid mesh position found in bearing data");
    endswitch
    
    idxnode = find([log_dat.nodes.label] == log_dat.bearings(i).cylindrical.nodes(imnode).label);
    R0 = log_dat.nodes(idxnode).R0;
    X0 = log_dat.nodes(idxnode).X0;
    Rb = log_dat.bearings(i).cylindrical.nodes(imnode).Rb;
    ob = log_dat.bearings(i).cylindrical.nodes(imnode).o;
    mesh.bearings(i).Xn = [(0.5 * dm) * cos(mesh.bearings(i).Phi);
			   (0.5 * dm) * sin(mesh.bearings(i).Phi);
			   mesh.bearings(i).z];

    mesh.bearings(i).A = [cos(mesh.bearings(i).Phi(:)), ...
			  sin(mesh.bearings(i).Phi(:)), ...
			  -mesh.bearings(i).z(:) .* sin(mesh.bearings(i).Phi(:)), ...
			  mesh.bearings(i).z(:) .* cos(mesh.bearings(i).Phi(:))];
    
    X = R0 * (Rb * mesh.bearings(i).Xn + ob) + X0;
    
    mesh.nodes(mesh.noffset(i) + (1:numel(mesh.bearings(i).nodeidx)), :) = X.';
  endfor

  mesh.elements.quad4 = zeros(eoffset, 4, "int32");

  for i=1:numel(log_dat.bearings)
    for j=1:numel(mesh.bearings(i).elements)
      mesh.elements.quad4(mesh.eoffset(i) + j, :) = log_dat.bearings(i).elements(mesh.bearings(i).elemidx(j)).nodes + mesh.noffset(i);
    endfor
  endfor
endfunction
