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
## @deftypefn {Function File} [@var{mesh}, @var{sol}] = mbdyn_post_export_model(@var{res}, @var{idx_t})
##
## Convert the output data from MBDyn to Finite Element model data suiteable for processing with mboct-fem-pkg
##
## @end deftypefn

function [mesh, sol] = mbdyn_post_export_model(res, idx_t)
  if (nargin < 1 || nargin > 2 || nargout > 2)
    print_usage();
  endif

  if (nargin < 2)
    idx_t = 1:numel(res.t);
  endif

  R = mbdyn_post_angles_to_rotation_mat([res.log_dat.nodes.label], res, res.log_dat, idx_t);

  elem_types_fe = {"beam2", "beam2", "beam3"};
  elem_types_mb = {"joints", "beams2", "beams3"};
  elem_node_num = int32([2, 2, 3]);
  elem_node_order = {[1:2], [1:2], [1,3,2]};
  
  num_aux_nodes = int32(0);

  for i=1:numel(elem_types_mb)
    num_aux_nodes += elem_node_num(i) * numel(getfield(res.log_dat, elem_types_mb{i}));
  endfor

  mesh.nodes = zeros(numel(res.log_dat.nodes) + num_aux_nodes, 6);

  mesh.nodes(1:numel(res.log_dat.nodes), :) = [[res.log_dat.nodes.X0].', [res.log_dat.nodes.Phi0].'];

  iaux_node = zeros(1, numel(elem_types_mb), "int32");
  
  mesh.elements = struct();
  
  for k=1:numel(elem_types_mb)
    elem_k = getfield(res.log_dat, elem_types_mb{k});
    elnodes = struct("nodes", cell(1, numel(elem_k)));
    
    for i=1:numel(elem_k)
      Xi = zeros(3, numel(elem_k(i).nodes));

      for j=1:numel(elem_k(i).nodes)
	node_idx = find([res.log_dat.nodes.label] == elem_k(i).nodes(j).label);

	if (numel(node_idx) ~= 1)
	  error("invalid node label %d", elem_k(i).nodes(j).label);
	endif
	
	Xi(:, j) = res.log_dat.nodes(node_idx).X0 + res.log_dat.nodes(node_idx).R0 * elem_k(i).nodes(j).offset(:);
      endfor

      node_id = numel(res.log_dat.nodes) + sum(iaux_node(1:k)) + (1:columns(Xi));
      mesh.nodes(node_id, 1:3) = Xi.';
      elnodes(i).nodes = node_id(elem_node_order{k});
      iaux_node(k) += columns(Xi);
    endfor

    if (isfield(mesh.elements, elem_types_fe{k}))
      elnodesprev = getfield(mesh.elements, elem_types_fe{k});
      elnodesprev(end + (1:numel(elnodes))) = elnodes;
    else
      elnodesprev = elnodes;
    endif
    
    mesh.elements = setfield(mesh.elements, elem_types_fe{k}, elnodesprev);
  endfor

  iaux_node = zeros(1, numel(elem_types_mb), "int32");
  
  for k=1:numel(elem_types_mb)
    elem_k = getfield(res.log_dat, elem_types_mb{k});
    elnodes = struct("nodes", cell(1, 2 * numel(elem_k)));
    ielem = int32(0);
    
    for i=1:numel(elem_k)
      node_id = numel(res.log_dat.nodes) + sum(iaux_node(1:k)) + (1:numel(elem_k(i).nodes));
      
      for j=[1, numel(elem_k(i).nodes)]
	node_idx = find([res.log_dat.nodes.label] == elem_k(i).nodes(j).label);

	if (numel(node_idx) ~= 1)
	  error("invalid node label %d", elem_k(i).nodes(j).label);
	endif
	
	elnodes(++ielem).nodes = int32([node_id(j), node_idx]);
      endfor      

      iaux_node(k) += numel(elem_k(i).nodes);
    endfor

    mesh.elements.beam2(end + (1:numel(elnodes))) = elnodes;
  endfor

  if (isfield(mesh.elements, "beam3"))
    beam3degen = struct("nodes", cell(1, 2 * numel(mesh.elements.beam3)));
    ielem = int32(0);
    for i=1:numel(mesh.elements.beam3)
      beam3degen(++ielem).nodes = mesh.elements.beam3(i).nodes([1, 3]);
      beam3degen(++ielem).nodes = mesh.elements.beam3(i).nodes([3, 2]);
    endfor
    mesh.elements.beam2(end + (1:numel(beam3degen))) = beam3degen;
    mesh.elements = rmfield(mesh.elements, "beam3");
  endif
  
  sol.t = res.t(idx_t);
  sol.def = zeros(rows(mesh.nodes), columns(mesh.nodes), numel(idx_t));

  for i=1:numel(res.log_dat.nodes)
    node_idx = find(res.node_id == res.log_dat.nodes(i).label);
    if (~isempty(node_idx))
      if (numel(node_idx) ~= 1)
	error("invalid node label %d", res.log_dat.nodes(i).label);
      endif
      
      for j=1:3
	sol.def(node_idx, j, :) = res.trajectory{node_idx}(idx_t, j) - res.log_dat.nodes(i).X0(j);
	sol.def(node_idx, j + 3, :) = res.trajectory{node_idx}(idx_t, j + 3) - res.log_dat.nodes(i).Phi0(j);
      endfor
    endif
  endfor
  
  iaux_node = zeros(1, numel(elem_types_mb), "int32");
  
  for k=1:numel(elem_types_mb)
    elem_k = getfield(res.log_dat, elem_types_mb{k});
    
    for i=1:numel(elem_k)
      node_id = numel(res.log_dat.nodes) + sum(iaux_node(1:k)) + (1:numel(elem_k(i).nodes));
      
      for j=1:numel(elem_k(i).nodes)
        node_idx = find(res.node_id == elem_k(i).nodes(j).label);

	if (~isempty(node_idx))
	  if (numel(node_idx) ~= 1)
	    error("invalid node label %d", res.log_dat.nodes(i).label);
	  endif
	  
	  Xj = res.trajectory{node_idx}(idx_t, 1:3).';
	  Phij = res.trajectory{node_idx}(idx_t, 4:6).';
	  Rj = R{node_idx};

	  for l=1:3
            sol.def(node_id(j), l, :) = Xj(l, :) - mesh.nodes(node_id(j), l);
	    sol.def(node_id(j), l + 3, :) = Phij(l, :) - mesh.nodes(node_id(j), l + 3);
	    for m=1:3
	      sol.def(node_id(j), l, :) += Rj(l, m, :) * elem_k(i).nodes(j).offset(m);
	    endfor
	  endfor
	endif
      endfor
      
      iaux_node(k) += numel(elem_k(i).nodes);
    endfor
  endfor
endfunction

%!demo
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 20;
%! l = 1000e-3;
%! h = 500e-3;
%! g = 9.81*1000;
%! D = 20e-3;
%! A = D^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = D^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! X = [linspace(0, l, N);
%!      zeros(2, N)];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! options.first_node_number = 3;
%! options.start_node = "1";
%! options.end_node = "2";
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_post_plot_model_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real g = %g;\n", g);
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fprintf(fd, " set: real l = %g;\n", l);
%!     fprintf(fd, " set: real h = %g;\n", h);
%!     fprintf(fd, " set: integer N = 10000;\n");
%!     fputs(fd, "set: real t1 = N;\n");
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: t1;\n");
%!     fputs(fd, "         time step: t1 / N;\n");
%!     fputs(fd, "         linear solver: umfpack, colamd, scale, row max column max, always, max iterations, 10;\n");
%!     fputs(fd, "         method: implicit euler;\n");
%!     fputs(fd, "         max iterations: 100;\n");
%!     fputs(fd, "         tolerance: 1.e-4, 1e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         output: iterations, solver condition number, stat, yes;\n");
%!     fputs(fd, "         nonlinear solver: nox, modified, 10,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 12,\n");
%!     fputs(fd, "             krylov subspace size, 12;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");
%!     fputs(fd, "       model: static;\n");
%!     fputs(fd, "       output meter: closest next, 0., forever, const, t1 / 10;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, " structural: 1, static, position, 0., 0., h, orientation, eye, velocity, null, angular velocity, null;\n");
%!     fputs(fd, " structural: 2, static, position, l, 0., h, orientation, eye, velocity, null, angular velocity, null;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, " joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., mult, time, g / t1;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   if (f_print_input_file)
%!     spawn_wait(spawn("cat", {fname}));
%!   endif
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   [res.t, res.trajectory, res.deformation, res.velocity, res.acceleration, res.node_id] = mbdyn_post_load_output_struct(options.output_file);
%!   res.log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   [mesh, sol] = mbdyn_post_export_model(res);
%!   if (~isempty(pkg("list", "mboct-fem-pkg")))
%!     pkg load mboct-fem-pkg;
%!     opts.print_and_exit = true;
%!     opts.print_to_file = fname;
%!     opts.rotation_angle = [pi/2, 0, 0];
%!     opts.skin_only = false;
%!     opts.show_element = false;
%!     opts.animation_delay = 0;
%!     opts.scale_def = 1;
%!     fem_post_sol_external(mesh, sol, opts);
%!     fn = dir([opts.print_to_file, "*.jpg"]);
%!     for i=1:numel(fn)
%!       [img, map, alpha] = imread(fullfile(fn(i).folder, fn(i).name));
%!       figure("visible", "off");
%!       imshow(img, map);
%!       title(sprintf("deformed mesh step %d", i));
%!     endfor
%!   endif
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
