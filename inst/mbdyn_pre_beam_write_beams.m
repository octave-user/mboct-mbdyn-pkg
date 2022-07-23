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
## @deftypefn {Function File} mbdyn_pre_beam_write_beams(@var{beam}, @var{output_file}, @var{options})
##
## Generate an MBDyn input file <@var{output_file}> containing all beam elements of beam model <@var{beam}>.
##
## @var{beam} @dots{} Return value from mbdyn_pre_beam_compute.
##
## @var{output_file} @dots{} If <@var{output_file}> is a string, open a new file using <@var{output_file}> as name.
## If <@var{output_file}> it is a file descriptor, write the output to that file descriptor.
##
## @var{open_mode} @dots{} If <@var{output_file}> is a string, pass <@var{open_mode}> to fopen (e.g. "wt", "at").
##
## @var{options}.first_reference_frame_number @dots{} String or integer number of the first reference frame of the beam.
##
## @var{options}.first_node_number @dots{} String or integer number of the first node of the beam.
##
## @var{options}.first_beam_number @dots{} String or integer number of the first beam element.
##
## @var{options}.constitutive_law_number @dots{} String name of an integer variable holding
## the id of the constitutive law for the complete beam.
##
## @var{options}.id_offset @dots{} Integer number to be added to reference frame-, node- and beam element- and drive caller-numbers.
##
## @var{options}.private_data_output @dots{} Cell array of string names for private data which should be recorded for each beam element.
##
## @var{options}.private_data_output_flag @dots{} String name of a boolean variable which enables/disables output of private data.
##
## @var{options}.private_data_first_drive_idx @dots{} String name of an integer variable holding the id of the first drive caller.
##
## @var{options}.output_flag @dots{} String name of a boolean variable which enables/disables output for the beam.
##
## @var{options}.output_Vi @dots{} Boolean value which enables output of the total strain energy for all beam elements.
##
## @var{options}.Vi_first_drv_idx @dots{} String name of an integer variable which holds the index of the drive caller
## for output of strain energy.
##
## @var{options}.reorient_cross_section @dots{} If this flag is enabled, the variation of the cross section orientation
## within one beam element will be minimized.
##
## @var{options}.start_node @dots{} If the beginning of the beam should be connected to an existing node, then <@var{start_node}> is
## the string name of an integer variable holding the existing node number.
##
## @var{options}.end_node @dots{} If the end of the beam should be connected to an existing node, then <@var{end_node}> is the
## string name of an integer variable holding the existing node number.
##
## @var{options}.weak_spring.stiffness @dots{} Fix all nodes of the beam to the ground by spring support elements with a total
## stiffness of <@var{stiffness}>.
##
## @var{options}.weak_spring.initial_time @dots{} Time when the spring support elements are activated.
##
## @var{options}.weak_spring.tau @dots{} Time duration the spring support elements exist. The stiffness will be reduced from
## <@var{stiffness}> at <@var{initial_time}> to zero at <@var{initial_time}> + <@var{tau}>.
##
## @end deftypefn

function mbdyn_pre_beam_write_beams(beam, output_file, options)
  if (nargin < 2)
    print_usage();
  endif

  if (nargin < 3)
    options = struct();
  endif

  if (~isfield(options, "first_reference_frame_number"))
    options.first_reference_frame_number = 1;
  endif

  if (~isfield(options, "first_node_number"))
    options.first_node_number = 1;
  endif

  if (~isfield(options, "first_beam_number"))
    options.first_beam_number = 1;
  endif

  if (~isfield(options, "constitutive_law_number"))
    error("constitutive_law_number is missing in options struct!");
  endif

  if (~isfield(options, "id_offset"))
    options.id_offset = 0;
  endif    

  if (~isfield(options, "open_mode"))
    options.open_mode = "wt";
  endif

  if (~isfield(options, "private_data_output"))
    options.private_data_output = {};
  endif

  if (~isfield(options, "private_data_output_flag"))
    options.private_data_output_flag = "yes";
  endif
  
  if (~isfield(options, "private_data_first_drv_idx"))
    options.private_data_first_drv_idx = {""};
  endif

  if (~isfield(options, "output_flag"))
    options.output_flag = "default";
  endif

  if (~ischar(options.first_reference_frame_number))
    options.first_reference_frame_number = sprintf("%d", options.first_reference_frame_number);
  endif

  if (~ischar(options.first_node_number))
    options.first_node_number = sprintf("%d", options.first_node_number);
  endif

  if (~ischar(options.first_beam_number))
    options.first_beam_number = sprintf("%d", options.first_beam_number);
  endif

  if (~isfield(options, "reorient_cross_sections"))
    options.reorient_cross_sections = true;
  endif

  if (~isfield(options, "output_Vi"))
    options.output_Vi = false;
  endif

  if (options.output_Vi && ~isfield(options, "Vi_first_drv_idx"))
    error("output of potential energy requested but parameter Vi_first_drv_idx is missing");
  endif

  fout = -1;
  owns_fd = false;

  unwind_protect
    if (ischar(output_file))
      owns_fd = true;

      [fout, msg] = fopen(output_file, options.open_mode);

      if (fout == -1)
        error("could not open file \"%s\": %s", output_file,msg);
      endif
    else
      fout = output_file;
    endif

    if (isfield(options, "weak_spring"))
      start_node = 1;
      end_node = columns(beam.Xn);

      if (isfield(options, "start_node"))
        ++start_node;
      endif

      if (isfield(options, "end_node"))
        --end_node;
      endif
      
      fprintf(fout, "\n# weak spring for curved beam node [%d:%d]\n", start_node, end_node);
      ispring = int32(0);
      
      for i=start_node:end_node
        for j=1:3
          fprintf(fout, "\n# curved beam: weak spring #%d,%d\n", i, j);
          fprintf(fout, "genel: %s + %d + %d - 1,\n", ...
                  options.weak_spring.element_number, ...
                  ++ispring, ...
                  options.id_offset);
          
          fprintf(fout, "\tspring support, %s + %d + %d - 1, structural, %d, algebraic,\n", ...
                  options.first_node_number, ...
                  i, ...
                  options.id_offset, ...
                  j);
          
          fprintf(fout, "\tlinear time variant viscoelastic generic, %s / %d,\n", ...
                  options.weak_spring.stiffness, ...
                  length(start_node:end_node));
          
          fprintf(fout, "\t\tramp, -1. / (%s), %s, %s + %s, 1.,\n", ...
                  options.weak_spring.tau, ...
                  options.weak_spring.initial_time, ...
                  options.weak_spring.initial_time, ...
                  options.weak_spring.tau);
          fprintf(fout, "\t\t0., null,\n");
          fprintf(fout, "\tposition, from node;\n");
        endfor
      endfor
    endif
    
    for i=1:length(beam.beams)
      fprintf(fout,"# curved beam: beam #%d\n", i);
      fprintf(fout,"beam3: %s + %d + %d - 1,\n", options.first_beam_number, i, options.id_offset);
      Rn = zeros(3, 3, 3);
      Rg = zeros(3, 3, 2);   
      
      if (options.reorient_cross_sections)
        ## required for mbdyn because rotation parameters could become singular
        ## if the relative angle is equal to or higher than 180 degrees
        
        for j=1:3
          Rn(:, :, j) = beam.Rn(:, :, beam.beams(i).nidx(j));
        endfor
        
        for j=2:3
          Rn(:, :, j) = Rn(:, :, j).' * beam_orientation_matrix(Rn(:, :, 1), Rn(:, :, j));
          mbdyn_pre_beam_check_rotation_matrix(Rn(:, :, j));
        endfor
        
        for j=1:2
          Rg(:, :, j) = beam.Rg(:, :, beam.beams(i).gidx(j));
        endfor
        
        for j=1:2
          Rg(:, :, j) = Rg(:, :, j).' * beam_orientation_matrix(Rn(:, :, 1), Rg(:, :, j));
          mbdyn_pre_beam_check_rotation_matrix(Rg(:, :, j));
        endfor
        
        Rn(:, :, 1) = eye(3);
      else
        for j=1:size(Rn, 3)
          Rn(:, :, j) = eye(3);
        endfor
        
        for j=1:size(Rg, 3)
          Rg(:, :, j) = eye(3);
        endfor
      endif
      
      for j=1:3
        fprintf(fout,"\t# node %d\n",j);

        if (j == 1 && i == 1 && isfield(options, "start_node"))
          fprintf(fout, "\t%s,\n", options.start_node);
        elseif (j == 3 && i == length(beam.beams) && isfield(options, "end_node"))
          fprintf(fout, "\t%s,\n", options.end_node);
        else
          fprintf(fout, "\t%s + %d + %d - 1,\n", options.first_node_number, beam.beams(i).nidx(j), options.id_offset);
        endif

        fprintf(fout, "\t\tposition, reference, %s + %d + %d - 1, null,\n", options.first_reference_frame_number, beam.beams(i).nidx(j), options.id_offset);
        fprintf(fout, "\t\torientation, reference, %s + %d + %d - 1, matr,\n%s\n", options.first_reference_frame_number, beam.beams(i).nidx(j), options.id_offset, sprintf("\t\t\t%.16e, %.16e, %.16e,\n", Rn(:, :, j).'));
      endfor

      for j=1:2
        fprintf(fout,"\t# orientation matrix section %d\n", j);
        fprintf(fout, "\treference, %s + %d + %d - 1, matr,\n%s\n", options.first_reference_frame_number, beam.beams(i).gidx(j), options.id_offset + columns(beam.Xn), sprintf("\t\t\t%.16e, %.16e, %.16e,\n", Rg(:, :, j).'));
        
        fprintf(fout,"\t# constitutive law section %d\n", j);
        fprintf(fout,"\treference, %s,\n", options.constitutive_law_number);
      endfor

      fprintf(fout, "\toutput, %s;\n", options.output_flag);

      for k=1:2
        for j=1:length(options.private_data_output)
          fprintf(fout, ...
                  "\n\tdrive caller: %s + (%d + %d - 1) * %d + %d,\n", ...
                  options.private_data_first_drv_idx, ...
                  i, ...
                  options.id_offset, ...
                  2 * length(options.private_data_output), ...
                  j - 1 + length(options.private_data_output) * (k - 1));
          fprintf(fout, "\telement, %s + %d + %d - 1, beam,\n", options.first_beam_number, i, options.id_offset);
          fprintf(fout, "\tstring, \"%s.%s\", direct, output, yes;\n", {"pI", "pII"}{k}, options.private_data_output{j});
        endfor
      endfor
    endfor

    if (options.output_Vi)
      force_name = {"Fx","Fy","Fz","Mx","My","Mz"};
      strain_name = {"ex", "ey", "ez", "kx", "ky", "kz"};

      fprintf(fout, "\n\tdrive caller: %s + %d", options.Vi_first_drv_idx, options.id_offset);
      fprintf(fout, ",\n\t\tmult, const, 0.5 * %.16e / %d", beam.sn(end), 2 * length(beam.beams));
      fprintf(fout, ",\n\t\t\tarray, %d", length(beam.beams) * 2 * length(force_name));

      for i=1:length(beam.beams)
        for k=1:2
          for j=1:length(force_name)
            fprintf(fout, ",\n\t\t\t\tmult");
            fprintf(fout, ",\n\t\t\t\t\telement, %s + %d + %d - 1, beam", options.first_beam_number, i, options.id_offset);
            fprintf(fout, ",\n\t\t\t\t\t\tstring, \"%s.%s\", direct", {"pI", "pII"}{k}, force_name{j});
            fprintf(fout, ",\n\t\t\t\t\telement, %s + %d + %d - 1, beam", options.first_beam_number, i, options.id_offset);
            fprintf(fout, ",\n\t\t\t\t\t\tstring, \"%s.%s\", direct", {"pI", "pII"}{k}, strain_name{j});
          endfor
        endfor
      endfor

      fprintf(fout, ",\n\toutput, %s;\n\n", options.private_data_output_flag);
    endif

  unwind_protect_cleanup
    if (owns_fd && fout ~= -1)
      fclose(fout);
    endif
  end_unwind_protect
endfunction

function R = beam_orientation_matrix(R1, R2)
  e1 = R2(:, 1);
  
  e3 = cross(e1, R1(:, 2));
  e2 = cross(R1(:, 3), e1);
  
  if (norm(e3) > norm(e2))
    e2 = cross(e3, e1);
  else
    e3 = cross(e1, e2);
  endif
  
  R = [ e1, e2, e3 ];
  
  for i=1:columns(R)
    R(:, i) /= norm(R(:, i));
  endfor

  mbdyn_pre_beam_check_rotation_matrix(R);
endfunction

%!test
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 10;
%! l = 1000e-3;
%! g = 9.81;
%! D = 10e-3;
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
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_write_beams_XXXXXX"));
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
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
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
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., mult, time, g;\n");
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
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert(bodies(i).node, int32(i));
%!   endfor
%!   w = zeros(1, numel(deformation));
%!   for i=1:numel(w)
%!     w(i) = deformation{i}(end, 3);
%!   endfor
%!   z = zeros(1, numel(log_dat.nodes));
%!   for i=1:numel(log_dat.nodes)
%!     z(i) = l - log_dat.nodes(i).X0(1);
%!   endfor
%!   wref = -rho * A * g * l^4 / (24 * E * Iy) * (3 - 4 * z / l + (z / l).^4);
%!   if (f_plot_deformation)
%!     figure("visible", "off");
%!     hold on;
%!     plot(z, wref, "-;wref;0");
%!     plot(z, w, "-;w;1");
%!     grid on;
%!     grid minor on;
%!     title("deformation of a cantilever beam");
%!     xlabel("z [m]");
%!     ylabel("w [m]");
%!   endif
%!   tol = 1e-2;
%!   assert(w, wref, tol * max(abs(wref)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 30;
%! R = 1000e-3;
%! F = 150;
%! D = 50e-3;
%! A = D^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = D^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! Phi = linspace(0, pi / 2, N);
%! X = [-R * cos(Phi);
%!      R * sin(Phi);
%!      zeros(1, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_write_beams_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-4;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
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
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "       forces: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fprintf(fd, " force: 1, absolute, %d, position, reference, node, null, 0., -1., 0., mult, time, F;\n", columns(beam.Xn));
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
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert(bodies(i).node, int32(i));
%!   endfor
%!   Xn = [log_dat.nodes.X0];
%!   Phin = atan2(Xn(2, :), -Xn(1, :));
%!   wref = F * R^3 / (2 * E * Iz) * Phin .* sin(Phin);
%!   uref = F * R^3 / (2 * E * Iz) * (sin(Phin) - Phin .* cos(Phin));
%!   Uref = [wref .* cos(Phin) + uref .* sin(Phin);
%!           -wref .* sin(Phin) + uref .* cos(Phin);
%!           zeros(1, numel(Phin))];
%!   U = zeros(3, numel(Phin));
%!   for i=1:columns(U)
%!     U(:, i) = deformation{i}(end, 1:3).';
%!   endfor
%!   tol = 1e-3;
%!   assert(U, Uref, tol * norm(Uref));
%!   wref_1 = F * R^3 / (2 * E * Iz);
%!   w_1 = deformation{end}(end, 1);
%!   tol = 1e-3;
%!   assert(w_1, wref_1, tol * max(abs(wref_1)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 10;
%! l = 1000e-3;
%! g = 9.81;
%! D = 10e-3;
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
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_write_beams_XXXXXX"));
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
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
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
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., mult, time, g;\n");
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
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert(bodies(i).node, int32(i));
%!   endfor
%!   w = zeros(1, numel(deformation));
%!   for i=1:numel(w)
%!     w(i) = deformation{i}(end, 3);
%!   endfor
%!   z = zeros(1, numel(log_dat.nodes));
%!   for i=1:numel(log_dat.nodes)
%!     z(i) = l - log_dat.nodes(i).X0(1);
%!   endfor
%!   wref = -rho * A * g * l^4 / (24 * E * Iy) * (3 - 4 * z / l + (z / l).^4);
%!   if (f_plot_deformation)
%!     figure("visible", "off");
%!     hold on;
%!     plot(z, wref, "-;wref;0");
%!     plot(z, w, "-;w;1");
%!     grid on;
%!     grid minor on;
%!     title("deformation of a cantilever beam");
%!     xlabel("z [m]");
%!     ylabel("w [m]");
%!   endif
%!   tol = 1e-2;
%!   assert(w, wref, tol * max(abs(wref)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 30;
%! R = 1000e-3;
%! F = 150;
%! D = 50e-3;
%! A = D^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = D^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! Phi = linspace(0, pi / 2, N);
%! X = [-R * cos(Phi);
%!      R * sin(Phi);
%!      zeros(1, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_write_beams_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-4;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
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
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "       forces: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fprintf(fd, " force: 1, absolute, %d, position, reference, node, null, 0., -1., 0., mult, time, F;\n", columns(beam.Xn));
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
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert(bodies(i).node, int32(i));
%!   endfor
%!   Xn = [log_dat.nodes.X0];
%!   Phin = atan2(Xn(2, :), -Xn(1, :));
%!   wref = F * R^3 / (2 * E * Iz) * Phin .* sin(Phin);
%!   uref = F * R^3 / (2 * E * Iz) * (sin(Phin) - Phin .* cos(Phin));
%!   Uref = [wref .* cos(Phin) + uref .* sin(Phin);
%!           -wref .* sin(Phin) + uref .* cos(Phin);
%!           zeros(1, numel(Phin))];
%!   U = zeros(3, numel(Phin));
%!   for i=1:columns(U)
%!     U(:, i) = deformation{i}(end, 1:3).';
%!   endfor
%!   tol = 1e-3;
%!   assert(U, Uref, tol * norm(Uref));
%!   wref_1 = F * R^3 / (2 * E * Iz);
%!   w_1 = deformation{end}(end, 1);
%!   tol = 1e-3;
%!   assert(w_1, wref_1, tol * max(abs(wref_1)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 10;
%! l = 1000e-3;
%! g = 9.81;
%! D = 10e-3;
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
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_write_beams_XXXXXX"));
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
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 20;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., mult, time, g;\n");
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
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert(bodies(i).node, int32(i));
%!   endfor
%!   w = zeros(1, numel(deformation));
%!   for i=1:numel(w)
%!     w(i) = deformation{i}(end, 3);
%!   endfor
%!   z = zeros(1, numel(log_dat.nodes));
%!   for i=1:numel(log_dat.nodes)
%!     z(i) = l - log_dat.nodes(i).X0(1);
%!   endfor
%!   wref = -rho * A * g * l^4 / (24 * E * Iy) * (3 - 4 * z / l + (z / l).^4);
%!   if (f_plot_deformation)
%!     figure("visible", "off");
%!     hold on;
%!     plot(z, wref, "-;wref;0");
%!     plot(z, w, "-;w;1");
%!     grid on;
%!     grid minor on;
%!     title("deformation of a cantilever beam");
%!     xlabel("z [m]");
%!     ylabel("w [m]");
%!   endif
%!   tol = 1e-2;
%!   assert(w, wref, tol * max(abs(wref)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! f_print_input_file = false;
%! f_plot_deformation = false;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 30;
%! R = 1000e-3;
%! F = 150;
%! D = 50e-3;
%! A = D^2 * pi / 4.;
%! Ay = 9. / 10. * A;
%! Az = Ay;
%! Iy = D^4 * pi / 64.;
%! Iz = Iy;
%! It = Iy + Iz;
%! E = 210000e6;
%! G = 81500e6;
%! rho = 7850;
%! Phi = linspace(0, pi / 2, N);
%! X = [-R * cos(Phi);
%!      R * sin(Phi);
%!      zeros(1, numel(Phi))];
%! options.A = "A";
%! options.Ip = "It";
%! options.Iy = "Iy";
%! options.Iz = "Iz";
%! options.constitutive_law_number = "const_law_id_beam";
%! options.rho = "rho";
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_write_beams_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: integer const_law_id_beam = 1;\n");
%!     fprintf(fd, " set: real F = %g;\n", F);
%!     fprintf(fd, " set: real A = %g;\n", A);
%!     fprintf(fd, " set: real Ay = %g;\n", Ay);
%!     fprintf(fd, " set: real Az = %g;\n", Az);
%!     fprintf(fd, " set: real Iy = %g;\n", Iy);
%!     fprintf(fd, " set: real Iz = %g;\n", Iz);
%!     fprintf(fd, " set: real It = %g;\n", It);
%!     fprintf(fd, " set: real E = %g;\n", E);
%!     fprintf(fd, " set: real G = %g;\n", G);
%!     fprintf(fd, " set: real rho = %g;\n", rho);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd, scale, row max column max, always;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-4;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "       model: static;\n");
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "       joints: 1;\n");
%!     fputs(fd, "       forces: 1;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fprintf(fd, " force: 1, absolute, %d, position, reference, node, null, 0., -1., 0., mult, time, F;\n", columns(beam.Xn));
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
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert(bodies(i).node, int32(i));
%!   endfor
%!   Xn = [log_dat.nodes.X0];
%!   Phin = atan2(Xn(2, :), -Xn(1, :));
%!   wref = F * R^3 / (2 * E * Iz) * Phin .* sin(Phin);
%!   uref = F * R^3 / (2 * E * Iz) * (sin(Phin) - Phin .* cos(Phin));
%!   Uref = [wref .* cos(Phin) + uref .* sin(Phin);
%!           -wref .* sin(Phin) + uref .* cos(Phin);
%!           zeros(1, numel(Phin))];
%!   U = zeros(3, numel(Phin));
%!   for i=1:columns(U)
%!     U(:, i) = deformation{i}(end, 1:3).';
%!   endfor
%!   tol = 1e-3;
%!   assert(U, Uref, tol * norm(Uref));
%!   wref_1 = F * R^3 / (2 * E * Iz);
%!   w_1 = deformation{end}(end, 1);
%!   tol = 1e-3;
%!   assert(w_1, wref_1, tol * max(abs(wref_1)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!demo
%! f_print_input_file = false;
%! f_plot_deformation = true;
%! if (f_plot_deformation)
%!  close("all");
%! endif
%! N = 10;
%! l = 1000e-3;
%! g = 9.81;
%! D = 10e-3;
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
%! beam = mbdyn_pre_beam_compute(X, N, 20);
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_write_beams_XXXXXX"));
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
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-2;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
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
%!     fprintf(fd, "     structural nodes: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     rigid bodies: %d;\n", columns(beam.Xn));
%!     fprintf(fd, "     beams: %d;\n", numel(beam.beams));
%!     fputs(fd, "     joints: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " constitutive law: 1, 6, linear elastic generic, diag, E * A, G * Ay, G * Az, G * It, E * Iy, E * Iz;\n");
%!     mbdyn_pre_beam_write_reference_frames(beam, fd, options);
%!     fputs(fd, " begin: nodes;\n");
%!     mbdyn_pre_beam_write_nodes(beam, fd, options);
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     mbdyn_pre_beam_write_bodies(beam, fd, options);
%!     mbdyn_pre_beam_write_beams(beam, fd, options);
%!     fputs(fd, "joint: 1, clamp, 1, node, node;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., mult, time, g;\n");
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
%!   if (~options.verbose)
%!     options.logfile = [fname, ".stdout"];
%!   endif
%!   mbdyn_solver_run(fname, options);
%!   [t, trajectory, deformation] = mbdyn_post_load_output_struct(options.output_file);
%!   log_dat = mbdyn_post_load_log(fname);
%!   bodies = mbdyn_post_load_log_body(fname);
%!   tol = 1e-5;
%!   for i=1:numel(log_dat.nodes)
%!     assert(log_dat.nodes(i).X0, beam.Xn(:, i), tol);
%!     assert(log_dat.nodes(i).R0, beam.Rn(:, :, i), tol);
%!   endfor
%!   for i=1:numel(log_dat.beams3)
%!     for j=1:3
%!       assert(log_dat.beams3(i).nodes(j).label, int32(beam.beams(i).nidx(j)));
%!     endfor
%!   endfor
%!   for i=1:numel(bodies)
%!     assert(bodies(i).node, int32(i));
%!   endfor
%!   w = zeros(1, numel(deformation));
%!   for i=1:numel(w)
%!     w(i) = deformation{i}(end, 3);
%!   endfor
%!   z = zeros(1, numel(log_dat.nodes));
%!   for i=1:numel(log_dat.nodes)
%!     z(i) = l - log_dat.nodes(i).X0(1);
%!   endfor
%!   wref = -rho * A * g * l^4 / (24 * E * Iy) * (3 - 4 * z / l + (z / l).^4);
%!   if (f_plot_deformation)
%!     figure("visible", "off");
%!     hold on;
%!     plot(z, wref, "-;wref;0");
%!     plot(z, w, "-;w;1");
%!     grid on;
%!     grid minor on;
%!     title("deformation of a cantilever beam");
%!     xlabel("z [m]");
%!     ylabel("w [m]");
%!   endif
%!   tol = 1e-2;
%!   assert(w, wref, tol * max(abs(wref)));
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, "*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
