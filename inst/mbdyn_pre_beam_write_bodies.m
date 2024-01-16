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
## @deftypefn {Function File} mbdyn_pre_beam_write_bodies(@var{beam}, @var{output_file}, @var{options})
##
## Generate an MBDyn input file <@var{output_file}> containing all bodies of beam model <@var{beam}>.
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
## @var{options}.first_body_number @dots{} String or integer number of the first beam element.
##
## @var{options}.id_offset @dots{} Integer offset for all body id's.
##
## @var{options}.start_node @dots{} If the beginning of the beam should be connected to an existing node, then <@var{start_node}> is
## the string name of an integer variable holding the existing node number.
##
## @var{options}.end_node @dots{} If the end of the beam should be connected to an existing node, then <@var{end_node}> is the
## string name of an integer variable holding the existing node number.
##
## @var{options}.rho @dots{} String name of a real variable holding the density of the beam.
##
## @var{options}.private_data_output_flag @dots{} String name of a boolean variable which enables/disables output of private data.
##
## @var{options}.node_type @dots{} Possible values are "static" and "dynamic".
##
## @var{options}.A @dots{} Beam cross section area.
##
## @var{options}.Ip @dots{} Beam cross section polar area moment of inertia.
##
## @var{options}.Iy @dots{} Beam cross section area moment of inertia along y.
##
## @var{options}.Iz @dots{} Beam cross section area moment of inertia along z.
##
## @end deftypefn

function mbdyn_pre_beam_write_bodies(beam, output_file, options)
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

  if (~isfield(options, "first_body_number"))
    options.first_body_number = 1;
  endif

  if (~isfield(options,"node_type"))
    options.node_type = "dynamic";
  endif

  if (~isfield(options, "id_offset"))
    options.id_offset = 0;
  endif

  if (~isfield(options, "open_mode"))
    options.open_mode = "wt";
  endif

  if (~isfield(options, "rho"))
    error("rho is not defined in options!");
  endif

  if (~ischar(options.first_node_number))
    options.first_node_number = sprintf("%d", options.first_node_number);
  endif

  if (~ischar(options.first_reference_frame_number))
    options.first_reference_frame_number = sprintf("%d", options.first_reference_frame_number);
  endif

  if (~ischar(options.first_body_number))
    options.first_body_number = sprintf("%d", options.first_body_number);
  endif

  if (~isfield(options,"output_flag"))
    options.output_flag = "default";
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

    for i=1:columns(beam.Xn)
      fprintf(fout, "\n# curved beam: body #%d\n", i);
      fprintf(fout, "body: %s + %d + %d - 1,", options.first_body_number, i, options.id_offset);

      node_label = [];

      switch (i)
        case 1
          if (isfield(options, "start_node"))
            node_label = options.start_node;
          endif
        case columns(beam.Xn)
          if (isfield(options, "end_node"))
            node_label = options.end_node;
          endif
        otherwise
      endswitch

      if (isempty(node_label))
        node_label = sprintf("%s + %d + %d - 1", options.first_node_number, i, options.id_offset);
      endif

      fprintf(fout, "%s,\n", node_label);

      fprintf(fout, "\t# mass\n");
      fprintf(fout, "\t%s * %s * %.16g,\n", options.rho, options.A, beam.bodies(i).ds);
      fprintf(fout, "\t# relative center of mass\n");
      fprintf(fout, "\treference, %s + %d + %d - 1, null,\n", options.first_reference_frame_number, i, options.id_offset);
      fprintf(fout, "\t# inertia matrix\n");
      fprintf(fout, "\tdiag,  %s * %s * %.16g,\n", options.rho, options.Ip, beam.bodies(i).ds);
      fprintf(fout, "\t       %s * %s * %.16g,\n", options.rho, options.Iy, beam.bodies(i).ds);
      fprintf(fout, "\t       %s * %s * %.16g,\n", options.rho, options.Iz, beam.bodies(i).ds);
      fprintf(fout, "\t# orientation of the inertia tensor\n");
      fprintf(fout, "\torientation, reference, %s + %d + %d - 1, eye,\n", options.first_reference_frame_number, i, options.id_offset);
      fprintf(fout, "\toutput, %s;\n", options.output_flag);
    endfor
  unwind_protect_cleanup
    if (owns_fd && fout ~= -1)
      fclose(fout);
    endif
  end_unwind_protect
endfunction
