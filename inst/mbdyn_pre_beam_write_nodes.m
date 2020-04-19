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
## @deftypefn {Function File} mbdyn_pre_beam_write_nodes(@var{beam}, @var{output_file}, @var{options})
##
## Generate an MBDyn input file <@var{output_file}> containing all nodes of beam model <@var{beam}>.
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
## @var{options}.id_offset @dots{} Integer offset for all body id's.
##
## @var{options}.start_node @dots{} If the beginning of the beam should be connected to an existing node, then <@var{start_node}> is
## the string name of an integer variable holding the existing node number.
##
## @var{options}.end_node @dots{} If the end of the beam should be connected to an existing node, then <@var{end_node}> is the
## string name of an integer variable holding the existing node number.
##
## @var{options}.node_type @dots{} Possible values are "static" and "dynamic".
##
## @end deftypefn

function mbdyn_pre_beam_write_nodes(beam, output_file, options)
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

  if (~isfield(options,"node_type"))
    options.node_type = "dynamic";
  endif

  if (~isfield(options, "id_offset"))
    options.id_offset = 0;
  endif

  if (~isfield(options, "open_mode"))
    options.open_mode = "wt";
  endif

  if (~ischar(options.first_node_number))
    options.first_node_number = sprintf("%d", options.first_node_number);
  endif

  if (~ischar(options.first_reference_frame_number))
    options.first_reference_frame_number = sprintf("%d", options.first_reference_frame_number);
  endif

  fout = -1;
  owns_fd = false;

  unwind_protect
    if (ischar(output_file))
      owns_fd = true;

      [fout,msg] = fopen(output_file, options.open_mode);

      if (fout == -1)
        error("could not open file \"%s\": %s", output_file,msg);
      endif
    else
      fout = output_file;
    endif

    start_node = 1;
    end_node = columns(beam.Xn);

    if (isfield(options, "start_node"))
      ++start_node;
    endif

    if (isfield(options, "end_node"))
      --end_node;
    endif

    fprintf(fout, "\n# curved beam node [%d:%d]\n", start_node, end_node);

    for i=start_node:end_node
      fprintf(fout, "\n# curved beam: node #%d\n", i);
      fprintf(fout, "structural: %s + %d + %d - 1, %s,\n", options.first_node_number, i, options.id_offset, options.node_type);
      fprintf(fout, "\tposition, reference, %s + %d + %d - 1, null,\n", options.first_reference_frame_number, i, options.id_offset);
      fprintf(fout, "\torientation, reference, %s + %d + %d - 1, eye,\n", options.first_reference_frame_number, i, options.id_offset);
      fprintf(fout, "\tvelocity, reference, %s + %d + %d - 1, null,\n", options.first_reference_frame_number, i, options.id_offset);
      fprintf(fout, "\tangular velocity, reference, %s + %d + %d - 1, null", options.first_reference_frame_number, i, options.id_offset);

      if (isfield(options,"output_flag"))
        fprintf(fout,",\n");
        fprintf(fout, "\toutput, %s", options.output_flag);
      endif

      fprintf(fout, ";\n");
    endfor

  unwind_protect_cleanup
    if (owns_fd && fout ~= -1)
      fclose(fout);
    endif
  end_unwind_protect
endfunction
