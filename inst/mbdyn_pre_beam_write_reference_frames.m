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
## @deftypefn {Function File} mbdyn_pre_beam_write_reference_frames(@var{beam}, @var{output_file}, @var{options})
##
## Generate an MBDyn input file <@var{output_file}> containing all reference frames of beam model <@var{beam}>.
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
## @end deftypefn

function mbdyn_pre_beam_write_reference_frames(beam, output_file, options)
  if (nargin < 2)
    print_usage();
  endif

  if (nargin < 3)
    options = struct();
  endif

  if (~isfield(options, "reference_frame"))
    options.reference_frame = "global";
  endif

  if (~isfield(options, "first_reference_frame_number"))
    options.first_reference_frame_number = 1;
  endif

  if (~isfield(options, "id_offset"))
    options.id_offset = 0;
  endif

  if (~isfield(options, "open_mode"))
    options.open_mode = "wt";
  endif

  if (~ischar(options.reference_frame))
    options.reference_frame = sprintf("%d", options.reference_frame);
  endif

  if (~ischar(options.first_reference_frame_number))
    options.first_reference_frame_number = sprintf("%d", options.first_reference_frame_number);
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

    mbdyn_curved_beam_write_frames(fout, beam.Xn, beam.Rn, options.id_offset, "\n# curved beam: node #", options);
    mbdyn_curved_beam_write_frames(fout, beam.Xg, beam.Rg, options.id_offset + columns(beam.Xn), "\n# curved beam: gauss point #", options);
  unwind_protect_cleanup
    if (owns_fd && fout ~= -1)
      fclose(fout);
    endif
  end_unwind_protect
endfunction

function mbdyn_curved_beam_write_frames(fout, X, R, offset, header, options)
  for i=1:columns(X)
    fprintf(fout, "%s%d\n", header, i);
    fprintf(fout, "reference: %s + %d + %d - 1,\n", options.first_reference_frame_number, i, offset);
    fprintf(fout, "\tposition, reference, %s,\n", options.reference_frame);

    for j=1:3
      fprintf(fout, "\t\t%.16g,\n", X(j,i));
    endfor

    fprintf(fout, "\torientation, reference, %s, matr,\n", options.reference_frame);

    for k=1:3
      fprintf(fout, "\t\t");

      for l=1:3
        fprintf(fout, "%.16g, ", R(k,l,i));
      endfor

      fprintf(fout, "\n");
    endfor

    fprintf(fout, "\tvelocity, reference, %s, null,\n", options.reference_frame);
    fprintf(fout, "\tangular velocity, reference, %s, null;\n\n", options.reference_frame);
  endfor
endfunction
