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
## @deftypefn {Function File} mbdyn_pre_beam_print_nodes(@var{beam}, @var{output_file}, @var{options})
##
## Print a list of nodes for beam model <@var{beam}> to file <@var{output_file}>.
##
## @var{beam} @dots{} Return value from mbdyn_pre_beam_compute.
##
## @var{output_file} @dots{} Output file name.
##
## @var{options}.open_mode @dots{} Mode passed to fopen (e.g. "wt", "at").
##
## @end deftypefn

function mbdyn_pre_beam_print_nodes(beam, output_file, options)
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
      [fout, msg] = fopen(output_file, options.open_mode);

      owns_fd = true;

      if (fout == -1)
        error("could not open file \"%s\": %s", output_file,msg);
      endif
    else
      fout = output_file;
    endif

    mbdyn_curved_beam_print_frames(fout, beam.Xn, beam.Rn, options.id_offset, "", options);
  unwind_protect_cleanup
    if (owns_fd && fout ~= -1)
      fclose(fout);
    endif
  end_unwind_protect
endfunction

function mbdyn_curved_beam_print_frames(fout, X, R, offset, header, options)
  fprintf(fout, "#\t");

  for j=1:3
    fprintf(fout, "X%d\t", j);
  endfor

  for k=1:3
    for l=1:3
      fprintf(fout, "R%d%d\t", l, k);
    endfor
  endfor

  fprintf(fout, "\n");

  for i=1:columns(X)
    fprintf(fout, "%s%d", header, i);

    for j=1:3
      fprintf(fout, "\t%.16g", X(j,i));
    endfor

    for k=1:3
      for l=1:3
        fprintf(fout, "\t%.16g", R(l, k, i));
      endfor
    endfor

    fprintf(fout, "\n");
  endfor
endfunction
