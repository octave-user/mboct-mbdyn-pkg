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
## @deftypefn {Function File} [@var{beam_id}, @var{F_I}, @var{M_I}, @var{F_II}, @var{M_II}] = mbdyn_post_load_output_act(@var{mbdyn_filename})
## @deftypefnx {} [@dots{}] = mbdyn_post_load_output_act(@var{mbdyn_filename}, @var{filter_beam_id})
## @deftypefnx {} [@dots{}] = mbdyn_post_load_output_act(@var{mbdyn_filename}, @var{filter_beam_id}, @var{append_rows})
##
## Loads internal forces and moments of beam2 and beam3 elements from MBDyn output file "<@var{mbdyn_filename}>.act".
##
## @var{filter_beam_id} @dots{} Array of beam labels to be loaded. If this array is empty, all beams are loaded.
##
## @var{append_rows} @dots{} Hint for memory reallocation.
##
## @var{beam_id} @dots{} The label of the beam.
##
## @var{F_I} @dots{} The three components of the force at the first evaluation point,
## oriented according to the reference frame of that beam section.
##
## @var{M_I} @dots{} The three components of the moment at the first evaluation point,
## oriented according to the reference frame of that beam section.
##
## The three-node beam element generates six additional columns:
##
## @var{F_II} @dots{} The three components of the force at the second evaluation point,
## oriented according to the reference frame of that beam section;
##
## @var{M_II} @dots{} The three components of the moment at the second evaluation point,
## oriented according to the reference frame of that beam section.
##
## @end deftypefn

function [beam_id, F_I, M_I, F_II, M_II] = mbdyn_post_load_output_act(mbdyn_filename, filter_beam_id, append_rows)
  if (nargin < 1 || nargin > 3 || nargout > 5)
    print_usage();
  endif

  if (nargin < 2)
    filter_beam_id = [];
  endif

  if (nargin < 3)
    append_rows = 1024;
  endif

  act_filename = mbdyn_post_output_filename(mbdyn_filename, ".act");

  column_count = (nargout - 1) * 3;

  [beam_id, data] = mbdyn_post_load_output(act_filename, column_count, filter_beam_id, append_rows);

  for i=1:length(data)
    if (nargout >= 2)
      F_I{i} = data{i}(:,1:3);
    endif
    if (nargout >= 3)
      M_I{i} = data{i}(:,4:6);
    endif
    if (nargout >= 4)
      F_II{i} = data{i}(:,7:9);
    endif
    if (nargout >= 5)
      M_II{i} = data{i}(:,10:12);
    endif
  endfor
endfunction
