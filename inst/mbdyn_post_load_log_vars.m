## Copyright (C) 2016 Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{vars},@var{node_idx},@var{force_idx}] = mbdyn_post_load_log_vars( @var{filename}, @var{node_id}, @var{force_id}, @var{options})
## Parse MBDyn's .log file "<@var{filename}.log>" and return all variable definitions.
##
## @var{vars} @dots{} Structure that receives the variables read from the .log file.
##
## @var{filename} @dots{} Name of the MBDyn .log file that contains the variables
##
## @end deftypefn

function [vars,node_idx,force_idx] = mbdyn_post_load_log_vars(mbdyn_filename, node_id, force_id, options)
  if (nargin < 1 || nargin > 4)
    print_usage();
  endif

  options.nodes = false;
  options.dof_info = true;
  options.vars = true;
  options.beams2 = false;
  options.beams3 = false;

  if (nargin >= 3)
    options.struct_force_id = force_id;
  endif

  data = mbdyn_post_load_log(mbdyn_filename, options);

  vars = data.vars;
endfunction
