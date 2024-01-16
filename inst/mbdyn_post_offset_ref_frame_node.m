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
## @deftypefn {Function File} @var{x} = mbdyn_post_offset_ref_frame_node(@var{reference_frames}, @var{nodes}, @var{vars}, @var{reference}, @var{node})
## @deftypefnx {} @var{x} = mbdyn_post_offset_ref_frame_node(@dots{}, @var{frame})
##
## Returns the offset between reference frame <@var{reference}> and structural node <@var{node}>
## measured in the global reference frame or in the reference frame of the node.
##
## @var{reference_frames} @dots{} Return value from mbdyn_post_load_output_rfm.
##
## @var{nodes} @dots{} Return value from mbdyn_post_load_log_nodes.
##
## @var{vars} @dots{} Return value from mbdyn_post_load_log_vars.
##
## @var{reference} @dots{} Character string name of the reference frame.
##
## @var{node} @dots{} Character string name of the node.
##
## @var{frame} @dots{} Reference frame to use for the offset. One of ("global", "node").
##
## @end deftypefn

function x = mbdyn_post_offset_ref_frame_node(reference_frames, nodes, vars, reference, node, frame)
  if (nargin < 5 || nargin > 6 || nargout > 1)
    print_usage();
  endif

  if (nargin < 6)
    frame = "global";
  endif

  switch(reference)
    case "node"
      x = zeros(3, 1);
    otherwise
      if (strncmp(reference, "ref_", length("ref_")))
        ref_id_i = getfield(vars, reference);
        ref_idx = find(reference_frames.ref_id == ref_id_i);

        if (length(ref_idx) == 0)
          error("reference id %d(%s) not found!", ref_id_i, reference);
        endif

        node_id_i = getfield(vars, node);

        node_idx = find([nodes.label] == node_id_i);

        if (length(node_idx) == 0)
          error("node %d(%s) not found!", node_id_i, node);
        endif

        x = reference_frames.position{ref_idx}(1, 1:3).' - nodes(node_idx).X0;

        switch (frame)
          case "global"
          case "node"
            x = nodes(node_idx).R0.' * x;
          otherwise
            error("invalid value for argument frame=\"%s\"", frame);
        endswitch
      else
        error("unsupported reference: %s!", reference);
      endif
  endswitch
endfunction
