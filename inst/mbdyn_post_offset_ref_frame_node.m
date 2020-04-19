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

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_offset_ref_frame_nodes_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-1;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, "     default output: reference frames;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         set: integer ref_id_0 = 1001;\n");
%!     fputs(fd, "         set: integer ref_id_1 = 1002;\n");
%!     fputs(fd, "         set: integer ref_id_2 = 1003;\n");
%!     fputs(fd, "         set: integer node_id_1 = 2001;\n");
%!     fputs(fd, "         reference: ref_id_0,\n");
%!     fputs(fd, "                 position, reference, global, 0.4, 0.2, 0.3,\n");
%!     fputs(fd, "                 orientation, reference, global, euler123, pi / 2, pi / 4, pi,\n");
%!     fputs(fd, "                 velocity, reference, global, null,\n");
%!     fputs(fd, "                 angular velocity, reference, global, null;\n");
%!     fputs(fd, "         reference: ref_id_1,\n");
%!     fputs(fd, "                 position, reference, ref_id_0, 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_0, euler123, 0., 0., pi / 2,\n");
%!     fputs(fd, "                 velocity, reference, ref_id_0, null,\n");
%!     fputs(fd, "                 angular velocity, reference, ref_id_0, null;\n");
%!     fputs(fd, "         reference: ref_id_2,\n");
%!     fputs(fd, "                 position, reference, ref_id_1, -0.1, 0.5, 2.2,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_1, euler123, pi / 3, -pi/2, -pi / 3,\n");
%!     fputs(fd, "                 velocity, reference, ref_id_1, null,\n");
%!     fputs(fd, "                 angular velocity, reference, ref_id_1, null;\n");
%!     fputs(fd, "         structural: node_id_1, dynamic,\n");
%!     fputs(fd, "                 position, reference, ref_id_0, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_0, eye,\n");
%!     fputs(fd, "                 velocity, reference, ref_id_0, null,\n");
%!     fputs(fd, "                 angular velocity, reference, ref_id_0, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, node_id_1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         force: 1, absolute, node_id_1, position, null, 1., 0., 0, 100.;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., 9.81;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [ref.ref_id, ref.position, ref.orientation, ref.velocity, ref.angular_velocity] = mbdyn_post_load_output_rfm(options.output_file);
%!   X = mbdyn_post_offset_ref_frame_node(ref, log_dat.nodes, log_dat.vars, "ref_id_2", "node_id_1", "node");
%!   assert(X, [0.6; 2.1; 5.5], 1e-6);
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
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_offset_ref_frame_nodes_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "         initial time: 0;\n");
%!     fputs(fd, "         final time: 1;\n");
%!     fputs(fd, "         time step: 1e-1;\n");
%!     fputs(fd, "         linear solver: naive, colamd;\n");
%!     fputs(fd, "         method: ms, 0.6;\n");
%!     fputs(fd, "         max iterations: 10;\n");
%!     fputs(fd, "         tolerance: 1.e-6;\n");
%!     fputs(fd, "         threads: assembly, 1;\n");
%!     fputs(fd, "         derivatives max iterations: 10;\n");
%!     fputs(fd, "         derivatives coefficient: auto;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     forces: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, "     default output: reference frames;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         set: integer ref_id_0 = 1001;\n");
%!     fputs(fd, "         set: integer ref_id_1 = 1002;\n");
%!     fputs(fd, "         set: integer ref_id_2 = 1003;\n");
%!     fputs(fd, "         set: integer node_id_1 = 2001;\n");
%!     fputs(fd, "         reference: ref_id_0,\n");
%!     fputs(fd, "                 position, reference, global, 0.4, 0.2, 0.3,\n");
%!     fputs(fd, "                 orientation, reference, global, euler123, pi / 2, pi / 4, pi,\n");
%!     fputs(fd, "                 velocity, reference, global, null,\n");
%!     fputs(fd, "                 angular velocity, reference, global, null;\n");
%!     fputs(fd, "         reference: ref_id_1,\n");
%!     fputs(fd, "                 position, reference, ref_id_0, 1.1, 2.2, 3.3,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_0, euler123, 0., 0., pi / 2,\n");
%!     fputs(fd, "                 velocity, reference, ref_id_0, null,\n");
%!     fputs(fd, "                 angular velocity, reference, ref_id_0, null;\n");
%!     fputs(fd, "         reference: ref_id_2,\n");
%!     fputs(fd, "                 position, reference, ref_id_1, -0.1, 0.5, 2.2,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_1, euler123, pi / 3, -pi/2, -pi / 3,\n");
%!     fputs(fd, "                 velocity, reference, ref_id_1, null,\n");
%!     fputs(fd, "                 angular velocity, reference, ref_id_1, null;\n");
%!     fputs(fd, "         structural: node_id_1, dynamic,\n");
%!     fputs(fd, "                 position, reference, ref_id_0, null,\n");
%!     fputs(fd, "                 orientation, reference, ref_id_0, eye,\n");
%!     fputs(fd, "                 velocity, reference, ref_id_0, null,\n");
%!     fputs(fd, "                 angular velocity, reference, ref_id_0, null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, node_id_1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, "         force: 1, absolute, node_id_1, position, null, 1., 0., 0, 100.;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., 9.81;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.verbose = false;
%!   options.logfile = [fname, ".stdout"];
%!   mbdyn_solver_run(fname, options);
%!   log_dat = mbdyn_post_load_log(options.output_file);
%!   [ref.ref_id, ref.position, ref.orientation, ref.velocity, ref.angular_velocity] = mbdyn_post_load_output_rfm(options.output_file);
%!   X = mbdyn_post_offset_ref_frame_node(ref, log_dat.nodes, log_dat.vars, "ref_id_2", "node_id_1", "node");
%!   assert(X, [0.6; 2.1; 5.5], 1e-6);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
