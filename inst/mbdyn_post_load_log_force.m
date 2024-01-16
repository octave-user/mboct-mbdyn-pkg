## Copyright (C) 2011(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{forces}] = mbdyn_post_load_log_force(@var{mbdyn_filename})
##
## Parses the MBDyn log-file "<@var{mbdyn_filename}>.log" and returns information about the structural forces.
##
## @var{forces}.label @dots{} The label of the structural force or couple.
##
## @var{forces}.node1 @dots{} The label of the first structural node the force is applied to.
##
## @var{forces}.arm1 @dots{}  The distance between the force and the first structural node.
##
## @var{forces}.node2 @dots{} If the force is an internal structural force node2 is the label of the second node the force is applied to.
##
## @var{forces}.arm2 @dots{}  If the force is an internal structural force arm2 is the distance between the force and the second node.
##
## @var{forces}.internal @dots{} Equal to one if the force is an internal structural force.
##
## @var{forces}.absolute @dots{} Equal to one if the force is a absolute force. Equal to zero if the force is a follower force.
##
## @var{force}.couple @dots{} Equal to zero if it is a force. Equal to one if it is a couple.
##
## @end deftypefn

function [forces] = mbdyn_post_load_log_force(mbdyn_filename)
  if (nargin < 1)
    print_usage();
  endif

  log_filename = mbdyn_post_output_filename(mbdyn_filename, ".log");

  fid = -1;

  unwind_protect
    [fid, msg] = fopen(log_filename, "rt");

    if (fid == -1)
      error("could not open file \"%s\": %s", log_filename, msg);
    endif

    i = 0;
    line_no = 0;

    forces = struct()([]);

    while (1)
      line = fgets(fid);

      if (~ischar(line) && line == -1)
        break;
      endif

      ++line_no;

      tag_end = find(line == ':');

      if (length(tag_end) >= 1)
        tag_end = tag_end(1);
        tag = line(1:tag_end);
        data = line(length(tag)+1:end);
        switch (tag)
          case {"structural internal follower force:",
                "structural internal follower couple:",
                "structural internal absolute force:",
                "structural internal absolute couple:",
                "structural follower force:",
                "structural follower couple:",
                "structural absolute force:",
                "structural absolute couple:"}

            ++i;

            [label, node1, arm1x, arm1y, arm1z, node2, arm2x, arm2y, arm2z, count] = sscanf(data,"%d %d %g %g %g %d %g %g %g", "C");

            if (count ~= 5 && count ~= 9)
              error("parse error in file \"%s\": line %d",log_filename,line_no);
            endif

            if (count >= 5)
              forces(i).label = int32(label);
              forces(i).node1 = int32(node1);
              forces(i).arm1 = [ arm1x; arm1y; arm1z ];
            endif

            if (count >= 9)
              forces(i).node2 = int32(node2);
              forces(i).arm2 = [ arm2x; arm2y; arm2z ];
            endif

            forces(i).internal = logical(length(strfind(tag, "internal")) > 0);
            forces(i).absolute = logical(length(strfind(tag, "absolute")) > 0);
            forces(i).couple   = logical(length(strfind(tag, "couple")) > 0);
        endswitch
      endif
    endwhile
  unwind_protect_cleanup
    if (fid ~= -1)
      fclose(fid);
    endif
  end_unwind_protect

  if (length(forces) > 0)
    [force_labels,force_idx] = sort([forces.label]);
    forces = forces(force_idx);
  endif
endfunction
