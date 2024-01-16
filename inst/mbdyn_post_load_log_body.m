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
## @deftypefn {Function File} [@var{bodies}] = mbdyn_post_load_log_body(@var{mbdyn_filename})
##
## Loads information about rigid and flexible bodies from MBDyn output file "<@var{mbdyn_filename}>.log".
##
## @var{bodies} @dots{} Struct array which contains the information about the bodies in the log file.
##
## @var{bodies}(@var{i}).label @dots{} Label of the body @var{i}.
##
## @var{bodies}(@var{i}).dm @dots{} Mass of the body.
##
## @var{bodies}(@var{i}).node @dots{} Label of the structural node the body is connected to.
##
## @var{bodies}(@var{i}).Xgc @dots{} Center of gravity of the body with respect to the node.
##
## @var{bodies}(@var{i}).J @dots{} Inertia matrix of the body in the reference frame of the node.
##
## @end deftypefn

function bodies = mbdyn_post_load_log_body(mbdyn_filename)
  if (nargin ~= 1 || nargout ~= 1)
    print_usage();
  endif

  bodies = struct();

  log_filename = mbdyn_post_output_filename(mbdyn_filename, ".log");

  fid = -1;

  unwind_protect
    [fid, msg] = fopen(log_filename, "rt");

    if (fid == -1)
      error("could not open file \"%s\": %s",log_filename,msg);
    endif

    lineno = 0;
    bodyno = 0;

    while (true)
      line = fgets(fid);

      if (feof(fid))
        break;
      endif

      lineno++;
      dm = nan;
      Xgc = nan(3,1);
      J = nan(3,3);

      tokens = {"body", "modal"};

      for i=1:length(tokens)
        token = cstrcat(tokens{i}, ":");

        if (strncmp(line, token, length(token)))
          [type, ...
           label, ...
           node, ...
           dm, ...
           Xgc(1), ...
           Xgc(2), ...
           Xgc(3), ...
           J(1,1), ...
           J(1,2), ...
           J(1,3), ...
           J(2,1), ...
           J(2,2), ...
           J(2,3), ...
           J(3,1), ...
           J(3,2), ...
           J(3,3), ...
           count] = sscanf(line, "%s %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g", "C");

          if (count >= 1)
            if (count ~= 16)
              error("invalid input at line %d of file \"%s\":\"%s\"", lineno, log_filename, line);
            endif

            bodyno++;
            bodies(bodyno).label = int32(label);
            bodies(bodyno).node = int32(node);
            bodies(bodyno).dm = dm;
            bodies(bodyno).Xgc = Xgc;
            bodies(bodyno).J = J;
            bodies(bodyno).type = tokens{i};
          endif
        endif
      endfor
    endwhile
  unwind_protect_cleanup
    if (fid ~= -1)
      fclose(fid);
    endif
  end_unwind_protect
endfunction
