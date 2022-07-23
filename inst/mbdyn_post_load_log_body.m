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

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_log_body_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 20;\n");
%!     fputs(fd, " set: real m1 = 1;\n");
%!     fputs(fd, " set: real m2 = 2;\n");
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
%!     fputs(fd, "         nonlinear solver: nox, modified, 10,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 12,\n");
%!     fputs(fd, "             krylov subspace size, 12;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     use automatic differentiation;\n");
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1., 2., 3.,\n");
%!     fputs(fd, "                 euler123, 0.1, 0.2, 0.3,\n");
%!     fputs(fd, "                 10., 20., 30.,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, null, 1., 0., 0, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, null, 0., 1., 0, F2;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
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
%!   bodies = mbdyn_post_load_log_body(fname);
%!   assert([bodies.label], int32([1,2]));
%!   assert([bodies.node], int32([1,2]));
%!   assert([bodies.dm], [1, 2]);
%!   assert([bodies.Xgc], [0.1, 0.01; 0.2, 0.02; 0.3, 0.03]);
%!   assert(bodies(1).J, diag([1.1, 1.2, 1.3]) - bodies(1).dm * skew(bodies(1).Xgc) * skew(bodies(1).Xgc), sqrt(eps) * bodies(1).dm);
%!   assert(bodies(2).J, diag([2.1, 2.2, 2.3]) - bodies(2).dm * skew(bodies(2).Xgc) * skew(bodies(2).Xgc), sqrt(eps) * bodies(2).dm);
%!   assert(bodies(1).type, "body");
%!   assert(bodies(2).type, "body");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_log_body_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 20;\n");
%!     fputs(fd, " set: real m1 = 1;\n");
%!     fputs(fd, " set: real m2 = 2;\n");
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
%!     fputs(fd, "         nonlinear solver: nox, modified, 10,\n");
%!     fputs(fd, "             keep jacobian matrix,\n");
%!     fputs(fd, "             inner iterations before assembly, 6,\n");
%!     fputs(fd, "             jacobian operator, newton krylov,\n");
%!     fputs(fd, "             solver, line search based,\n");
%!     fputs(fd, "             forcing term, type 2,\n");
%!     fputs(fd, "             direction, newton,\n");
%!     fputs(fd, "             weighted rms absolute tolerance, 0.,\n");
%!     fputs(fd, "             weighted rms relative tolerance, 0.,\n");
%!     fputs(fd, "             linear solver, gmres,\n");
%!     fputs(fd, "             linear solver max iterations, 12,\n");
%!     fputs(fd, "             krylov subspace size, 12;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1., 2., 3.,\n");
%!     fputs(fd, "                 euler123, 0.1, 0.2, 0.3,\n");
%!     fputs(fd, "                 10., 20., 30.,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, null, 1., 0., 0, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, null, 0., 1., 0, F2;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
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
%!   bodies = mbdyn_post_load_log_body(fname);
%!   assert([bodies.label], int32([1,2]));
%!   assert([bodies.node], int32([1,2]));
%!   assert([bodies.dm], [1, 2]);
%!   assert([bodies.Xgc], [0.1, 0.01; 0.2, 0.02; 0.3, 0.03]);
%!   assert(bodies(1).J, diag([1.1, 1.2, 1.3]) - bodies(1).dm * skew(bodies(1).Xgc) * skew(bodies(1).Xgc), sqrt(eps) * bodies(1).dm);
%!   assert(bodies(2).J, diag([2.1, 2.2, 2.3]) - bodies(2).dm * skew(bodies(2).Xgc) * skew(bodies(2).Xgc), sqrt(eps) * bodies(2).dm);
%!   assert(bodies(1).type, "body");
%!   assert(bodies(2).type, "body");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_log_body_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 20;\n");
%!     fputs(fd, " set: real m1 = 1;\n");
%!     fputs(fd, " set: real m2 = 2;\n");
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
%!     fputs(fd, "         nonlinear solver: line search, default solver options, heavy nonlinear, divergence check, no;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1., 2., 3.,\n");
%!     fputs(fd, "                 euler123, 0.1, 0.2, 0.3,\n");
%!     fputs(fd, "                 10., 20., 30.,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, null, 1., 0., 0, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, null, 0., 1., 0, F2;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
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
%!   bodies = mbdyn_post_load_log_body(fname);
%!   assert([bodies.label], int32([1,2]));
%!   assert([bodies.node], int32([1,2]));
%!   assert([bodies.dm], [1, 2]);
%!   assert([bodies.Xgc], [0.1, 0.01; 0.2, 0.02; 0.3, 0.03]);
%!   assert(bodies(1).J, diag([1.1, 1.2, 1.3]) - bodies(1).dm * skew(bodies(1).Xgc) * skew(bodies(1).Xgc), sqrt(eps) * bodies(1).dm);
%!   assert(bodies(2).J, diag([2.1, 2.2, 2.3]) - bodies(2).dm * skew(bodies(2).Xgc) * skew(bodies(2).Xgc), sqrt(eps) * bodies(2).dm);
%!   assert(bodies(1).type, "body");
%!   assert(bodies(2).type, "body");
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_post_load_log_body_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     fputs(fd, " set: real g = 9.81;\n");
%!     fputs(fd, " set: real F1 = 100;\n");
%!     fputs(fd, " set: real F2 = 20;\n");
%!     fputs(fd, " set: real m1 = 1;\n");
%!     fputs(fd, " set: real m2 = 2;\n");
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
%!     fputs(fd, "     use automatic differentiation;\n");

%!     fputs(fd, "     structural nodes: 2;\n");
%!     fputs(fd, "     rigid bodies: 2;\n");
%!     fputs(fd, "     forces: 2;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, "         structural: 2, dynamic,\n");
%!     fputs(fd, "                 1., 2., 3.,\n");
%!     fputs(fd, "                 euler123, 0.1, 0.2, 0.3,\n");
%!     fputs(fd, "                 10., 20., 30.,\n");
%!     fputs(fd, "                 null;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, m1, 0.1, 0.2, 0.3, diag, 1.1, 1.2, 1.3;\n");
%!     fputs(fd, "         body: 2, 2, m2, 0.01, 0.02, 0.03, diag, 2.1, 2.2, 2.3;\n");
%!     fputs(fd, "         force: 1, absolute, 1, position, null, 1., 0., 0, F1;\n");
%!     fputs(fd, "         force: 2, absolute, 2, position, null, 0., 1., 0, F2;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., g;\n");
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
%!   bodies = mbdyn_post_load_log_body(fname);
%!   assert([bodies.label], int32([1,2]));
%!   assert([bodies.node], int32([1,2]));
%!   assert([bodies.dm], [1, 2]);
%!   assert([bodies.Xgc], [0.1, 0.01; 0.2, 0.02; 0.3, 0.03]);
%!   assert(bodies(1).J, diag([1.1, 1.2, 1.3]) - bodies(1).dm * skew(bodies(1).Xgc) * skew(bodies(1).Xgc), sqrt(eps) * bodies(1).dm);
%!   assert(bodies(2).J, diag([2.1, 2.2, 2.3]) - bodies(2).dm * skew(bodies(2).Xgc) * skew(bodies(2).Xgc), sqrt(eps) * bodies(2).dm);
%!   assert(bodies(1).type, "body");
%!   assert(bodies(2).type, "body");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
