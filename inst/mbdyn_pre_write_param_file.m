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
## @deftypefn {Function File} mbdyn_pre_write_param_file(@var{output_file}, @var{param},  @var{options})
##
## Write all scalar variables in struct <@var{param}> to file "<@var{output_file}>".
##
## @var{output_file} @dots{} Output file which can be included in an MBDyn input file.
##
## @var{param} @dots{} Data structure with scalar variables.
##
## @var{options}.open_mode @dots{} Argument to be passed to fopen.
##
## @end deftypefn

function mbdyn_pre_write_param_file(output_file, param, options)
  if (nargin < 2 || nargin > 3)
    print_usage();
  endif

  if (~ischar(output_file) && ~isscalar(output_file))
    error("output_file must be a filename or a file descriptor!");
  endif

  if (~isstruct(param))
    error("param must be a struct!");
  endif

  if (nargin < 3)
    options = struct();
  endif

  if (~isstruct(options))
    error("options must be a struct!");
  endif

  if (~isfield(options, "open_mode"))
    options.open_mode = "wt";
  endif

  fout = -1;
  owns_fd = false;

  unwind_protect
    if (ischar(output_file))
      owns_fd = true;

      [fout, msg] = fopen(output_file, options.open_mode);

      if (fout == -1)
        error("could not open file \"%s\": %s",  output_file, msg);
      endif
    else
      fout = output_file;
    endif

    param_names = fieldnames(param);

    for i=1:length(param_names)
      param_val = getfield(param, param_names{i});
      if (isscalar(param_val))
        if (isinteger(param_val))
          fprintf(fout, "set: integer %s = %d;\n", param_names{i}, param_val);
        elseif (isbool(param_val))
          fprintf(fout, "set: bool %s = %s;\n",  param_names{i},  {"FALSE", "TRUE"}{param_val + 1});
        elseif (isreal(param_val))
          if (param_val > realmax)
            param_val = realmax;
          elseif (param_val < -realmax)
            param_val = -realmax;
          endif
          fprintf(fout, "set: real %s = %.16e;\n", param_names{i}, param_val);
        endif
      elseif (ischar(param_val))
        fprintf(fout, "set: string %s = \"%s\";\n", param_names{i}, param_val);
      endif
    endfor
  unwind_protect_cleanup
    if (owns_fd && fout ~= -1)
      fclose(fout);
    endif
  end_unwind_protect
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_write_param_file_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     param.x = -1.5;
%!     param.i = int32(2);
%!     param.b1 = true;
%!     param.b2 = false;
%!     param.s = "123 456 789";
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "         problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "   initial time: 0;\n");
%!     fputs(fd, "   final time: 1;\n");
%!     fputs(fd, "   time step: 0.1;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "   abstract nodes: 1;\n");
%!     fputs(fd, "   genels: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, " abstract: 1, differential;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, " genel: 1, mass, 1, abstract, algebraic, 1.;\n");
%!     fputs(fd, " genel: 2, spring support, 1, abstract, algebraic, 1.;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!     opts.output_file = fname;
%!     opts.verbose = false;
%!     opts.logfile = [fname, ".stdout"];
%!     mbdyn_solver_run(fname, opts);
%!     log_dat = mbdyn_post_load_log(opts.output_file);
%!     assert(log_dat.vars.x, param.x);
%!     assert(log_dat.vars.i, param.i);
%!     assert(log_dat.vars.b1, param.b1);
%!     assert(log_dat.vars.b2, param.b2);
%!     assert(log_dat.vars.s, param.s);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!   endif
%!   files = dir([fname, ".*"]);
%!   for i=1:numel(files)
%!     unlink(fullfile(files(i).folder, files(i).name));
%!   endfor
%! end_unwind_protect

%!demo
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_write_param_file_XXXXXX"));
%!     if (fd == -1)
%!       error("failed to open temporary file");
%!     endif
%!     param.x = -1.5;
%!     param.i = int32(2);
%!     param.b1 = true;
%!     param.b2 = false;
%!     param.s = "123 456 789";
%!     mbdyn_pre_write_param_file(fd, param);
%!     fputs(fd, " begin: data;\n");
%!     fputs(fd, "   problem: initial value;\n");
%!     fputs(fd, " end: data;\n");
%!     fputs(fd, " begin: initial value;\n");
%!     fputs(fd, "   initial time: 0;\n");
%!     fputs(fd, "   final time: 1;\n");
%!     fputs(fd, "   time step: 0.1;\n");
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "   abstract nodes: 1;\n");
%!     fputs(fd, "   genels: 2;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, " abstract: 1, differential;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, " genel: 1, mass, 1, abstract, algebraic, 1.;\n");
%!     fputs(fd, " genel: 2, spring support, 1, abstract, algebraic, 1.;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!     opts.output_file = fname;
%!     opts.verbose = false;
%!     opts.logfile = [fname, ".stdout"];
%!     mbdyn_solver_run(fname, opts);
%!     log_dat = mbdyn_post_load_log(opts.output_file);
%!     assert(log_dat.vars.x, param.x);
%!     assert(log_dat.vars.i, param.i);
%!     assert(log_dat.vars.b1, param.b1);
%!     assert(log_dat.vars.b2, param.b2);
%!     assert(log_dat.vars.s, param.s);
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!   endif
%!   files = dir([fname, ".*"]);
%!   for i=1:numel(files)
%!     unlink(fullfile(files(i).folder, files(i).name));
%!   endfor
%! end_unwind_protect
