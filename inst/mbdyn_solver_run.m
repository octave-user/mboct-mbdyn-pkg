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
## @deftypefn {Function File} @var{info} = mbdyn_solver_run (@var{mbdyn_filename}, @var{every}, @var{showAll}, @var{f_run_mbdyn}, @var{f_run_mbdyn2easyanim}, @var{f_runEasyAnim}, @var{working_directory}, @var{logfile}, @var{verbose}, @var{output_file}, @var{f_pedantic}, @var{mbdyn_command})
## @deftypefnx {} @var{info} = mbdyn_solver_run (@var{mbdyn_filename}, @var{options})
##
## Runs the multibody dynamics software MBDyn.
## For the first form all parameters are passed as separate arguments.
## For the second form, all parameters except the input file are passed as
## members of the struct argument <@var{options}>.
##
## @var{mbdyn_filename} @dots{} MBDyn input file without extension .mbdyn.
##
## @var{every} @dots{} Timestep interval which is processed by mbdyn2easyanim.
##
## @var{showAll} @dots{} If true all items (beams, joints) are displayed by EasyAnim.
##
## @var{f_run_mbdyn} @dots{} If true the command mbdyn is executed.
##
## @var{f_run_mbdyn2easyanim} @dots{} If true the command mbdyn2easyanim is executed after mbdyn.
##
## @var{f_runEasyAnim} @dots{} If true the command EasyAnim is executed.
##
## @var{working_directory} @dots{} Working directory where all commands are executed.
##
## @var{logfile} @dots{} If numel(logfile) > 0 all output of MBDyn is redirected to that file.
##
## @var{verbose} @dots{} If true the command is displayed before it is executed.
##
## @var{mbdyn_command} @dots{} The filename of mbdyn's binary.
##
## @var{output_file} @dots{} File name for output data.
##
## @var{options} @dots{} Data structure containing all of the previous arguments as fields.
##
## @end deftypefn

function info = mbdyn_solver_run(mbdyn_filename, every_or_options = 1, showAll = 1, f_run_mbdyn = true, f_run_mbdyn2easyanim = true, ...
                                 f_runEasyAnim = false, working_directory = pwd(), logfile = "", verbose = true, output_file = "", ...
                                 f_pedantic = false, mbdyn_command = "mbdyn", f_silent = false)

  info.total_steps = int32(-1);
  info.total_iter = int32(-1);
  info.total_jac = int32(-1);
  info.total_err = -1;
  info.total_cpu = -1;

  if (~exist("mbdyn_filename", "var"))
    print_usage();
  endif

  every = 1;

  if (isstruct(every_or_options))
    options = every_or_options;

    if (isfield(options, "f_silent"))
      f_silent = options.f_silent;
    else
      f_silent = false;
    endif

    if (isfield(options, "every"))
      every = options.every;
    endif

    if (isfield(options, "showAll"))
      showAll = options.showAll;
    endif

    if (isfield(options, "f_run_mbdyn"))
      f_run_mbdyn = options.f_run_mbdyn;
    endif

    if (isfield(options, "f_run_mbdyn2easyanim"))
      f_run_mbdyn2easyanim = options.f_run_mbdyn2easyanim;
    endif

    if (isfield(options, "f_runEasyAnim"))
      f_runEasyAnim = options.f_runEasyAnim;
    endif

    if (isfield(options, "working_directory"))
      working_directory = options.working_directory;
    endif

    if (isfield(options, "logfile"))
      logfile = options.logfile;
    endif

    if (isfield(options, "verbose"))
      verbose = options.verbose;
    endif

    if (isfield(options, "output_file"))
      output_file = options.output_file;
    endif

    if (isfield(options, "f_pedantic"))
      f_pedantic = options.f_pedantic;
    else
      f_pedantic = false;
    endif

    if (isfield(options, "mbdyn_command"))
      mbdyn_command = options.mbdyn_command;
    endif
  else
    every = every_or_options;
  endif

  logfile = mbdyn_convert_path(logfile);
  output_file = mbdyn_convert_path(output_file);
  mbdyn_filename = mbdyn_convert_path(mbdyn_filename);
  working_directory = mbdyn_convert_path(working_directory);

  cwd = pwd();

  unwind_protect
    cd(working_directory);

    if (numel(logfile))
      redirect = sprintf("> '%s' 2>&1", mbdyn_convert_path(logfile));
    else
      redirect = "";
    endif

    if (length(output_file) > 0)
      output_flag = sprintf("-o '%s'", mbdyn_convert_path(output_file));
    else
      output_flag = "";
      output_file = mbdyn_convert_path(mbdyn_post_output_filename(mbdyn_filename));
    endif

    if (f_pedantic)
      pedantic_flag = "--pedantic";
    else
      pedantic_flag = "";
    endif

    if (f_silent)
      silent_flag = "--silent";
    else
      silent_flag = "";
    endif

    mbdyn_path_init();
    
    if (f_run_mbdyn)
      run_command(sprintf("%s %s %s -f \"%s\" %s %s", ...
                          mbdyn_command, ...
                          pedantic_flag, ...
                          silent_flag, ...
                          mbdyn_filename, ...
                          output_flag, ...
                          redirect), verbose);
    endif

    if (numel(logfile) && exist(logfile, "file"))
      fd = -1;
      persistent cmd = "\
BEGIN {\
total_steps=0;\
total_iter=0;\
total_jac=0;\
total_err=0.;\
total_cpu=0.;\
};\
/^End of simulation at time/{ total_steps = $8; };\
/^total iterations:/{ total_iter = $3; };\
/^total Jacobian matrices:/{ total_jac = $4; };\
/^total error:/{  total_err = $3; };\
/^Elapsed time/{ total_cpu = $3; };\
/^The simulation required/{ total_cpu = $4; };\
END {\
printf(\"%d\t%d\t%d\t%g\t%g\\n\", total_steps, total_iter, total_jac, total_err, total_cpu);\
};";
      
      [status, output] = shell(sprintf("exec awk -F ' ' '%s' '%s'", cmd, mbdyn_convert_path(logfile)), true, "sync");

      if (status == 0 && numel(output))
        [info.total_steps, ...
         info.total_iter, ...
         info.total_jac, ...
         info.total_err, ...
         info.total_cpu, ...
         count, ...
         err_msg] = sscanf(output, "%d\t%d\t%d\t%g\t%g", "C");

        if (verbose && nargout == 0)
          fprintf(stderr, "MBDyn terminated\n");
          fprintf(stderr, "total steps: %d\n", info.total_steps);
          fprintf(stderr, "total iterations: %d\n", info.total_iter);
          fprintf(stderr, "total Jacobian matrices: %d\n", info.total_jac);
          fprintf(stderr, "total error: %g\n", info.total_err);
          fprintf(stderr, "total CPU time: %d\n", info.total_cpu);
        endif
      endif
    endif

    [info_mov, err_mov] = stat([output_file, ".mov"]);
    [info_log, err_log] = stat([output_file, ".log"]);

    if (err_mov == 0 && err_log == 0)
      if (f_run_mbdyn2easyanim)
        if (numel(logfile) && f_run_mbdyn)
          redirect = sprintf(">> \"%s\" 2>&1", logfile);
        else
          redirect = "";
        endif
        run_command(sprintf("mbdyn2easyanim.sh -v every=%d -v showAll=%d \"%s\" %s", every, showAll, output_file, redirect), verbose);
      endif

      if (f_runEasyAnim)
        run_command(sprintf("EasyAnimm \"%s.vol\"", output_file), verbose);
      endif
    endif
  unwind_protect_cleanup
    cd(cwd);
  end_unwind_protect
endfunction

function run_command(command, verbose)
  if (verbose)
    fprintf(stderr, "command: \"%s\" is executed in %s ...\n", command, pwd());
  endif

  command = ["exec ", command];

  ## On Windows systems the system function always starts cmd.exe
  ## which cannot be used to execute mbdyn-launch.sh
  status = shell(command, false, "sync");

  if (0 ~= status)
    error("command '%s' failed with status %d", command, status);
  endif
endfunction

%!test
%! fd = -1;
%! unwind_protect
%!   unwind_protect
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_solver_run_XXXXXX"));
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
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., 9.81;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.logfile = [fname, ".stdout"];
%!   options.verbose = false;
%!   mbdyn_solver_run(fname, options);
%!   [t, dt, niter, reserr, solerr, solconv, output_flag] = mbdyn_post_load_output_out(options.output_file, 1024, false);
%!   [node_id, node_data] = mbdyn_post_load_output([options.output_file, ".mov"], 18, [1], numel(t), 18, 1, false);
%!   g = -9.81;
%!   assert(node_id, int32(1));
%!   assert(t, (0:0.1:1.1).', 1e-5);
%!   assert(dt, [0; repmat(0.1, numel(t) - 1, 1)]);
%!   assert(node_data{1}, ...
%!          [zeros(numel(t), 2), ...
%!          0.5 * g * t.^2, ...
%!          zeros(numel(t), 5), ...
%!          g * t, ...
%!          zeros(numel(t), 5), ...
%!          repmat(g, numel(t), 1), ...
%!          zeros(numel(t), 3)], ...
%!          1e-5);
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
%!     [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_solver_run_XXXXXX"));
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
%!     fputs(fd, " end: initial value;\n");
%!     fputs(fd, " begin: control data;\n");
%!     fputs(fd, "     structural nodes: 1;\n");
%!     fputs(fd, "     rigid bodies: 1;\n");
%!     fputs(fd, "     gravity;\n");
%!     fputs(fd, " end: control data;\n");
%!     fputs(fd, " begin: nodes;\n");
%!     fputs(fd, "         structural: 1, dynamic,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 eye,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 null,\n");
%!     fputs(fd, "                 accelerations, yes;\n");
%!     fputs(fd, " end: nodes;\n");
%!     fputs(fd, " begin: elements;\n");
%!     fputs(fd, "         body: 1, 1, 1., null, diag, 1., 1., 1.;\n");
%!     fputs(fd, " gravity: uniform, 0., 0., -1., 9.81;\n");
%!     fputs(fd, " end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   options.output_file = fname;
%!   options.logfile = [fname, ".stdout"];
%!   options.verbose = false;
%!   mbdyn_solver_run(fname, options);
%!   [t, dt, niter, reserr, solerr, solconv, output_flag] = mbdyn_post_load_output_out(options.output_file, 1024, false);
%!   [node_id, node_data] = mbdyn_post_load_output([options.output_file, ".mov"], 18, [1], numel(t), 18, 1, false);
%!   g = -9.81;
%!   assert(node_id, int32(1));
%!   assert(t, (0:0.1:1.1).', 1e-5);
%!   assert(dt, [0; repmat(0.1, numel(t) - 1, 1)]);
%!   assert(node_data{1}, ...
%!          [zeros(numel(t), 2), ...
%!          0.5 * g * t.^2, ...
%!          zeros(numel(t), 5), ...
%!          g * t, ...
%!          zeros(numel(t), 5), ...
%!          repmat(g, numel(t), 1), ...
%!          zeros(numel(t), 3)], ...
%!          1e-5);
%!   figure("visible", "off");
%!   subplot(3, 1, 1);
%!   plot(t, node_data{1}(:, 3), "-;x(t);1");
%!   xlabel("t [s]");
%!   ylabel("x [m]");
%!   grid on;
%!   grid minor on;
%!   title("trajectory");
%!   subplot(3, 1, 2);
%!   plot(t, node_data{1}(:, 9), "-;v(t);1");
%!   xlabel("t [s]");
%!   ylabel("v [m/s]");
%!   grid on;
%!   grid minor on;
%!   title("velocity");
%!   subplot(3, 1, 3);
%!   plot(t, node_data{1}(:, 15), "-;a(t);1");
%!   xlabel("t [s]");
%!   ylabel("a [m/s^2]");
%!   grid on;
%!   grid minor on;
%!   title("acceleration");
%! unwind_protect_cleanup
%!   if (fd ~= -1)
%!     unlink(fname);
%!     files = dir([fname, ".*"]);
%!     for i=1:numel(files)
%!       unlink(fullfile(files(i).folder, files(i).name));
%!     endfor
%!   endif
%! end_unwind_protect
