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
## @deftypefn {Function File} mbdyn_path_init()
##
## Set the environment variables needed to run mbdyn, mbdyn2easyanim.sh and abs2rel.awk.
##
## @end deftypefn

function mbdyn_path_init()  
  [status, output] = shell("exec which mbdyn 2>&1", true, "sync");

  if (status ~= 0)
    env_path = getenv("PATH");

    install_paths = {"/usr/local/mbdyn/bin", "/usr/local/bin", "/opt/mbdyn/bin"};

    if (ispc())
      exe_suffix = ".exe";
    else
      exe_suffix = "";
    endif

    for i=1:numel(install_paths)
      idx = strfind(env_path, install_paths{i});

      if (numel(idx))
        continue;
      endif

      [info, err, msg] = stat(fullfile(install_paths{i}, ["mbdyn", exe_suffix]));

      if (err == 0 && S_ISREG(info.mode))
        putenv("PATH", [env_path, ":", install_paths{i}]);
        break;
      endif
    endfor
  endif

  [status, output] = shell("exec awk -g -f abs2rel.awk 2>&1", true, "sync");

  if (status ~= 0)
    install_paths = {"/usr/local/mbdyn/share/awk", "/usr/local/share/awk", "/opt/mbdyn/share/awk"};

    for i=1:numel(install_paths)
      awk_path = getenv("AWKPATH");
      
      idx = strfind(awk_path, install_paths{i});

      if (numel(idx))
        continue;
      endif

      [info, err, msg] = stat(fullfile(install_paths{i}, "abs2rel.awk"));

      if (ispc() || (err == 0 && S_ISREG(info.mode)))
        putenv("AWKPATH", [awk_path, ":", install_paths{i}]);
	if (~ispc())
          break;
	endif
      endif
    endfor
  endif
endfunction
