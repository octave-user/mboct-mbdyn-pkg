function path_unix = mbdyn_convert_path(path_input)
  path_unix = path_input;

  if (ispc())
    path_unix(find(path_unix == '\')) = '/';
  endif
endfunction
