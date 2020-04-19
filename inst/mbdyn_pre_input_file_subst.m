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
## @deftypefn {Function File} mbdyn_pre_input_file_subst(@var{mbdyn_filename}, @var{log_filename}, @var{output_filename})
##
## This function parses the MBDyn input file <@var{mbdyn_filename}> and replaces all occurrences of variable names with the values found in the file <@var{log_filename}>.
##
## @var{log_filename} @dots{} MBDyn .log file where all variable values are stored.
##
## @var{output_filename} @dots{} New MBDyn file where all variables are substituted.
##
## @end deftypefn

function mbdyn_pre_input_file_subst(mbdyn_filename, log_filename, output_filename, verbose)
  if (nargin < 3 || nargin > 4)
    print_usage();
  endif

  if (nargin < 4)
    verbose = false;
  endif

  vars = mbdyn_post_load_log_vars(log_filename);

  fout = -1;

  unwind_protect
    [fout, msg] = fopen(output_filename,"wt");

    if (fout == -1)
      error("could not open file \"%s\": %s",output_filename,msg);
    endif

    parse_mbdyn_input_file(mbdyn_filename,vars,fout,verbose);
  unwind_protect_cleanup
    if (fout ~= -1)
      fclose(fout);
    endif
  end_unwind_protect
endfunction

function parse_mbdyn_input_file(mbdyn_filename,vars,fout,verbose)
  persistent delim = "+-*/() \t\",;^[]{}:=#~\r\n";

  fin = -1;

  unwind_protect
    [fin, msg] = fopen(mbdyn_filename, "rt");

    if (fin == -1)
      error("could not open file \"%s\": %s", mbdyn_filename, msg);
    endif

    line_number = 0;
    empty_line = false;

    while (1)
      line = fgets(fin);

      if (~ischar(line) && line == -1)
        break;
      endif

      ++line_number;
      output = "";
      tok = "";
      in_tok = false;
      num_tok = 0;
      is_set = false;
      is_node = false;
      is_dof = false;
      is_include = false;
      include_file = "";
      in_quote = false;
      in_number = false;
      is_blank = false;
      number_end = [];

      i = 0;

      while (++i <= length(line))
        switch (line(i))
          case "\""
            in_quote = ~in_quote;
        endswitch

        switch (line(i))
          case {"\t", " "}
            is_blank = true;
          otherwise
            is_blank = false;
        endswitch

        is_delim = false;

        if (in_quote || is_blank)
          in_number = false;
        else
          [val, count, errmsg, pos] = sscanf(line(i:end), "%g");
          if (count > 0)
            in_number = true;
            number_end = i + pos - 1;
          endif
        endif

        if (in_number)
          if (i >= number_end)
            in_number = false;
            number_end = [];
          endif
        endif

        if (~in_number)
          for j=1:length(delim)
            if(line(i) == delim(j))
              is_delim = true;
              break;
            endif
          endfor
        endif

        if (is_delim)
          if (in_tok)
            if (length(tok) > 0)
              if (verbose)
                fprintf(stderr,"%s:%d:%d>%s\n",mbdyn_filename,line_number,num_tok,tok);
              endif

              if (num_tok == 1 && strcmp(tok, "set"))
                is_set = true;
              endif

              if (num_tok == 2 && strcmp(tok, "node"))
                is_node = true;
              endif

              if (num_tok == 2 && strcmp(tok, "dof"))
                is_dof = true;
              endif

              if (num_tok == 1 && strcmp(tok, "include"))
                is_include = true;
              endif

              if (num_tok == 2 && is_include)
                include_file = tok;
              endif

              if (~is_builtin(tok))
                if (isfield(vars, tok) && ~in_number)
                  val = getfield(vars,tok);

                  if (fix(val) == val)
                    format = "(%d)";
                  else
                    format = "(%g)";
                  endif

                  tok = sprintf(format, val);
                elseif (verbose && ~iskeyword(tok))
                  fprintf(stderr, "warning: variable=\"%s\" not found!", tok);
                endif
              endif
              output = cstrcat(output, tok);
            endif
            tok = "";
            in_tok = false;
          endif

          if (strcmp(line(i), "#"))
            output = cstrcat(output, "\n");
            break;
          endif
          output = cstrcat(output, delim(j));
        else
          if (~in_tok)
            in_tok = true;
            ++num_tok;
          endif
          tok(end + 1) = line(i);
        endif
      endwhile

      if (is_include)
        parse_mbdyn_input_file(include_file, vars, fout, verbose);
      else
        if (~is_set || is_dof || is_node)
          if (length(output) == 1 && output(1) == "\n")
            if (empty_line)
              output = "";
            endif
            empty_line = true;
          else
            empty_line = false;
          endif
          if (length(output) > 0)
            fputs(fout, output);
          endif
        endif
      endif
    endwhile
  unwind_protect_cleanup
    if (fin ~= -1)
      fclose(fin);
    endif
  end_unwind_protect
endfunction

function f_builtin = is_builtin(tok)
  switch(tok)
    case {"Var",
          "Time"}
      f_builtin = true;
    otherwise
      f_builtin = false;
  endswitch
endfunction

function f_keyword = is_keyword(tok)
  switch(tok)
    case {"set",
          "real",
          "integer",
          "begin",
          "end",
          "include",
          "reference",
          "null",
          "control",
          "data",
          "elements",
          "nodes",
          "sin",
          "cos",
          "tan",
          "atan",
          "exp",
          "log"}
      f_keyword = true;
    otherwise
      f_keyword = false;
  endswitch
endfunction
