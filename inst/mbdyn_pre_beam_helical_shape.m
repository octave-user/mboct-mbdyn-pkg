## Copyright (C) 2015(-2020) Reinhard <octave-user@a1.net>
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
## @deftypefn {Function File} [@var{X}, @var{shape}] = mbdyn_pre_beam_helical_shape(@var{options})
## Compute coordinates at the center-line of a helical spring.
##
## @var{options}.na @dots{} Number of active coils.
##
## @var{options}.ni @dots{} Number of inactive coils.
##
## @var{options}.L @dots{} Free length.
##
## @var{options}.d @dots{} Wire diameter.
##
## @var{options}.D @dots{} Mean coil diameter.
##
## @var{options}.ground_surface @dots{} Angle which defines the flat area of ground end coils.
##
## @seealso{mbdyn_pre_beam_helical_compute}
## @end deftypefn

function [X, shape] = mbdyn_pre_beam_helical_shape(options)
  if (nargin ~= 1 || nargout > 2)
    print_usage();
  endif

  if (~isfield(options, "type"))
    options.type = "variable diameter, variable pitch";
  endif

  if (~isfield(options, "interpolation"))
    options.interpolation = 'spline';
  endif

  if (~isfield(options, "interpolation_z"))
    options.interpolation_z = options.interpolation;
  endif

  if (~isfield(options, "interpolation_D"))
    options.interpolation_D = options.interpolation;
  endif

  if (~isfield(options, "D"))
    if (isfield(options, "Di") && isfield(options, "d"))
      options.D = options.Di + options.d;
    endif
  endif

  if (~isfield(options, "helix_direction"))
    options.helix_direction = 1;
  endif

  if (~isfield(options, "inactive_coil_alignment"))
    options.inactive_coil_alignment = true;
  endif

  if (isfield(options, "ground_surface"))
    switch (options.type)
      case {"inactive coils, const pitch, const diameter", ...
            "inactive coils, const pitch, variable diameter", ...
            "const pitch, const diameter", ...
            "const pitch, variable diameter"}
        ## Ok
      otherwise
        error("option \"ground_surface\" not implemented for type \"%s\"", options.type);
    endswitch
  else
    options.ground_surface = 0;
  endif

  switch(options.type )
    case "const pitch, const diameter"
      if (~isfield(options, "D"))
        error("options.D is required for type %s", options.type);
      endif

      if (~isfield(options, "d"))
        error("options.d is required for type %s", options.type);
      endif

      if (~isfield(options, "L"))
        error("options.L is required for type %s", options.type);
      endif

      if (~isfield(options, "na"))
        error("options.na is required for type %s", options.type);
      endif

      if (~isfield(options, "ni"))
        if (~isfield(options, "ni1") || ~isfield(options, "ni2"))
          error("options.ni is required for type %s", options.type);
        endif
      endif

      if (~isfield(options, "ni1"))
        options.ni1 = options.ni;
      endif

      if (~isfield(options, "ni2"))
        options.ni2 = options.ni;
      endif

      shape.L = options.L;

      dhg = options.d * options.ground_surface / (2 * pi);

      options.z = [ options.d/2 - dhg + options.d * options.ni1, ...
                    options.L - options.d/2 + dhg - options.d * options.ni2 ];

      options.D = repmat(options.D, 1, 2);

      options.Phi = options.helix_direction * 2 * pi * [ options.inactive_coil_alignment * options.ni1, ...
                                                         options.inactive_coil_alignment * options.ni1 + options.na ];

      options.type = "variable diameter, variable pitch";
      options.interpolation_z = 'linear';
      options.interpolation_D = 'linear';

    case {"inactive coils, const pitch, const diameter", "inactive coils, const pitch, variable diameter"}
      if (~isfield(options, "D"))
        error("options.D is required for type %s", options.type);
      endif

      if (~isfield(options, "d"))
        error("options.d is required for type %s", options.type);
      endif

      if (~isfield(options, "L"))
        error("options.L is required for type %s", options.type);
      endif

      if (~isfield(options, "na"))
        error("options.na is required for type %s", options.type);
      endif

      if (~isfield(options, "ni"))
        if (~isfield(options, "ni1") || ~isfield(options, "ni2"))
          error("options.ni is required for type %s", options.type);
        endif
      endif

      if (~isfield(options, "ni1"))
        options.ni1 = options.ni;
      endif

      if (~isfield(options, "ni2"))
        options.ni2 = options.ni;
      endif

      shape.L = options.L;

      dhg = options.d * options.ground_surface / (2 * pi);

      options.z = [ options.d/2 - dhg, ...
                    options.d/2 - dhg + options.d * options.ni1, ...
                    options.L - options.d/2 + dhg - options.d * options.ni2, ...
                    options.L - options.d/2 + dhg];

      switch (options.type)
        case "inactive coils, const pitch, variable diameter"
          if (isfield(options, "z_rel_D"))
            options.z_D = options.z_rel_D * options.L;
          elseif (~isfield(options, "z_D"))
            error("missing field z_D required for type \"%s\"", options.type);
          endif
        otherwise
          options.D = repmat(options.D,1,4);
      endswitch

      options.Phi = options.helix_direction * 2 * pi * [ 0, ...
                                                         options.ni1, ...
                                                         options.ni1 + options.na, ...
                                                         options.ni1 + options.na + options.ni2 ];

      options.type = "variable diameter, variable pitch";
      options.interpolation_z = 'linear';

      if (~isfield(options, "interpolation_D"))
        options.interpolation_D = 'linear';
      endif

    case {"const pitch, variable diameter"}
      if (~isfield(options, "D"))
        error("options.D is required for type %s", options.type);
      endif

      if (~isfield(options, "d"))
        error("options.d is required for type %s", options.type);
      endif

      if (~isfield(options, "L"))
        error("options.L is required for type %s", options.type);
      endif

      if (~isfield(options, "na"))
        error("options.na is required for type %s", options.type);
      endif

      if (~isfield(options, "ni"))
        if (~isfield(options, "ni1") || ~isfield(options, "ni2"))
          error("options.ni is required for type %s", options.type);
        endif
      endif

      if (~isfield(options, "ni1"))
        options.ni1 = options.ni;
      endif

      if (~isfield(options, "ni2"))
        options.ni2 = options.ni;
      endif

      shape.L = options.L;

      dhg = options.d * options.ground_surface / (2 * pi);

      options.z = [options.d/2 - dhg + options.d * options.ni1, ...
                   options.L - options.d/2 + dhg - options.d * options.ni2];

      if (isfield(options, "z_rel_D"))
        options.z_D = options.z_rel_D * options.L;
      elseif (~isfield(options, "z_D"))
        error("missing field z_D required for type \"%s\"", options.type);
      endif

      options.Phi = options.helix_direction * 2 * pi * [ options.ni1, ...
                                                         options.ni1 + options.na ];

      options.type = "variable diameter, variable pitch";
      options.interpolation_z = 'linear';

      if (~isfield(options, "interpolation_D"))
        options.interpolation_D = 'linear';
      endif
    case "const pitch, 3 diameters"

      if (~isfield(options, "D"))
        error("options.D is required for type %s", options.type);
      endif

      if (length(options.D) ~= 3 )
        error("length(D) must be three for type %d", options.type);
      endif

      if (~isfield(options, "d"))
        error("options.d is required for type %s", options.type);
      endif

      if (~isfield(options, "L"))
        error("options.L is required for type %s", options.type);
      endif

      if (~isfield(options, "na"))
        error("options.na is required for type %s", options.type);
      endif

      if (~isfield(options, "ni"))
        if (~isfield(options, "ni1") || ~isfield(options, "ni2"))
          error("options.ni is required for type %s", options.type);
        endif
      endif

      if (~isfield(options, "ni1"))
        options.ni1 = options.ni;
      endif

      if (~isfield(options, "ni2"))
        options.ni2 = options.ni;
      endif

      shape.L = options.L;

      options.z = [ options.d/2 + options.d * options.ni1, ...
                    (options.d * (options.ni1 - options.ni2) + options.L)/2, ...
                    options.L - options.d/2 - options.d * options.ni2 ];

      options.Phi = options.helix_direction * 2 * pi * [ options.inactive_coil_alignment * options.ni1, ...
                                                         options.inactive_coil_alignment * options.ni1 + options.na / 2., ...
                                                         options.inactive_coil_alignment * options.ni1 + options.na ];

      options.type = "variable diameter, variable pitch";

      options.interpolation_z = 'linear';
      options.interpolation_D = 'pchip';

    case { "relative variable pitch, variable diameter", "relative variable pitch, const diameter" }

      if (~isfield(options, "D"))
        error("options.D is required for type %s", options.type);
      endif

      if (~isfield(options, "z_rel"))
        error("options.z_rel is required for type %s", options.type);
      endif

      switch (options.type )
        case "relative variable pitch, const diameter"
          options.D = repmat(options.D, 1, length(options.z_rel));
      endswitch

      if (~isfield(options, "z_D"))
        if (~isfield(options, "z_rel_D"))
          options.z_rel_D = linspace(0, 1, length(options.D));
        endif
      endif

      if (~isfield(options, "d"))
        error("options.d is required for type %s", options.type);
      endif

      if (~isfield(options, "L"))
        error("options.L is required for type %s", options.type);
      endif

      if (~isfield(options, "na"))
        error("options.na is required for type %s", options.type);
      endif

      if (~isfield(options, "ni"))
        if (~isfield(options, "ni1") || ~isfield(options, "ni2"))
          error("options.ni is required for type %s", options.type);
        endif
      endif

      if (~isfield(options, "ni1"))
        options.ni1 = options.ni;
      endif

      if (~isfield(options, "ni2"))
        options.ni2 = options.ni;
      endif

      if (max(options.z_rel) > 1 || min(options.z_rel) < 0 )
        error("options.z_rel must be between zero and one for type %s", options.type);
      endif

      shape.L = options.L;

      options.z = options.d/2 + options.d * options.ni1 ...
                  + options.z_rel * (options.L - options.d * (1 + options.ni1 + options.ni2 ));

      if (~isfield(options, "z_D"))
        options.z_D = options.d/2 + options.d * options.ni1 ...
                      + options.z_rel_D * (options.L - options.d * (1 + options.ni1 + options.ni2 ));
      endif

      options.Phi = options.helix_direction * 2 * pi ...
                    * linspace(options.inactive_coil_alignment * options.ni1, ...
                               options.inactive_coil_alignment * options.ni1 + options.na, length(options.z));

      options.type = "variable diameter, variable pitch";
  endswitch

  if (~isfield(options, "Phi"))
    error("options.Phi is required for type %s", options.type);
  endif

  if (~isfield(options, "N"))
    if (~isfield(options, "N_coil"))
      options.N_coil = 3600;
    endif

    options.N = options.N_coil * abs(options.Phi(end) - options.Phi(1)) / (2 * pi );
  endif

  if (~isfield(options, "N_begin" ))
    options.N_begin = 20;
  endif

  if (~isfield(options, "N_end"))
    options.N_end = options.N_begin;
  endif

  shape.Phi = linspace(options.Phi(1), options.Phi(end), options.N);

  if (length(shape.Phi >= 4))
    shape.Phi = [ linspace(shape.Phi(1), shape.Phi(2), options.N_begin), ...
                  shape.Phi(3:end-2), ...
                  linspace(shape.Phi(end-1), shape.Phi(end), options.N_end) ];
  endif

  switch(options.type )
    case "variable diameter, variable pitch"
      if (~isfield(options, "z"))
        error("options.z is required for type %s", options.type);
      endif

      if (~isfield(options, "D"))
        error("options.D is required for type %s", options.type);
      endif

      if (~isfield(options, "Phi"))
        error("options.Phi is required for type %s", options.type);
      endif

      shape.z = interp1(options.Phi, options.z, shape.Phi, options.interpolation_z, "extrap");
      shape.z(1) = options.z(1);
      shape.z(end) = options.z(end);

      if (~isfield(options, "z_D"))
        options.z_D = options.z;
      endif

      shape.D = interp1(options.z_D, options.D, shape.z, options.interpolation_D, "extrap");
    otherwise
      error("type \"%s\" not implemented");
  endswitch

  X = [ shape.D/2 .* cos(shape.Phi);
        shape.D/2 .* sin(shape.Phi);
        shape.z ];

  if (max(max(isnan(X))))
    error("interpolation failed: interp1 returned NaN");
  endif

  if (isfield(options, "L"))
    shape.L = options.L;
  endif

  if (~isfield(shape, "L"))
    error("parameter L must be provided for custom typ %s", options.type);
  endif
endfunction
