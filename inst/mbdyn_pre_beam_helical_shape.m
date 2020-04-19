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
      options.N_coil = 200;
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

%!test
%! close all;
%! Di = [ 10.5e-3; 3e-3; 10.5e-3 ];
%! d = 1.5e-3;
%! n = 7.5;
%! options.L = 30e-3;
%! options.D = Di + d;
%! options.z = options.L * [ 0; 0.05; 0.95; 1 ];
%! options.z_D = options.L * [ 0; 0.5; 1 ];
%! options.Phi = 2 * pi * n * [ 0; 0.2; 0.8; 1 ];
%! options.N = 20;
%! options.interpolation = 'linear';
%! options.interpolation_D = 'spline';
%! [X, shape] = mbdyn_pre_beam_helical_shape(options);
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! plot(options.Phi, options.z, 'x;z_i(Phi_i);3');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! plot(options.z_D, options.D, 'x;D_i(z_i);3');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! set(gca(), "DataAspectRatio",[1,1,1]);

%!test
%! close all;
%! Di = 10.5e-3;
%! options.d = 1.5e-3;
%! options.na = 2;
%! options.ni = 2;
%! options.L = 30e-3;
%! options.D = Di + options.d;
%! options.N = 20;
%! options.type = "inactive coils, const pitch, const diameter";
%! number_of_elements_per_coil = 7;
%! N = ceil(number_of_elements_per_coil*options.na + 2 * options.ni);
%! interpolation_points = 1;
%! [X, shape] = mbdyn_pre_beam_helical_shape(options);
%! beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
%! mbdyn_pre_beam_plot(beam,struct("Rn",true, "Rg",false, "s",options.d));
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! set(gca(), "DataAspectRatio",[1,1,1]);

%!test
%! close all;
%! Di = 10.5e-3;
%! options.d = 1.5e-3;
%! options.na = 2;
%! options.ni = 3-270/360;
%! options.L = 30e-3;
%! options.D = Di + options.d;
%! options.N = 20;
%! options.N_begin = 10;
%! options.N_end = 10;
%! options.type = "const pitch, const diameter";
%! number_of_elements_per_coil = 7;
%! N = ceil(number_of_elements_per_coil*options.na);
%! interpolation_points = 1;
%! [X, shape] = mbdyn_pre_beam_helical_shape(options);
%! beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! nfig = figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! mbdyn_pre_beam_plot(beam,struct("Rn",true, "Rg",false, "s",options.d, "figure",nfig));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! set(gca(), "DataAspectRatio",[1,1,1]);

%!test
%! close all;
%! options.d = 1.5e-3;
%! options.na = 7.25;
%! options.ni = 0.25;
%! options.L = 30e-3;
%! options.Di = [12.45e-3, 10.45e-3, 12.45e-3 ];
%! options.type = "const pitch, 3 diameters";
%! number_of_elements_per_coil = 7;
%! N = ceil(number_of_elements_per_coil*options.na);
%! interpolation_points = 1;
%! [X, shape] = mbdyn_pre_beam_helical_shape(options);
%! beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
%! mbdyn_pre_beam_plot(beam,struct("Rn",true, "Rg",false, "s",options.d));
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! plot(shape.Phi, shape.z, 'x;z_i(Phi_i);3');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! set(gca(), "DataAspectRatio",[1,1,1]);

%!test
%! close all;
%! options.d = 1.6e-3;
%! options.na = 6;
%! options.ni = 3;
%! options.L = 34.7e-3;
%! options.Di = 12.4e-3;
%! options.z_rel = [ 0
%! 0.125106383
%! 0.301702128
%! 0.513191489
%! 0.726382979
%! 0.892340426
%! 1];
%! options.interpolation = "linear";
%! options.type = "relative variable pitch, const diameter";
%! number_of_elements_per_coil = 7;
%! N = ceil(number_of_elements_per_coil*options.na);
%! interpolation_points = 1;
%! [X, shape] = mbdyn_pre_beam_helical_shape(options);
%! beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
%! mbdyn_pre_beam_plot(beam,struct("Rn",true, "Rg",false, "s",options.d));
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! % plot(shape.Phi, shape.z, 'x;z_i(Phi_i);3');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! set(gca(), "DataAspectRatio",[1,1,1]);

%!test
%! close all;
%! param.shape.d = 1.6e-3;
%! param.shape.na = 5.5;
%! param.shape.ni = 3;
%! param.shape.L = 34.48e-3;
%! param.shape.Di = [13.45e-3, 12.45e-3, 13.45e-3];
%! param.shape.z_rel_D = [0, 0.5, 1];
%! param.shape.z_rel = [ 0, 0.38, 0.62, 1 ];
%! param.shape.interpolation_z = "linear";
%! param.shape.interpolation_D = "linear";
%! param.shape.type = "relative variable pitch, variable diameter";
%! param.shape.inactive_coil_alignment = false;
%! number_of_elements_per_coil = 7;
%! N = ceil(number_of_elements_per_coil*param.shape.na);
%! interpolation_points = 1;
%! [X, shape] = mbdyn_pre_beam_helical_shape(param.shape);
%! beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
%! mbdyn_pre_beam_plot(beam,struct("Rn",true, "Rg",false, "s",param.shape.d));
%! fd = -1;
%! unwind_protect
%! [fd, fname] = mkstemp(fullfile(tempdir(), "mbdyn_pre_beam_helical_shape_XXXXXX"), true);
%! mbdyn_pre_beam_print_nodes(beam, fd);
%! unwind_protect_cleanup
%! fclose(fd);
%! assert(unlink(fname), 0);
%! end_unwind_protect
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! set(gca(), "DataAspectRatio",[1,1,1]);

%!test
%! close all;
%! spring.d = 2.6e-3;
%! spring.L = 31.4e-3;
%! spring.Di = [ 13.55e-3; 13.55e-3; 13.9e-3; 13.9e-3; 13.55e-3; 13.55e-3];
%! spring.na = 4.6;
%! spring.ni = 2.7;
%! spring.z_D = [ -spring.d; 0; 6.5e-3; spring.L - 6.5e-3; spring.L; spring.L + spring.d ];
%! spring.N = int32(50);
%! spring.interpolation_D = 'linear';
%! spring.type = "inactive coils, const pitch, variable diameter";
%! spring.ground_surface = 340 * pi / 180;  # angle of ground area at each side
%! [X, shape] = mbdyn_pre_beam_helical_shape(spring);
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot3(X(1,:),X(2,:),X(3,:));
%! xlabel('x [m]');
%! ylabel('y [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! set(gca(), "DataAspectRatio",[1,1,1]);

%!demo
%! close all;
%! spring.d = 2.6e-3;
%! spring.L = 31.4e-3;
%! spring.Di = [ 13.55e-3; 13.55e-3; 13.9e-3; 13.9e-3; 13.55e-3; 13.55e-3];
%! spring.na = 4.6;
%! spring.ni = 2.7;
%! spring.z_D = [ -spring.d; 0; 6.5e-3; spring.L - 6.5e-3; spring.L; spring.L + spring.d ];
%! spring.N = int32(50);
%! spring.interpolation_D = 'linear';
%! spring.type = "inactive coils, const pitch, variable diameter";
%! spring.ground_surface = 340 * pi / 180;  # angle of ground area at each side
%! [X, shape] = mbdyn_pre_beam_helical_shape(spring);
%! figure("visible","off");
%! hold on;
%! plot(shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [rad]');
%! ylabel('z [m]');
%! title('pitch versus angle');
%! figure("visible","off");
%! hold on;
%! plot(shape.z, shape.D, '-;D(z);1');
%! grid on;
%! grid minor on;
%! xlabel('z [m]');
%! ylabel('D [m]');
%! title('coil diameter versus height');
%! figure("visible","off");
%! plot(X(1,:), X(3,:));
%! xlabel('x [m]');
%! zlabel('z [m]');
%! grid on;
%! grid minor on;
%! title('helical beam shape');
%! set(gca(), "DataAspectRatio", [1, 1, 1]);
