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
## @deftypefn {Function File} [@var{beam}, @var{X}, @var{shape}] = mbdyn_pre_beam_helical_compute(@var{shape_options}, @var{number_of_elements_per_coil}, @var{interpolation_points})
## Build a beam model for a helical spring.
## @seealso{mbdyn_pre_beam_helical_shape}
## @end deftypefn

function [beam, X, shape] = mbdyn_pre_beam_helical_compute(shape_options, number_of_elements_per_coil, interpolation_points)
    [X, shape] = mbdyn_pre_beam_helical_shape(shape_options);
    number_of_coils = abs( shape.Phi(end) - shape.Phi(1) ) / ( 2 * pi );
    N = ceil(double(number_of_elements_per_coil) * number_of_coils);
    beam = mbdyn_pre_beam_compute(X,N,interpolation_points);
endfunction

%!test
%! close all;
%! options.d = 1.5e-3;
%! options.na = 1.2;
%! options.ni = 3;
%! options.L = 30e-3;
%! options.Di = 10.6e-3;
%! options.N_coil = 20;
%! options.type = "const pitch, const diameter";
%! number_of_elements_per_coil = 7;
%! interpolation_points = 1;
%! [beam, X, shape] = mbdyn_pre_beam_helical_compute(options, number_of_elements_per_coil, interpolation_points);
%! nfig = figure("visible", "off");
%! hold on;
%! mbdyn_pre_beam_plot(beam,struct("Rn",true,"Rg",false,"s",options.d,"figure",nfig));
%! plot3(X(1,:),X(2,:),X(3,:),'-;X;3');
%! figure("visible","off");
%! hold on;
%! plot(180 / pi * shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [deg]');
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
%! set(gca(),"DataAspectRatio",[1,1,1]);

%!demo
%! close all;
%! options.d = 1.5e-3;
%! options.na = 1.2;
%! options.ni = 3;
%! options.L = 30e-3;
%! options.Di = 10.6e-3;
%! options.N_coil = 20;
%! options.type = "const pitch, const diameter";
%! number_of_elements_per_coil = 7;
%! interpolation_points = 1;
%! [beam, X, shape] = mbdyn_pre_beam_helical_compute(options, number_of_elements_per_coil, interpolation_points);
%! figure("visible","off");
%! hold on;
%! plot(180 / pi * shape.Phi, shape.z,'-;z(Phi);1');
%! grid on;
%! grid minor on;
%! xlabel('Phi [deg]');
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
%! figure_list();
