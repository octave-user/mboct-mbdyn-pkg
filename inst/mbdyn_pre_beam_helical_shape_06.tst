## mbdyn_pre_beam_helical_shape.tst:06
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
%! [fd, fname] = mkstemp(fullfile(tempdir(), "oct-mbdyn_pre_beam_helical_shape_XXXXXX"), true);
%! mbdyn_pre_beam_print_nodes(beam, fd);
%! unwind_protect_cleanup
%! fclose(fd);
%! assert_simple(unlink(fname), 0);
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
%! daspect(ones(1, 3));
