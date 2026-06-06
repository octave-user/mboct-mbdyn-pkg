## mbdyn_post_ehd_load_output.tst:10
%!test
%! try
%!   ####################################################################################################################################
%!   ## EHD TEST CASE according D.E. Sander, H. Allmaier, H.H. Priebsch, M. Witt and A. Skiadas
%!   ## Simulation of journal bearing friction in severe mixed lubrication - Validation and effect of surface smoothing due to running-in.
%!   ####################################################################################################################################
%!   close all;
%!   pkg load mboct-fem-pkg;
%!   output_file = "";
%!   output_file = tempname();
%!   if (ispc())
%!     output_file(output_file == "\\") = "/";
%!   endif
%!                                 %unwind_protect
%!   SI_unit_meter = 1e-3;
%!   SI_unit_kilogram = 1e-3;
%!   SI_unit_second = 1e-3;
%!   SI_unit_kelvin = 1;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   param.enable_Patir_Cheng = true;
%!   param.E = 210000e6 / SI_unit_pascal; ## Young's modulus of the shaft
%!   param.nu = 0.3; ## Poisson's ratio of the shaft
%!   param.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3); ## density of the shaft
%!   param.E2 = 70000e6 / SI_unit_pascal; ## Young's modulus of the conrod
%!   param.nu2 = 0.3; ## Poisson's ratio of the conrod
%!   param.rho2 = 2700 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.E3 = 210000e6 / SI_unit_pascal; ## Young's modulus of the support bearings
%!   param.nu3 = 0.3; ## Poisson's ratio of the support bearings
%!   param.rho3 = 7850 / (SI_unit_kilogram / SI_unit_meter^3); ## density of the support bearings
%!   param.d = 47.8e-3 / SI_unit_meter; ## journal diameter of the conrod big end bearing (test bearing)
%!   param.w = 17.2e-3 / SI_unit_meter; ## width of the conrod big end bearing (test bearing)
%!   param.ds = 54e-3 / SI_unit_meter; ## journal diameter of the support bearings
%!   param.wg = 5e-3 / SI_unit_meter; ## width of 180deg oil supply groove of the support bearings (assumption)
%!   param.dg = 5e-3 / SI_unit_meter; ## diameter of oil supply groove for the big end bearing (assumption)
%!   param.Psi = 1e-3; ## relative clearance of the conrod big end bearing (assumption)
%!   param.Psis = 1e-3; ## relative clearance of the support bearings (assumption)
%!   param.ws = 25e-3 / SI_unit_meter; ## width of the support bearings
%!   param.h1 = 2e-3 / SI_unit_meter; ## axial gap between conrod big end bearing and support bearing (assumption)
%!   param.sigma1 = 0.13e-6 / SI_unit_meter;
%!   param.sigma2 = 0.28e-6 / SI_unit_meter;
%!   param.delta1 = 0.21e-6 / SI_unit_meter;
%!   param.delta2 = 0.39e-6 / SI_unit_meter;
%!   param.Rq1 = 0.20e-6 / SI_unit_meter;
%!   param.Rq2 = 0.34e-6 / SI_unit_meter;
%!   param.gamma1 = 4;
%!   param.gamma2 = 12;
%!   param.sigma = sqrt(param.sigma1^2 + param.sigma2^2);
%!   param.delta = param.delta1 + param.delta2;
%!   param.K = 0.001;
%!   param.beta = param.sigma / (param.K / (16 * sqrt(2) * pi / 15 * 0.05^2))^2;
%!   param.eta = 0.05 / (param.sigma * param.beta);
%!   assert(param.K, 16 * sqrt(2) * pi / 15 * (param.eta * param.beta * param.sigma)^2 * sqrt(param.sigma / param.beta));
%!   param.Ered = 53.3e9 / SI_unit_pascal;
%!   param.mu = 0.02;
%!   param.sigma0 = 1e-6 / SI_unit_meter;
%!   param.l1 = param.w / 2 + param.ws / 2 + param.h1;
%!   param.l2 = param.w / 2 + param.ws + param.h1;
%!   param.l3 = 182e-3 / SI_unit_meter;
%!   param.dl = 50e-3 / SI_unit_meter;
%!   param.l4 = 50e-3 / SI_unit_meter;
%!   param.wr = 90e-3 / SI_unit_meter;
%!   param.d2 = 140e-3 / SI_unit_meter;
%!   param.l5 = 162e-3 / SI_unit_meter;
%!   param.l6 = 57e-3 / SI_unit_meter;
%!   param.l7 = 54e-3 / SI_unit_meter;
%!   param.l8 = 150e-3 / SI_unit_meter;
%!   param.l9 = 246e-3 / SI_unit_meter;
%!   param.l10 = 57e-3 / SI_unit_meter;
%!   param.l11 = 10e-3 / SI_unit_meter;
%!   param.l12 = 50e-3 / SI_unit_meter;
%!   param.l13 = 70e-3 / SI_unit_meter;
%!   param.l14 = 132e-3 / SI_unit_meter;
%!   param.l15 = 20e-3 / SI_unit_meter;
%!   param.l16 = 20e-3 / SI_unit_meter;
%!   param.pref = 100e6 / SI_unit_pascal; ## reference pressure
%!   param.rhol = 900 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.etal = 6.8e-3 / (SI_unit_pascal * SI_unit_second);
%!   param.betal = 2.4e9 / SI_unit_pascal;
%!   param.Tl = 100 / SI_unit_kelvin;
%!   param.pcl = 100 / SI_unit_pascal;
%!   param.fact_etav = 1e-3;
%!   param.pside = 1e5 / SI_unit_pascal;
%!   param.pin = 1e5 / SI_unit_pascal;
%!   param.hm = param.d * pi / 10; ## general mesh size
%!   param.hb = param.d * pi / 40; ## mesh size at the bearing surface
%!   param.number_of_nodes_x = 200;
%!   param.number_of_nodes_z = 25;
%!   param.num_modes_cms = int32(10); ## number of dynamic Craig Bampton modes
%!   param.num_modes_bearing = int32(70); ## number of bearing modes
%!   param.omega = [150,200,250,300,400,500,1000,3000] * pi / 30 / SI_unit_second^-1;
%!   param.F1 = [4e3, 8e3] / SI_unit_newton;
%!   empty_cell = cell(1, 3);
%!   ## steel shaft
%!   group_defs = struct("id", empty_cell, ...
%!                       "name", empty_cell, ...
%!                       "R", empty_cell, ...
%!                       "X0", empty_cell, ...
%!                       "Xi", empty_cell, ...
%!                       "type", empty_cell, ...
%!                       "geometry", empty_cell, ...
%!                       "compliance_matrix", empty_cell);
%!   group_defs(1).id = 1;
%!   group_defs(1).name = "node_id_drive";
%!   group_defs(1).R = [0, 0,-1;
%!                      0, 1, 0;
%!                      1, 0, 0];
%!   group_defs(1).X0 = [-param.l2; 0; 0];
%!   group_defs(1).type = "cylinder";
%!   group_defs(1).geometry.rmin = 0;
%!   group_defs(1).geometry.rmax = 0.5 * param.ds;
%!   group_defs(1).geometry.zmin = 0;
%!   group_defs(1).geometry.zmax = 0;
%!   group_defs(1).compliance_matrix.matrix_type = "none";
%!   group_defs(2).id = 2;
%!   group_defs(2).name = "node_id_support_bearing_journal_1";
%!   group_defs(2).R = [0, 0, -1;
%!                      0, 1,  0;
%!                      1, 0,  0];
%!   group_defs(2).X0 = [-param.l1; 0; 0];
%!   group_defs(2).type = "cylinder";
%!   group_defs(2).geometry.rmin = 0.5 * param.ds;
%!   group_defs(2).geometry.rmax = 0.5 * param.ds;
%!   group_defs(2).geometry.zmin = -0.5 * param.ws;
%!   group_defs(2).geometry.zmax = 0.5 * param.ws;
%!   group_defs(2).compliance_matrix.matrix_type = "modal substruct total";
%!   group_defs(2).compliance_matrix.bearing_type = "journal";
%!   group_defs(2).compliance_matrix.bearing_model = "EHD/FD";
%!   group_defs(2).compliance_matrix.reference_pressure = param.pref;
%!   group_defs(2).compliance_matrix.number_of_nodes_x = 50;
%!   group_defs(2).compliance_matrix.number_of_nodes_z = 5;
%!   group_defs(2).compliance_matrix.number_of_modes = param.num_modes_bearing;
%!   group_defs(2).compliance_matrix.include_rigid_body_modes = true;
%!   group_defs(2).bearing = "elem_id_support_bearing_1";
%!   group_defs(3).id = 3;
%!   group_defs(3).name = "node_id_big_end_bearing_journal";
%!   group_defs(3).R = [0, 0, -1;
%!                      0, 1,  0;
%!                      1, 0,  0];
%!   group_defs(3).X0 = [0; 0; 0];
%!   group_defs(3).type = "cylinder";
%!   group_defs(3).geometry.rmin = 0.5 * param.d;
%!   group_defs(3).geometry.rmax = 0.5 * param.d;
%!   group_defs(3).geometry.zmin = -0.5 * param.w;
%!   group_defs(3).geometry.zmax = 0.5 * param.w;
%!   group_defs(3).compliance_matrix.matrix_type = "modal substruct total";
%!   group_defs(3).compliance_matrix.bearing_type = "journal";
%!   group_defs(3).compliance_matrix.bearing_model = "EHD/FD";
%!   group_defs(3).compliance_matrix.reference_pressure = param.pref;
%!   group_defs(3).compliance_matrix.number_of_nodes_x = param.number_of_nodes_x;
%!   group_defs(3).compliance_matrix.number_of_nodes_z = param.number_of_nodes_z;
%!   group_defs(3).compliance_matrix.number_of_modes = param.num_modes_bearing;
%!   group_defs(3).compliance_matrix.include_rigid_body_modes = true;
%!   group_defs(3).bearing = "elem_id_big_end_bearing";
%!   group_defs(4).id = 4;
%!   group_defs(4).name = "node_id_support_bearing_journal_2";
%!   group_defs(4).R = [0, 0, -1;
%!                      0, 1,  0;
%!                      1, 0,  0];
%!   group_defs(4).X0 = [param.l1; 0; 0];
%!   group_defs(4).type = "cylinder";
%!   group_defs(4).geometry.rmin = 0.5 * param.ds;
%!   group_defs(4).geometry.rmax = 0.5 * param.ds;
%!   group_defs(4).geometry.zmin = -0.5 * param.ws;
%!   group_defs(4).geometry.zmax = 0.5 * param.ws;
%!   group_defs(4).compliance_matrix.matrix_type = "modal substruct total";
%!   group_defs(4).compliance_matrix.bearing_type = "journal";
%!   group_defs(4).compliance_matrix.bearing_model = "EHD/FD";
%!   group_defs(4).compliance_matrix.reference_pressure = param.pref;
%!   group_defs(4).compliance_matrix.number_of_nodes_x = param.number_of_nodes_x;
%!   group_defs(4).compliance_matrix.number_of_nodes_z = param.number_of_nodes_z;
%!   group_defs(4).compliance_matrix.number_of_modes = param.num_modes_bearing;
%!   group_defs(4).compliance_matrix.include_rigid_body_modes = true;
%!   group_defs(4).bearing = "elem_id_support_bearing_2";
%!   empty_cell = cell(1, 2);
%!   ## conrod
%!   group_defs2 = struct("id", empty_cell, ...
%!                        "name", empty_cell, ...
%!                        "R", empty_cell, ...
%!                        "X0", empty_cell, ...
%!                        "Xi", empty_cell, ...
%!                        "type", empty_cell, ...
%!                        "geometry", empty_cell, ...
%!                        "compliance_matrix", empty_cell);
%!   group_defs2(1).id = 1;
%!   group_defs2(1).name = "node_id_small_end_bearing_shell";
%!   group_defs2(1).R = [0, 0,-1;
%!                       0, 1, 0;
%!                       1, 0, 0];
%!   group_defs2(1).X0 = [0; 0; param.l3];
%!   group_defs2(1).type = "cylinder";
%!   group_defs2(1).geometry.rmin = 0.5 * param.dl;
%!   group_defs2(1).geometry.rmax = 0.5 * param.dl;
%!   group_defs2(1).geometry.zmin = -0.5 * param.w;
%!   group_defs2(1).geometry.zmax = 0.5 * param.w;
%!   group_defs2(1).compliance_matrix.matrix_type = "none";
%!   group_defs2(2).id = 2;
%!   group_defs2(2).name = "node_id_big_end_bearing_shell";
%!   group_defs2(2).R = [0, 0, -1;
%!                       0, 1,  0;
%!                       1, 0,  0];
%!   group_defs2(2).X0 = [0; 0; 0];
%!   group_defs2(2).type = "cylinder";
%!   group_defs2(2).geometry.rmin = 0.5 * param.d;
%!   group_defs2(2).geometry.rmax = 0.5 * param.d;
%!   group_defs2(2).geometry.zmin = -0.5 * param.w;
%!   group_defs2(2).geometry.zmax = 0.5 * param.w;
%!   group_defs2(2).compliance_matrix.matrix_type = "modal substruct total";
%!   group_defs2(2).compliance_matrix.bearing_type = "shell";
%!   group_defs2(2).compliance_matrix.bearing_model = "EHD/FD";
%!   group_defs2(2).compliance_matrix.reference_pressure = param.pref;
%!   group_defs2(2).compliance_matrix.number_of_nodes_x = param.number_of_nodes_x;
%!   group_defs2(2).compliance_matrix.number_of_nodes_z = param.number_of_nodes_z;
%!   group_defs2(2).compliance_matrix.number_of_modes = param.num_modes_bearing;
%!   group_defs2(2).compliance_matrix.include_rigid_body_modes = true;
%!   group_defs2(2).bearing = "elem_id_big_end_bearing";
%!   empty_cell = cell(1, 2);
%!   ## support bracket
%!   group_defs3 = struct("id", empty_cell, ...
%!                        "name", empty_cell, ...
%!                        "R", empty_cell, ...
%!                        "X0", empty_cell, ...
%!                        "Xi", empty_cell, ...
%!                        "type", empty_cell, ...
%!                        "geometry", empty_cell, ...
%!                        "compliance_matrix", empty_cell);
%!   group_defs3(1).id = 1;
%!   group_defs3(1).name = "node_id_supporting_area";
%!   group_defs3(1).R = eye(3);
%!   group_defs3(1).X0 = [0; 0; -param.l14];
%!   group_defs3(1).type = "box";
%!   group_defs3(1).geometry.xmin = -0.5 * param.ws;
%!   group_defs3(1).geometry.xmax = 0.5 * param.ws;
%!   group_defs3(1).geometry.ymin = -0.5 * param.l9;
%!   group_defs3(1).geometry.ymax = 0.5 * param.l9;
%!   group_defs3(1).geometry.zmin = 0;
%!   group_defs3(1).geometry.zmax = 0;
%!   group_defs3(1).compliance_matrix.matrix_type = "none";
%!   group_defs3(2).id = 2;
%!   group_defs3(2).name = "node_id_support_bearing_shell";
%!   group_defs3(2).R = [0, 0, -1;
%!                       0, 1,  0;
%!                       1, 0,  0];
%!   group_defs3(2).X0 = [0; 0; 0];
%!   group_defs3(2).type = "cylinder";
%!   group_defs3(2).geometry.rmin = 0.5 * param.ds;
%!   group_defs3(2).geometry.rmax = 0.5 * param.ds;
%!   group_defs3(2).geometry.zmin = -0.5 * param.ws;
%!   group_defs3(2).geometry.zmax = 0.5 * param.ws;
%!   group_defs3(2).compliance_matrix.matrix_type = "modal substruct total";
%!   group_defs3(2).compliance_matrix.bearing_type = "shell";
%!   group_defs3(2).compliance_matrix.bearing_model = "EHD/FD";
%!   group_defs3(2).compliance_matrix.reference_pressure = param.pref;
%!   group_defs3(2).compliance_matrix.number_of_nodes_x = param.number_of_nodes_x;
%!   group_defs3(2).compliance_matrix.number_of_nodes_z = param.number_of_nodes_z;
%!   group_defs3(2).compliance_matrix.number_of_modes = param.num_modes_bearing;
%!   group_defs3(2).compliance_matrix.include_rigid_body_modes = true;
%!   group_defs3(2).bearing = "elem_id_support_bearing";

%!   fd = -1;
%!   unwind_protect
%!     fd = fopen([output_file, "_support.geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", [output_file, "_support.geo"]);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Geometry.OCCUnionUnify = 0;\n");
%!     fprintf(fd, "ds = %g;\n", param.ds);
%!     fprintf(fd, "ws = %g;\n", param.ws);
%!     fprintf(fd, "l7 = %g;\n", param.l7);
%!     fprintf(fd, "l8 = %g;\n", param.l8);
%!     fprintf(fd, "l9 = %g;\n", param.l9);
%!     fprintf(fd, "l10 = %g;\n", param.l10);
%!     fprintf(fd, "l11 = %g;\n", param.l11);
%!     fprintf(fd, "l12 = %g;\n", param.l12);
%!     fprintf(fd, "l13 = %g;\n", param.l13);
%!     fprintf(fd, "l14 = %g;\n", param.l14);
%!     fprintf(fd, "l15 = %g;\n", param.l15);
%!     fprintf(fd, "l16 = %g;\n", param.l16);
%!     fprintf(fd, "hm = %g;\n", param.hm);
%!     fprintf(fd, "hb = %g;\n", param.hb);
%!     fputs(fd, "Point(1) = {-0.5 * ws, -0.5 * l12, l7 + l11};\n");
%!     fputs(fd, "Point(2) = {-0.5 * ws, 0.5 * l12, l7 + l11};\n");
%!     fputs(fd, "Point(3) = {-0.5 * ws, 0.5 * l13, l7};\n");
%!     fputs(fd, "Point(4) = {-0.5 * ws, 0.5 * l8, l7};\n");
%!     fputs(fd, "Point(5) = {-0.5 * ws, 0.5 * l8, -l14 + l10};\n");
%!     fputs(fd, "Point(6) = {-0.5 * ws, 0.5 * l9, -l14 + l10};\n");
%!     fputs(fd, "Point(7) = {-0.5 * ws, 0.5 * l9, -l14};\n");
%!     fputs(fd, "Point(8) = {-0.5 * ws, 0.5 * l16, -l14};\n");
%!     fputs(fd, "Point(9) = {-0.5 * ws, 0.5 * l16, -l14 + l15};\n");
%!     fputs(fd, "Point(10) = {-0.5 * ws, -0.5 * l16, -l14 + l15};\n");
%!     fputs(fd, "Point(11) = {-0.5 * ws, -0.5 * l16, -l14};\n");
%!     fputs(fd, "Point(12) = {-0.5 * ws, -0.5 * l9, -l14};\n");
%!     fputs(fd, "Point(13) = {-0.5 * ws, -0.5 * l9, -l14 + l10};\n");
%!     fputs(fd, "Point(14) = {-0.5 * ws, -0.5 * l8, -l14 + l10};\n");
%!     fputs(fd, "Point(15) = {-0.5 * ws, -0.5 * l8, l7};\n");
%!     fputs(fd, "Point(16) = {-0.5 * ws, -0.5 * l13, l7};\n");
%!     fputs(fd, "Point(17) = {-0.5 * ws, 0, 0};\n");
%!     fputs(fd, "Point(18) = {-0.5 * ws, 0.5 * ds, 0};\n");
%!     fputs(fd, "Point(19) = {-0.5 * ws, 0, -0.5 * ds};\n");
%!     fputs(fd, "Point(20) = {-0.5 * ws, -0.5 * ds, 0};\n");
%!     fputs(fd, "Point(21) = {-0.5 * ws, 0, 0.5 * ds};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 9};\n");
%!     fputs(fd, "Line(9) = {9, 10};\n");
%!     fputs(fd, "Line(10) = {10, 11};\n");
%!     fputs(fd, "Line(11) = {11, 12};\n");
%!     fputs(fd, "Line(12) = {12, 13};\n");
%!     fputs(fd, "Line(13) = {13, 14};\n");
%!     fputs(fd, "Line(14) = {14, 15};\n");
%!     fputs(fd, "Line(15) = {15, 16};\n");
%!     fputs(fd, "Line(16) = {16, 1};\n");
%!     fputs(fd, "Circle(17) = {18, 17, 19};\n");
%!     fputs(fd, "Circle(18) = {19, 17, 20};\n");
%!     fputs(fd, "Circle(19) = {20, 17, 21};\n");
%!     fputs(fd, "Circle(20) = {21, 17, 18};\n");
%!     fputs(fd, "Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};\n");
%!     fputs(fd, "Curve Loop(2) = {17,18,19,20};\n");
%!     fputs(fd, "Plane Surface(1) = {1, 2};\n");
%!     fputs(fd, "v = Extrude{ws, 0, 0}{Surface{1};};\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{v[1]}; } } = hm;\n");
%!     fputs(fd, "MeshSize{PointsOf{Surface{19, 18, 21, 20}; } } = hb;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!     fputs(fd, "Physical Surface(\"support\", 1) = {12, 8};\n");
%!     fputs(fd, "Physical Surface(\"support_bearing\", 2) = {19, 18, 21, 20};\n");
%!     fputs(fd, "Physical Volume(\"support\", 1) = {v[1]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   ##spawn_wait(spawn("gmsh", {[output_file, "_support.geo"]}));
%!   pid = spawn("gmsh", {"-format", "msh2", "-order", "2", "-3", [output_file, "_support.geo"]});
%!   status = spawn_wait(pid);
%!   mesh3 = fem_pre_mesh_import([output_file, "_support.msh"], "gmsh");
%!   mesh3 = fem_pre_mesh_reorder(mesh3);
%!   mesh3.groups.tria6 = fem_pre_mesh_groups_create(mesh3, group_defs3, 1e-3 * param.d).tria6;
%!   mesh3.materials.tet10 = zeros(rows(mesh3.elements.tet10), 1, "int32");
%!   mesh3.materials.tet10(mesh3.groups.tet10(find([mesh3.groups.tet10.id == 1])).elements) = 1;
%!   mesh3.material_data(1).E = param.E3;
%!   mesh3.material_data(1).nu = param.nu3;
%!   mesh3.material_data(1).rho = param.rho3;
%!   cms_opt3.invariants = true;
%!   cms_opt3.refine_max_iter = int32(10);
%!   cms_opt3.number_of_threads = mbdyn_solver_num_threads_default();
%!   cms_opt3.verbose = false;
%!   cms_opt3.modes.number = param.num_modes_cms;
%!   cms_opt3.element.name = "elem_id_support";
%!   node_set3 = int32(rows(mesh3.nodes) + (1:numel(group_defs3)));
%!   node_names3 = {group_defs3.name};
%!   cms_opt3.nodes.modal.number = node_set3(1);
%!   cms_opt3.nodes.modal.name = node_names3{1};
%!   for j=1:numel(node_set3) - 1
%!     cms_opt3.nodes.interfaces(j).number = node_set3(j + 1);
%!     cms_opt3.nodes.interfaces(j).name = node_names3{j + 1};
%!   endfor
%!   idx_grp_itf3 = find([group_defs3.id] > 1);
%!   idx_grp_modal3 = find([group_defs3.id] == 1);
%!   mesh3.nodes([cms_opt3.nodes.interfaces.number], 1:3) = [group_defs3(idx_grp_itf3).X0].';
%!   mesh3.nodes(cms_opt3.nodes.modal.number, 1:3) = group_defs3(idx_grp_modal3).X0.';
%!   mesh3.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh3, [group_defs3(idx_grp_itf3).id], [node_set3(idx_grp_itf3)]);
%!   mesh3.elements.rbe2 = fem_pre_mesh_rbe2_from_surf(mesh3, [group_defs3(idx_grp_modal3).id], [node_set3(idx_grp_modal3)]);
%!   load_case3.locked_dof = false(rows(mesh3.nodes), columns(mesh3.nodes));
%!   bearing_surf3 = repmat(struct("group_idx", [], "X0", [], "R", [], "options", [], "name", [], "bearing", []), 1, numel(group_defs3));
%!   num_comp_mat3 = int32(0);
%!   for j=1:numel(group_defs3)
%!     switch group_defs3(j).compliance_matrix.matrix_type
%!       case "none"
%!       otherwise
%!         ++num_comp_mat3;
%!         bearing_surf3(num_comp_mat3).group_idx = find([mesh3.groups.tria6.id] == group_defs3(j).id);
%!         bearing_surf3(num_comp_mat3).group_id_interface = group_defs3(j).id + int32(1000);
%!         bearing_surf3(num_comp_mat3).absolute_tolerance = 1e-3 * param.d;
%!         bearing_surf3(num_comp_mat3).relative_tolerance = 0;
%!         bearing_surf3(num_comp_mat3).name = group_defs3(j).name;
%!         bearing_surf3(num_comp_mat3).bearing = group_defs3(j).bearing;
%!         bearing_surf3(num_comp_mat3).X0 = group_defs3(j).X0;
%!         bearing_surf3(num_comp_mat3).R = group_defs3(j).R;
%!         bearing_surf3(num_comp_mat3).options = group_defs3(j).compliance_matrix;
%!         bearing_surf3(num_comp_mat3).master_node_no = node_set3(j);
%!         switch group_defs3(j).type
%!           case "cylinder"
%!             bearing_surf3(num_comp_mat3).r = mean([group_defs3(j).geometry.rmax, group_defs3(j).geometry.rmin]);
%!             bearing_surf3(num_comp_mat3).w = group_defs3(j).geometry.zmax - group_defs3(j).geometry.zmin;
%!           otherwise
%!             error("bearing geometry type \"%s\" not implemented", group_defs3(j).type);
%!         endswitch
%!         bearing_surf3(num_comp_mat3).nodes = mesh3.groups.tria6(find([[mesh3.groups.tria6].id] == group_defs3(j).id)).nodes;
%!     endswitch
%!   endfor
%!   bearing_surf3 = bearing_surf3(1:num_comp_mat3);
%!   [mesh3, mat_ass3, dof_map3, cms_opt3, comp_mat3, bearing_surf3, sol_eig3] = fem_ehd_pre_comp_mat_linear_mesh(mesh3, load_case3, cms_opt3, bearing_surf3);
%!   for j=1:2
%!     cms_opt3_tmp = cms_opt3;
%!     cms_opt3_tmp.element.name = sprintf("%s_%d", cms_opt3.element.name, j);
%!     cms_opt3_tmp.nodes.modal.name = sprintf("%s_%d", cms_opt3.nodes.modal.name, j);
%!     cms_opt3_tmp.nodes.interfaces.name = sprintf("%s_%d", cms_opt3.nodes.interfaces.name, j);
%!     fem_cms_export([output_file, sprintf("_support_cms_%d", j)], mesh3, dof_map3, mat_ass3, cms_opt3_tmp);
%!   endfor
%!   for j=1:numel(comp_mat3)
%!     comp_mat_file = [output_file, "_", bearing_surf3(j).bearing, "_", bearing_surf3(j).options.bearing_type, ".dat"];
%!     fem_ehd_pre_comp_mat_export(comp_mat3(j), bearing_surf3(j).options, comp_mat_file);
%!   endfor
%!   fd = -1;
%!   unwind_protect
%!     fd = fopen([output_file, "_conrod.geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", [output_file, "_conrod.geo"]);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Geometry.OCCUnionUnify = 0;\n");
%!     fprintf(fd, "d = %g;\n", param.d);
%!     fprintf(fd, "d2 = %g;\n", param.d2);
%!     fprintf(fd, "w = %g;\n", param.w);
%!     fprintf(fd, "dl = %g;\n", param.dl);
%!     fprintf(fd, "wr = %g;\n", param.wr);
%!     fprintf(fd, "l3 = %g;\n", param.l3);
%!     fprintf(fd, "l4 = %g;\n", param.l4);
%!     fprintf(fd, "l5 = %g;\n", param.l5);
%!     fprintf(fd, "l6 = %g;\n", param.l6);
%!     fprintf(fd, "hm = %g;\n", param.hm);
%!     fprintf(fd, "hb = %g;\n", param.hb);
%!     fputs(fd, "Point(1) = {-0.5 * w, -0.5 * wr, l3 + l4};\n");
%!     fputs(fd, "Point(2) = {-0.5 * w, 0.5 * wr, l3 + l4};\n");
%!     fputs(fd, "Point(3) = {-0.5 * w, 0.5 * wr, Sqrt((d2/2)^2-(wr/2)^2)};\n");
%!     fputs(fd, "Point(4) = {-0.5 * w, Sqrt((d2/2)^2 - (l6/2)^2), 0.5 * l6};\n");
%!     fputs(fd, "Point(5) = {-0.5 * w, 0.5 * l5, 0.5 * l6};\n");
%!     fputs(fd, "Point(6) = {-0.5 * w, 0.5 * l5, -0.5 * l6};\n");
%!     fputs(fd, "Point(7) = {-0.5 * w, Sqrt((d2/2)^2 - (l6/2)^2), -0.5 * l6};\n");
%!     fputs(fd, "Point(8) = {-0.5 * w, -Sqrt((d2/2)^2 - (l6/2)^2), -0.5 * l6};\n");
%!     fputs(fd, "Point(9) = {-0.5 * w, -0.5 * l5, -0.5 * l6};\n");
%!     fputs(fd, "Point(10) = {-0.5 * w, -0.5 * l5, 0.5 * l6};\n");
%!     fputs(fd, "Point(11) = {-0.5 * w, -Sqrt((d2/2)^2 - (l6/2)^2), 0.5 * l6};\n");
%!     fputs(fd, "Point(12) = {-0.5 * w, -0.5 * wr, Sqrt((d2/2)^2-(wr/2)^2)};\n");
%!     fputs(fd, "Point(13) = {-0.5 * w, 0, 0};\n");
%!     fputs(fd, "Point(14) = {-0.5 * w, 0.5 * d, 0};\n");
%!     fputs(fd, "Point(15) = {-0.5 * w, 0, -0.5 * d};\n");
%!     fputs(fd, "Point(16) = {-0.5 * w, -0.5 * d, 0};\n");
%!     fputs(fd, "Point(17) = {-0.5 * w, 0, 0.5 * d};\n");
%!     fputs(fd, "Point(18) = {-0.5 * w, 0, l3};\n");
%!     fputs(fd, "Point(19) = {-0.5 * w, 0.5 * dl, l3};\n");
%!     fputs(fd, "Point(20) = {-0.5 * w, 0, l3 - 0.5 * dl};\n");
%!     fputs(fd, "Point(21) = {-0.5 * w, -0.5 * dl, l3};\n");
%!     fputs(fd, "Point(22) = {-0.5 * w, 0, l3 + 0.5 * dl};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Circle(3) = {3, 13, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Circle(7) = {7, 13, 8};\n");
%!     fputs(fd, "Line(8) = {8, 9};\n");
%!     fputs(fd, "Line(9) = {9, 10};\n");
%!     fputs(fd, "Line(10) = {10, 11};\n");
%!     fputs(fd, "Circle(11) = {11, 13, 12};\n");
%!     fputs(fd, "Line(12) = {12, 1};\n");
%!     fputs(fd, "Circle(13) = {14, 13, 15};\n");
%!     fputs(fd, "Circle(14) = {15, 13, 16};\n");
%!     fputs(fd, "Circle(15) = {16, 13, 17};\n");
%!     fputs(fd, "Circle(16) = {17, 13, 14};\n");
%!     fputs(fd, "Circle(17) = {19, 18, 20};\n");
%!     fputs(fd, "Circle(18) = {20, 18, 21};\n");
%!     fputs(fd, "Circle(19) = {21, 18, 22};\n");
%!     fputs(fd, "Circle(20) = {22, 18, 19};\n");
%!     fputs(fd, "Curve Loop(1) = {1,2,3,4,5,6,7,8,9,10,11,12};\n");
%!     fputs(fd, "Curve Loop(2) = {13, 14, 15, 16};\n");
%!     fputs(fd, "Curve Loop(3) = {17, 18, 19, 20};\n");
%!     fputs(fd, "Plane Surface(1) = {1, 2, 3};\n");
%!     fputs(fd, "v = Extrude{w, 0, 0}{Surface{1};};\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{v[1]}; } } = hm;\n");
%!     fputs(fd, "MeshSize{PointsOf{Surface{14, 15, 16, 17}; } } = hb;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!     fputs(fd, "Physical Surface(\"big_end\", 1) = {14, 15, 16, 17};\n");
%!     fputs(fd, "Physical Surface(\"small_end\", 2) = {18, 19, 20, 21};\n");
%!     fputs(fd, "Physical Volume(\"conrod\", 1) = {v[1]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-order", "2", "-3", [output_file, "_conrod.geo"]});
%!   status = spawn_wait(pid);
%!   mesh2 = fem_pre_mesh_import([output_file, "_conrod.msh"], "gmsh");
%!   mesh2 = fem_pre_mesh_reorder(mesh2);
%!   mesh2.groups.tria6 = fem_pre_mesh_groups_create(mesh2, group_defs2, 1e-3 * param.d).tria6;
%!   mesh2.materials.tet10 = zeros(rows(mesh2.elements.tet10), 1, "int32");
%!   mesh2.materials.tet10(mesh2.groups.tet10(find([mesh2.groups.tet10.id == 1])).elements) = 1;
%!   mesh2.material_data(1).E = param.E2;
%!   mesh2.material_data(1).nu = param.nu2;
%!   mesh2.material_data(1).rho = param.rho2;
%!   cms_opt2.invariants = true;
%!   cms_opt2.refine_max_iter = int32(10);
%!   cms_opt2.number_of_threads = mbdyn_solver_num_threads_default();
%!   cms_opt2.verbose = false;
%!   cms_opt2.modes.number = param.num_modes_cms;
%!   cms_opt2.element.name = "elem_id_conrod";
%!   node_set2 = int32(rows(mesh2.nodes) + (1:numel(group_defs2)));
%!   node_names2 = {group_defs2.name};
%!   cms_opt2.nodes.modal.number = node_set2(1);
%!   cms_opt2.nodes.modal.name = node_names2{1};
%!   for j=1:numel(node_set2) - 1
%!     cms_opt2.nodes.interfaces(j).number = node_set2(j + 1);
%!     cms_opt2.nodes.interfaces(j).name = node_names2{j + 1};
%!   endfor
%!   idx_grp_itf2 = find([group_defs2.id] > 1);
%!   idx_grp_modal2 = find([group_defs2.id] == 1);
%!   mesh2.nodes([cms_opt2.nodes.interfaces.number], 1:3) = [group_defs2(idx_grp_itf2).X0].';
%!   mesh2.nodes(cms_opt2.nodes.modal.number, 1:3) = group_defs2(idx_grp_modal2).X0.';
%!   mesh2.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh2, [group_defs2(idx_grp_modal2).id, group_defs2(idx_grp_itf2).id], [node_set2(idx_grp_modal2), node_set2(idx_grp_itf2)]);
%!   load_case2.locked_dof = false(rows(mesh2.nodes), columns(mesh2.nodes));
%!   bearing_surf2 = repmat(struct("group_idx", [], "X0", [], "R", [], "options", [], "name", [], "bearing", []), 1, numel(group_defs2));
%!   num_comp_mat2 = int32(0);
%!   for j=1:numel(group_defs2)
%!     switch group_defs2(j).compliance_matrix.matrix_type
%!       case "none"
%!       otherwise
%!         ++num_comp_mat2;
%!         bearing_surf2(num_comp_mat2).group_idx = find([mesh2.groups.tria6.id] == group_defs2(j).id);
%!         bearing_surf2(num_comp_mat2).group_id_interface = group_defs2(j).id + int32(1000);
%!         bearing_surf2(num_comp_mat2).absolute_tolerance = 1e-3 * param.d;
%!         bearing_surf2(num_comp_mat2).relative_tolerance = 0;
%!         bearing_surf2(num_comp_mat2).name = group_defs2(j).name;
%!         bearing_surf2(num_comp_mat2).bearing = group_defs2(j).bearing;
%!         bearing_surf2(num_comp_mat2).X0 = group_defs2(j).X0;
%!         bearing_surf2(num_comp_mat2).R = group_defs2(j).R;
%!         bearing_surf2(num_comp_mat2).options = group_defs2(j).compliance_matrix;
%!         bearing_surf2(num_comp_mat2).master_node_no = node_set2(j);
%!         switch group_defs2(j).type
%!           case "cylinder"
%!             bearing_surf2(num_comp_mat2).r = mean([group_defs2(j).geometry.rmax, group_defs2(j).geometry.rmin]);
%!             bearing_surf2(num_comp_mat2).w = group_defs2(j).geometry.zmax - group_defs2(j).geometry.zmin;
%!           otherwise
%!             error("bearing geometry type \"%s\" not implemented", group_defs2(j).type);
%!         endswitch
%!         bearing_surf2(num_comp_mat2).nodes = mesh2.groups.tria6(find([[mesh2.groups.tria6].id] == group_defs2(j).id)).nodes;
%!     endswitch
%!   endfor
%!   bearing_surf2 = bearing_surf2(1:num_comp_mat2);
%!   [mesh2, mat_ass2, dof_map2, cms_opt2, comp_mat2, bearing_surf2, sol_eig2] = fem_ehd_pre_comp_mat_linear_mesh(mesh2, load_case2, cms_opt2, bearing_surf2);
%!   fem_cms_export([output_file, "_conrod_cms"], mesh2, dof_map2, mat_ass2, cms_opt2);
%!   for j=1:numel(comp_mat2)
%!     comp_mat_file = [output_file, "_", bearing_surf2(j).bearing, "_", bearing_surf2(j).options.bearing_type, ".dat"];
%!     fem_ehd_pre_comp_mat_export(comp_mat2(j), bearing_surf2(j).options, comp_mat_file);
%!   endfor
%!   fd = -1;
%!   unwind_protect
%!     fd = fopen([output_file, "_shaft.geo"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", [output_file, "_shaft.geo"]);
%!     endif
%!     fputs(fd, "SetFactory(\"OpenCASCADE\");\n");
%!     fputs(fd, "Geometry.OCCUnionUnify = 0;\n");
%!     fprintf(fd, "d = %g;\n", param.d);
%!     fprintf(fd, "w = %g;\n", param.w);
%!     fprintf(fd, "ds = %g;\n", param.ds);
%!     fprintf(fd, "ws = %g;\n", param.ws);
%!     fprintf(fd, "h = %g;\n", param.h1);
%!     fprintf(fd, "hm = %g;\n", param.hm);
%!     fprintf(fd, "hb = %g;\n", param.hb);
%!     fputs(fd, "Point(1) = {-ws - w / 2 - h, 0, 0};\n");
%!     fputs(fd, "Point(2) = {-ws - w / 2 - h, 0, ds / 2};\n");
%!     fputs(fd, "Point(3) = {-w / 2 - h, 0, ds / 2};\n");
%!     fputs(fd, "Point(4) = {-w / 2 - h, 0, d / 2};\n");
%!     fputs(fd, "Point(5) = {-w / 2, 0, d / 2};\n");
%!     fputs(fd, "Point(6) = {w / 2, 0, d / 2};\n");
%!     fputs(fd, "Point(7) = {w / 2 + h, 0, d / 2};\n");
%!     fputs(fd, "Point(8) = {w / 2 + h, 0, ds / 2};\n");
%!     fputs(fd, "Point(9) = {ws + w / 2 + h, 0, ds / 2};\n");
%!     fputs(fd, "Point(10) = {ws + w / 2 + h, 0, 0};\n");
%!     fputs(fd, "Line(1) = {1, 2};\n");
%!     fputs(fd, "Line(2) = {2, 3};\n");
%!     fputs(fd, "Line(3) = {3, 4};\n");
%!     fputs(fd, "Line(4) = {4, 5};\n");
%!     fputs(fd, "Line(5) = {5, 6};\n");
%!     fputs(fd, "Line(6) = {6, 7};\n");
%!     fputs(fd, "Line(7) = {7, 8};\n");
%!     fputs(fd, "Line(8) = {8, 9};\n");
%!     fputs(fd, "Line(9) = {9, 10};\n");
%!     fputs(fd, "Line(10) = {10, 1};\n");
%!     fputs(fd, "Line Loop(1) = {1,2,3,4,5,6,7,8,9,10};\n");
%!     fputs(fd, "Plane Surface(1) = {1};\n");
%!     fputs(fd, "v[] = Extrude {{1, 0, 0}, {0, 0, 0}, 2*Pi} { Surface{1}; };\n");
%!     fputs(fd, "MeshSize{PointsOf{Volume{v[1]}; } } = hm;\n");
%!     fputs(fd, "MeshSize{PointsOf{Surface{6,3,9,5,4,7,8}; } } = hb;\n");
%!     fputs(fd, "Mesh.HighOrderOptimize = 2;\n");
%!     fputs(fd, "Physical Surface(\"big_end_bearing\", 1) = {6};\n");
%!     fputs(fd, "Physical Surface(\"support_bearing1\", 2) = {3};\n");
%!     fputs(fd, "Physical Surface(\"support_bearing2\", 3) = {9};\n");
%!     fputs(fd, "Physical Surface(\"torque_input\", 4) = {2};\n");
%!     #fputs(fd, "Physical Surface(\"transition\", 5) = {5, 4, 7, 8};\n");
%!     fputs(fd, "Physical Volume(\"shaft\",1) = {v[1]};\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   pid = spawn("gmsh", {"-format", "msh2", "-order", "2", "-3", [output_file, "_shaft.geo"]});
%!   status = spawn_wait(pid);
%!   mesh = fem_pre_mesh_import([output_file, "_shaft.msh"], "gmsh");
%!   mesh = fem_pre_mesh_reorder(mesh);
%!   mesh.groups.tria6 = fem_pre_mesh_groups_create(mesh, group_defs, 1e-3 * param.d).tria6;
%!   mesh.materials.tet10 = zeros(rows(mesh.elements.tet10), 1, "int32");
%!   mesh.materials.tet10(mesh.groups.tet10(find([mesh.groups.tet10.id == 1])).elements) = 1;
%!   mesh.material_data(1).E = param.E;
%!   mesh.material_data(1).nu = param.nu;
%!   mesh.material_data(1).rho = param.rho;
%!   cms_opt.invariants = true;
%!   cms_opt.refine_max_iter = int32(10);
%!   cms_opt.number_of_threads = mbdyn_solver_num_threads_default();
%!   cms_opt.verbose = false;
%!   cms_opt.modes.number = param.num_modes_cms;
%!   cms_opt.element.name = "elem_id_shaft";
%!   node_set = int32(rows(mesh.nodes) + (1:numel(group_defs)));
%!   node_names = {group_defs.name};
%!   cms_opt.nodes.modal.number = node_set(1);
%!   cms_opt.nodes.modal.name = node_names{1};
%!   for j=1:numel(node_set) - 1
%!     cms_opt.nodes.interfaces(j).number = node_set(j + 1);
%!     cms_opt.nodes.interfaces(j).name = node_names{j + 1};
%!   endfor
%!   idx_grp_itf = find([group_defs.id] > 1);
%!   idx_grp_modal = find([group_defs.id] == 1);
%!   mesh.nodes([cms_opt.nodes.interfaces.number], 1:3) = [group_defs(idx_grp_itf).X0].';
%!   mesh.nodes(cms_opt.nodes.modal.number, 1:3) = group_defs(idx_grp_modal).X0.';
%!   mesh.elements.rbe3 = fem_pre_mesh_rbe3_from_surf(mesh, [group_defs(idx_grp_modal).id, group_defs(idx_grp_itf).id], [node_set(idx_grp_modal), node_set(idx_grp_itf)]);
%!   load_case.locked_dof = false(rows(mesh.nodes), columns(mesh.nodes));
%!   bearing_surf = repmat(struct("group_idx", [], "X0", [], "R", [], "options", [], "name", [], "bearing", []), 1, numel(group_defs));
%!   num_comp_mat = int32(0);
%!   for j=1:numel(group_defs)
%!     switch group_defs(j).compliance_matrix.matrix_type
%!       case "none"
%!       otherwise
%!         ++num_comp_mat;
%!         bearing_surf(num_comp_mat).group_idx = find([mesh.groups.tria6.id] == group_defs(j).id);
%!         bearing_surf(num_comp_mat).group_id_interface = group_defs(j).id + int32(1000);
%!         bearing_surf(num_comp_mat).absolute_tolerance = 1e-3 * param.d;
%!         bearing_surf(num_comp_mat).relative_tolerance = 0;
%!         bearing_surf(num_comp_mat).name = group_defs(j).name;
%!         bearing_surf(num_comp_mat).bearing = group_defs(j).bearing;
%!         bearing_surf(num_comp_mat).X0 = group_defs(j).X0;
%!         bearing_surf(num_comp_mat).R = group_defs(j).R;
%!         bearing_surf(num_comp_mat).options = group_defs(j).compliance_matrix;
%!         bearing_surf(num_comp_mat).master_node_no = node_set(j);
%!         switch group_defs(j).type
%!           case "cylinder"
%!             bearing_surf(num_comp_mat).r = mean([group_defs(j).geometry.rmax, group_defs(j).geometry.rmin]);
%!             bearing_surf(num_comp_mat).w = group_defs(j).geometry.zmax - group_defs(j).geometry.zmin;
%!           otherwise
%!             error("bearing geometry type \"%s\" not implemented", group_defs(j).type);
%!         endswitch
%!         bearing_surf(num_comp_mat).nodes = mesh.groups.tria6(find([[mesh.groups.tria6].id] == group_defs(j).id)).nodes;
%!     endswitch
%!   endfor
%!   bearing_surf = bearing_surf(1:num_comp_mat);
%!   [mesh, mat_ass, dof_map, cms_opt, comp_mat, bearing_surf, sol_eig] = fem_ehd_pre_comp_mat_linear_mesh(mesh, load_case, cms_opt, bearing_surf);
%!   fem_cms_export([output_file, "_shaft_cms"], mesh, dof_map, mat_ass, cms_opt);
%!   for j=1:numel(comp_mat)
%!     comp_mat_file = [output_file, "_", bearing_surf(j).bearing, "_", bearing_surf(j).options.bearing_type, ".dat"];
%!     fem_ehd_pre_comp_mat_export(comp_mat(j), bearing_surf(j).options, comp_mat_file);
%!   endfor
%!   empty_cell = cell(1, 4);
%!   cms_data = struct("mesh", empty_cell, "dof_map", empty_cell, "cms_opt", empty_cell, "mat_ass", empty_cell, "load_case", empty_cell);
%!   cms_data(1).mesh = mesh;
%!   cms_data(1).dof_map = dof_map;
%!   cms_data(1).cms_opt = cms_opt;
%!   cms_data(1).mat_ass = mat_ass;
%!   cms_data(1).load_case = load_case;
%!   cms_data(2).mesh = mesh2;
%!   cms_data(2).dof_map = dof_map2;
%!   cms_data(2).cms_opt = cms_opt2;
%!   cms_data(2).mat_ass = mat_ass2;
%!   cms_data(2).load_case = load_case2;
%!   cms_data(3).mesh = mesh3;
%!   cms_data(3).dof_map = dof_map3;
%!   cms_data(3).cms_opt = cms_opt3;
%!   cms_data(3).mat_ass = mat_ass3;
%!   cms_data(3).load_case = load_case3;
%!   cms_data(3).cms_opt.element.name = [cms_data(3).cms_opt.element.name, "_1"];
%!   cms_data(3).cms_opt.nodes.modal.name = [cms_data(3).cms_opt.nodes.modal.name, "_1"];
%!   cms_data(3).cms_opt.nodes.interfaces.name = [cms_data(3).cms_opt.nodes.interfaces.name, "_1"];
%!   cms_data(4).mesh = mesh3;
%!   cms_data(4).dof_map = dof_map3;
%!   cms_data(4).cms_opt = cms_opt3;
%!   cms_data(4).mat_ass = mat_ass3;
%!   cms_data(4).load_case = load_case3;
%!   cms_data(4).cms_opt.element.name = [cms_data(4).cms_opt.element.name, "_2"];
%!   cms_data(4).cms_opt.nodes.modal.name = [cms_data(4).cms_opt.nodes.modal.name, "_2"];
%!   cms_data(4).cms_opt.nodes.interfaces.name = [cms_data(4).cms_opt.nodes.interfaces.name, "_2"];
%!   options_mbdyn.f_run_mbdyn2easyanim = false;
%!   empty_cell = cell(numel(param.omega), numel(param.F1));
%!   data = struct("res", empty_cell);
%!   for j=1:numel(param.F1)
%!     for i=1:numel(param.omega)
%!       options_mbdyn.output_file = sprintf("%s_%d_%d_mbd", output_file, i, j);
%!       fd = -1;
%!       unwind_protect
%!         fd = fopen([options_mbdyn.output_file, ".mbd"], "w");
%!         if (fd == -1)
%!           error("failed to open file \"%s\"", [options_mbdyn.output_file, ".mbd"]);
%!         endif
%!         fprintf(fd, "set: real omega = %g;\n", param.omega(i));
%!         fputs(fd, "set: real t2 = pi / omega;\n");
%!         fputs(fd, "set: real t1 = 0.5 * t2;\n");
%!         fputs(fd, "set: real dt = 2. * pi / omega / 100;\n");
%!         fprintf(fd, "set: real l1 = %g;\n", param.l1);
%!         fprintf(fd, "set: real l2 = %g;\n", param.l2);
%!         fprintf(fd, "set: real l3 = %g;\n", param.l3);
%!         fprintf(fd, "set: real l14 = %g;\n", param.l14);
%!         fprintf(fd, "set: real pcl = %g;\n", param.pcl);
%!         fprintf(fd, "set: real etal = %g;\n", param.etal);
%!         fprintf(fd, "set: real betal = %g;\n", param.betal);
%!         fprintf(fd, "set: real rhol = %g;\n", param.rhol);
%!         fprintf(fd, "set: real Tl = %g;\n", param.Tl);
%!         fprintf(fd, "set: real fact_etav = %g;\n", param.fact_etav);
%!         fprintf(fd, "set: real d = %g;\n", param.d);
%!         fprintf(fd, "set: real w = %g;\n", param.w);
%!         fprintf(fd, "set: real Psi = %g;\n", param.Psi);
%!         fprintf(fd, "set: real ds = %g;\n", param.ds);
%!         fprintf(fd, "set: real ws = %g;\n", param.ws);
%!         fprintf(fd, "set: real Psis = %g;\n", param.Psis);
%!         fprintf(fd, "set: real wg = %g;\n", param.wg);
%!         fprintf(fd, "set: real dg = %g;\n", param.dg);
%!         fprintf(fd, "set: real pside = %g;\n", param.pside);
%!         fprintf(fd, "set: real pin = %g;\n", param.pin);
%!         fprintf(fd, "set: real sigma = %g;\n", param.sigma);
%!         fprintf(fd, "set: real eta = %g;\n", param.eta);
%!         fprintf(fd, "set: real beta = %g;\n", param.beta);
%!         fprintf(fd, "set: real delta = %g;\n", param.delta);
%!         fprintf(fd, "set: real Rq1 = %g;\n", param.Rq1);
%!         fprintf(fd, "set: real Rq2 = %g;\n", param.Rq2);
%!         fprintf(fd, "set: real gamma1 = %g;\n", param.gamma1);
%!         fprintf(fd, "set: real gamma2 = %g;\n", param.gamma2);
%!         fprintf(fd, "set: real Ered = %g;\n", param.Ered);
%!         fprintf(fd, "set: real E = %g;\n", param.E);
%!         fprintf(fd, "set: real nu = %g;\n", param.nu);
%!         fprintf(fd, "set: real mu = %g;\n", param.mu);
%!         fprintf(fd, "set: real sigma0 = %g;\n", param.sigma0);
%!         fprintf(fd, "set: real pref = %g;\n", param.pref);
%!         fprintf(fd, "set: real F1 = %g;\n", param.F1(j));
%!         fputs(fd, "set: integer ref_id_assembly = 1001;\n");
%!         fputs(fd, "set: integer ref_id_shaft = 1002;\n");
%!         fputs(fd, "set: integer ref_id_drive = 1003;\n");
%!         fputs(fd, "set: integer ref_id_big_end_bearing_journal = 1004;\n");
%!         fputs(fd, "set: integer ref_id_small_end_bearing_shell = 1005;\n");
%!         fputs(fd, "set: integer ref_id_support_bearing_journal_1 = 1006;\n");
%!         fputs(fd, "set: integer ref_id_support_bearing_journal_2 = 1007;\n");
%!         fputs(fd, "set: integer ref_id_supporting_area_1 = 1008;\n");
%!         fputs(fd, "set: integer ref_id_supporting_area_2 = 1009;\n");
%!         fputs(fd, "set: integer node_id_drive = 2001;\n");
%!         fputs(fd, "set: integer node_id_big_end_bearing_journal = 2002;\n");
%!         fputs(fd, "set: integer node_id_support_bearing_journal_1 = 2003;\n");
%!         fputs(fd, "set: integer node_id_support_bearing_journal_2 = 2004;\n");
%!         fputs(fd, "set: integer node_id_small_end_bearing_shell = 2005;\n");
%!         fputs(fd, "set: integer node_id_big_end_bearing_shell = 2006;\n");
%!         fputs(fd, "set: integer node_id_supporting_area_1 = 2007;\n");
%!         fputs(fd, "set: integer node_id_support_bearing_shell_1 = 2008;\n");
%!         fputs(fd, "set: integer node_id_supporting_area_2 = 2009;\n");
%!         fputs(fd, "set: integer node_id_support_bearing_shell_2 = 2010;\n");
%!         fputs(fd, "set: integer joint_id_drive = 3001;\n");
%!         fputs(fd, "set: integer joint_id_small_end_bearing = 3002;\n");
%!         fputs(fd, "set: integer joint_id_supporting_area1 = 3003;\n");
%!         fputs(fd, "set: integer joint_id_supporting_area2 = 3004;\n");
%!         fputs(fd, "set: integer elem_id_shaft = 3005;\n");
%!         fputs(fd, "set: integer elem_id_conrod = 3006;\n");
%!         fputs(fd, "set: integer elem_id_support_1 = 3007;\n");
%!         fputs(fd, "set: integer elem_id_support_2 = 3008;\n");
%!         fputs(fd, "set: integer elem_id_support_bearing_1 = 4001;\n");
%!         fputs(fd, "set: integer elem_id_support_bearing_2 = 4002;\n");
%!         fputs(fd, "set: integer elem_id_big_end_bearing = 4003;\n");
%!         fputs(fd, "set: integer force_id_load = 5001;\n");
%!         fputs(fd, "begin: data;\n");
%!         fputs(fd, "        problem: initial value;\n");
%!         fputs(fd, "end: data;\n");
%!         fputs(fd, "begin: initial value;\n");
%!         fputs(fd, "        initial time: 0;\n");
%!         fputs(fd, "        final time: t2;\n");
%!         fputs(fd, "        time step: dt;\n");
%!         fputs(fd, "        method: hybrid, ms, 0.6;\n");
%!         fputs(fd, "        linear solver: pardiso, grad, scale, row max column max, always, max iterations, 250;\n");
%!         fprintf(fd, "      threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!         fprintf(fd, "      threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!         fputs(fd, "          nonlinear solver: linesearch, modified, 0, default solver options, heavy nonlinear, divergence check, no, lambda min, 0., tolerance x, 1e-4, verbose, yes, print convergence info, yes;\n");
%!         fputs(fd, "        tolerance: 1e-4, test, minmax, 1e-4, test, norm;\n");
%!         fputs(fd, "        max iterations: 100;\n");
%!         fputs(fd, "        derivatives tolerance: 1e-3;\n");
%!         fputs(fd, "        derivatives max iterations: 5;\n");
%!         fputs(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!         fputs(fd, "        output: iterations, cpu time, solver condition number, stat, yes;\n");
%!         fputs(fd, "end: initial value;\n");
%!         fputs(fd, "begin: control data;\n");
%!         fputs(fd, "       ## output meter: closest next, t2, forever, dt;\n");
%!         fputs(fd, "       use automatic differentiation;\n");
%!         fputs(fd, "       output precision: 16;\n");
%!         fputs(fd, "       max iterations: 0;\n");
%!         fputs(fd, "       print: dof stats, to file;\n");
%!         fputs(fd, "       print: dof description, to file;\n");
%!         fputs(fd, "       print: equation description, to file;\n");
%!         fputs(fd, "       structural nodes: 10;\n");
%!         fputs(fd, "       joints: 8;\n");
%!         fputs(fd, "       loadable elements: 3;\n");
%!         fputs(fd, "       forces: 1;\n");
%!         fputs(fd, "end: control data;\n");
%!         fputs(fd, "reference: ref_id_assembly,\n");
%!         fputs(fd, "        position, reference, global, null,\n");
%!         fputs(fd, "        orientation, reference, global, eye,\n");
%!         fputs(fd, "        velocity, reference, global, null,\n");
%!         fputs(fd, "        angular velocity, reference, global, null;\n");
%!         fputs(fd, "reference: ref_id_shaft,\n");
%!         fputs(fd, "        position, reference, ref_id_assembly, null,\n");
%!         fputs(fd, "        orientation, reference, ref_id_assembly, eye,\n");
%!         fputs(fd, "        velocity, reference, ref_id_assembly, null,\n");
%!         fputs(fd, "        angular velocity, reference, ref_id_assembly,\n");
%!         fputs(fd, "                omega,\n");
%!         fputs(fd, "                0.,\n");
%!         fputs(fd, "                0.;\n");
%!         fputs(fd, "reference: ref_id_drive,\n");
%!         fputs(fd, "        position, reference, ref_id_shaft, -l2, 0., 0.,\n");
%!         fputs(fd, "        orientation, reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        velocity, reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        angular velocity, reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_big_end_bearing_journal,\n");
%!         fputs(fd, "        position, reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        orientation, reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        velocity, reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        angular velocity, reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_small_end_bearing_shell,\n");
%!         fputs(fd, "        position, reference, ref_id_shaft, 0., 0., l3,\n");
%!         fputs(fd, "        orientation, reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        velocity, reference, ref_id_assembly, null,\n");
%!         fputs(fd, "        angular velocity, reference, ref_id_assembly, null;\n");
%!         fputs(fd, "reference: ref_id_support_bearing_journal_1,\n");
%!         fputs(fd, "        position, reference, ref_id_shaft, -l1, 0., 0.,\n");
%!         fputs(fd, "        orientation, reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        velocity, reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        angular velocity, reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_support_bearing_journal_2,\n");
%!         fputs(fd, "        position, reference, ref_id_shaft, l1, 0., 0.,\n");
%!         fputs(fd, "        orientation, reference, ref_id_shaft, eye,\n");
%!         fputs(fd, "        velocity, reference, ref_id_shaft, null,\n");
%!         fputs(fd, "        angular velocity, reference, ref_id_shaft, null;\n");
%!         fputs(fd, "reference: ref_id_supporting_area_1,\n");
%!         fputs(fd, "        position, reference, ref_id_support_bearing_journal_1, 0., 0., -l14,\n");
%!         fputs(fd, "        orientation, reference, ref_id_support_bearing_journal_1, eye,\n");
%!         fputs(fd, "        velocity, reference, ref_id_assembly, null,\n");
%!         fputs(fd, "        angular velocity, reference, ref_id_assembly, null;\n");
%!         fputs(fd, "reference: ref_id_supporting_area_2,\n");
%!         fputs(fd, "        position, reference, ref_id_support_bearing_journal_2, 0., 0., -l14,\n");
%!         fputs(fd, "        orientation, reference, ref_id_support_bearing_journal_2, eye,\n");
%!         fputs(fd, "        velocity, reference, ref_id_assembly, null,\n");
%!         fputs(fd, "        angular velocity, reference, ref_id_assembly, null;\n");
%!         fputs(fd, "begin: nodes;\n");
%!         fputs(fd, "        structural: node_id_drive, modal,\n");
%!         fputs(fd, "                position, reference, ref_id_drive, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_drive, eye,\n");
%!         fputs(fd, "                velocity, reference, ref_id_drive, null,\n");
%!         fputs(fd, "                angular velocity, reference, ref_id_drive, null;\n");
%!         fputs(fd, "        structural: node_id_big_end_bearing_journal, static,\n");
%!         fputs(fd, "                position, reference, ref_id_big_end_bearing_journal, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_big_end_bearing_journal, eye,\n");
%!         fputs(fd, "                velocity, reference, ref_id_big_end_bearing_journal, null,\n");
%!         fputs(fd, "                angular velocity, reference, ref_id_big_end_bearing_journal, null;\n");
%!         fputs(fd, "        structural: node_id_support_bearing_journal_1, static,\n");
%!         fputs(fd, "                position, reference, ref_id_support_bearing_journal_1, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_support_bearing_journal_1, eye,\n");
%!         fputs(fd, "                velocity, reference, ref_id_support_bearing_journal_1, null,\n");
%!         fputs(fd, "                angular velocity, reference, ref_id_support_bearing_journal_1, null;\n");
%!         fputs(fd, "        structural: node_id_support_bearing_journal_2, static,\n");
%!         fputs(fd, "                position, reference, ref_id_support_bearing_journal_2, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_support_bearing_journal_2, eye,\n");
%!         fputs(fd, "                velocity, reference, ref_id_support_bearing_journal_2, null,\n");
%!         fputs(fd, "                angular velocity, reference, ref_id_support_bearing_journal_2, null;\n");
%!         fputs(fd, "        structural: node_id_small_end_bearing_shell, modal,\n");
%!         fputs(fd, "                position, reference, ref_id_small_end_bearing_shell, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_small_end_bearing_shell, eye,\n");
%!         fputs(fd, "                velocity, reference, ref_id_small_end_bearing_shell, null,\n");
%!         fputs(fd, "                angular velocity, reference, ref_id_small_end_bearing_shell, null;\n");
%!         fputs(fd, "        structural: node_id_big_end_bearing_shell, static,\n");
%!         fputs(fd, "                position, reference, ref_id_big_end_bearing_journal, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_big_end_bearing_journal, eye,\n");
%!         fputs(fd, "                velocity, reference, ref_id_assembly, null,\n");
%!         fputs(fd, "                angular velocity, reference, ref_id_assembly, null;\n");
%!         fputs(fd, "        structural: node_id_supporting_area_1, modal,\n");
%!         fputs(fd, "                position, reference, ref_id_supporting_area_1, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_supporting_area_1, eye,\n");
%!         fputs(fd, "                velocity, reference, ref_id_supporting_area_1, null,\n");
%!         fputs(fd, "                angular velocity, reference, ref_id_supporting_area_1, null;\n");
%!         fputs(fd, "        structural: node_id_support_bearing_shell_1, static,\n");
%!         fputs(fd, "                position, reference, ref_id_support_bearing_journal_1, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_support_bearing_journal_1, eye,\n");
%!         fputs(fd, "                velocity, reference, ref_id_assembly, null,\n");
%!         fputs(fd, "                angular velocity, reference, ref_id_assembly, null;\n");
%!         fputs(fd, "        structural: node_id_supporting_area_2, modal,\n");
%!         fputs(fd, "                position, reference, ref_id_supporting_area_2, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_supporting_area_2, eye,\n");
%!         fputs(fd, "                velocity, reference, ref_id_supporting_area_2, null,\n");
%!         fputs(fd, "                angular velocity, reference, ref_id_supporting_area_2, null;\n");
%!         fputs(fd, "        structural: node_id_support_bearing_shell_2, static,\n");
%!         fputs(fd, "                position, reference, ref_id_support_bearing_journal_2, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_support_bearing_journal_2, eye,\n");
%!         fputs(fd, "                velocity, reference, ref_id_assembly, null,\n");
%!         fputs(fd, "                angular velocity, reference, ref_id_assembly, null;\n");
%!         fputs(fd, "end: nodes;\n");
%!         fputs(fd, "begin: elements;\n");
%!         fputs(fd, "        joint: joint_id_drive, total pin joint,\n");
%!         fputs(fd, "                node_id_drive,\n");
%!         fputs(fd, "                position, reference, ref_id_drive, null,\n");
%!         fputs(fd, "                position orientation, reference, ref_id_drive, eye,\n");
%!         fputs(fd, "                rotation orientation, reference, ref_id_drive, eye,\n");
%!         fputs(fd, "                position, reference, ref_id_drive, null,\n");
%!         fputs(fd, "                position orientation, reference, ref_id_drive, eye,\n");
%!         fputs(fd, "                rotation orientation, reference, ref_id_drive, eye,\n");
%!         fputs(fd, "                position constraint,\n");
%!         fputs(fd, "                        active,\n");
%!         fputs(fd, "                        inactive,\n");
%!         fputs(fd, "                        inactive,\n");
%!         fputs(fd, "                        component,\n");
%!         fputs(fd, "                                null,\n");
%!         fputs(fd, "                                null,\n");
%!         fputs(fd, "                                null,\n");
%!         fputs(fd, "                orientation constraint,\n");
%!         fputs(fd, "                        angular velocity,\n");
%!         fputs(fd, "                        inactive,\n");
%!         fputs(fd, "                        inactive,\n");
%!         fputs(fd, "                        component,\n");
%!         fputs(fd, "                                omega,\n");
%!         fputs(fd, "                                null,\n");
%!         fputs(fd, "                                null;\n");
%!         fputs(fd, "        joint: joint_id_small_end_bearing, total pin joint,\n");
%!         fputs(fd, "                node_id_small_end_bearing_shell,\n");
%!         fputs(fd, "                position, reference, ref_id_small_end_bearing_shell, null,\n");
%!         fputs(fd, "                position orientation, reference, ref_id_small_end_bearing_shell, eye,\n");
%!         fputs(fd, "                rotation orientation, reference, ref_id_small_end_bearing_shell, eye,\n");
%!         fputs(fd, "                position, reference, ref_id_small_end_bearing_shell, null,\n");
%!         fputs(fd, "                position orientation, reference, ref_id_small_end_bearing_shell, eye,\n");
%!         fputs(fd, "                rotation orientation, reference, ref_id_small_end_bearing_shell, eye,\n");
%!         fputs(fd, "                position constraint,\n");
%!         fputs(fd, "                        active,\n");
%!         fputs(fd, "                        active,\n");
%!         fputs(fd, "                        inactive,\n");
%!         fputs(fd, "                        null,\n");
%!         fputs(fd, "                orientation constraint,\n");
%!         fputs(fd, "                        inactive,\n");
%!         fputs(fd, "                        active,\n");
%!         fputs(fd, "                        active,\n");
%!         fputs(fd, "                        component,\n");
%!         fputs(fd, "                                null,\n");
%!         fputs(fd, "                                null,\n");
%!         fputs(fd, "                                null;\n");
%!         fputs(fd, "\n");
%!         fputs(fd, "        force: force_id_load, absolute,\n");
%!         fputs(fd, "          node_id_small_end_bearing_shell,\n");
%!         fputs(fd, "          position,reference,node, null,\n");
%!         fputs(fd, "          0.,0.,-1., string, \"F1 * (sin(pi / 2. * Time / t1)^2 * (Time <= t1) + (Time > t1))\";\n");
%!         fputs(fd, "\n");
%!         fputs(fd, "        joint: joint_id_supporting_area1, clamp, node_id_supporting_area_1, node, node;\n");
%!         fputs(fd, "        joint: joint_id_supporting_area2, clamp, node_id_supporting_area_2, node, node;\n");
%!         fprintf(fd, "        include: \"%s_shaft_cms.elm\";\n", output_file);
%!         fprintf(fd, "        include: \"%s_conrod_cms.elm\";\n", output_file);
%!         fprintf(fd, "        include: \"%s_support_cms_1.elm\";\n", output_file);
%!         fprintf(fd, "        include: \"%s_support_cms_2.elm\";\n", output_file);
%!         fputs(fd, "\n");
%!         fputs(fd, "        user defined: elem_id_support_bearing_1, hydrodynamic plain bearing2,\n");
%!         fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!         fputs(fd, "                density, rhol,\n");
%!         fputs(fd, "                betal,\n");
%!         fputs(fd, "                pressure, pcl,\n");
%!         fputs(fd, "                viscosity, etal,\n");
%!         fputs(fd, "                temperature, Tl,\n");
%!         fputs(fd, "            viscosity vapor, factor, fact_etav,\n");
%!         fputs(fd, "                mesh, linear finite difference,\n");
%!         fputs(fd, "                enable mcp, no,\n");
%!         if (param.enable_Patir_Cheng)
%!           fputs(fd, "                flow factors, patir cheng, sigma, Rq1, Rq2, lambdax, gamma1, gamma2, lambdaz, 1., 1.,\n");
%!         endif
%!         fputs(fd, "                geometry, cylindrical,\n");
%!         fputs(fd, "                        mesh position, at bearing,\n");
%!         fputs(fd, "                        bearing width, ws,\n");
%!         fputs(fd, "                        shaft diameter, ds,\n");
%!         fputs(fd, "                        bearing diameter, ds * (1. + Psis),\n");
%!         fputs(fd, "                shaft node, node_id_drive,\n");
%!         fputs(fd, "                offset, reference, ref_id_support_bearing_journal_1, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_support_bearing_journal_1, 3, -1., 0., 0.,\n");
%!         fputs(fd, "                                                                          2, 0., 1., 0.,\n");
%!         fputs(fd, "                bearing node, node_id_supporting_area_1,\n");
%!         fputs(fd, "                offset, reference, ref_id_support_bearing_journal_1, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_support_bearing_journal_1, 3, -1., 0., 0.,\n");
%!         fputs(fd, "                                                                          2, 0., 1., 0.,\n");
%!         fprintf(fd, "                        number of nodes z, %d,\n", numel(comp_mat3.bearing_surf.grid_z));
%!         fprintf(fd, "                        number of nodes Phi, %d,\n", numel(comp_mat3.bearing_surf.grid_x) + 1);
%!         fputs(fd, "                        boundary conditions,\n");
%!         fputs(fd, "                                        pressure, pside,\n");
%!         fputs(fd, "                                        pressure, pside,\n");
%!         fputs(fd, "                        lubrication grooves, 1,\n");
%!         fputs(fd, "                                at bearing,\n");
%!         fputs(fd, "                                        pressure, pin,\n");
%!         fputs(fd, "                                        position, 0., 0.,\n");
%!         fputs(fd, "                                        rectangle, width, ds * pi / 2., height, wg,\n");
%!         fputs(fd, "            contact model,\n");
%!         fputs(fd, "              greenwood tripp,\n");
%!         fputs(fd, "                 E1, 2 * (1 - nu^2) * Ered,\n");
%!         fputs(fd, "                 nu1, nu,\n");
%!         fputs(fd, "                 E2, 2 * (1 - nu^2) * Ered,\n");
%!         fputs(fd, "                 nu2, nu,\n");
%!         fputs(fd, "                      unit system, consistent,\n");
%!         fputs(fd, "                      standard deviation, sigma,\n");
%!         fputs(fd, "                      asperity density, eta,\n");
%!         fputs(fd, "                      asperity curvature, beta,\n");
%!         fputs(fd, "                      asperity mean height, delta,\n");
%!         fputs(fd, "              friction model, lugre,\n");
%!         fputs(fd, "                 method, implicit euler,\n");
%!         fputs(fd, "                 coulomb friction coefficient, mu,\n");
%!         fputs(fd, "                 micro slip stiffness, sigma0,\n");
%!         fputs(fd, "                compliance model,\n");
%!         fputs(fd, "                  double nodal,\n");
%!         fprintf(fd, "                  matrix at shaft, from file, \"%s_elem_id_support_bearing_1_journal.dat\",\n", output_file);
%!         fputs(fd, "                  modal element, elem_id_shaft,\n");
%!         fprintf(fd, "                  matrix at bearing, from file, \"%s_elem_id_support_bearing_shell.dat\",\n", output_file);
%!         fputs(fd, "                  modal element, elem_id_support_1,\n");
%!         fputs(fd, "                  axial displacement, small,\n");
%!         fputs(fd, "                pressure dof scale, pref,\n");
%!         fputs(fd, "                reynolds equation scale, dt / (rhol * Psis * ds),\n");
%!         fputs(fd, "                elasticity equation scale, dt / (Psis * ds),\n");
%!         fputs(fd, "                output pressure, yes,\n");
%!         fputs(fd, "                output contact pressure, yes,\n");
%!         fputs(fd, "                output stress, no,\n");
%!         fputs(fd, "                output density, yes,\n");
%!         fputs(fd, "                output friction loss, yes,\n");
%!         fputs(fd, "                output clearance, yes,\n");
%!         fputs(fd, "                output reaction force, yes,\n");
%!         fputs(fd, "                output mesh, yes,\n");
%!         fputs(fd, "                output total deformation, yes,\n");
%!         fputs(fd, "                output, yes;\n");
%!         fputs(fd, "\n");
%!         fputs(fd, "        user defined: elem_id_support_bearing_2, hydrodynamic plain bearing2,\n");
%!         fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!         fputs(fd, "                density, rhol,\n");
%!         fputs(fd, "                betal,\n");
%!         fputs(fd, "                pressure, pcl,\n");
%!         fputs(fd, "                viscosity, etal,\n");
%!         fputs(fd, "                temperature, Tl,\n");
%!         fputs(fd, "            viscosity vapor, factor, fact_etav,\n");
%!         fputs(fd, "                mesh, linear finite difference,\n");
%!         fputs(fd, "                enable mcp, no,\n");
%!         if (param.enable_Patir_Cheng)
%!           fputs(fd, "                flow factors, patir cheng, sigma, Rq1, Rq2, lambdax, gamma1, gamma2, lambdaz, 1., 1.,\n");
%!         endif
%!         fputs(fd, "                geometry, cylindrical,\n");
%!         fputs(fd, "                        mesh position, at bearing,\n");
%!         fputs(fd, "                        bearing width, ws,\n");
%!         fputs(fd, "                        shaft diameter, ds,\n");
%!         fputs(fd, "                        bearing diameter, ds * (1. + Psis),\n");
%!         fputs(fd, "                shaft node, node_id_drive,\n");
%!         fputs(fd, "                offset, reference, ref_id_support_bearing_journal_2, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_support_bearing_journal_2, 3, -1., 0., 0.,\n");
%!         fputs(fd, "                                                                          2, 0., 1., 0.,\n");
%!         fputs(fd, "                bearing node, node_id_supporting_area_2,\n");
%!         fputs(fd, "                offset, reference, ref_id_support_bearing_journal_2, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_support_bearing_journal_2, 3, -1., 0., 0.,\n");
%!         fputs(fd, "                                                                          2, 0., 1., 0.,\n");
%!         fprintf(fd, "                        number of nodes z, %d,\n", numel(comp_mat3.bearing_surf.grid_z));
%!         fprintf(fd, "                        number of nodes Phi, %d,\n", numel(comp_mat3.bearing_surf.grid_x) + 1);
%!         fputs(fd, "                        boundary conditions,\n");
%!         fputs(fd, "                                        pressure, pside,\n");
%!         fputs(fd, "                                        pressure, pside,\n");
%!         fputs(fd, "                        lubrication grooves, 1,\n");
%!         fputs(fd, "                                at bearing,\n");
%!         fputs(fd, "                                        pressure, pin,\n");
%!         fputs(fd, "                                        position, 0., 0.,\n");
%!         fputs(fd, "                                        rectangle, width, ds * pi / 2., height, wg,\n");
%!         fputs(fd, "            contact model,\n");
%!         fputs(fd, "              greenwood tripp,\n");
%!         fputs(fd, "                 E1, 2 * (1 - nu^2) * Ered,\n");
%!         fputs(fd, "                 nu1, nu,\n");
%!         fputs(fd, "                 E2, 2 * (1 - nu^2) * Ered,\n");
%!         fputs(fd, "                 nu2, nu,\n");
%!         fputs(fd, "                      unit system, consistent,\n");
%!         fputs(fd, "                      standard deviation, sigma,\n");
%!         fputs(fd, "                      asperity density, eta,\n");
%!         fputs(fd, "                      asperity curvature, beta,\n");
%!         fputs(fd, "                      asperity mean height, delta,\n");
%!         fputs(fd, "              friction model, lugre,\n");
%!         fputs(fd, "                 method, implicit euler,\n");
%!         fputs(fd, "                 coulomb friction coefficient, mu,\n");
%!         fputs(fd, "                 micro slip stiffness, sigma0,\n");
%!         fputs(fd, "                compliance model,\n");
%!         fputs(fd, "                  double nodal,\n");
%!         fprintf(fd, "                  matrix at shaft, from file, \"%s_elem_id_support_bearing_2_journal.dat\",\n", output_file);
%!         fputs(fd, "                  modal element, elem_id_shaft,\n");
%!         fprintf(fd, "                  matrix at bearing, from file, \"%s_elem_id_support_bearing_shell.dat\",\n", output_file);
%!         fputs(fd, "                  modal element, elem_id_support_2,\n");
%!         fputs(fd, "                  axial displacement, small,\n");
%!         fputs(fd, "                pressure dof scale, pref,\n");
%!         fputs(fd, "                reynolds equation scale, dt / (rhol * Psis * ds),\n");
%!         fputs(fd, "                elasticity equation scale, dt / (Psis * ds),\n");
%!         fputs(fd, "                output pressure, yes,\n");
%!         fputs(fd, "                output contact pressure, yes,\n");
%!         fputs(fd, "                output stress, no,\n");
%!         fputs(fd, "                output density, yes,\n");
%!         fputs(fd, "                output friction loss, yes,\n");
%!         fputs(fd, "                output clearance, yes,\n");
%!         fputs(fd, "                output reaction force, yes,\n");
%!         fputs(fd, "                output mesh, yes,\n");
%!         fputs(fd, "                output total deformation, yes,\n");
%!         fputs(fd, "                output, yes;\n");
%!         fputs(fd, "\n");
%!         fputs(fd, "        user defined: elem_id_big_end_bearing, hydrodynamic plain bearing2,\n");
%!         fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!         fputs(fd, "                density, rhol,\n");
%!         fputs(fd, "                betal,\n");
%!         fputs(fd, "                pressure, pcl,\n");
%!         fputs(fd, "                viscosity, etal,\n");
%!         fputs(fd, "                temperature, Tl,\n");
%!         fputs(fd, "            viscosity vapor, factor, fact_etav,\n");
%!         fputs(fd, "                mesh, linear finite difference,\n");
%!         fputs(fd, "                enable mcp, no,\n");
%!         if (param.enable_Patir_Cheng)
%!           fputs(fd, "                flow factors, patir cheng, sigma, Rq1, Rq2, lambdax, gamma1, gamma2, lambdaz, 1., 1.,\n");
%!         endif
%!         fputs(fd, "                geometry, cylindrical,\n");
%!         fputs(fd, "                        mesh position, at bearing,\n");
%!         fputs(fd, "                        bearing width, w,\n");
%!         fputs(fd, "                        shaft diameter, d,\n");
%!         fputs(fd, "                        bearing diameter, d * (1. + Psi),\n");
%!         fputs(fd, "                shaft node, node_id_drive,\n");
%!         fputs(fd, "                offset, reference, ref_id_big_end_bearing_journal, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_big_end_bearing_journal, 3, -1., 0., 0.,\n");
%!         fputs(fd, "                                                                          2, 0., 1., 0.,\n");
%!         fputs(fd, "                bearing node, node_id_small_end_bearing_shell,\n");
%!         fputs(fd, "                offset, reference, ref_id_big_end_bearing_journal, null,\n");
%!         fputs(fd, "                orientation, reference, ref_id_big_end_bearing_journal, 3, -1., 0., 0.,\n");
%!         fputs(fd, "                                                                          2, 0., 1., 0.,\n");
%!         fprintf(fd, "                        number of nodes z, %d,\n", numel(comp_mat2.bearing_surf.grid_z));
%!         fprintf(fd, "                        number of nodes Phi, %d,\n", numel(comp_mat2.bearing_surf.grid_x) + 1);
%!         fputs(fd, "                        boundary conditions,\n");
%!         fputs(fd, "                                        pressure, pside,\n");
%!         fputs(fd, "                                        pressure, pside,\n");
%!         fputs(fd, "                        lubrication grooves, 1,\n");
%!         fputs(fd, "                                at bearing,\n");
%!         fputs(fd, "                                        pressure, pin,\n");
%!         fputs(fd, "                                        position, d * (1. + Psi) * pi / 2, 0.,\n");
%!         fputs(fd, "                                        circle, radius, dg / 2.,\n");
%!         fputs(fd, "            contact model,\n");
%!         fputs(fd, "              greenwood tripp,\n");
%!         fputs(fd, "                 E1, 2 * (1 - nu^2) * Ered,\n");
%!         fputs(fd, "                 nu1, nu,\n");
%!         fputs(fd, "                 E2, 2 * (1 - nu^2) * Ered,\n");
%!         fputs(fd, "                 nu2, nu,\n");
%!         fputs(fd, "                      unit system, consistent,\n");
%!         fputs(fd, "                      standard deviation, sigma,\n");
%!         fputs(fd, "                      asperity density, eta,\n");
%!         fputs(fd, "                      asperity curvature, beta,\n");
%!         fputs(fd, "                      asperity mean height, delta,\n");
%!         fputs(fd, "              friction model, lugre,\n");
%!         fputs(fd, "                 method, implicit euler,\n");
%!         fputs(fd, "                 coulomb friction coefficient, mu,\n");
%!         fputs(fd, "                 micro slip stiffness, sigma0,\n");
%!         fputs(fd, "                compliance model,\n");
%!         fputs(fd, "                  double nodal,\n");
%!         fprintf(fd, "                  matrix at shaft, from file, \"%s_elem_id_big_end_bearing_journal.dat\",\n", output_file);
%!         fputs(fd, "                  modal element, elem_id_shaft,\n");
%!         fprintf(fd, "                  matrix at bearing, from file, \"%s_elem_id_big_end_bearing_shell.dat\",\n", output_file);
%!         fputs(fd, "                  modal element, elem_id_conrod,\n");
%!         fputs(fd, "                  axial displacement, small,\n");
%!         fputs(fd, "                pressure dof scale, pref,\n");
%!         fputs(fd, "                reynolds equation scale, dt / (rhol * Psi * d),\n");
%!         fputs(fd, "                elasticity equation scale, dt / (Psi * d),\n");
%!         fputs(fd, "                output pressure, yes,\n");
%!         fputs(fd, "                output contact pressure, yes,\n");
%!         fputs(fd, "                output stress, no,\n");
%!         fputs(fd, "                output density, yes,\n");
%!         fputs(fd, "                output friction loss, yes,\n");
%!         fputs(fd, "                output clearance, yes,\n");
%!         fputs(fd, "                output reaction force, yes,\n");
%!         fputs(fd, "                output mesh, yes,\n");
%!         fputs(fd, "                output total deformation, yes,\n");
%!         fputs(fd, "                output, yes;\n");
%!         fputs(fd, "end: elements;\n");
%!       unwind_protect_cleanup
%!         if (fd ~= -1)
%!           fclose(fd);
%!         endif
%!       end_unwind_protect
%!       mbdyn_solver_run([options_mbdyn.output_file, ".mbd"], options_mbdyn);
%!       data(i, j).res.log_dat = mbdyn_post_load_log(options_mbdyn.output_file);
%!       [data(i, j).res.t, ...
%!        data(i, j).res.trajectory, ...
%!        data(i, j).res.deformation, ...
%!        data(i, j).res.velocity, ...
%!        data(i, j).res.acceleration, ...
%!        data(i, j).res.node_id, ...
%!        data(i, j).res.force, ...
%!        data(i, j).res.force_id, ...
%!        data(i, j).res.force_node_id, ...
%!        data(i, j).res.orientation_description] = mbdyn_post_load_output_struct(options_mbdyn.output_file);
%!       [data(i, j).res.elem_id, data(i, j).res.q, data(i, j).res.qdot, data(i, j).res.qddot] = mbdyn_post_load_output_mod(options_mbdyn.output_file, numel(data(i, j).res.t));
%!       [data(i, j).res.joint_id, data(i, j).res.local_reaction, data(i, j).res.global_reaction] = mbdyn_post_load_output_jnt(options_mbdyn.output_file);
%!       data(i, j).res.bearings = mbdyn_post_ehd_load_output(options_mbdyn.output_file, data(i, j).res.log_dat);
%!       opt_scale.scale_type = "least square";
%!       opt_scale.scale = 5000;
%!       opt_post.print_and_exit = false;
%!       opt_scale.output_stress = FEM_SCA_STRESS_VMIS;
%!       opt_post.step_start = -1;
%!       opt_post.step_inc = -1;
%!       opt_post.step_end = intmax;
%!       opt_post.elem_types = {"tet10"};
%!       fem_post_cms_sol_export(cms_data, options_mbdyn.output_file, options_mbdyn.output_file, opt_scale, opt_post);
%!       opt_post_h.deformation_scale = 10000;
%!       opt_post_h.start = 1;
%!       opt_post_h.step = 1;
%!       opt_post_h.step_end = intmax;
%!       data(i, j).res.mesh = mbdyn_post_ehd_create_mesh(data(i, j).res.log_dat);
%!       mbdyn_post_ehd_export_mesh(data(i, j).res.mesh, [options_mbdyn.output_file, "_hydro.msh"]);
%!       output_files = mbdyn_post_ehd_export_data(data(i, j).res.mesh, data(i, j).res, [options_mbdyn.output_file, "_hydro.msh"], 1:numel(data(i, j).res.t), opt_post_h);
%!     endfor
%!   endfor
%!   omega_ref_h_10MPa = [3000, 150] * pi / 30 / SI_unit_second^-1;
%!   h_ref_10MPa = [2e-6, 0.7e-6] / SI_unit_meter;
%!   omega_ref_10MPa = [3000,  1000,  500,  400,  300,   250,   200,  150] * pi / 30 / SI_unit_second^-1;
%!   M_ref_10MPa     = [1.06, 0.544, 0.43, 0.46, 0.58, 0.747, 1.076, 1.71] / (SI_unit_newton * SI_unit_meter);
%!   omega_ref_5MPa = [3000,  1000,   500,   400,   300,  250, 150] * pi / 30 / SI_unit_second^-1;
%!   M_ref_5MPa     = [0.96,  0.46, 0.297, 0.257, 0.243, 0.34, 0.5] / (SI_unit_newton * SI_unit_meter);
%!   omega_ref = {omega_ref_5MPa, omega_ref_10MPa};
%!   M_ref = {M_ref_5MPa, M_ref_10MPa};
%!   M_shaft = zeros(numel(param.omega), numel(param.F1));
%!   min_h = zeros(numel(param.omega), numel(param.F1));
%!   for j=1:numel(param.F1)
%!     for i=1:numel(param.omega)
%!       joint_idx = data(i, j).res.log_dat.vars.joint_id_drive == data(i, j).res.joint_id;
%!       bearing_idx = data(i, j).res.log_dat.vars.elem_id_big_end_bearing == [data(i, j).res.bearings.label];
%!       M_shaft(i, j) = -data(i, j).res.global_reaction{joint_idx}(end, 4);
%!       min_h(i, j) = min(data(i, j).res.bearings(bearing_idx).columns.h_n(:));
%!     endfor
%!   endfor
%!   for j=1:numel(param.F1)
%!     figure("visible", "off");
%!     hold on;
%!     plot(omega_ref{j} * 30 / pi * SI_unit_second^-1, M_ref{j} * SI_unit_newton * SI_unit_meter, "-;Sander New;k");
%!     plot(param.omega * 30 / pi * SI_unit_second^-1, M_shaft(:, j) * SI_unit_newton * SI_unit_meter, "-;MBDyn;r");
%!     xlabel("n [rpm]");
%!     ylabel("M [Nm]");
%!     title(sprintf("shaft input torque versus speed F1=%gN", param.F1(j) * SI_unit_newton));
%!     grid on;
%!     grid minor on;
%!   endfor
%!   for j=1:numel(param.F1)
%!     M_bearing = zeros(numel(param.omega), 3);
%!     for i=1:numel(param.omega)
%!       for k=1:3
%!         M_bearing(i, k) = -data(i, j).res.bearings(k).columns.M1(end, 1);
%!       endfor
%!     endfor
%!     figure("visible", "off");
%!     area(param.omega * 30 / pi * SI_unit_second^-1, M_bearing * SI_unit_newton * SI_unit_meter);
%!     xlabel("n [rpm]");
%!     ylabel("M [Nm]");
%!     legend("support bearing 1", "support bearing 2", "big end bearing");
%!     title(sprintf("frictional torque versus speed F1=%gN", param.F1(j) * SI_unit_newton));
%!     grid minor on;
%!   endfor
%!   for j=1:numel(param.F1)
%!     M_big_end_bearing = zeros(numel(param.omega), 2);
%!     for i=1:numel(param.omega)
%!       M_big_end_bearing(i, 1) = -data(i, j).res.bearings(3).columns.Pff(end) / param.omega(i);
%!       M_big_end_bearing(i, 2) = -data(i, j).res.bearings(3).columns.Pfc(end) / param.omega(i);
%!     endfor
%!     M_big_end_bearing(:, 3) = M_big_end_bearing(:, 1) + M_big_end_bearing(:, 2);
%!     figure("visible", "off");
%!     hold on;
%!     plot(param.omega * 30 / pi * SI_unit_second^-1, M_big_end_bearing(:, 1) * SI_unit_newton * SI_unit_meter, "-;hydrodynamic friction torque;r");
%!     plot(param.omega * 30 / pi * SI_unit_second^-1, M_big_end_bearing(:, 2) * SI_unit_newton * SI_unit_meter, "-;asperity friction torque;b");
%!     plot(param.omega * 30 / pi * SI_unit_second^-1, M_big_end_bearing(:, 3) * SI_unit_newton * SI_unit_meter, "-;total friction torque;k");
%!     xlabel("n [rpm]");
%!     ylabel("M [Nm]");
%!     title(sprintf("frictional torque versus speed F1=%gN", param.F1(j) * SI_unit_newton));
%!     grid minor on;
%!   endfor
%!   for j=1:numel(param.F1)
%!     max_p = max_pc = zeros(numel(param.omega), 1);
%!     for i=1:numel(param.omega)
%!       max_p(i) = max(data(i, j).res.bearings(3).columns.p(:));
%!       max_pc(i) = max(data(i, j).res.bearings(3).columns.pc(:));
%!     endfor
%!     figure("visible", "off");
%!     hold on;
%!     plot(param.omega * 30 / pi * SI_unit_second^-1, 1e-6 * max_p * SI_unit_pascal, "-;max(p);r");
%!     plot(param.omega * 30 / pi * SI_unit_second^-1, 1e-6 * max_pc * SI_unit_pascal, "-;max(pc);b");
%!     xlabel("n [rpm]");
%!     ylabel("p [MPa]");
%!     grid minor on;
%!     title(sprintf("maximum pressure versus speed F1=%gN", param.F1(j) * SI_unit_newton));
%!   endfor
%!   figure("visible", "off");
%!   hold on;
%!   plot(omega_ref_h_10MPa * 30 / pi * SI_unit_second^-1, 1e6 * h_ref_10MPa * SI_unit_meter, "-;Sander New;k");
%!   plot(param.omega * 30 / pi * SI_unit_second^-1, 1e6 * min_h(:, 2) * SI_unit_meter, "-;MBDyn F1=8kN;r");
%!   xlabel("n [rpm]");
%!   ylabel("h [\\mum]");
%!   grid on;
%!   grid minor on;
%!   title("minimum film thickness versus speed");
%!   figure_list();
%!                                 %unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       #status = unlink(fullfile(fn(i).folder, fn(i).name));
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%!                                 %end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
