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
%!   SI_unit_meter = 1;
%!   SI_unit_kilogram = 1;
%!   SI_unit_second = 1;
%!   SI_unit_kelvin = 1;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   param.E = 210000e6 / SI_unit_pascal;
%!   param.nu = 0.3;
%!   param.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.E2 = 70000e6 / SI_unit_pascal;
%!   param.nu2 = 0.3;
%!   param.rho2 = 2700 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.E3 = 210000e6 / SI_unit_pascal;
%!   param.nu3 = 0.3;
%!   param.rho3 = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   param.d = 47.8e-3 / SI_unit_meter;
%!   param.w = 17.2e-3 / SI_unit_meter;
%!   param.ds = 54e-3 / SI_unit_meter;
%!   param.ws = 25e-3 / SI_unit_meter;
%!   param.h1 = 2e-3 / SI_unit_meter;
%!   param.l1 = param.w / 2 + param.ws / 2 + param.h1;
%!   param.l2 = param.w / 2 + param.ws + param.h1;
%!   param.l3 = 150e-3 / SI_unit_meter;
%!   param.dl = 50e-3 / SI_unit_meter;
%!   param.l4 = 50e-3 / SI_unit_meter;
%!   param.wr = 100e-3 / SI_unit_meter;
%!   param.d2 = 140e-3 / SI_unit_meter;
%!   param.l5 = 150e-3 / SI_unit_meter;
%!   param.l6 = 70e-3 / SI_unit_meter;
%!   param.l7 = 60e-3 / SI_unit_meter;
%!   param.l8 = 150e-3 / SI_unit_meter;
%!   param.l9 = 200e-3 / SI_unit_meter;
%!   param.l10 = 50e-3 / SI_unit_meter;
%!   param.l11 = 10e-3 / SI_unit_meter;
%!   param.l12 = 50e-3 / SI_unit_meter;
%!   param.l13 = 70e-3 / SI_unit_meter;
%!   param.l14 = 100e-3 / SI_unit_meter;
%!   param.l15 = 20e-3 / SI_unit_meter;
%!   param.l16 = 20e-3 / SI_unit_meter;
%!   param.pref = 100e6 / SI_unit_pascal;
%!   param.h = param.d * pi / 10;
%!   param.num_modes_cms = 10;
%!   param.num_modes_bearing = 10;
%!   opt_mesh.mesh.jacobian_range = [0.5, 1.5];
%!   opt_mesh.verbose = false;
%!   opt_mesh.output_file = [output_file, "_msh"];
%!   empty_cell = cell(1, 3);
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
%!   group_defs(2).name = "node_id_support_journal1";
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
%!   group_defs(2).compliance_matrix.mesh_size = param.h;
%!   group_defs(2).compliance_matrix.number_of_modes = param.num_modes_bearing;
%!   group_defs(2).compliance_matrix.include_rigid_body_modes = true;
%!   group_defs(2).bearing = "elem_id_support_bearing1";
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
%!   group_defs(3).compliance_matrix.mesh_size = param.h;
%!   group_defs(3).compliance_matrix.number_of_modes = param.num_modes_bearing;
%!   group_defs(3).compliance_matrix.include_rigid_body_modes = true;
%!   group_defs(3).bearing = "elem_id_big_end_bearing";
%!   group_defs(4).id = 4;
%!   group_defs(4).name = "node_id_support_journal2";
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
%!   group_defs(4).compliance_matrix.mesh_size = param.h;
%!   group_defs(4).compliance_matrix.number_of_modes = param.num_modes_bearing;
%!   group_defs(4).compliance_matrix.include_rigid_body_modes = true;
%!   group_defs(4).bearing = "elem_id_support_bearing2";
%!   empty_cell = cell(1, 2);
%!   group_defs2 = struct("id", empty_cell, ...
%!                        "name", empty_cell, ...
%!                        "R", empty_cell, ...
%!                        "X0", empty_cell, ...
%!                        "Xi", empty_cell, ...
%!                        "type", empty_cell, ...
%!                        "geometry", empty_cell, ...
%!                        "compliance_matrix", empty_cell);
%!   group_defs2(1).id = 1;
%!   group_defs2(1).name = "node_id_small_end";
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
%!   group_defs2(2).compliance_matrix.bearing_type = "journal";
%!   group_defs2(2).compliance_matrix.bearing_model = "EHD/FD";
%!   group_defs2(2).compliance_matrix.reference_pressure = param.pref;
%!   group_defs2(2).compliance_matrix.mesh_size = param.h;
%!   group_defs2(2).compliance_matrix.number_of_modes = param.num_modes_bearing;
%!   group_defs2(2).compliance_matrix.include_rigid_body_modes = true;
%!   group_defs2(2).bearing = "elem_id_big_end_bearing";
%!   empty_cell = cell(1, 2);
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
%!   group_defs3(2).name = "node_id_support_shell";
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
%!   group_defs3(2).compliance_matrix.bearing_type = "journal";
%!   group_defs3(2).compliance_matrix.bearing_model = "EHD/FD";
%!   group_defs3(2).compliance_matrix.reference_pressure = param.pref;
%!   group_defs3(2).compliance_matrix.mesh_size = param.h;
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
%!     fputs(fd, "Physical Surface(\"big_end_bearing\", 1) = {6};\n");
%!     fputs(fd, "Physical Surface(\"support_bearing1\", 2) = {3};\n");
%!     fputs(fd, "Physical Surface(\"support_bearing2\", 3) = {9};\n");
%!     fputs(fd, "Physical Surface(\"torque_input\", 4) = {2};\n");
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
%!   return
%!   options_mbdyn.output_file = [output_file, "_mbd"];
%!   options_mbdyn.f_run_mbdyn2easyanim = false;
%!   fd = -1;
%!   unwind_protect
%!     fd = fopen([output_file, "_mbd.mbdyn"], "w");
%!     if (fd == -1)
%!       error("failed to open file \"%s\"", [output_file, "_mbd.mbdyn"]);
%!     endif
%!     fputs(fd, "set: integer ref_id_assembly = 1001;\n");
%!     fputs(fd, "set: integer ref_id_shell_bearing = 1002;\n");
%!     fputs(fd, "set: integer ref_id_shell_support = 1003;\n");
%!     fputs(fd, "set: integer ref_id_journal_bearing = 1004;\n");
%!     fputs(fd, "set: integer ref_id_journal_support = 1005;\n");
%!     fputs(fd, "set: integer node_id_shell_support = 2001;\n");
%!     fputs(fd, "set: integer node_id_shell_bearing = 2002;\n");
%!     fputs(fd, "set: integer node_id_journal_bearing = 2004;\n");
%!     fputs(fd, "set: integer joint_id_shell_support = 3001;\n");
%!     fputs(fd, "set: integer joint_id_journal_support = 3002;\n");
%!     fputs(fd, "set: integer elem_id_diaphragm_cms = 3003;\n");
%!     fputs(fd, "set: integer elem_id_bearing = 3005;\n");
%!     fputs(fd, "set: integer number_of_nodes_x;\n");
%!     fputs(fd, "set: integer number_of_nodes_z;\n");
%!     fprintf(fd, "include: \"%s.set\";\n", output_file);
%!     fprintf(fd, "include: \"%s_shell.set\";\n", output_file);
%!     fputs(fd, "set: real omega1 = U1 / (0.5 * Di - cr);\n");
%!     fputs(fd, "set: real omega2 = U2 / (0.5 * Di);\n");
%!     fputs(fd, "set: real omega = omega1 - omega2;\n");
%!     fputs(fd, "begin: data;\n");
%!     fputs(fd, "        problem: initial value;\n");
%!     fputs(fd, "end: data;\n");
%!     fputs(fd, "begin: initial value;\n");
%!     fputs(fd, "        initial time: 0;\n");
%!     fputs(fd, "        final time: t1;\n");
%!     fputs(fd, "        time step: dt;\n");
%!     fputs(fd, "        method: implicit euler;\n");
%!     fputs(fd, "        linear solver: pardiso, grad, scale, row max column max, always;\n");
%!     fprintf(fd, "      threads: assembly, %d;\n", mbdyn_solver_num_threads_default());
%!     fprintf(fd, "      threads: solver, %d;\n", mbdyn_solver_num_threads_default());
%!     fputs(fd, "        nonlinear solver: mcp newton min fb;\n");
%!     fputs(fd, "        tolerance: 1e-4, test, minmax, 1e-4, test, norm;\n");
%!     fputs(fd, "        max iterations: 100;\n");
%!     fputs(fd, "        derivatives tolerance: 1e-3;\n");
%!     fputs(fd, "        derivatives max iterations: 5;\n");
%!     fputs(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!     fputs(fd, "        output: iterations, cpu time, solver condition number, stat, yes;\n");
%!     fputs(fd, "        threads: assembly, 1;\n");
%!     fputs(fd, "end: initial value;\n");
%!     fputs(fd, "begin: control data;\n");
%!     fputs(fd, "       use automatic differentiation;\n");
%!     fputs(fd, "       output precision: 16;\n");
%!     fputs(fd, "       max iterations: 0;\n");
%!     fputs(fd, "       print: dof stats, to file;\n");
%!     fputs(fd, "       print: dof description, to file;\n");
%!     fputs(fd, "       print: equation description, to file;\n");
%!     fputs(fd, "       structural nodes: 3;\n");
%!     fputs(fd, "       joints: 3;\n");
%!     fputs(fd, "       loadable elements: 1;\n");
%!     fputs(fd, "end: control data;\n");
%!     fputs(fd, "reference: ref_id_assembly,\n");
%!     fputs(fd, "        position, reference, global, null,\n");
%!     fputs(fd, "        orientation, reference, global, eye,\n");
%!     fputs(fd, "        velocity, reference, global, null,\n");
%!     fputs(fd, "        angular velocity, reference, global, null;\n");
%!     fputs(fd, "reference: ref_id_shell_bearing,\n");
%!     fputs(fd, "        position, reference, ref_id_assembly, null,\n");
%!     fputs(fd, "        orientation, reference, ref_id_assembly, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_assembly, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_assembly,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                omega2;\n");
%!     fputs(fd, "reference: ref_id_shell_support,\n");
%!     fputs(fd, "        position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "        orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_shell_bearing, null;\n");
%!     fputs(fd, "reference: ref_id_journal_bearing,\n");
%!     fputs(fd, "        position, reference, ref_id_assembly,\n");
%!     fputs(fd, "                  epsilon * cr * cos(delta),\n");
%!     fputs(fd, "                  epsilon * cr * sin(delta),\n");
%!     fputs(fd, "                  0.,\n");
%!     fputs(fd, "        orientation, reference, ref_id_assembly, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_assembly, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_assembly,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                0.,\n");
%!     fputs(fd, "                omega1;\n");
%!     fputs(fd, "reference: ref_id_journal_support,\n");
%!     fputs(fd, "        position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "        orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "        velocity, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "        angular velocity, reference, ref_id_journal_bearing, null;\n");
%!     fputs(fd, "begin: nodes;\n");
%!     fputs(fd, "        structural: node_id_shell_support, modal,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_support, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shell_support, eye,\n");
%!     fputs(fd, "                velocity, reference, ref_id_shell_support, null,\n");
%!     fputs(fd, "                angular velocity, reference, ref_id_shell_support, null;\n");
%!     fputs(fd, "        structural: node_id_shell_bearing, static,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                velocity, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                angular velocity, reference, ref_id_shell_bearing, null;\n");
%!     fputs(fd, "        structural: node_id_journal_bearing, static,\n");
%!     fputs(fd, "                position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                velocity, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                angular velocity, reference, ref_id_journal_bearing, null;\n");
%!     fputs(fd, "end: nodes;\n");
%!     fputs(fd, "begin: elements;\n");
%!     fputs(fd, "        joint: joint_id_shell_support, total pin joint,\n");
%!     fputs(fd, "                node_id_shell_support,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                position constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        component,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                orientation constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        component,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                mult, time, omega2;\n");
%!     fputs(fd, "        joint: joint_id_journal_support, total pin joint,\n");
%!     fputs(fd, "                node_id_journal_bearing,\n");
%!     fputs(fd, "                position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                position, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                position orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                rotation orientation, reference, ref_id_journal_bearing, eye,\n");
%!     fputs(fd, "                position constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        null,\n");
%!     fputs(fd, "                orientation constraint,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        active,\n");
%!     fputs(fd, "                        component,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                null,\n");
%!     fputs(fd, "                                mult, time, omega1;\n");
%!     fprintf(fd, "        include: \"%s_cms.elm\";\n", output_file);
%!     fputs(fd, "        user defined: elem_id_bearing, hydrodynamic plain bearing2,\n");
%!     fputs(fd, "            hydraulic fluid, linear compressible,\n");
%!     fputs(fd, "                density, rhol,\n");
%!     fputs(fd, "                betal,\n");
%!     fputs(fd, "                pressure, pcl,\n");
%!     fputs(fd, "                viscosity, etal,\n");
%!     fputs(fd, "                temperature, Tl,\n");
%!     fputs(fd, "            viscosity vapor, factor, fact_etav,\n");
%!     fputs(fd, "                mesh, linear finite difference,\n");
%!     fputs(fd, "                enable mcp, yes,\n");
%!     fputs(fd, "                geometry, cylindrical,\n");
%!     fputs(fd, "                        mesh position, at bearing,\n");
%!     fputs(fd, "                        bearing width, Wo,\n");
%!     fputs(fd, "                        shaft diameter, Di - 2 * cr,\n");
%!     fputs(fd, "                        bearing diameter, Di,\n");
%!     fputs(fd, "                shaft node, node_id_journal_bearing,\n");
%!     fputs(fd, "                offset, reference, ref_id_journal_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_journal_bearing, eye,\n");
%!     switch (group_defs(2).compliance_matrix.matrix_type)
%!       case "modal substruct total"
%!         fputs(fd, "                bearing node, node_id_shell_support,\n");
%!       otherwise
%!         fputs(fd, "                bearing node, node_id_shell_bearing,\n");
%!     endswitch
%!     fputs(fd, "                offset, reference, ref_id_shell_bearing, null,\n");
%!     fputs(fd, "                orientation, reference, ref_id_shell_bearing, eye,\n");
%!     fputs(fd, "                        number of nodes z, number_of_nodes_z,\n");
%!     fputs(fd, "                        number of nodes Phi, number_of_nodes_x,\n");
%!     fputs(fd, "                        boundary conditions,\n");
%!     fputs(fd, "                                        pressure, pside,\n");
%!     fputs(fd, "                                        pressure, pside,\n");
%!     fputs(fd, "                        lubrication grooves, 1,\n");
%!     fputs(fd, "                                at bearing,\n");
%!     fputs(fd, "                                        pressure, pin,\n");
%!     fputs(fd, "                                        position, 0., 0.,\n");
%!     fputs(fd, "                                        rectangle, width, d7, height, Wo,\n");
%!     fputs(fd, "                compliance model,\n");
%!     fprintf(fd, "                        matrix, 1, from file, \"%s_elem_id_bearing.dat\",\n", output_file);
%!     fputs(fd, "                        E1, E,\n");
%!     fputs(fd, "                        nu1, nu,\n");
%!     switch (group_defs(2).compliance_matrix.matrix_type)
%!       case "modal substruct total"
%!         fputs(fd, "                        modal element, elem_id_diaphragm_cms,\n");
%!     endswitch
%!     fputs(fd, "                pressure dof scale, pref,\n");
%!     fputs(fd, "                reynolds equation scale, dt / (rhol * cr),\n");
%!     fputs(fd, "                elasticity equation scale, dt / cr,\n");
%!     fputs(fd, "                output pressure, yes,\n");
%!     fputs(fd, "                output stress, yes,\n");
%!     fputs(fd, "                output density, yes,\n");
%!     fputs(fd, "                output friction loss, yes,\n");
%!     fputs(fd, "                output clearance, yes,\n");
%!     fputs(fd, "                output reaction force, yes,\n");
%!     fputs(fd, "                output mesh, yes,\n");
%!     fputs(fd, "                output total deformation, yes,\n");
%!     fputs(fd, "                output, yes;\n");
%!     fputs(fd, "end: elements;\n");
%!   unwind_protect_cleanup
%!     if (fd ~= -1)
%!       fclose(fd);
%!     endif
%!   end_unwind_protect
%!   mbdyn_solver_run([output_file, "_mbd.mbdyn"], options_mbdyn);
%!   res.log_dat = mbdyn_post_load_log(options_mbdyn.output_file);
%!   [res.t, ...
%!    res.trajectory, ...
%!    res.deformation, ...
%!    res.velocity, ...
%!    res.acceleration, ...
%!    res.node_id, ...
%!    res.force, ...
%!    res.force_id, ...
%!    res.force_node_id, ...
%!    res.orientation_description] = mbdyn_post_load_output_struct(options_mbdyn.output_file);
%!   [res.elem_id, res.q, res.qdot, res.qddot] = mbdyn_post_load_output_mod(options_mbdyn.output_file, numel(res.t));
%!   res.bearings = mbdyn_post_ehd_load_output(options_mbdyn.output_file, res.log_dat);
%!   res.log_dat.vars = mbdyn_post_id_to_index(res, res.log_dat.vars);
%!   res.elem_idx_cms = find(res.elem_id == res.log_dat.vars.elem_id_diaphragm_cms);
%!   res.elem_idx_bearing = find(res.log_dat.vars.elem_id_bearing == [res.log_dat.bearings.label]);
%!   cms_data.mesh = mesh;
%!   cms_data.dof_map = dof_map;
%!   cms_data.cms_opt = cms_opt;
%!   res.sol_dyn = fem_post_cms_sol_import(options_mbdyn.output_file, cms_data);
%!   Phi = res.bearings.xi(1, :) / (0.5 * param.Di);
%!   p = interp2(res.bearings.xi, res.bearings.zi, res.bearings.columns.p(:, :, end), res.bearings.xi(1, :), 0);
%!   w = interp2(res.bearings.xi, res.bearings.zi, res.bearings.columns.wtot(:, :, end), res.bearings.xi(1, :), 0);
%!   h = interp2(res.bearings.xi, res.bearings.zi, res.bearings.columns.h(:, :, end), res.bearings.xi(1, :), 0);
%!   if (f_plot)
%!     figure("visible", "off");
%!     ax = plotyy(180 / pi * Phi, 1e-3 * p * SI_unit_pascal, 180 / pi * Phi, 1e6 * w * SI_unit_meter);
%!     figure("visible", "off");
%!     plot(res.t * SI_unit_second, 1e-3 * max(res.bearings.columns.p_n, [], 2) * SI_unit_pascal, "-;max(p(t));r");
%!     xlabel("t [s]");
%!     ylabel("p [kPa]");
%!     grid on;
%!     grid minor on;
%!     title("convergence history of peak pressure");
%!     figure("visible", "off");
%!     plot(res.t * SI_unit_second, 1e6 * max(res.bearings.columns.wtot_n, [], 2) * SI_unit_pascal, "-;max(w(t));r");
%!     xlabel("t [s]");
%!     ylabel("w [um]");
%!     grid on;
%!     grid minor on;
%!     title("convergence history of peak deformation");
%!     figure("visible", "off");
%!     xlabel("Phi [deg]");
%!     ylabel(ax(1), "p [kPa]");
%!     ylabel(ax(2), "w [um]");
%!     grid on;
%!     grid minor on;
%!     for i=1:2
%!       xlim(ax(i), [0, 360]);
%!     endfor
%!     xticks(0:30:360);
%!     title("midplane pressure and deformation");
%!     figure("visible","off");
%!     hold on;
%!     set(plot(180 / pi * Phi, 1e6 * w * SI_unit_meter, "-;w(Phi) [um];r"), "linewidth", 5);
%!     set(plot(ref_data_w(:, 1), ref_data_w(:, 2), "-;reference w(Phi) [um];k"), "linewidth", 3);
%!     xlabel("Phi [deg]");
%!     ylabel("w [um]");
%!     grid on;
%!     grid minor on;
%!     xlim([0,360]);
%!     xticks(0:30:360);
%!     title("midplane deformation");
%!     figure("visible","off");
%!     hold on;
%!     set(plot(180 / pi * Phi, 1e-3 * p * SI_unit_pascal, "-;p(Phi) [kPa];r"), "linewidth", 5);
%!     set(plot(ref_p(:, 1), ref_p(:, 2), "-;reference p(Phi) [kPa];k"), "linewidth", 3);
%!     xlabel("Phi [deg]");
%!     ylabel("p* [kPa]");
%!     grid on;
%!     grid minor on;
%!     xlim([0,360]);
%!     xticks(0:30:360);
%!     title("midplane pressure");
%!     figure("visible","off");
%!     set(plot(180 / pi * Phi, p * param.cr^2 / (0.5 * param.Di * param.U1 * param.etal), "-;p*(Phi) [1];r"), "linewidth", 5);
%!     xlabel("Phi [deg]");
%!     ylabel("p* [1]");
%!     grid on;
%!     grid minor on;
%!     xlim([0,360]);
%!     xticks(0:30:360);
%!     title("midplane nondimensional pressure");
%!     figure("visible", "off");
%!     hold on;
%!     set(plot(180 / pi * Phi, 1e6 * h * SI_unit_meter, "-;h(Phi) [um];r"), "linewidth", 5);
%!     set(plot(180 / pi * Phi, 1e6 * w * SI_unit_meter, "-;w(Phi) [um];b"), "linewidth", 5);
%!     xlabel("Phi [deg]");
%!     ylabel(ax(2), "h, w [um]");
%!     grid on;
%!     grid minor on;
%!     xlim([0, 360]);
%!     title("midplane clearance and deformation");
%!     figure("visible","off");
%!     contourf(180/pi * res.bearings.xi / (0.5 * param.Di), 1e3 * res.bearings.zi * SI_unit_meter, 1e-3 * res.bearings.columns.p(:, :, end) * SI_unit_pascal);
%!     colormap jet;
%!     colorbar;
%!     xlabel("Phi [deg]");
%!     ylabel("z [mm]");
%!     grid on;
%!     grid minor on;
%!     title("pressure distribution p [kPa]");
%!     figure("visible","off");
%!     contourf(180/pi * res.bearings.xi / (0.5 * param.Di), 1e3 * res.bearings.zi * SI_unit_meter, 1e6 * res.bearings.columns.wtot(:, :, end) * SI_unit_meter);
%!     colormap jet;
%!     colorbar;
%!     xlabel("Phi [deg]");
%!     ylabel("z [mm]");
%!     grid on;
%!     grid minor on;
%!     title("radial deformation [um]");
%!     figure("visible","off");
%!     contourf(180/pi * res.bearings.xi / (0.5 * param.Di), 1e3 * res.bearings.zi * SI_unit_meter, res.bearings.columns.rho(:, :, end) / param.rhol);
%!     colormap jet;
%!     colorbar;
%!     xlabel("Phi [deg]");
%!     ylabel("z [mm]");
%!     grid on;
%!     grid minor on;
%!     title("fractional film content [1]");
%!     figure_list();
%!   endif
%!   p_int = interp1(180 / pi * Phi, 1e-3 * p * SI_unit_pascal, ref_p(:, 1), "linear");
%!   w_int = interp1(180 / pi * Phi, 1e6 * w * SI_unit_meter, ref_data_w(:, 1), "linear");
%!   assert_simple(p_int, ref_p(:, 2), 0.15 * max(abs(ref_p(:, 2))));
%!   assert_simple(mean(abs(p_int - ref_p(:, 2))) < 0.05 * max(abs(ref_p(:, 2))));
%!   ## Don't check the first point because of issues related to the interpolation of the compliance matrix near to the slot
%!   assert_simple(w_int(2:end), ref_data_w(2:end, 2), 0.07 * max(abs(ref_data_w(:, 2))));
%!   assert_simple(mean(abs(w_int(2:end) - ref_data_w(2:end, 2))) < 0.04 * max(abs(ref_data_w(:, 2))));
%! #unwind_protect_cleanup
%!   if (numel(output_file))
%!     fn = dir([output_file, "*"]);
%!     for i=1:numel(fn)
%!       status = unlink(fullfile(fn(i).folder, fn(i).name));
%!       if (status ~= 0)
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! #end_unwind_protect
%! catch
%!   gtest_error = lasterror();
%!   gtest_fail(gtest_error, evalin("caller", "__file"));
%!   rethrow(gtest_error);
%! end_try_catch
