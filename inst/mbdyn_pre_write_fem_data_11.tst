## mbdyn_pre_write_fem_data.tst:11
%!test
%! ## TEST11
%! pkg load mboct-fem-pkg;
%! filename = "";
%! unwind_protect
%!   filename = tempname();
%!   if (ispc())
%!     filename(filename == "\\") = "/";
%!   endif
%!   SI_unit_meter = 1e-3;
%!   SI_unit_second = 1;
%!   SI_unit_kilogram = 1e3;
%!   SI_unit_newton = SI_unit_kilogram * SI_unit_meter / SI_unit_second^2;
%!   SI_unit_pascal = SI_unit_newton / SI_unit_meter^2;
%!   SI_unit_rad = 1;
%!   geometry.l = 4e-3 / SI_unit_meter;
%!   geometry.w = 6e-3 / SI_unit_meter;
%!   geometry.h = 8e-3 / SI_unit_meter;
%!   h = [geometry.l; geometry.w; geometry.h];
%!   t1 = 1;
%!   dt = t1 / 40;
%!   model = "static";
%!   method = "implicit euler";
%!   options.verbose = false;
%!   material.E = 210000e6 / SI_unit_pascal;
%!   material.ET = 2100e6 / SI_unit_pascal;
%!   material.sigmayv = 235000e6 / SI_unit_pascal;
%!   material.nu = 0.3;
%!   material.G = material.E / (2 * (1 + material.nu));
%!   material.kappa = material.E / (3 * (1 - 2 * material.nu));
%!   material.rho = 7850 / (SI_unit_kilogram / SI_unit_meter^3);
%!   material.beta = 1 / SI_unit_second^-1;
%!   material.delta = 0.25;
%!   material.theta = 0.5;
%!   material.tau = 0.5 / SI_unit_second;
%!   elem_types = {"tet10h", ...
%!                 "tet10upc", ...
%!                 "iso8", ...
%!                 "iso8upc", ...
%!                 "iso20", ...
%!                 "iso20r", ...
%!                 "iso20upc", ...
%!                 "iso20upcr", ...
%!                 "iso27", ...
%!                 "penta15", ...
%!                 "penta15upc", ...
%!                };
%!   mat_types = {"linear elastic generic", ...
%!                "hookean linear elastic isotropic", ...
%!                "neo hookean elastic", ...
%!                "bilinear isotropic hardening", ...
%!                "mooney rivlin elastic", ...
%!                "linear viscoelastic generic", ...
%!                "hookean linear viscoelastic isotropic", ...
%!                "neo hookean viscoelastic", ...
%!                "linear viscoelastic maxwell1", ...
%!                "linear viscoelastic maxwelln", ...
%!               };
%!   f_transfinite_mesh = [true, ...
%!                         false, ...
%!                        ];
%!   boundary_cond = { #"symmetry", ...
%!                     "three point", ...
%!                     #"two surfaces one line", ...
%!                   };
%!   load_type = {"traction", "pressure", "prestrain"};
%!   sigma = material.sigmayv * diag([1.5, 0.3, 0.9, 0.2, 0.2, 0.2]);
%!   epsilon0 = [1; 2; 3; 0.4; 0.5; 0.6];
%!   for idx_sigma=1:columns(sigma)
%!     sigmav = sqrt(sum(sigma(1:3, idx_sigma).^2) - (sigma(1, idx_sigma) * sigma(2, idx_sigma) + sigma(2, idx_sigma) * sigma(3, idx_sigma) + sigma(1, idx_sigma) * sigma(3, idx_sigma)) + 3 * sum(sigma(4:6, idx_sigma).^2));
%!     for idx_load_type=1:numel(load_type)
%!       switch (load_type{idx_load_type})
%!         case "pressure"
%!           if (norm(sigma(1:3, idx_sigma)) == 0)
%!             continue; ## Uniform pressure cannot cause shear stress
%!           endif
%!         case "prestrain"
%!           if (idx_sigma > columns(epsilon0))
%!             continue;
%!           endif
%!       endswitch
%!       for idx_boundary_cond=1:numel(boundary_cond)
%!         switch (boundary_cond{idx_boundary_cond})
%!           case "symmetry"
%!             switch (load_type{idx_load_type})
%!             case "prestrain"
%!               sym_cond = norm(epsilon0(4:6, idx_sigma));
%!             otherwise
%!               sym_cond = norm(sigma(4:6, idx_sigma));
%!             endswitch
%!             if (sym_cond)
%!               continue; ## Deformation cannot be symmetric because of shear stress
%!             endif
%!         endswitch
%!         for idx_transfinite=1:numel(f_transfinite_mesh)
%!           for idx_mat_type=1:numel(mat_types)
%!             material.type = mat_types{idx_mat_type};
%!             for idx_elem_type=1:numel(elem_types)
%!               t2 = 0;
%!               elem_type = elem_types{idx_elem_type};
%!               switch (elem_type)
%!               case {"iso8upc", "iso20upc", "iso20upcr", "penta15upc", "tet10upc"}
%!                 switch (material.type)
%!                 case {"hookean linear elastic isotropic", "mooney rivlin elastic", "bilinear isotropic hardening"}
%!                 otherwise
%!                   ## incompressible version of constitutive law not implemented yet
%!                   continue;
%!                 endswitch
%!               endswitch
%!               switch (load_type{idx_load_type})
%!               case "prestrain"
%!                 switch (mat_types{idx_mat_type})
%!                 case {"linear elastic generic", "hookean linear elastic isotropic", "neo hookean elastic", "mooney rivlin elastic", "linear viscoelastic generic", "hookean linear viscoelastic isotropic", "neo hookean viscoelastic", "bilinear isotropic hardening"}
%!                 otherwise
%!                   ## prestrain is not implemented yet
%!                   continue;
%!                 endswitch
%!                 switch (mat_types{idx_mat_type})
%!                 case {"linear viscoelastic generic", "hookean linear viscoelastic isotropic", "neo hookean viscoelastic"}
%!                   t2 = 1000;
%!                 endswitch
%!                 switch (mat_types{idx_mat_type})
%!                 case "neo hookean viscoelastic"
%!                   continue; ## FIXME: test not passed yet
%!                 endswitch
%!               endswitch
%!               switch (elem_type)
%!                 case {"iso20upcr", "iso20r"}
%!                   if (f_transfinite_mesh(idx_transfinite))
%!                     elem_factor_h = [0.5; 2; 2]; ## avoid hourglass instability
%!                   else
%!                     elem_factor_h = [0.5; 0.5; 0.5];
%!                   endif
%!                 otherwise
%!                   elem_factor_h = [2; 2; 2];
%!               endswitch
%!               switch (mat_types{idx_mat_type})
%!                 case {"neo hookean elastic", "neo hookean viscoelastic", "mooney rivlin elastic", "linear viscoelastic generic", "hookean linear viscoelastic isotropic"}
%!                   switch(boundary_cond{idx_boundary_cond})
%!                   case {"three point", "two surfaces one line"}
%!                     switch (idx_sigma)
%!                     case {4, 5, 6}
%!                       ## shear deformation with those materials and elements not passed yet because the Jacobian may become singular
%!                       switch (elem_type)
%!                       case {"tet10h", "tet10upc", "penta15", "penta15upc"}
%!                         continue;
%!                       otherwise
%!                         if (~f_transfinite_mesh(idx_transfinite))
%!                           continue;
%!                         endif
%!                       endswitch
%!                     endswitch
%!                   endswitch
%!               endswitch
%!               file_prefix = sprintf("%s_%d_%d_%d_%d_%d_%d", filename, idx_sigma, idx_load_type, idx_boundary_cond, idx_transfinite, idx_mat_type, idx_elem_type);
%!               geo_file = [file_prefix, "_gmsh.geo"];
%!               mesh_file = [file_prefix, "_gmsh.msh"];
%!               nodes_file = [file_prefix, "_mbd.nod"];
%!               elem_file = [file_prefix, "_mbd.elm"];
%!               set_file = [file_prefix, "_mbd.set"];
%!               csl_file = [file_prefix, "_mbd.csl"];
%!               control_file = [file_prefix, "_mbd.con"];
%!               initial_value_file = [file_prefix, "_mbd.inv"];
%!               input_file = [file_prefix, "_mbd_inp.mbdyn"];
%!               output_file = [file_prefix, "_mbd_out"];
%!               opt_mbd.output_file = output_file;
%!               if (~options.verbose)
%!                 opt_mbd.logfile = [opt_mbd.output_file, ".stdout"];
%!               endif
%!               opt_mbd.mbdyn_command = "mbdyn -C";
%!               opt_mbd.f_run_mbdyn = true;
%!               switch (elem_type)
%!                 case {"iso8", "iso8upc"}
%!                   mesh_order = 1;
%!                   elem_type_solid = {elem_type};
%!                   elem_type_surf = {"iso4"};
%!                 case "iso20"
%!                   mesh_order = 2;
%!                   elem_type_solid = {"iso20", "penta15"};
%!                   elem_type_surf = {"quad8", "tria6h"};
%!                 case "iso20upc"
%!                   mesh_order = 2;
%!                   elem_type_solid = {elem_type};
%!                   elem_type_surf = {"quad8"};
%!                 case "iso20upcr"
%!                   mesh_order = 2;
%!                   elem_type_solid = {elem_type};
%!                   elem_type_surf = {"quad8r"};
%!                 case "iso27"
%!                   mesh_order = 2;
%!                   elem_type_solid = {"iso27"};
%!                   elem_type_surf = {"quad9"};
%!                 case "iso20r"
%!                   mesh_order = 2;
%!                   elem_type_solid = {"iso20r", "penta15"};
%!                   elem_type_surf = {"quad8r", "tria6h"};
%!                 case {"penta15", "penta15upc"}
%!                   mesh_order = 2;
%!                   elem_type_solid = {elem_type};
%!                   elem_type_surf = {"quad8", "tria6h"};
%!                 case {"tet10h", "tet10upc"}
%!                   mesh_order = 2;
%!                   elem_type_solid = {elem_type};
%!                   elem_type_surf = {"tria6h"};
%!                 otherwise
%!                   error("unknown element type \"%s\"", elem_type);
%!               endswitch
%!               fd = -1;
%!               unwind_protect
%!                 [fd, msg] = fopen(geo_file, "w");
%!                 if (fd == -1)
%!                   error("failed to open file \"%s.geo\"", geo_file);
%!                 endif
%!                 fprintf(fd, "SetFactory(\"OpenCASCADE\");\n");
%!                 fprintf(fd, "a=%g;\n", geometry.l);
%!                 fprintf(fd, "b=%g;\n", geometry.w);
%!                 fprintf(fd, "c=%g;\n", geometry.h);
%!                 for i=1:3
%!                   fprintf(fd, "h%s = %g;\n", {"x","y","z"}{i}, h(i) * elem_factor_h(i));
%!                 endfor
%!                 switch (elem_type)
%!                   case {"iso20", "iso20upc", "iso20upcr", "iso20r", "penta15", "penta15upc"}
%!                     fputs(fd, "Mesh.SecondOrderIncomplete=1;\n");
%!                   otherwise
%!                     fputs(fd, "Mesh.SecondOrderIncomplete=0;\n");
%!                 endswitch
%!                 fprintf(fd, "Mesh.ElementOrder = %d;\n", mesh_order);
%!                 fputs(fd, "Point(1) = { 0, -0.5 * b, -0.5 * c};\n");
%!                 fputs(fd, "Point(2) = { a, -0.5 * b, -0.5 * c};\n");
%!                 fputs(fd, "Point(3) = { a,  0.5 * b, -0.5 * c};\n");
%!                 fputs(fd, "Point(4) = { 0,  0.5 * b, -0.5 * c};\n");
%!                 fputs(fd, "Line(1) = {4,3};\n");
%!                 fputs(fd, "Line(2) = {3,2};\n");
%!                 fputs(fd, "Line(3) = {2,1};\n");
%!                 fputs(fd, "Line(4) = {1,4};\n");
%!                 switch (elem_type)
%!                 case {"iso20upcr", "iso20r"}
%!                   num_layers = 2; ## because of hourglass instability
%!                 otherwise
%!                   num_layers = 1;
%!                 endswitch
%!                 fprintf(fd, "num_layers = %d;\n", num_layers);
%!                 if (f_transfinite_mesh(idx_transfinite))
%!                   fprintf(fd, "Transfinite Curve(1) = Max(num_layers + 1, Round(a / hx) + 1);\n");
%!                   fprintf(fd, "Transfinite Curve(2) = Max(num_layers + 1, Round(b / hy) + 1);\n");
%!                   fprintf(fd, "Transfinite Curve(3) = Max(num_layers + 1, Round(a / hx) + 1);\n");
%!                   fprintf(fd, "Transfinite Curve(4) = Max(num_layers + 1, Round(b / hy) + 1);\n");
%!                 endif
%!                 fputs(fd, "Line Loop(5) = {2,3,4,1};\n");
%!                 fputs(fd, "Plane Surface(6) = {5};\n");
%!                 if (f_transfinite_mesh(idx_transfinite))
%!                   fputs(fd, "Transfinite Surface(6) = {2,3,4,1};\n");
%!                 endif
%!                 fputs(fd, "tmp[] = Extrude {0,0.0,c} {\n");
%!                 switch (elem_type)
%!                   case {"tet10h", "tet10upc"}
%!                     fputs(fd, "  Surface{6};\n");
%!                   otherwise
%!                     fprintf(fd, "  Surface{6}; Layers{Max(num_layers, Round(c/hz))}; Recombine;\n");
%!                 endswitch
%!                 fputs(fd, "};\n");
%!                 f_unstruct_mesh_size = false;
%!                 switch (elem_type)
%!                   case {"iso8", "iso8upc", "iso20", "iso20r", "iso20upc", "iso20upcr", "iso27"}
%!                     f_unstruct_mesh_size = ~f_transfinite_mesh(idx_transfinite);
%!                     fputs(fd, "Recombine Surface{6, tmp[0]};\n");
%!                   otherwise
%!                     f_unstruct_mesh_size = true;
%!                 endswitch
%!                 if (f_unstruct_mesh_size)
%!                   fprintf(fd, "MeshSize{PointsOf{Volume{tmp[1]};}} = %.16e;\n", 2 * mean(h .* elem_factor_h));
%!                 endif
%!                 switch (elem_type)
%!                   case {"tet10h", "tet10upc"}
%!                     if (~f_transfinite_mesh(idx_transfinite))
%!                       fputs(fd, "Mesh.HighOrderOptimize=2;\n");
%!                       fputs(fd, "Mesh.OptimizeThreshold=0.99;\n");
%!                     endif
%!                 endswitch
%!                 fputs(fd, "Physical Volume(\"volume\",1) = {tmp[1]};\n");
%!                 fputs(fd, "Physical Surface(\"clamp\",2) = {tmp[4]};\n");
%!                 fputs(fd, "Physical Surface(\"load+x\",3) = {tmp[2]};\n");
%!                 fputs(fd, "Physical Surface(\"load+y\",6) = {tmp[5]};\n");
%!                 fputs(fd, "Physical Surface(\"load+z\",7) = {tmp[0]};\n");
%!                 fputs(fd, "Physical Surface(\"load-x\",8) = {tmp[4]};\n");
%!                 fputs(fd, "Physical Surface(\"load-y\",9) = {tmp[3]};\n");
%!                 fputs(fd, "Physical Surface(\"load-z\",10) = {6};\n");
%!                 fputs(fd, "Physical Surface(\"symmetry-xy\",4) = {6};\n");
%!                 fputs(fd, "Physical Surface(\"symmetry-xz\",5) = {tmp[3]};\n");
%!               unwind_protect_cleanup
%!                 if (fd ~= -1)
%!                   fclose(fd);
%!                 endif
%!               end_unwind_protect
%!               pid = spawn("gmsh", {"-format", "msh2", "-3", geo_file});
%!               status = spawn_wait(pid);
%!               if (status ~= 0)
%!                 warning("gmsh failed with status %d", status);
%!               endif
%!               opt_msh.elem_type = {elem_type_solid{:}, elem_type_surf{:}};
%!               mesh = fem_pre_mesh_reorder(fem_pre_mesh_import(mesh_file, "gmsh", opt_msh));
%!               opt_mbd_mesh = struct();
%!               switch (model)
%!                 case "dynamic"
%!                   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_DYNAMIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!                 case "static"
%!                   opt_mbd_mesh.struct_nodes.type = repmat(MBDYN_NODE_TYPE_STATIC_STRUCT_DISP, rows(mesh.nodes), 1);
%!               endswitch
%!               grp_idx_volume = zeros(1, numel(elem_type_solid), "int32");
%!               for i=1:numel(elem_type_solid)
%!                 if (~isfield(mesh.groups, elem_type_solid{i}))
%!                   continue;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_solid{i}).id] == 1);
%!                 if (isempty(idx))
%!                   continue;
%!                 endif
%!                 grp_idx_volume(i) = idx;
%!               endfor
%!               grp_idx_load_px = grp_idx_load_py = grp_idx_load_pz = grp_idx_clamp = grp_idx_symmetry_xy = grp_idx_symmetry_xz = zeros(1, numel(elem_type_surf), "int32");
%!               grp_idx_load_mx = grp_idx_load_my = grp_idx_load_mz = zeros(1, numel(elem_type_surf), "int32");
%!               for i=1:numel(elem_type_surf)
%!                 if (~isfield(mesh.groups, elem_type_surf{i}))
%!                   continue;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 2);
%!                 if (~isempty(idx))
%!                   grp_idx_clamp(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 3);
%!                 if (~isempty(idx))
%!                   grp_idx_load_px(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 4);
%!                 if (~isempty(idx))
%!                   grp_idx_symmetry_xy(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 5);
%!                 if (~isempty(idx))
%!                   grp_idx_symmetry_xz(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 6);
%!                 if (~isempty(idx))
%!                   grp_idx_load_py(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 7);
%!                 if (~isempty(idx))
%!                   grp_idx_load_pz(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 8);
%!                 if (~isempty(idx))
%!                   grp_idx_load_mx(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 9);
%!                 if (~isempty(idx))
%!                   grp_idx_load_my(i) = idx;
%!                 endif
%!                 idx = find([getfield(mesh.groups, elem_type_surf{i}).id] == 10);
%!                 if (~isempty(idx))
%!                   grp_idx_load_mz(i) = idx;
%!                 endif
%!               endfor
%!               load_case_dof.locked_dof = false(rows(mesh.nodes), 6);
%!               grp_idx_p1 = find((mesh.nodes(:, 1) == 0) & (mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 3) == -0.5 * geometry.h));
%!               grp_idx_p2 = find((mesh.nodes(:, 1) == geometry.l) & (mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 3) == -0.5 * geometry.h));
%!               grp_idx_p3 = find((mesh.nodes(:, 1) == geometry.l) & (mesh.nodes(:, 2) == 0.5 * geometry.w) & (mesh.nodes(:, 3) == -0.5 * geometry.h));
%!               grp_idx_p4 = find((mesh.nodes(:, 1) == 0) & (mesh.nodes(:, 2) == 0.5 * geometry.w) & (mesh.nodes(:, 3) == -0.5 * geometry.h));
%!               grp_idx_p5 = find((mesh.nodes(:, 1) == 0) & (mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 3) == 0.5 * geometry.h));
%!               grp_idx_p6 = find((mesh.nodes(:, 1) == geometry.l) & (mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 3) == 0.5 * geometry.h));
%!               grp_idx_p7 = find((mesh.nodes(:, 1) == geometry.l) & (mesh.nodes(:, 2) == 0.5 * geometry.w) & (mesh.nodes(:, 3) == 0.5 * geometry.h));
%!               grp_idx_p8 = find((mesh.nodes(:, 1) == 0) & (mesh.nodes(:, 2) == 0.5 * geometry.w) & (mesh.nodes(:, 3) == 0.5 * geometry.h));
%!               switch (boundary_cond{idx_boundary_cond})
%!                 case "symmetry"
%!                   for i=1:numel(elem_type_surf)
%!                     if (grp_idx_clamp(i))
%!                       load_case_dof.locked_dof(getfield(mesh.groups, elem_type_surf{i})(grp_idx_clamp(i)).nodes, 1) = true;
%!                     endif
%!                     if (grp_idx_symmetry_xz(i))
%!                       load_case_dof.locked_dof(getfield(mesh.groups, elem_type_surf{i})(grp_idx_symmetry_xz(i)).nodes, 2) = true;
%!                     endif
%!                     if (grp_idx_symmetry_xy(i))
%!                       load_case_dof.locked_dof(getfield(mesh.groups, elem_type_surf{i})(grp_idx_symmetry_xy(i)).nodes, 3) = true;
%!                     endif
%!                   endfor
%!                 case "three point"
%!                   load_case_dof.locked_dof(grp_idx_p1, 1:3) = true;
%!                   load_case_dof.locked_dof(grp_idx_p2, 2:3) = true;
%!                   load_case_dof.locked_dof(grp_idx_p4, 3) = true;
%!                 case "two surfaces one line"
%!                     switch (idx_sigma)
%!                     case {1, 2, 3, 4}
%!                       load_case_dof.locked_dof(mesh.nodes(:, 2) == -0.5 * geometry.w, 2) = true;
%!                       load_case_dof.locked_dof(mesh.nodes(:, 3) == -0.5 * geometry.h, 3) = true;
%!                       load_case_dof.locked_dof((mesh.nodes(:, 2) == -0.5 * geometry.w) & (mesh.nodes(:, 1) == 0), 1) = true;
%!                     case 5
%!                       load_case_dof.locked_dof(mesh.nodes(:, 3) == -0.5 * geometry.h, 3) = true;
%!                       load_case_dof.locked_dof(mesh.nodes(:, 1) == 0, 1) = true;
%!                       load_case_dof.locked_dof((mesh.nodes(:, 3) == -0.5 * geometry.h) & (mesh.nodes(:, 2) == -0.5 * geometry.w), 2) = true;
%!                     case 6
%!                       load_case_dof.locked_dof(mesh.nodes(:, 3) == -0.5 * geometry.h, 3) = true;
%!                       load_case_dof.locked_dof(mesh.nodes(:, 2) == -0.5 * geometry.w, 2) = true;
%!                       load_case_dof.locked_dof((mesh.nodes(:, 3) == -0.5 * geometry.h) & (mesh.nodes(:, 1) == 0), 1) = true;
%!                     otherwise
%!                       error("invalid boundary condition");
%!                     endswitch
%!                 otherwise
%!                   error("unkown boundary condition");
%!               endswitch
%!               mesh.material_data = material;
%!               switch (mesh.material_data.type)
%!               case "linear viscoelastic maxwelln"
%!                 mesh.material_data.tau = repmat(mesh.material_data.tau, 1, 2);
%!                 mesh.material_data.theta = repmat(mesh.material_data.theta / 2, 1, 2);
%!               endswitch
%!               mesh.materials = struct();
%!               for i=1:numel(elem_type_solid)
%!                 if (~isfield(mesh.elements, elem_type_solid{i}))
%!                   continue;
%!                 endif
%!                 elem_mat = zeros(rows(getfield(mesh.elements, elem_type_solid{i})), 1, "int32");
%!                 elem_mat(getfield(mesh.groups, elem_type_solid{i})(grp_idx_volume(i)).elements) = 1;
%!                 mesh.materials = setfield(mesh.materials, elem_type_solid{i}, elem_mat);
%!               endfor
%!               opt_mbd_mesh.forces.time_function = "time";
%!               opt_mbd_mesh.surface_loads.time_function = opt_mbd_mesh.forces.time_function;
%!               load_case.pressure = struct();
%!               load_case.traction = struct();
%!               load_case.traction_abs = struct();
%!               for i=1:numel(grp_idx_load_px)
%!                 if (~isfield(mesh.elements, elem_type_surf{i}))
%!                   continue;
%!                 endif
%!                 if (grp_idx_load_px(i) > 0)
%!                   elem_px = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_px(i)).elements;
%!                 else
%!                   elem_px = [];
%!                 endif
%!                 if (grp_idx_load_py(i) > 0)
%!                   elem_py = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_py(i)).elements;
%!                 else
%!                   elem_py = [];
%!                 endif
%!                 if (grp_idx_load_pz(i) > 0)
%!                   elem_pz = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_pz(i)).elements;
%!                 else
%!                   elem_pz = [];
%!                 endif
%!                 if (grp_idx_load_mx(i) > 0)
%!                   elem_mx = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_mx(i)).elements;
%!                 else
%!                   elem_mx = [];
%!                 endif
%!                 if (grp_idx_load_my(i) > 0)
%!                   elem_my = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_my(i)).elements;
%!                 else
%!                   elem_my = [];
%!                 endif
%!                 if (grp_idx_load_mz(i) > 0)
%!                   elem_mz = getfield(mesh.groups, elem_type_surf{i})(grp_idx_load_mz(i)).elements;
%!                 else
%!                   elem_mz = [];
%!                 endif
%!                 elem_nodes = [getfield(mesh.elements, elem_type_surf{i})(elem_px, :);
%!                               getfield(mesh.elements, elem_type_surf{i})(elem_py, :);
%!                               getfield(mesh.elements, elem_type_surf{i})(elem_pz, :);
%!                               getfield(mesh.elements, elem_type_surf{i})(elem_mx, :);
%!                               getfield(mesh.elements, elem_type_surf{i})(elem_my, :);
%!                               getfield(mesh.elements, elem_type_surf{i})(elem_mz, :)];
%!                 switch (load_type{idx_load_type})
%!                   case "pressure"
%!                     elem_press = [repmat(-sigma(1, idx_sigma), numel(elem_px), columns(elem_nodes));
%!                                   repmat(-sigma(2, idx_sigma), numel(elem_py), columns(elem_nodes));
%!                                   repmat(-sigma(3, idx_sigma), numel(elem_pz), columns(elem_nodes));
%!                                   repmat(-sigma(1, idx_sigma), numel(elem_mx), columns(elem_nodes));
%!                                   repmat(-sigma(2, idx_sigma), numel(elem_my), columns(elem_nodes));
%!                                   repmat(-sigma(3, idx_sigma), numel(elem_mz), columns(elem_nodes))];
%!                     load_case.pressure = setfield(load_case.pressure, ...
%!                                                   elem_type_surf{i}, ...
%!                                                   struct("elements", elem_nodes, ...
%!                                                          "p", elem_press));
%!                   case {"traction", "traction_abs"}
%!                     elem_trac = zeros(rows(elem_nodes), columns(elem_nodes), 3);
%!                     ioffset = 0;               elem_trac(ioffset + (1:numel(elem_px)), :, 1) += sigma(1, idx_sigma);
%!                     ioffset += numel(elem_px); elem_trac(ioffset + (1:numel(elem_py)), :, 2) += sigma(2, idx_sigma);
%!                     ioffset += numel(elem_py); elem_trac(ioffset + (1:numel(elem_pz)), :, 3) += sigma(3, idx_sigma);
%!                     ioffset += numel(elem_pz); elem_trac(ioffset + (1:numel(elem_mx)), :, 1) -= sigma(1, idx_sigma);
%!                     ioffset += numel(elem_mx); elem_trac(ioffset + (1:numel(elem_my)), :, 2) -= sigma(2, idx_sigma);
%!                     ioffset += numel(elem_my); elem_trac(ioffset + (1:numel(elem_mz)), :, 3) -= sigma(3, idx_sigma);
%!                     ioffset += numel(elem_mz);
%!                     ioffset = 0;               elem_trac(ioffset + (1:numel(elem_px)), :, 3) += sigma(6, idx_sigma);
%!                     ioffset += numel(elem_px); elem_trac(ioffset + (1:numel(elem_py)), :, 2) += 0;
%!                     ioffset += numel(elem_py); elem_trac(ioffset + (1:numel(elem_pz)), :, 1) += sigma(6, idx_sigma);
%!                     ioffset += numel(elem_pz); elem_trac(ioffset + (1:numel(elem_mx)), :, 3) -= sigma(6, idx_sigma);
%!                     ioffset += numel(elem_mx); elem_trac(ioffset + (1:numel(elem_my)), :, 2) -= 0;
%!                     ioffset += numel(elem_my); elem_trac(ioffset + (1:numel(elem_mz)), :, 1) -= sigma(6, idx_sigma);
%!                     ioffset += numel(elem_mz);
%!                     ioffset = 0;               elem_trac(ioffset + (1:numel(elem_px)), :, 2) += sigma(4, idx_sigma);
%!                     ioffset += numel(elem_px); elem_trac(ioffset + (1:numel(elem_py)), :, 1) += sigma(4, idx_sigma);
%!                     ioffset += numel(elem_py); elem_trac(ioffset + (1:numel(elem_pz)), :, 1) += 0;
%!                     ioffset += numel(elem_pz); elem_trac(ioffset + (1:numel(elem_mx)), :, 2) -= sigma(4, idx_sigma);
%!                     ioffset += numel(elem_mx); elem_trac(ioffset + (1:numel(elem_my)), :, 1) -= sigma(4, idx_sigma);
%!                     ioffset += numel(elem_my); elem_trac(ioffset + (1:numel(elem_mz)), :, 1) -= 0;
%!                     ioffset += numel(elem_mz);
%!                     ioffset = 0;               elem_trac(ioffset + (1:numel(elem_px)), :, 3) += 0;
%!                     ioffset += numel(elem_px); elem_trac(ioffset + (1:numel(elem_py)), :, 3) += sigma(5, idx_sigma);
%!                     ioffset += numel(elem_py); elem_trac(ioffset + (1:numel(elem_pz)), :, 2) += sigma(5, idx_sigma);
%!                     ioffset += numel(elem_pz); elem_trac(ioffset + (1:numel(elem_mx)), :, 2) -= 0;
%!                     ioffset += numel(elem_mx); elem_trac(ioffset + (1:numel(elem_my)), :, 3) -= sigma(5, idx_sigma);
%!                     ioffset += numel(elem_my); elem_trac(ioffset + (1:numel(elem_mz)), :, 2) -= sigma(5, idx_sigma);
%!                     ioffset += numel(elem_mz);
%!                     load_case = setfield(load_case, ...
%!                                          load_type{idx_load_type}, ...
%!                                          setfield(getfield(load_case, load_type{idx_load_type}), ...
%!                                                   elem_type_surf{i}, ...
%!                                                   struct("elements", elem_nodes, ...
%!                                                          "f", elem_trac)));
%!                 case "prestrain"
%!                     mesh.material_data.extra_data = ", prestrain, reference, tpl_drive_id_epsilon0";
%!                 endswitch
%!               endfor
%!               opt_mbd_mesh = mbdyn_pre_solid_write_nodes(mesh, nodes_file, opt_mbd_mesh);
%!               opt_mbd_mesh = mbdyn_pre_solid_write_const_laws(mesh, csl_file, opt_mbd_mesh);
%!               opt_mbd_mesh = mbdyn_pre_solid_write_elements(mesh, load_case_dof, load_case, elem_file, opt_mbd_mesh);
%!               idx_joint = int32(0);
%!               unwind_protect
%!                 [fd, msg] = fopen(set_file, "wt");
%!                 if (fd == -1)
%!                   error("failed to open file \"%s\": %s", set_file, msg);
%!                 endif
%!                 fprintf(fd, "set: integer number_of_nodes = %d;\n", opt_mbd_mesh.struct_nodes.number);
%!                 fprintf(fd, "set: integer number_of_nodes_hydraulic = %d;\n", opt_mbd_mesh.hydraulic_nodes.number);
%!                 fprintf(fd, "set: integer number_of_solids = %d;\n", opt_mbd_mesh.solids.number);
%!                 fprintf(fd, "set: integer number_of_genels = %d;\n", opt_mbd_mesh.genels.number);
%!                 fprintf(fd, "set: integer number_of_forces = %d;\n", opt_mbd_mesh.forces.number);
%!                 fprintf(fd, "set: integer number_of_joints = %d;\n", idx_joint);
%!                 fprintf(fd, "set: integer number_of_surface_loads = %d;\n", opt_mbd_mesh.surface_loads.number);
%!                 fprintf(fd, "set: real t1 = %.16e;\n", t1);
%!                 fprintf(fd, "set: real t2 = %.16e;\n", t2);
%!                 fprintf(fd, "set: real dt = %.16e;\n", dt);
%!                 switch (load_type{idx_load_type})
%!                 case "prestrain"
%!                   fprintf(fd, "set: integer tpl_drive_id_epsilon0 = 2001;\n");
%!                   for i=1:6
%!                     fprintf(fd, "set: real epsilon0_%d = %.16e;\n", i, epsilon0(i, idx_sigma));
%!                   endfor
%!                 endswitch
%!               unwind_protect_cleanup
%!                 if (fd ~= -1)
%!                   fclose(fd);
%!                 endif
%!                 fd = -1;
%!               end_unwind_protect
%!               unwind_protect
%!                 [fd, msg] = fopen(control_file, "wt");
%!                 if (fd == -1)
%!                   error("failed to open file \"%s\": %s", control_file, msg);
%!                 endif
%!                 switch (model)
%!                   case "static"
%!                     fprintf(fd, "model: %s;\n", model);
%!                   case "dynamic"
%!                     fprintf(fd, "# model: %s;\n", model);
%!                 endswitch
%!               unwind_protect_cleanup
%!                 if (fd ~= -1)
%!                   fclose(fd);
%!                 endif
%!                 fd = -1;
%!               end_unwind_protect
%!               unwind_protect
%!                 [fd, msg] = fopen(initial_value_file, "wt");
%!                 if (fd == -1)
%!                   error("failed to open file \"%s\": %s", initial_value_file, msg);
%!                 endif
%!                 fprintf(fd, "method: %s;\n", method);
%!               unwind_protect_cleanup
%!                 if (fd ~= -1)
%!                   fclose(fd);
%!                 endif
%!                 fd = -1;
%!               end_unwind_protect
%!               fd = -1;
%!               unwind_protect
%!                 [fd, msg] = fopen(input_file, "wt");
%!                 if (fd == -1)
%!                   error("failed to open file \"%s\": %s", input_file, msg);
%!                 endif
%!                 fprintf(fd, "include: \"%s\";\n", set_file);
%!                 fprintf(fd, "begin: data;\n");
%!                 fprintf(fd, "        problem: initial value; # the default\n");
%!                 fprintf(fd, "end: data;\n");
%!                 fprintf(fd, "begin: initial value;\n");
%!                 fprintf(fd, "        initial time: 0;\n");
%!                 fprintf(fd, "        final time: t1 + t2;\n");
%!                 fprintf(fd, "        time step: dt;\n");
%!                 if (t2 > 0)
%!                   fprintf(fd, "        strategy: change, piecewise linear, 4, 0., dt, t1, dt, t1 + 0.1 * (t2 - t1), (t2 - t1) / t1 * dt, t2, (t2 - t1) / t1 * dt;\n");
%!                 endif
%!                 fprintf(fd, "        max iterations: 100;\n");
%!                 tolerance_type = "sepnorm";
%!                 switch (elem_type)
%!                 case {"tet10upc", "iso8upc", "iso20upc", "iso20upcr", "penta15upc"}
%!                   ## FIXME: test "sepnorm" does not work well with u/p-c elements; use norm instead
%!                   tolerance_type = "norm";
%!                 endswitch
%!                 switch (load_type{idx_load_type})
%!                 case "prestrain"
%!                   ## because there is no load
%!                   tolerance_type = "norm";
%!                 endswitch
%!                 switch (tolerance_type)
%!                 case "norm"
%!                   fprintf(fd, "        tolerance: 1e-6, test, norm, 1e-12, test, norm;\n");
%!                 case "sepnorm"
%!                   fprintf(fd, "        tolerance: 1e-12, test, sepnorm, 1e-12, test, norm;\n");
%!                 endswitch
%!                 fprintf(fd, "        output: messages;\n");
%!                 fprintf(fd, "        output: iterations, solver condition number, stat, yes;\n");
%!                 fprintf(fd, "        nonlinear solver: nox,\n");
%!                 fprintf(fd, "                          modified, 30,\n");
%!                 fprintf(fd, "                          keep jacobian matrix,\n");
%!                 fprintf(fd, "                          use preconditioner as solver, no,\n");
%!                 fprintf(fd, "                          linesearch method, backtrack,\n");
%!                 fprintf(fd, "                          direction, newton,\n");
%!                 fprintf(fd, "                          jacobian operator, newton krylov,\n");
%!                 fprintf(fd, "                          forcing term, constant,\n");
%!                 fprintf(fd, "                          linear solver tolerance, 1e-12,\n");
%!                 fprintf(fd, "                          inner iterations before assembly, 15,\n");
%!                 fprintf(fd, "                          linear solver max iterations, 300,\n");
%!                 fprintf(fd, "                          krylov subspace size, 300,\n");
%!                 fprintf(fd, "                          minimum step, 1e-12,\n");
%!                 fprintf(fd, "                          recovery step type, constant,\n");
%!                 fprintf(fd, "                          recovery step, 1e-12,\n");
%!                 fprintf(fd, "                          verbose, 3,\n");
%!                 fprintf(fd, "                          print convergence info, no;\n");
%!                 fprintf(fd, "        linear solver: umfpack, grad, scale, iterative, always, max iterations, 3;\n");
%!                 fprintf(fd, "        derivatives coefficient: 1e-6, auto;\n");
%!                 fprintf(fd, "        derivatives tolerance: 1e-5, 1e-5;\n");
%!                 fprintf(fd, "        derivatives max iterations: 10;\n");
%!                 fprintf(fd, "        threads: assembly, 1;\n");
%!                 fprintf(fd, "        threads: solver, 1;\n");
%!                 fprintf(fd, "        output: cpu time;\n");
%!                 fprintf(fd, "        include: \"%s\";\n", initial_value_file);
%!                 fprintf(fd, "end: initial value;\n");
%!                 fprintf(fd, "begin: control data;\n");
%!                 fprintf(fd, "       output meter: closest next, t1 + t2, forever, dt;\n");
%!                 fprintf(fd, "       skip initial joint assembly;\n");
%!                 fprintf(fd, "       output precision: 16;\n");
%!                 fprintf(fd, "       include: \"%s\";\n", control_file);
%!                 fprintf(fd, "       default output: all, solids, accelerations;\n");
%!                 fprintf(fd, "       structural nodes: number_of_nodes;\n");
%!                 fprintf(fd, "       hydraulic nodes: number_of_nodes_hydraulic;\n");
%!                 fprintf(fd, "       solids: number_of_solids;\n");
%!                 fprintf(fd, "       genels: number_of_genels;\n");
%!                 fprintf(fd, "       forces: number_of_forces;\n");
%!                 fprintf(fd, "       surface loads: number_of_surface_loads;\n");
%!                 fprintf(fd, "       joints: number_of_joints;\n");
%!                 fprintf(fd, "       use automatic differentiation;\n");
%!                 fprintf(fd, "end: control data;\n");
%!                 switch (load_type{idx_load_type})
%!                 case "prestrain"
%!                    fprintf(fd, "template drive caller: tpl_drive_id_epsilon0, 6,\n");
%!                    fprintf(fd, " green lagrange strain, single,\n");
%!                    for i=1:6
%!                      fprintf(fd, "  epsilon0_%d,\n", i);
%!                    endfor
%!                    fprintf(fd, "  min, 2, const, 1., mult, time, const, 1. / t1;\n");
%!                 endswitch
%!                 fprintf(fd, "include: \"%s\";\n", csl_file);
%!                 fprintf(fd, "begin: nodes;\n");
%!                 fprintf(fd, "       include: \"%s\";\n", nodes_file);
%!                 fprintf(fd, "end: nodes;\n");
%!                 fprintf(fd, "begin: elements;\n");
%!                 fprintf(fd, "       include: \"%s\";\n", elem_file);
%!                 fprintf(fd, "end: elements;\n");
%!               unwind_protect_cleanup
%!                 if (fd ~= -1)
%!                   fclose(fd);
%!                 endif
%!                 fd = -1;
%!               end_unwind_protect
%!               if (options.verbose)
%!                 shell(sprintf("cat \"%s\" | nl", set_file));
%!                 shell(sprintf("cat \"%s\" | nl", input_file));
%!                 shell(sprintf("cat \"%s\" | nl", nodes_file));
%!                 shell(sprintf("cat \"%s\" | nl", csl_file));
%!                 shell(sprintf("cat \"%s\" | nl", elem_file));
%!               endif
%!               fprintf(stderr, "element type: %s\n", elem_type);
%!               info = mbdyn_solver_run(input_file, opt_mbd);
%!               [mesh_sol, sol] = mbdyn_post_load_output_sol(output_file);
%!               [genel_id, genel_data] = mbdyn_post_load_output([output_file, ".gen"], 1, [], numel(sol.t), 1);
%!               switch (load_type{idx_load_type})
%!               case {"pressure", "traction", "traction_abs"}
%!                 [surfl_id, surfl_data] = mbdyn_post_load_output([output_file, ".prl"], 3, [], numel(sol.t), 3);
%!               endswitch
%!               F11 = (sol.def(grp_idx_p2, 1, :)(:) - sol.def(grp_idx_p1, 1, :)(:)) ./ (mesh.nodes(grp_idx_p2, 1, :)(:) - mesh.nodes(grp_idx_p1, 1, end)(:)) + 1;
%!               F22 = (sol.def(grp_idx_p4, 2, :)(:) - sol.def(grp_idx_p1, 2, :)(:)) ./ (mesh.nodes(grp_idx_p4, 2, :)(:) - mesh.nodes(grp_idx_p1, 2, end)(:)) + 1;
%!               F33 = (sol.def(grp_idx_p5, 3, :)(:) - sol.def(grp_idx_p1, 3, :)(:)) ./ (mesh.nodes(grp_idx_p5, 3, :)(:) - mesh.nodes(grp_idx_p1, 3, end)(:)) + 1;
%!               F12 = (sol.def(grp_idx_p4, 1, :)(:) - sol.def(grp_idx_p1, 1, :)(:)) ./ (mesh.nodes(grp_idx_p4, 2, :)(:) - mesh.nodes(grp_idx_p1, 2, end)(:));
%!               F31 = (sol.def(grp_idx_p2, 3, :)(:) - sol.def(grp_idx_p1, 3, :)(:)) ./ (mesh.nodes(grp_idx_p2, 1, :)(:) - mesh.nodes(grp_idx_p1, 1, end)(:));
%!               F32 = (sol.def(grp_idx_p4, 3, :)(:) - sol.def(grp_idx_p1, 3, :)(:)) ./ (mesh.nodes(grp_idx_p4, 2, :)(:) - mesh.nodes(grp_idx_p1, 2, end)(:));
%!               F21 = (sol.def(grp_idx_p2, 2, :)(:) - sol.def(grp_idx_p1, 2, :)(:)) ./ (mesh.nodes(grp_idx_p2, 1, :)(:) - mesh.nodes(grp_idx_p1, 1, end)(:));
%!               F13 = (sol.def(grp_idx_p5, 1, :)(:) - sol.def(grp_idx_p1, 1, :)(:)) ./ (mesh.nodes(grp_idx_p5, 3, :)(:) - mesh.nodes(grp_idx_p1, 3, end)(:));
%!               F23 = (sol.def(grp_idx_p5, 2, :)(:) - sol.def(grp_idx_p1, 2, :)(:)) ./ (mesh.nodes(grp_idx_p5, 3, :)(:) - mesh.nodes(grp_idx_p1, 3, end)(:));
%!               F = zeros(3, 3, numel(sol.t));
%!               for i=1:numel(sol.t)
%!                 F(1, 1, :) = F11;
%!                 F(1, 2, :) = F12;
%!                 F(1, 3, :) = F13;
%!                 F(2, 1, :) = F21;
%!                 F(2, 2, :) = F22;
%!                 F(2, 3, :) = F23;
%!                 F(3, 1, :) = F31;
%!                 F(3, 2, :) = F32;
%!                 F(3, 3, :) = F33;
%!               endfor
%!               G = C = zeros(size(F));
%!               for i=1:numel(sol.t)
%!                 G(:, :, i) = 0.5 * (F(:, :, i).' * F(:, :, i) - eye(3));
%!                 C(:, :, i) = F(:, :, i).' * F(:, :, i);
%!               endfor
%!               Epsilon = epsilon = zeros(6, size(G, 3));
%!               for i=1:3
%!                 epsilon(i, :) = sqrt(1 + 2 * G(i, i, :)(:).') - 1;
%!                 Epsilon(i, :) = G(i, i, :)(:).';
%!               endfor
%!               epsilon(4, :) = 2 * G(1, 2, :)(:).' ./ ((1 + epsilon(1, :)) .* (1 + epsilon(2, :)));
%!               epsilon(5, :) = 2 * G(2, 3, :)(:).' ./ ((1 + epsilon(2, :)) .* (1 + epsilon(3, :)));
%!               epsilon(6, :) = 2 * G(3, 1, :)(:).' ./ ((1 + epsilon(3, :)) .* (1 + epsilon(1, :)));
%!               Epsilon(4, :) = 2 * G(1, 2, :)(:).';
%!               Epsilon(5, :) = 2 * G(2, 3, :)(:).';
%!               Epsilon(6, :) = 2 * G(3, 1, :)(:).';
%!               tau_ref = sigma(:, idx_sigma);
%!               sin_gamma = zeros(3, 1);
%!               sin_gamma_cnt = 0;
%!               for i=1:numel(elem_type_solid)
%!                 if (~isfield(sol.strain.epsilon, elem_type_solid{i}))
%!                   continue;
%!                 endif
%!                 sin_gamma += mean(mean(getfield(sol.strain.epsilon, elem_type_solid{i})(:, :, 4:6, end), 1), 2)(:);
%!                 ++sin_gamma_cnt;
%!               endfor
%!               sin_gamma /= sin_gamma_cnt;
%!               cos_gamma = sqrt(1 - sin_gamma.^2);
%!               tan_gamma = sin_gamma ./ cos_gamma;
%!               tau_ref(1) += 2 * sigma(6, idx_sigma) * tan_gamma(3) + 2 * sigma(4, idx_sigma) * tan_gamma(1);
%!               tau_ref(2) += 2 * sigma(5, idx_sigma) * tan_gamma(2);
%!               tol_epsilon = 1e-8;
%!               tol_sigma = 1e-10;
%!               tol = 1e-9;
%!               tol_F = 1e-8;
%!               switch (load_type{idx_load_type})
%!               case "prestrain"
%!                 for j=1:numel(elem_type_solid)
%!                   if (~isfield(sol.strain.epsilon, elem_type_solid{j}))
%!                     continue;
%!                   endif
%!                   epsilon_res = getfield(sol.strain.epsilon, elem_type_solid{j});
%!                   sigma_res = getfield(sol.stress.tau, elem_type_solid{j});
%!                   assert_simple(max(max(max(max(abs(sigma_res))))) < tol_sigma * material.E);
%!                   for i=1:numel(sol.t)
%!                     for k=1:size(epsilon_res, 1)
%!                       for l=1:size(epsilon_res, 2)
%!                         assert_simple(epsilon_res(k, l, :, i)(:), epsilon0(:, idx_sigma), tol_epsilon * norm(epsilon0(:, idx_sigma)));
%!                       endfor
%!                     endfor
%!                   endfor
%!                 endfor
%!                 continue;
%!               endswitch
%!               Fsurfl_sum = zeros(3, 1);
%!               for i=1:numel(surfl_data)
%!                 Fsurfl_sum += surfl_data{i}(end, :)(:);
%!               endfor
%!               for i=1:numel(elem_type_solid)
%!                 if (~isfield(mesh.elements, elem_type_solid{i}))
%!                   continue;
%!                 endif
%!                 tau_res = getfield(sol.stress.tau, elem_type_solid{i});
%!                 for j=1:size(tau_res, 1)
%!                   for k=1:size(tau_res, 2)
%!                     assert_simple(tau_res(j, k, :, end)(:), tau_ref, tol * norm(tau_ref));
%!                   endfor
%!                 endfor
%!               endfor
%!               Fref = max(abs([geometry.w * geometry.h * tau_ref(1), ...
%!                               geometry.l * geometry.h * tau_ref(2), ...
%!                               geometry.l * geometry.w * tau_ref(3)]));
%!               for i=1:numel(genel_data)
%!                 assert_simple(all(all(abs(genel_data{i}) < tol_F * Fref)));
%!               endfor
%!               assert_simple(norm(Fsurfl_sum) < tol_F * Fref);
%!               for j=1:numel(elem_type_solid)
%!                 if (~isfield(sol.strain.epsilon, elem_type_solid{j}))
%!                   continue;
%!                 endif
%!                 epsilon_res = getfield(sol.strain.epsilon, elem_type_solid{j});
%!                 for i=1:numel(sol.t)
%!                   for k=1:size(epsilon_res, 1)
%!                     for l=1:size(epsilon_res, 2)
%!                       assert_simple(epsilon_res(k, l, :, i)(:), epsilon(:, i), tol_epsilon);
%!                     endfor
%!                   endfor
%!                 endfor
%!               endfor
%!               S_res = zeros(3, 3, numel(sol.t));
%!               switch (mesh.material_data.type)
%!                 case {"linear elastic generic", "hookean linear elastic isotropic", "bilinear isotropic hardening"}
%!                   switch (mesh.material_data.type)
%!                     case {"bilinear isotropic hardening"}
%!                       if (sigmav > mesh.material_data.sigmayv)
%!                         continue;
%!                       endif
%!                   endswitch
%!                   H = fem_pre_mat_isotropic(mesh.material_data.E, mesh.material_data.nu);
%!                   sigma_res = H * Epsilon;
%!                   for i=1:3
%!                     S_res(i, i, :) = sigma_res(i, :);
%!                   endfor
%!                   S_res(1, 2, :) = S_res(2, 1, :) = sigma_res(4, :);
%!                   S_res(2, 3, :) = S_res(3, 2, :) = sigma_res(5, :);
%!                   S_res(3, 1, :) = S_res(1, 3, :) = sigma_res(6, :);
%!                 case {"neo hookean elastic", "mooney rivlin elastic"}
%!                   mu = mesh.material_data.E / (2 * (1 + mesh.material_data.nu));
%!                   lambda = mesh.material_data.E * mesh.material_data.nu / ((1 + mesh.material_data.nu ) * (1 - 2 * mesh.material_data.nu));
%!                   for i=1:numel(sol.t)
%!                     IC = trace(C(:, :, i));
%!                     IIC = 1/2 * (trace(C(:, :, i))^2 - trace(C(:, :, i)^2));
%!                     IIIC = det(C(:, :, i));
%!                     invC = inv(C(:, :, i));
%!                     for k=1:3
%!                       for l=1:3
%!                         switch (mesh.material_data.type)
%!                         case "neo hookean elastic"
%!                           S_res(k, l, i) = mu * (k == l) + (lambda * (IIIC - sqrt(IIIC)) - mu) * invC(k, l);
%!                         case {"mooney rivlin elastic"}
%!                           C1 = mesh.material_data.G / (2 * (1 + mesh.material_data.delta));
%!                           C2 = mesh.material_data.delta * C1;
%!                           S_res(k, l, i) = 2 * (C1 * IIIC^(-1/3) * (k == l) + C2 * IIIC^(-2/3) * (IC * (k == l) - C(k, l, i)) + (1/2 * mesh.material_data.kappa * (IIIC - sqrt(IIIC)) - 1/3 * C1 * IC * IIIC^(-1/3) - 2/3 * C2 * IIC * IIIC^(-2/3)) * invC(k, l));
%!                         endswitch
%!                       endfor
%!                     endfor
%!                   endfor
%!               endswitch
%!               switch (mesh.material_data.type)
%!               case {"linear viscoelastic generic", "hookean linear viscoelastic isotropic", "neo hookean viscoelastic", "linear viscoelastic maxwelln", "linear viscoelastic maxwell1"}
%!                 ## TODO: viscoelastic case is not handled yet
%!               otherwise
%!                 tau_res = zeros(3, 3, numel(sol.t));
%!                 for i=1:numel(sol.t)
%!                   tau_res(:, :, i) = F(:, :, i) * S_res(:, :, i) * F(:, :, i).' / det(F(:, :, i));
%!                 endfor
%!                 Tau_res = zeros(6, numel(sol.t));
%!                 for i=1:3
%!                   Tau_res(i, :) = tau_res(i, i, :);
%!                 endfor
%!                 Tau_res(4, :) = tau_res(1, 2, :);
%!                 Tau_res(5, :) = tau_res(2, 3, :);
%!                 Tau_res(6, :) = tau_res(3, 1, :);
%!                 assert_simple(Tau_res(:, end), tau_ref, tol * norm(tau_ref));
%!               endswitch
%!             endfor
%!           endfor
%!         endfor
%!       endfor
%!     endfor
%!   endfor
%! unwind_protect_cleanup
%!   if (~isempty(filename))
%!     fn = dir([filename, "*"]);
%!     for i=1:numel(fn)
%!       if (0 ~= unlink(fullfile(fn(i).folder, fn(i).name)))
%!         warning("failed to remove file \"%s\"", fn(i).name);
%!       endif
%!     endfor
%!   endif
%! end_unwind_protect
