function lf = compute_lf_femfuns(dp_elec, rad_elec)

%% Create a spherical volume conductor with two spheres of radius 7 and 10 cm at the origin
csvol.o = [0, 0,0];
csvol.r = [0.07 0.1];
cfg = [];
cfg.numvertices = 1000;
csbnd = ft_prepare_mesh(cfg, csvol);

%% pick a few electrode positions on top of the sphere
sel = find(csbnd(1).pos(:,3)>0); sel = sel(1:100:end);
elec = [];
elec.elecpos = csbnd(1).pos(sel,:);
for i=1:length(sel)
  elec.label{i} = sprintf('elec%d', i);
end
elec.unit = 'm';
% update the electrode sets to the latest standards
elec = ft_datatype_sens(elec);

% combine the brain surface with electrode surfaces and get inner points of the electrodes (here 'elecmarkers')
% dp_elec = 0.005; %height  of the electrode cylinder
% rad_elec = 0.002; %radius of the electrode cylinder
[dented_elsurf,elecmarkers] = add_electrodes(csbnd(1), elec.elecpos, rad_elec, dp_elec);
merged_surfs = add_surf(dented_elsurf,csbnd(2)); %combine with the scalp

%create volumetric tetrahedral mesh
 [tet_node,tet_elem] = s2m(merged_surfs.pos,merged_surfs.tri, 1, 1, 'tetgen', [point_in_surf(csbnd(1));point_in_surf(csbnd(2));elecmarkers]);
%label the electrode surfaces where they make contact with the brain
el_faces = label_surf(tet_elem, 3:7, 1);

%% Construct the FT mesh structure
mesh.unit = 'm';
mesh.pos = tet_node;
mesh.tet = tet_elem(:,1:4);
mesh.tri = el_faces(:,1:3);
mesh.boundary = el_faces(:,4);
mesh.boundarylabel = elec.label;
mesh.tissue = tet_elem(:,5);
mesh.tissuelabel = [{'brain'}, {'skull'},elec.label(:)'];

%% construct a vol to create the FT sourcemodel
vol.pos = mesh.pos;
vol.tet = mesh.tet;
vol.tissue = mesh.tissue;
vol.tissuelabel = mesh.tissuelabel;
vol.unit = mesh.unit;
vol.type = 'simbio';

cfg                 = [];
cfg.resolution      = 0.01; %in the same unit as the mesh &co
cfg.headmodel       = vol;
cfg.inwardshift     = 1; %shifts dipoles away from surfaces
sourcemodel         = ft_prepare_sourcemodel(cfg);
sourcemodel.inside(:) = false;

elec.label{1} = ['stim_' elec.label{1}];
elec.label{2} = ['ground_' elec.label{2}];
mesh.boundarylabel = elec.label;
elec.params{1} = {1, 100, 'int'};
elec.params{2} = {0, 100, 'int'};
conductivities = [0.33 0.01 0 0 1e10 1e10 1e10];
lf = femfuns_leadfield(mesh,conductivities,sourcemodel,elec);