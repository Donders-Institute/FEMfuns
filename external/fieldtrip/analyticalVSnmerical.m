%% validate the numerical result from FEMfuns with analytical solution in FT
% run the part of tutorial about dipole in a sphere here: https://www.fieldtriptoolbox.org/development/project/femfuns/#running-a-simulation-with-fieldtrip-and-femfuns-combined
clear; clc;
addpath PATH2FIELDTRIP
ft_defaults

addpath PATH2fieldtrip/external/iso2mesh/
cd PATH2FEMfuns/external/fieldtrip/
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
%%

% combine the brain surface with electrode surfaces and get inner points of the electrodes (here 'elecmarkers')
dp_elec = 0.005; %height  of the electrode cylinder
rad_elec = 0.002; %radius of the electrode cylinder
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
cfg.resolution      = 0.05; %in the same unit as the mesh &co - 10 dipoles inside, 17 dipoles outside brain
cfg.headmodel       = vol;
cfg.inwardshift     = 1; %shifts dipoles away from surfaces
sourcemodel         = ft_prepare_sourcemodel(cfg);

%% conductivities for brain, SKULL and metal electrodes are set
conductivities = [0.33 0.01 1e10 1e10 1e10 1e10 1e10];
lf_rec = femfuns_leadfield(mesh,conductivities,sourcemodel,elec);

disp(lf_rec)
%       dim: [3 3 3]
%       pos: [27×3 double]
%      unit: 'cm'
%    inside: [27×1 logical]
%       cfg: [1×1 struct]
% leadfield: {1×27 cell}
%     label: {'elec1'  'elec2'  'elec3'  'elec4'  'elec5'}

%% compute analytical solution in fieldtrip

headmodel = ft_headmodel_singlesphere(csbnd(1), 'conductivity', 0.33);
cfg = [];
cfg.sourcemodel = sourcemodel;    %% where are the sources?
cfg.headmodel   = headmodel;      %% how do currents spread?
cfg.elec        = elec; %% where are the sensors?

% how do sources and sensors connect?
lf_analytical = ft_prepare_leadfield(cfg);
disp(lf_analytical)


%% error computation

lf_an_matrix = cell2mat(lf_analytical.leadfield);
lf_femfuns_matrix = cell2mat(lf_rec.leadfield)*3e10; %scaling factor?

ratioRMS = 100*rms(lf_an_matrix-lf_femfuns_matrix,2)./rms(lf_an_matrix,2);

% RRMSE: relative root mean square error in percentage (see e.g., https://doi.org/10.1016/j.rser.2015.11.058)
rrmse = 100*sqrt((1/numel(lf_femfuns_matrix)) * sum(sum((lf_an_matrix-lf_femfuns_matrix).^2)))/sum(sum(lf_an_matrix));







