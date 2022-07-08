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
rad_elec = 0.003; %radius of the electrode cylinder
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
% mesh.boundarylabel = elec.label; %%%% this is important
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
cfg.resolution      = 0.05; %in the same unit as the mesh &co - 11 dipoles inside, 114 dipoles outside brain
cfg.headmodel       = vol;
cfg.inwardshift     = 1; %shifts dipoles away from surfaces
sourcemodel         = ft_prepare_sourcemodel(cfg);

%% compute analytical solution in fieldtrip
% headmodel = ft_headmodel_concentricspheres(csbnd, 'conductivity', [0.33
% 0.01]);  % here I pass csbnd that is the surface mesh sphere
% here I pass the surf mesh with the electrodes

compl_mesh = csbnd;
compl_mesh(1).pos = dented_elsurf.pos;
compl_mesh(1).tri = dented_elsurf.tri(:,1:3);


headmodel = ft_headmodel_concentricspheres(compl_mesh, 'conductivity', [0.33 0.01]);  

cfg = [];
cfg.sourcemodel = sourcemodel;    %% where are the sources?
cfg.headmodel   = headmodel;      %% how do currents spread?
cfg.elec        = elec; %% where are the sensors?
lf_analytical_useless = ft_prepare_leadfield(cfg);%Warning: electrodes do not lie on skin surface -> using radial projection

%in debugging mode I saved the electrodes projected into the skin surface.
% from these:
% sens.elecpos
% 
% ans =
% 
%     0.0000         0    0.0700
%    -0.0097    0.0426    0.0547
%     0.0154   -0.0525    0.0436
%     0.0460   -0.0432    0.0304
%    -0.0019   -0.0682    0.0156
% to these:
% sens.elecpos
% 
% ans =
% 
%    -0.0000   -0.0000    0.1000
%    -0.0139    0.0608    0.0782
%     0.0220   -0.0750    0.0623
%     0.0657   -0.0617    0.0434
%    -0.0028   -0.0975    0.0223

%% run analytical again in the skin electrodes
% how do sources and sensors connect?

elec_scalp = [];
elec_scalp.elecpos = [-0.0000   -0.0000    0.1000;    -0.0139    0.0608    0.0782;     0.0220   -0.0750    0.0623;     0.0657   -0.0617    0.0434;    -0.0028   -0.0975    0.0223];
for i=1:length(sel)
  elec_scalp.label{i} = sprintf('elec%d', i);
end
elec_scalp.unit = 'm';
% update the electrode sets to the latest standards
elec_scalp = ft_datatype_sens(elec_scalp);
cfg.elec = elec_scalp;
lf_analytical = ft_prepare_leadfield(cfg); %% why it is projecting again into skin?????

   

disp(lf_analytical)
%                 dim: [3 3 3]
%                 pos: [27×3 double]
%                unit: 'm'
%              inside: [27×1 logical]
%                 cfg: [1×1 struct]
%           leadfield: {1×27 cell}
%               label: {'elec1'  'elec2'  'elec3'  'elec4'  'elec5'}
%     leadfielddimord: '{pos}_chan_ori'



%% compute leadfield with femfuns 

conductivities_bio = [0.33 0.01];
elec_cond = repmat(0.01,size(elec_scalp.chanpos,1),1); % I give here the skull conductivity to the electrodes
conductivities = [conductivities_bio elec_cond'];

elec_scalp.params = cell(5,1);

lf_rec = femfuns_leadfield(mesh,conductivities,sourcemodel,elec_scalp);

disp(lf_rec)
%       dim: [3 3 3]
%       pos: [27×3 double]
%      unit: 'cm'
%    inside: [27×1 logical]
%       cfg: [1×1 struct]
% leadfield: {1×27 cell}
%     label: {'elec1'  'elec2'  'elec3'  'elec4'  'elec5'}


%% are they all CAR?

lf_an_matrix = cell2mat(lf_analytical.leadfield);
lf_femfuns_matrix = cell2mat(lf_rec.leadfield);


car_an = mean(lf_an_matrix,1);
car_femfuns = mean(lf_femfuns_matrix,1);

figure, plot(car_an), hold on, plot(car_femfuns)

%% NO! rereferencing

lf_femfuns_matrix = lf_femfuns_matrix - repmat(car_femfuns,5,1);


%% rdm and mag ana vs femfuns
close all
for i=1:5
    num_norm = norm(lf_femfuns_matrix(i,:));
    ana_norm = norm(lf_an_matrix(i,:));
    rdm(i) = 50*norm(lf_femfuns_matrix(i,:)./num_norm - lf_an_matrix(i,:)./ana_norm);
    mag(i) = (num_norm/ana_norm);
end

figure, plot(rdm), title('rdm ana vs femfuns') %values in percentage below 10%
figure, plot(mag), title('mag ana vs femfuns') %values should be around 1
