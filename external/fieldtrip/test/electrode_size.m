%% compare numerical results in FEMfuns with different size electrodes

clear; clc;
addpath PATH2FIELDTRIP
ft_defaults

%addpath PATH2fieldtrip/external/iso2mesh/
cd PATH2FEMfuns/external/fieldtrip/
addpath PATH2compute_lf_femfuns_folder
%%

dp_elec = 0.005;
rad_elec = 0.002;
lf_small = compute_lf_femfuns(dp_elec, rad_elec);


dp_elec = 0.01;
rad_elec = 0.004;
lf_big = compute_lf_femfuns(dp_elec, rad_elec);

%%
figure, plot(lf_small.leadfield(3:end)*1e9)
hold on
plot(lf_big.leadfield(3:end)*1e9)

