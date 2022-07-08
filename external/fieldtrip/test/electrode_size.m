%% compare numerical results in FEMfuns with different size electrodes

clear; clc;
addpath PATH2FIELDTRIP
ft_defaults

%addpath PATH2fieldtrip/external/iso2mesh/
cd PATH2FEMfuns/external/fieldtrip/
addpath PATH2compute_lf_femfuns_folder

%% small electrode
dp_elec = 0.005;
rad_elec = 0.002;
lf_small = compute_lf_femfuns(dp_elec, rad_elec);

%% visualize results
% a vtu (brain_solution_example000000.vtu) and a psd (brain_solution_example.psd) file are created by the femfunscaller.py and saved in the same folder where the script is launched from (subfolder pvds_dir)
% you can manually modify the filename before running the code for the bigger electrode. you can subsequently visualize the results in paraview.

%% big electrode
dp_elec = 0.01;
rad_elec = 0.004;
lf_big = compute_lf_femfuns(dp_elec, rad_elec);

%%
figure, plot(lf_small.leadfield(3:end)*1e9)
hold on
plot(lf_big.leadfield(3:end)*1e9)

%% compare potentials in paraview with Figure1.png
