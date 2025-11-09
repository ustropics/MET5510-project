%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% FILENAME: plt.m

% DESCRIPTION: Master plotting script for the Hoskins-West modified Eady model. 
% Loads precomputed diagnostic fields from calc.m, generates all diagnostic 
% plots (geopotential, winds, temperature, PV, vertical velocity, Hovmöller, 
% background state), and saves figures to disk. No computation is performed.

% INPUT:
% - Precomputed params: from calc_wave-m0_eMode-n.mat
% - Model parameters: from cfg()

% OUTPUT:
% - Saves all plots to: figures/m0/

% MATH/FUNCTIONS:
%   All fields already computed in calc.m
%   Plotting uses standard contourf, pcolor, quiver, etc.

% VARIABLES:
% - gpt_h, temp, ug, vg, pvfield, wfield: 3D diagnostic fields
% - *_hovmoler*: Time-longitude evolution at fixed latitude/level
% - Ubar, BPVy: Background state (zonal wind, PV gradient)
% - xx, yy, zz: Coordinate vectors (longitude, latitude, height)
% - time: Time vector for Hovmöller (days)
% - hlat: Latitude for Hovmöller (degrees)
% - m0, n_mode: Wavenumber and mode index
% - growth_rate, omega: Stability parameters
% - fig_path: Output directory for plots

% - Use "plot_one('list')" to get a full list of plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
addpath(['functions', filesep])
addpath(['plots', filesep])
addpath(['cmaps', filesep])

%% Load config and computed params
params = cfg();
load(fullfile(params.data_dir, params.calc_filename));

jj = params.jj;
kk = params.kk;
ll = params.ll;

% Unpack params from calculated date file
xx = params.xx; yy = params.yy; zz = params.zz; time = params.time;
m0 = params.m0; n_mode = params.n_mode;
hlat = params.hlat;

fig_path = params.plot_dir;
fig_path2 = fullfile(fig_path, 'combined');

%% Ensure plot directory exists
if ~exist(fig_path2, 'dir')
    mkdir(fig_path2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot fluxes
plot_momentum_flux(zz, yy, vg, ug, m0, n_mode, fig_path);
plot_meridional_hflux(zz, yy, vg, temp, m0, n_mode, fig_path);
plot_vertical_hflux(zz, yy, wfield, temp, m0, n_mode, fig_path);

% Hovmoller: Geopotential height
plot_gph_hovmoller(xx, time, gpt_h_hovmoler1, m0, n_mode, fig_path, 1);
plot_gph_hovmoller(xx, time, gpt_h_hovmoler51, m0, n_mode, fig_path, 51);

% Plot eigenvector amplitude
eVec_amp = zeros(jj + 1, kk + 1);

for l = 1 : ll
    [j, k] = l2jk(l);
    eVec_amp(j, k) = XV(l) .* conj(XV(l));
end

plot_evec_amp(yy, zz, eVec_amp, m0, n_mode, growth_rate, omega, fig_path);

% Combined plots
combined_momentum_vertical_hflux(vg, ug, wfield, temp, m0, n_mode, fig_path2)
combined_momentum_meridional_hflux(vg, ug, temp, m0, n_mode, fig_path2)
combined_ubar_evec_amp(yy, zz, Ubar, eVec_amp, m0, n_mode, ...
                                 growth_rate, omega, fig_path2)

% Zonal wind at levels 1, 26, 51
hlevels = [1, 26, 51];
for h = hlevels
    plot_zonal_wind(xx, yy, ug, h, m0, n_mode, fig_path);
    plot_meridional_wind(xx, yy, vg, h, m0, n_mode, fig_path)
    plot_temp(xx, yy, temp, h, m0, n_mode, fig_path);
    plot_pvfield(xx, yy, pvfield, h, m0, n_mode, fig_path);
    plot_gph(xx, yy, gpt_h, h, m0, n_mode, fig_path);
    combined_meridional_wind_temp(xx, yy, vg, temp, h, m0, n_mode, fig_path2);
    combined_zonal_wind_temp(xx, yy, ug, temp, h, m0, n_mode, fig_path2)
    combined_gph_temp(xx, yy, gpt_h, temp, h, m0, n_mode, fig_path2)
    combined_pvfield_temp(xx, yy, pvfield, temp, h, m0, n_mode, fig_path2)
end

% Top boundary geopotential
plot_gph_xsec(xx, zz, gpt_h, jj, m0, n_mode, fig_path);


% Plot cross-sections
plot_meridional_xsec(xx, zz, vg, jj, m0, n_mode, fig_path);
combined_gph_meridional_xsec(xx, zz, gpt_h, vg, jj, m0, n_mode, fig_path2)

% Zonal wind Hovmoller
plot_zonal_hovmoller(xx, time, ug_hovmoler1, hlat, 1, m0, n_mode, fig_path);
plot_zonal_hovmoller(xx, time, ug_hovmoler51, hlat, 51, m0, n_mode, fig_path);

% Background state
plot_ubar(yy, zz, Ubar, m0, n_mode, fig_path);
plot_dpvdym_int(yy, zz, BPVy, m0, n_mode, fig_path);
plot_dpvdym_bnd(yy, BPVy, params.beta, params.kk, m0, n_mode, fig_path);
combined_bg_flow(yy, zz, jj, params.kk, Ubar, BPVy, m0, n_mode, fig_path);

disp('All plots generated successfully in:');
disp(fig_path);
