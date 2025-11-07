%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plt.m
%
% DESCRIPTION: Master plotting script for the Hoskins-West modified Eady model. 
% Loads precomputed diagnostic fields from calc.m, generates all diagnostic 
% plots (geopotential, winds, temperature, PV, vertical velocity, Hovmöller, 
% background state), and saves figures to disk. No computation is performed.
%
% INPUT:
% - Precomputed results: from calc_wave-m0_eMode-n.mat
% - Model parameters: from cfg()
%
% OUTPUT:
% - Saves all plots to: figures/m0/
%
% MATH/FUNCTIONS:
%   All fields already computed in calc.m
%   Plotting uses standard contourf, pcolor, quiver, etc.
%
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
% - Functions used: plot_gph, plot_zvu, plot_zvt, plot_zwt, plot_gph_hovmoller,
%   plot_evec_amp, plot_zonal_wind, plot_meridional_wind, plot_temperature,
%   plot_gph_top, plot_pvfield, plot_meridional_xsec, plot_ug_hovmoller,
%   plot_ubar_contour, plot_dpvdym_int, plot_dpvdym_boundaries, plot_background_flow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCRIPT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;
addpath(['functions', filesep])
addpath(['plots', filesep])
addpath(['config', filesep])

%% Load config and computed results
params = cfg();
load(fullfile(params.data_dir, params.calc_filename), 'results');

jj = params.jj;
kk = params.kk;
ll = params.ll;

% Unpack results from calculated date file
xx = results.xx; yy = results.yy; zz = results.zz; time = results.time;
gpt_h = results.gpt_h; temp = results.temp; ug = results.ug; vg = results.vg;
pvfield = results.pvfield; wfield = results.wfield; XV = results.XV;
gpt_h_hovmoler1 = results.gpt_h_hovmoler1;
gpt_h_hovmoler51 = results.gpt_h_hovmoler51;
ug_hovmoler1 = results.ug_hovmoler1;
ug_hovmoler51 = results.ug_hovmoler51;
Ubar = results.Ubar; BPVy = results.BPVy;
m0 = results.m0; n_mode = results.n_mode;
growth_rate = results.growth_rate; omega = results.omega;
hlat = results.hlat;

fig_path = params.plot_dir;

%% Ensure plot directory exists
if ~exist(fig_path, 'dir')
    mkdir(fig_path);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Geopotential height (full field)
plot_gph(xx, zz, gpt_h, jj, m0, n_mode, fig_path);
plot_momentum_flux(vg, ug, m0, n_mode, fig_path);
plot_meridional_hflux(vg, temp, m0, n_mode, fig_path);
plot_vertical_hflux(wfield, temp, m0, n_mode, fig_path);

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

% Zonal wind at levels 1, 25, 51
hlevels = [1, 25, 51];
for h = hlevels
    fprintf('Plotting zonal wind at hlevel = %d\n', h);
    plot_zonal_wind(xx, yy, ug, h, m0, n_mode, fig_path);
end

% Meridional wind at levels 1, 25, 51
for h = hlevels
    fprintf('Plotting meridional wind at hlevel = %d\n', h);
    plot_meridional_wind(xx, yy, vg, h, m0, n_mode, fig_path);
end

% Temperature at levels 1, 25, 51
for h = hlevels
    fprintf('Plotting temperature at hlevel = %d\n', h);
    plot_temp(xx, yy, temp, h, m0, n_mode, fig_path);
end

% Top boundary geopotential
plot_gph_top(xx, yy, gpt_h, 51, m0, n_mode, fig_path);

% Potential vorticity
plot_pvfield(xx, yy, pvfield, m0, n_mode, fig_path);

% Meridional wind cross-section
plot_meridional_xsec(xx, zz, vg, jj, m0, n_mode, fig_path);

% Zonal wind Hovmoller
plot_ug_hovmoller(xx, time, ug_hovmoler1, hlat, 1, m0, n_mode, fig_path);
plot_ug_hovmoller(xx, time, ug_hovmoler51, hlat, 51, m0, n_mode, fig_path);

% Background state
plot_ubar(yy, zz, Ubar, m0, n_mode, fig_path);
plot_dpvdym_int(yy, zz, BPVy, m0, n_mode, fig_path);
plot_dpvdym_bnd(yy, BPVy, params.beta, params.kk, m0, n_mode, fig_path);
plot_bg_flow(yy, zz, jj, params.kk, Ubar, BPVy, m0, n_mode, fig_path);

disp('All plots generated successfully in:');
disp(fig_path);
