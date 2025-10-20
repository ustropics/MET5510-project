%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: eady_plot.m

% Description: Script for plotting results from the Eady wave model, including
% eigenvector amplitude, geopotential height, Hovmoller diagram, zonal wind,
% meridional wind, temperature, potential vorticity, additional cross-sections,
% Hovmoller diagrams, and background state fields. Loads data from 'eady_wave_#.mat',
% computes necessary fields, and saves plots to 'output/plots/'.

% Functions used:
% - eady_config: Loads model parameters
% - plot_evec_amp: Plots eigenvector amplitude contour
% - plot_geopotential_height: Plots geopotential height contour
% - plot_hovmoller: Plots Hovmoller diagram
% - plot_zonal_wind: Plots zonal wind contour
% - plot_meridional_wind: Plots meridional wind contour
% - plot_temperature: Plots temperature contour
% - plot_gph_top: Plots geopotential height at top boundary
% - plot_pvfield: Plots potential vorticity contour
% - plot_vg_cross_section: Plots meridional wind vertical cross-section
% - plot_ug_hovmoller: Plots Hovmoller diagram for zonal wind
% - plot_ubar_contour: Plots Ubar contour
% - plot_dpvdym_int: Plots d(PVbar)/dy interior contour
% - plot_dpvdym_boundaries: Plots d(PVbar)/dy at boundaries with beta
% - plot_background_flow: Plots combined background flow (Ubar, interior PV, boundaries)
% - XV2field: Computes 3D field from streamfunction vector
% - XV2streamxtime: Computes Hovmoller data for streamfunction
% - XV2ugxtime: Computes Hovmoller data for zonal wind
% - XV2XVy: Computes meridional derivative of streamfunction
% - XV2XVx: Computes zonal derivative of streamfunction
% - XV2XVz: Computes vertical derivative of streamfunction
% - jk2l: Converts 2D indices to linear index
% - l2jk: Converts linear index to 2D indices

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global jj kk ll BPVy NN2 f0 dy m0 dz Lx Ubar beta cplx

addpath(['functions', filesep])
addpath(['config', filesep])

%% Load parameters from eady_config.m
params = eady_config();

% Assign global variables
jj = params.jj;
kk = params.kk;
ll = params.ll;
Lx = params.Lx;
f0 = params.f0;
beta = params.beta;
gg = params.gg;
Theta0 = params.Theta0;
HH = params.HH;
dy = params.dy;
dz = params.dz;
NN2 = params.NN2;
cplx = params.cplx;
m0 = params.m0;

% Assign plotting parameters
ii = params.ii;
dx = params.dx;
xx = params.xx;
yy = params.yy;
zz = params.zz;
hlat = params.hlat;
hlevel = params.hlevel;
time = params.time;
n_mode = params.n_mode;

%% Load data from eady_wave_#.mat
eady_data = fullfile(params.eady_data_dir, params.eady_data_filename);
load(eady_data)

%% Select mode and compute properties
XV = zeros(ll, 1);
XV(:) = eigVec3(:, n_mode);
omega = imag(eigVal3(n_mode));
phase_speed = -omega / (2 * pi * m0 / Lx);
growth_rate = real(eigVal3(n_mode));
eFolding = (1 / growth_rate) / 86400; % in days

%% Normalize XV
valuemax = max(XV2field(XV, ii, dx) * f0 / gg, [], 'all');
XV = (10 / valuemax) * XV;

%% Compute fields
QV = B * XV;
XVy = XV2XVy(XV);
XVx = XV2XVx(XV);
XVz = XV2XVz(XV);

gpt_h = XV2field(XV, ii, dx) * f0 / gg;
temp = (f0 * HH / 287) * XV2field(XVz, ii, dx);
ug = -XV2field(XVy, ii, dx);
vg = XV2field(XVx, ii, dx);
pvfield = XV2field(QV, ii, dx);

%% Compute Hovmoller diagrams
gpt_h_hovmoler = XV2streamxtime(XV, ii, dx, omega, hlat, hlevel) * f0 / gg;
ug_hovmoler = XV2ugxtime(XVy, ii, dx, omega, hlat, hlevel);

%% Plot eigenvector amplitude
eVec_amp = zeros(jj + 1, kk + 1);
for l = 1 : ll
    [j, k] = l2jk(l);
    eVec_amp(j, k) = XV(l) .* conj(XV(l));
end
plot_evec_amp(yy, zz, eVec_amp, m0, n_mode, growth_rate, omega, model);

%% Plot geopotential height
plot_gph(xx, zz, gpt_h, jj, model, m0);

%% Plot Hovmoller diagram
plot_hovmoller(xx, time, gpt_h_hovmoler, model, m0);

%% Plot zonal wind
plot_zonal_wind(xx, yy, ug, model, m0);

%% Plot meridional wind
plot_meridional_wind(xx, yy, vg, model, m0);

%% Plot temperature at mid-level
plot_temperature(xx, yy, temp, kk, model, m0);

%% Plot geopotential height at top boundary
plot_gph_top(xx, yy, gpt_h, model, m0);

%% Plot potential vorticity at surface
plot_pvfield(xx, yy, pvfield, model, m0);

%% Plot meridional wind vertical cross-section
plot_vg_cross_section(xx, zz, vg, jj, model, m0);

%% Plot Hovmoller diagram for zonal wind
plot_ug_hovmoller(xx, time, ug_hovmoler, hlat, hlevel, model, m0);

%% Plot Ubar contour
plot_ubar_contour(yy, zz, Ubar, model, m0);

%% Plot d(PVbar)/dy interior
plot_dpvdym_int(yy, zz, BPVy, model, m0);

%% Plot d(PVbar)/dy at boundaries with beta
plot_dpvdym_boundaries(yy, BPVy, beta, kk, model, m0);

%% Plot combined background flow
plot_background_flow(yy, zz, Ubar, params.qy_surf, params.qy_trop, model, m0);