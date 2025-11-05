%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: hwme_plot.m

% DESCRIPTION: Script for plotting results from the Hoskins-West modified 
% Eady-type model, including eigenvector amplitude, geopotential height, 
% Hovmoller diagrams for geopotential height and zonal wind at hlevel = 1 and 51, 
% zonal wind at hlevel = 1, 25, and 51, meridional wind at hlevel = 1, 25, and 51, 
% temperature at hlevel = 1, 25, and 51, potential vorticity, additional 
% cross-sections, and background state fields. Loads data from 'hwme_wave_#.mat', 
% computes necessary fields, and saves plots to 'output/plots/'.

% SCRIPTS:
% - hwme_config.m: Loads model parameters

% PLOTS:
% - plot_evec_amp: Plots eigenvector amplitude contour
% - plot_geopotential_height: Plots geopotential height contour
% - plot_hovmoller: Plots Hovmoller diagram for geopotential height (hlevel = 1 and 51)
% - plot_zonal_wind: Plots zonal wind contour (hlevel = 1, 25, 51)
% - plot_meridional_wind: Plots meridional wind contour (hlevel = 1, 25, 51)
% - plot_temperature: Plots temperature contour (hlevel = 1, 25, 51)
% - plot_gph_top: Plots geopotential height at top boundary
% - plot_pvfield: Plots potential vorticity contour
% - plot_vg_cross_section: Plots meridional wind vertical cross-section
% - plot_ug_hovmoller: Plots Hovmoller diagram for zonal wind (hlevel = 1 and 51)
% - plot_ubar_contour: Plots Ubar contour
% - plot_dpvdym_int: Plots d(PVbar)/dy interior contour
% - plot_dpvdym_boundaries: Plots d(PVbar)/dy at boundaries with beta

% FUNCTIONS:
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

global jj kk ll xx yy zz HH ii BPVy NN2 f0 gg dx dy m0 dz Lx Theta0 ... 
    Ubar beta cplx hlat hlevel time n_mode hwme_data fig_path

addpath(['functions', filesep])
addpath(['config', filesep])
addpath(['plots', filesep])

%% Load parameters from hwme_config.m
params = hwme_config();

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
fig_path = params.hwme_plot_dir;

% Assign plotting parameters
ii = params.ii;
dx = params.dx;
xx = params.xx;
yy = params.yy;
zz = params.zz;
hlat = params.hlat;
time = params.time;
n_mode = params.n_mode;

%% Load data from hwme_wave_#.mat
hwme_data = fullfile(params.hwme_data_dir, params.hwme_data_filename);
load(hwme_data)

%% Select mode and compute properties
XV = zeros(ll, 1);
XV(:) = eigVec2(:, n_mode);
omega = imag(eigVal2(n_mode));
phase_speed = -omega / (2 * pi * m0 / Lx);
growth_rate = real(eigVal2(n_mode));
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

[F1, F2, F3] = eady_F123(params, XV, Ubar);
w = (f0/NN2)*(gg^-1)*(F1+F2+F3);
wfield = w2wfield(w,ii,dx);

%% Compute Hovmoller diagrams for hlevel = 1 and hlevel = 51
gpt_h_hovmoler1 = XV2streamxtime(XV, ii, dx, omega, hlat, 1) * f0 / gg;
gpt_h_hovmoler51 = XV2streamxtime(XV, ii, dx, omega, hlat, 51) * f0 / gg;
ug_hovmoler1 = XV2ugxtime(XVy, ii, dx, omega, hlat, 1);
ug_hovmoler51 = XV2ugxtime(XVy, ii, dx, omega, hlat, 51);

% [max1, ind1] = max(eVec_amp(:));
% [jmax, kmax] = ind2sub(size(eVec_amp), ind1);

if ~exist(params.hwme_plot_dir, 'dir')
    mkdir(params.hwme_plot_dir);
end

% ckplot_zvt(vg, temp, m0, n_mode, fig_path);

% %% Plot eigenvector amplitude
% eVec_amp = zeros(jj + 1, kk + 1);
% for l = 1 : ll
%     [j, k] = l2jk(l);
%     eVec_amp(j, k) = XV(l) .* conj(XV(l));
% end
% plot_evec_amp(yy, zz, eVec_amp, m0, n_mode, growth_rate, omega, model, fig_path);

% Plot combo plot for hwme
% plot_background_flow(yy, zz, jj, kk, Ubar, BPVy, model, m0, n_mode, fig_path);

% Plot geopotential height
% plot_gph(xx, zz, gpt_h, jj, model, m0, n_mode, fig_path);


% Plot Hovmoller diagrams for geopotential height (hlevel = 1 and 51)
plot_gph_hovmoller(xx, time, gpt_h_hovmoler1, model, m0, n_mode, fig_path, 1);
plot_gph_hovmoller(xx, time, gpt_h_hovmoler51, model, m0, n_mode, fig_path, 51);

% Plot zonal wind at hlevel = 1, 25, and 51
% fprintf('Generating zonal wind plot for hlevel = 1\n'); % Debug output
% plot_zonal_wind(xx, yy, ug, 1, model, m0, n_mode, fig_path);
% fprintf('Generating zonal wind plot for hlevel = 25\n'); % Debug output
% plot_zonal_wind(xx, yy, ug, 25, model, m0, n_mode, fig_path);
% fprintf('Generating zonal wind plot for hlevel = 51\n'); % Debug output
% plot_zonal_wind(xx, yy, ug, 51, model, m0, n_mode, fig_path);

% Plot meridional wind at hlevel = 1, 25, and 51
fprintf('Generating meridional wind plot for hlevel = 1\n'); % Debug output
plot_meridional_wind(xx, yy, vg, 1, model, m0, n_mode, fig_path);
fprintf('Generating meridional wind plot for hlevel = 25\n'); % Debug output
plot_meridional_wind(xx, yy, vg, 25, model, m0, n_mode, fig_path);
fprintf('Generating meridional wind plot for hlevel = 51\n'); % Debug output
plot_meridional_wind(xx, yy, vg, 51, model, m0, n_mode, fig_path);

% % Plot temperature at hlevel = 1, 25, and 51
% fprintf('Generating temperature plot for hlevel = 1\n'); % Debug output
% plot_temperature(xx, yy, temp, 1, model, m0, n_mode, fig_path);
% fprintf('Generating temperature plot for hlevel = 25\n'); % Debug output
% plot_temperature(xx, yy, temp, 25, model, m0, n_mode, fig_path);
% fprintf('Generating temperature plot for hlevel = 51\n'); % Debug output
% plot_temperature(xx, yy, temp, 51, model, m0, n_mode, fig_path);
% 
% % Plot geopotential height at top boundary
% plot_gph_top(xx, yy, gpt_h, model, m0, n_mode, fig_path);
% 
% % Plot potential vorticity at surface
% plot_pvfield(xx, yy, pvfield, model, m0, n_mode, fig_path);
% 
% % Plot meridional wind vertical cross-section
% plot_vg_cross_section(xx, zz, vg, jj, model, m0, n_mode, fig_path);
% 
% % Plot Hovmoller diagrams for zonal wind (hlevel = 1 and 51)
% plot_ug_hovmoller(xx, time, ug_hovmoler1, hlat, 1, model, m0, n_mode, fig_path);
% plot_ug_hovmoller(xx, time, ug_hovmoler51, hlat, 51, model, m0, n_mode, fig_path);
% 
% % Plot Ubar contour
% plot_ubar_contour(yy, zz, Ubar, model, m0, n_mode, fig_path);
% 
% % Plot d(PVbar)/dy interior
% plot_dpvdym_int(yy, zz, BPVy, model, m0, n_mode, fig_path);
% 
% % Plot d(PVbar)/dy at boundaries with beta
% plot_dpvdym_boundaries(yy, BPVy, beta, kk, model, m0, n_mode, fig_path);
% 
% % Plot combo plot for hwme
% plot_hwme_bg_flow(yy, zz, jj, kk, Ubar, BPVy, model, m0, n_mode, fig_path);
% 
% disp('All plots generated successfully.');