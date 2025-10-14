%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: hwe_plot.m

% Description: Script for loading data from 'hwe_wave_#.mat' and 
% generating plots of the Hoskins-West Eady-type Model's background flow, 
% including Ubar, d(PVbar)/dy interior, and boundary PV gradients, along 
% with perturbation fields (meridional wind, temperature, and geopotential 
% height) for a specified mode, saving results to the 'plots' folder.

% Functions used: 
% - XV2field: Transforms eigenvector to 3D field.
% - XV2XVz: Computes z-derivative of eigenvector.
% - XVx2field: Converts x-derivative of eigenvector to 3D field.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('functions'); % add functions folder
addpath('config'); % add config folder

% Load constants from hwe_config.m and assign global variables
params = hwe_config();

global jj kk ll BPVy NN2 f0 dy dz m0 Lx Ubar beta cplx gg Theta0 n_mode HH

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS/VARIABLES %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jj = params.jj;
kk = params.kk;
ll = params.ll;
NN2 = params.NN2;
f0 = params.f0;
dy = params.dy;
dz = params.dz;
m0 = params.m0;
Lx = params.Lx;
beta = params.beta;
cplx = params.cplx;
gg = params.gg;
Theta0 = params.Theta0;
n_mode = params.n_mode;
HH = params.HH;
Ly = params.Ly;

%% Load data from hwe_wave_#.mat
hwe_data_dir = params.hwe_data_dir;
hwe_plot_dir = params.hwe_plot_dir;

hwe_data_filename = params.hwe_data_filename;
hwe_data = fullfile(hwe_data_dir, hwe_data_filename);
load(hwe_data);

%% Define coordinates
ii = 360; % longitude grid points
dx = Lx/ii; % longitude grid spacing
xx = 0.0:360/ii:360; % longitude grid
yy = linspace(45-25, 45+25, jj+1); % latitude grid
zz = linspace(0.0, 10, kk+1); % height grid

%% Parameters from the model 
DeltaT_hor = params.DeltaT_hor; % K
mu = params.mu; % mu parameter (set to 0 for Eady-type)
U0 = params.U0; % mean zonal wind offset
y0 = params.y0; % center y

% Create plots folder if it doesn't exist
if ~exist(hwe_plot_dir, 'dir')
    mkdir(hwe_plot_dir);
end

%% Combined Plot: Background flow for Modified Hoskins-West Eady-type model
fig = figure('units', 'inch', 'position', [1,1,24,8], 'Visible', 'off');

% Subplot 1: Ubar
subplot(1,3,1);
contourf(yy, zz, Ubar', 20, 'linestyle', 'none');
colorbar;
xlabel('Latitude (degrees)');
ylabel('Height (km)');
title('Ubar (m/s)');
set(gca, 'FontSize', 12, 'FontWeight', 'Bold');

% Subplot 2: d(PVbar)/dy interior (analytical computation)
dPVdy_interior = hwe_PV2intgrad(params);

subplot(1,3,2);
contourf(yy, zz, dPVdy_interior', 20, 'linestyle', 'none');
colorbar;
xlabel('Latitude (degrees)');
ylabel('Height (km)');
title('d[PVbar]/dy (interior)');
set(gca, 'FontSize', 12, 'FontWeight', 'Bold');
hold on;
caxis([min(dPVdy_interior(:)) max(dPVdy_interior(:))]); % Dynamic range

% Subplot 3: Boundary PV gradients
[dPVdy_surf, dPVdy_trop] = hwe_PV2bndgrad(params);
dPVdy_beta = beta * ones(1, jj+1);

subplot(1,3,3);
yy_plot = yy;
plot(yy_plot, dPVdy_surf*1e11, 'r-', 'LineWidth', 2); hold on;
plot(yy_plot, dPVdy_trop*1e11, 'b-', 'LineWidth', 2); hold on;
plot(yy_plot, dPVdy_beta*1e11, 'k-', 'LineWidth', 2);
xlabel('Latitude (degrees)');
ylabel('d[PVbar]/dy at surf./trop./beta × 10^{11} (s^{-1})');
legend('(\partial q / \partial y)_{Surf.}','(\partial q / \partial y)_{Trop.}','\beta', 'Location', 'best');
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'Bold');

% Overall title
sgtitle(sprintf('Modified Hoskins-West Eady-type Model''s background flow Delta T = %.0f; U_0 = %.0f', DeltaT_hor, U0));

HW_BG_FILE = fullfile(hwe_plot_dir, 'hwe_background_flow_Eady.png');
saveas(fig, HW_BG_FILE);

%% Plot perturbations (optional, similar to Eady)
n_mode = params.n_mode;
XV = zeros(ll, 1);
XV(:) = eigVec3(:, n_mode);

% Normalize for visualization
gpt_h = XV2field(XV, ii, dx) * f0/gg;
[valuemax, ~] = max(abs(gpt_h(:)));
XV = (10/valuemax) * XV;

% Recompute fields after normalization
gpt_h = XV2field(XV, ii, dx) * f0/gg;
XVz = XV2XVz(XV);
temp = (f0 * HH / 287) * XV2field(XVz, ii, dx);
vg_field = XVx2field(XV, ii, dx); % Meridional wind

%% Cross-section plots at mid-latitude (jj/2+1 ≈ 45°)
% Plot: Meridional Wind Perturbation (v') 
fig4 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
contourf(xx, zz, squeeze(vg_field(:, jj/2+1, :))', -0.0005:0.0001:0.0005, 'linestyle', 'none');
colorbar;
xlabel('Longitude (degrees)');
ylabel('Height (km)');
title(sprintf('Meridional Wind Perturbation (m/s) - Mode %d', n_mode));
set(gca, 'FontSize', 16, 'FontWeight', 'Bold');

HW_MERIDIONAL_CROSS = fullfile(hwe_plot_dir, 'hwe_meridional_wind_perturbation_cross_section_Eady.png');
saveas(fig4, HW_MERIDIONAL_CROSS);

% Plot: Temperature Perturbation (T') 
fig5 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
contourf(xx, zz, squeeze(temp(:, jj/2+1, :))', -0.05:0.01:0.05, 'linestyle', 'none');
colorbar;
xlabel('Longitude (degrees)');
ylabel('Height (km)');
title(sprintf('Temperature Perturbation (K) - Mode %d', n_mode));
set(gca, 'FontSize', 16, 'FontWeight', 'Bold');

HW_TEMP_CROSS = fullfile(hwe_plot_dir, 'hwe_temperature_perturbation_cross_section_Eady.png');
saveas(fig5, HW_TEMP_CROSS);

% Plot: Geopotential height perturbation (for completeness)
fig6 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
contourf(xx, zz, squeeze(gpt_h(:, jj/2+1, :))', 20, 'linestyle', 'none');
colorbar;
xlabel('Longitude (degrees)');
ylabel('Height (km)');
title(sprintf('Perturbation Geopotential Height - Mode %d', n_mode));
set(gca, 'FontSize', 16, 'FontWeight', 'Bold');

HW_GPH_PERT = fullfile(hwe_plot_dir, 'hwe_geopotential_perturbation_Eady.png');
saveas(fig6, HW_GPH_PERT);