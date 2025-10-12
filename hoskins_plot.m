%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: hoskins_plot.m

% Description: Script for loading data from 'HoskinsWest_wave.mat' and 
% generating plots of the Modified Hoskins-West Model's background flow, 
% including Ubar, d(PVbar)/dy interior, and boundary PV gradients, along 
% with perturbation fields (meridional wind, temperature, and geopotential 
% height) for a specified mode, saving results to the 'plots' folder.

% Functions used: 
% - XV2field: Transforms eigenvector to 3D field.
% - XV2XVz: Computes z-derivative of eigenvector.
% - XVx2field: Converts x-derivative of eigenvector to 3D field.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data from HoskinsWest_wave.mat
load('HoskinsWest_wave.mat');

global jj kk ll BPVy NN2 f0 dy dz m0 Lx Ubar beta cplx gg Theta0 HH

%% Define coordinates
ii = 360; % longitude grid points
dx = Lx/ii; % longitude grid spacing
xx = 0.0:360/ii:360; % longitude grid
yy = linspace(45-25, 45+25, jj+1); % latitude grid
zz = linspace(0.0, 10, kk+1); % height grid

%% Parameters from the model (consistent with rossby_main.m)
DeltaT_hor = 60; % K (from PDF page 2)
Lambda = gg * DeltaT_hor / (f0 * Theta0 * Ly); % shear rate
mu = 0.5; % mu parameter (from Modified Hoskins-West, adjust to 1 for page 21)
U0 = 0; % mean zonal wind offset
L_d = sqrt(NN2) * HH / f0; % rossby deformation radius
m_y = 2 * pi / Ly; % meridional wavenumber
gamma = m_y * L_d; % gamma for sinh term
y0 = Ly / 2; % center y

% Create plots folder if it doesn't exist
if ~exist('plots', 'dir')
    mkdir('plots');
end

%% Combined Plot: Background flow for Modified Hoskins-West Model
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
dPVdy_interior = zeros(jj+1, kk+1);
for j = 1:jj+1
    y = (j-1) * dy;
    for k = 1:kk+1
        z = (k-1) * dz;
        if abs(gamma) > 1e-10
            sinh_term = sinh(gamma * z / HH) / sinh(gamma);
            cosh_term = cosh(gamma * z / HH) / sinh(gamma);
        else
            sinh_term = z / HH; % Approximation for small gamma
            cosh_term = 1;
        end

        % analytical d(PVbar)/dy from Modified Hoskins-West Model
        % d(PVbar)/dy = β - (∂²U/∂y²) + (f₀²/N²) (∂²U/∂z²)
        d2U_dy2 = -Lambda * HH * m_y^2 * sinh_term * cos(m_y * (y - y0)); % Second y derivative
        d2U_dz2 = Lambda * HH * (1/HH - mu * 2 * z / HH^2 + gamma * cosh_term * cos(m_y * (y - y0))); % Second z derivative
        dPVdy_interior(j,k) = beta + d2U_dy2 - (f0^2 / NN2) * d2U_dz2;
    end
end

subplot(1,3,2);
contourf(yy, zz, dPVdy_interior', 20, 'linestyle', 'none');
colorbar;
xlabel('Latitude (degrees)');
ylabel('Height (km)');
title('d[PVbar]/dy (interior)');
set(gca, 'FontSize', 12, 'FontWeight', 'Bold');
hold on;
% Add text annotations similar to professor's plot
caxis([min(dPVdy_interior(:)) max(dPVdy_interior(:))]); % Dynamic range

% Subplot 3: Boundary PV gradients
subplot(1,3,3);
dPVdy_surf = zeros(1, jj+1);
dPVdy_trop = zeros(1, jj+1);

for j = 1:jj+1
    y = (j-1) * dy;
    % Surface (z=0) and Tropopause (z=HH) PV gradients
    if abs(gamma) > 1e-10
        sinh_0 = sinh(0) / sinh(gamma); % 0 at surface
        sinh_H = sinh(gamma) / sinh(gamma); % 1 at tropopause
        dU_dz_surf = Lambda * HH * (1 + gamma * cosh(0) * cos(m_y * (y - y0)) / sinh(gamma));
        dU_dz_trop = Lambda * HH * (1 - mu + gamma * cosh(gamma) * cos(m_y * (y - y0)) / sinh(gamma));
    else
        dU_dz_surf = Lambda * (1 + z/HH - mu * (z/HH)^2); % Approx at z=0
        dU_dz_trop = Lambda * (1 - mu); % Approx at z=HH
    end
    f02_over_N2 = f0^2 / NN2;
    dPVdy_surf(j) = -f02_over_N2 * dU_dz_surf;
    dPVdy_trop(j) = f02_over_N2 * dU_dz_trop;
end

dPVdy_beta = beta * ones(1, jj+1);

yy_plot = yy;
plot(yy_plot, dPVdy_surf*1e10, 'r-', 'LineWidth', 2); hold on;
plot(yy_plot, dPVdy_trop*1e10, 'b-', 'LineWidth', 2);
plot(yy_plot, dPVdy_beta*1e10, 'k-', 'LineWidth', 2);
xlabel('Latitude (degrees)');
ylabel('d[PVbar]/dy at surf./trop./beta × 10^{10} (s^{-1})');
legend('(\partial q / \partial y)_{Surf.}','(\partial q / \partial y)_{Trop.}','\beta', 'Location', 'best');
grid on;
set(gca, 'FontSize', 12, 'FontWeight', 'Bold');

% Overall title
sgtitle(sprintf('Modified Hoskins-West Model''s background flow \Delta T = 60; \mu = %.1f; U_0 = 0', mu));

saveas(fig, 'plots/background_flow_HW.png');

%% Plot perturbations (optional, similar to Eady)
n_mode = 7;
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

% Cross-section plots at mid-latitude (jj/2+1 ≈ 45°)
% Plot: Meridional Wind Perturbation (v') 
fig4 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
contourf(xx, zz, squeeze(vg_field(:, jj/2+1, :))', -0.0005:0.0001:0.0005, 'linestyle', 'none');
colorbar;
xlabel('Longitude (degrees)');
ylabel('Height (km)');
title(sprintf('Meridional Wind Perturbation (m/s) - Mode %d', n_mode));
set(gca, 'FontSize', 16, 'FontWeight', 'Bold');
saveas(fig4, 'plots/meridional_wind_perturbation_cross_section_HW.png');

% Plot: Temperature Perturbation (T') 
fig5 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
contourf(xx, zz, squeeze(temp(:, jj/2+1, :))', -0.05:0.01:0.05, 'linestyle', 'none');
colorbar;
xlabel('Longitude (degrees)');
ylabel('Height (km)');
title(sprintf('Temperature Perturbation (K) - Mode %d', n_mode));
set(gca, 'FontSize', 16, 'FontWeight', 'Bold');
saveas(fig5, 'plots/temperature_perturbation_cross_section_HW.png');

% Optional: Geopotential height perturbation (for completeness)
fig6 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
contourf(xx, zz, squeeze(gpt_h(:, jj/2+1, :))', 20, 'linestyle', 'none');
colorbar;
xlabel('Longitude (degrees)');
ylabel('Height (km)');
title(sprintf('Perturbation Geopotential Height - Mode %d', n_mode));
set(gca, 'FontSize', 16, 'FontWeight', 'Bold');
saveas(fig6, 'plots/geopotential_perturbation_mode7_HW.png');