% plot_Eady_results.m
% Loads Eady_wave.mat and plots background flow (Ubar, d(PVbar)/dy), 
% boundary PV gradients (surface/tropopause/beta), and cross-section perturbations
% (meridional wind and temperature) as in PDF page 3
% Saves plots to 'plots' folder with black background

clear all;
close all;

% Load data from Eady_wave.mat
load('Eady_wave.mat');

global jj kk ll BPVy NN2 f0 dy dz m0 Lx Ubar beta cplx gg Theta0 HH

% Define coordinates
ii = 360; % longitude grid points
dx = Lx/ii; % longitude grid spacing
xx = 0.0:360/ii:360; % longitude grid
yy = linspace(45-25, 45+25, jj+1); % latitude grid
zz = linspace(0.0, 10, kk+1); % height grid

% Update Ubar to match ΔT = 60 K, U₀ = 0 (from PDF page 2)
DeltaT_hor = 60; % K (as per PDF page 2)
Lambda = gg * DeltaT_hor / (f0 * Theta0 * Ly); % Recalculate shear
Ubar = zeros(jj+1, kk+1); % Reset Ubar
for j = 1:jj+1
    for k = 1:kk+1
        Ubar(j,k) = Lambda * ((k-1) * dz); % u = Λz, U₀ = 0
    end
end

% Create plots folder if it doesn't exist
if ~exist('plots', 'dir')
    mkdir('plots');
end

%% Plot 1: Ubar (zonal wind) vs. latitude and height (background flow)
fig1 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
contourf(yy, zz, Ubar', 20, 'linestyle', 'none');
colorbar;
xlabel('Latitude (degrees)');
ylabel('Height (km)');
title('Ubar (m/s) - Eady Model Background Flow');
set(gca, 'FontSize', 16, 'FontWeight', 'Bold');
saveas(fig1, 'plots/Ubar_Eady_background.png');

%% Plot 2: Interior d(PVbar)/dy (should be ~0 due to beta=0)
dPVdy_interior = zeros(jj+1, kk+1);
for j = 2:jj
    for k = 1:kk+1
        if j == 1
            dPVdy_interior(j,k) = (Ubar(j+1,k) - Ubar(j,k)) / dy;
        elseif j == jj
            dPVdy_interior(j,k) = (Ubar(j,k) - Ubar(j-1,k)) / dy;
        else
            dPVdy_interior(j,k) = (Ubar(j+1,k) - Ubar(j-1,k)) / (2*dy);
        end
    end
end

fig2 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
contourf(yy, zz, dPVdy_interior', 20, 'linestyle', 'none');
colorbar;
xlabel('Latitude (degrees)');
ylabel('Height (km)');
title('Interior d(PVbar)/dy = 0 - Eady Model');
set(gca, 'FontSize', 16, 'FontWeight', 'Bold');
caxis([-1e-10 1e-10]);
saveas(fig2, 'plots/dPVdy_interior_Eady.png');

%% Plot 3: Boundary PV gradients (surface, tropopause, and beta)
dPVdy_surf = zeros(1, jj+1);
dPVdy_trop = zeros(1, jj+1);

Lambda_surf = (Ubar(1,2) - Ubar(1,1)) / dz; % Surface shear
Lambda_trop = (Ubar(1,kk+1) - Ubar(1,kk)) / dz; % Tropopause shear

f02_over_N2 = f0^2 / NN2;
for j = 1:jj+1
    dPVdy_surf(j) = -f02_over_N2 * Lambda_surf;
    dPVdy_trop(j) = +f02_over_N2 * Lambda_trop;
end

dPVdy_beta = zeros(1, jj+1);

fig3 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
yy_plot = yy;
plot(yy_plot, dPVdy_surf*1e10, 'b-', 'LineWidth', 2); hold on;
plot(yy_plot, dPVdy_trop*1e10, 'r-', 'LineWidth', 2);
plot(yy_plot, dPVdy_beta*1e10, 'k-', 'LineWidth', 2);
xlabel('Latitude (degrees)');
ylabel('d(PVbar)/dy × 10^{-10} (s^{-1})');
title('Eady Model Boundary PV Gradients');
legend('Surface (∂q/∂y)_{surf}', 'Tropopause (∂q/∂y)_{trop}', 'β = 0 (f-plane)', ...
       'Location', 'best');
grid on;
set(gca, 'FontSize', 14, 'FontWeight', 'Bold');
saveas(fig3, 'plots/boundary_PV_gradients_Eady.png');

%% Print boundary values for verification
fprintf('Eady Model Boundary PV Gradients (ΔT = 60 K, U₀ = 0):\n');
fprintf('Surface PV gradient: %.2e s^-1\n', dPVdy_surf(1));
fprintf('Tropopause PV gradient: %.2e s^-1\n', dPVdy_trop(1));
fprintf('Beta (f-plane): %.2e s^-1 m^-1\n', beta);
fprintf('f₀²/N² = %.2e s^-2\n', f02_over_N2);
fprintf('Vertical shear Λ = %.2f s^-1\n', Lambda_surf);

%% Plot perturbations (including cross-sections from PDF page 3)
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
% Plot 4: Meridional Wind Perturbation (v') - Matches Figure 7.7(a)
fig4 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
contourf(xx, zz, squeeze(vg_field(:, jj/2+1, :))', -0.0005:0.0001:0.0005, 'linestyle', 'none');
colorbar;
xlabel('Longitude (degrees)');
ylabel('Height (km)');
title(sprintf('Meridional Wind Perturbation (m/s) - Mode %d', n_mode));
set(gca, 'FontSize', 16, 'FontWeight', 'Bold');
saveas(fig4, 'plots/meridional_wind_perturbation_cross_section.png');

% Plot 5: Temperature Perturbation (T') - Matches Figure 7.7(b)
fig5 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
contourf(xx, zz, squeeze(temp(:, jj/2+1, :))', -0.05:0.01:0.05, 'linestyle', 'none');
colorbar;
xlabel('Longitude (degrees)');
ylabel('Height (km)');
title(sprintf('Temperature Perturbation (K) - Mode %d', n_mode));
set(gca, 'FontSize', 16, 'FontWeight', 'Bold');
saveas(fig5, 'plots/temperature_perturbation_cross_section.png');

% Optional: Geopotential height perturbation (for completeness)
fig6 = figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
contourf(xx, zz, squeeze(gpt_h(:, jj/2+1, :))', 20, 'linestyle', 'none');
colorbar;
xlabel('Longitude (degrees)');
ylabel('Height (km)');
title(sprintf('Perturbation Geopotential Height - Mode %d', n_mode));
set(gca, 'FontSize', 16, 'FontWeight', 'Bold');
saveas(fig6, 'plots/geopotential_perturbation_mode7.png');