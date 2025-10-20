%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: rossby_dpvdym_boundaries.m

% Description: Plots d(PVbar)/dy at the surface, troposphere mid-level, and with
% beta effect at boundaries, with latitude on the x-axis and d(PVbar)/dy values
% on the y-axis, and saves it as an image.

% Input:
% - yy: Latitude coordinates (degrees)
% - BPVy: Meridional gradient of background PV (jj+1 x kk+1 array, s^-1)
% - beta: Beta parameter (m^-1 s^-1)
% - kk: Number of height grid points

% Output:
% - Saves plot to 'output/plots/dpvdym_boundaries.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_dpvdym_boundaries(yy, BPVy, beta, kk, model, m0)
figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')

% Extract d(PVbar)/dy at surface (k=1), mid-level (k=kk/2+1), and top (k=kk+1)
dpvdy_surface = BPVy(:, 1) * 1e12; % Scale to 10^-12 s^-1
dpvdy_mid = BPVy(:, floor(kk/2)+1) * 1e12;
dpvdy_top = BPVy(:, end) * 1e12;

% Convert beta to 10^-12 s^-1 for consistency (assuming Lx as length scale)
beta_scaled = beta * 1e12 * (6.37e6); % Approximate Earth radius in meters

% Plot
plot(yy, dpvdy_surface, 'b-', 'LineWidth', 2, 'DisplayName', 'Surface');
hold on;
plot(yy, dpvdy_mid, 'r-', 'LineWidth', 2, 'DisplayName', 'Troposphere Mid-Level');
plot(yy, dpvdy_top, 'g-', 'LineWidth', 2, 'DisplayName', 'Top Boundary');
plot(yy, beta_scaled * ones(size(yy)), 'k--', 'LineWidth', 2, 'DisplayName', 'Beta Effect');
hold off;

xlabel('Latitude (degrees)')
ylabel('d(PVbar)/dy (10^-12 s^-1)')
title('d(PVbar)/dy at Boundaries with Beta');
legend('Location', 'best');
grid on;

% Save plot
saveas(gcf, ['output', filesep, 'figures', filesep, 'dpvdym_boundaries_', model, '_m0_', num2str(m0), '.png']);
close(gcf);
end