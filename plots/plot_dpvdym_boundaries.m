%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_dpvdym_boundaries.m

% DESCRIPTION: Plots d(PVbar)/dy at the surface, troposphere mid-level, and with
% beta effect at boundaries, with latitude on the x-axis and d(PVbar)/dy values
% on the y-axis, and saves it as an image.

% INPUT:
% - yy: Latitude coordinates (degrees)
% - BPVy: Meridional gradient of background PV (jj+1 x kk+1 array, s^-1)
% - beta: Beta parameter (m^-1 s^-1)
% - kk: Number of height grid points

% OUTPUT:
% - Saves plot to 'output/figures/dpvdym_boundaries.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_dpvdym_boundaries(yy, BPVy, beta, kk, model, m0, n_mode, fig_path)

    %% Extract d(PVbar)/dy at surface (k=1), mid-level (k=kk/2+1), and top (k=kk+1)
    dpvdy_surface = BPVy(:, 1) * 1e12; % Scale to 10^-12 s^-1
    dpvdy_mid = BPVy(:, floor(kk/2)+1) * 1e12;
    dpvdy_top = BPVy(:, end) * 1e12;

    % Convert beta to 10^-12 s^-1 for consistency (assuming Lx as length scale)
    beta_scaled = beta * 1e12 * (6.37e6); % Approximate Earth radius in meters

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    plot(yy, dpvdy_surface, 'b-', 'LineWidth', 2, 'DisplayName', 'Surface');
    hold on;
    plot(yy, dpvdy_mid, 'r-', 'LineWidth', 2, 'DisplayName', 'Troposphere Mid-Level');
    plot(yy, dpvdy_top, 'g-', 'LineWidth', 2, 'DisplayName', 'Top Boundary');
    plot(yy, beta_scaled * ones(size(yy)), 'k--', 'LineWidth', 2, 'DisplayName', 'Beta Effect');
    hold off;

    xlabel('Latitude (degrees)')
    ylabel('d(PVbar)/dy (10^-12 s^-1)')

    title_str = ['d(PVbar)/dy at Boundaries with Beta', ...
    ', wave number = ', num2str(m0), ...
    ', eMode = ', num2str(n_mode)];

    legend('Location', 'best');
    grid on;

    %% Save figure
    outFile = fullfile(fig_path, [model, '_dpvdym_boundaries_', '_nmode-', num2str(n_mode), '_m0-', num2str(m0), '.png']);
    saveas(gcf, outFile);
    close(gcf);

end