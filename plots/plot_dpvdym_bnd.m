%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_background_flow.m
%
% DESCRIPTION: Generates a three-panel figure visualizing the background state 
% of the Hoskins-West Modified Eady (HWME) model: mean zonal wind Ubar(y,z), 
% interior meridional PV gradient d(PVbar)/dy = 0, and boundary d(PVbar)/dy 
% with beta. Useful for verifying thermal wind balance and PV structure.
%
% INPUT:
% - yy: 1D array of latitude coordinates (degrees), length jj+1
% - zz: 1D array of height coordinates (km), length kk+1
% - jj: Number of meridional grid points (interior-adjacent)
% - kk: Number of vertical grid points (including boundaries)
% - Ubar: 2D mean zonal wind (m/s), size (jj+1, kk+1)
% - BPVy: 2D meridional PV gradient (s⁻¹), size (jj+1, kk+1)
% - m0: Zonal wavenumber (integer)
% - n_mode: Selected eigenmode index
% - fig_path: Directory path for saving the figure
%
% OUTPUT:
% - Saves figure to: fig_path/background-flow_eMode-n_m0-m.png
%
% MATH/FUNCTIONS:
% - qy_surf = BPVy(:,1) / H,   qy_trop = -BPVy(:,end) / H
% - Scaled for plotting: ×10¹¹ s⁻¹ m⁻¹
%
% VARIABLES:
% - qy_surf, qy_trop: Boundary PV gradients (s⁻¹ m⁻¹)
% - qy_surf_plot, qy_trop_plot: Scaled versions (×10¹¹)
% - beta_plot: Scaled beta (constant line)
% - contourf: Filled contours for Ubar and interior d(PVbar)/dy
% - plot: Line plots for boundaries
% - annotation: Adds boundary equations
% - Figure size: 20×8 inches, font size 15

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_dpvdym_bnd(yy, BPVy, beta, kk, m0, n_mode, fig_path)

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

    title_str = ['d(PVbar)/dy at Boundaries with Beta (', ...
    'zonal wave # = ', num2str(m0), ...
    ', eMode # = ', num2str(n_mode)];

    title(title_str)
    legend('Location', 'best');
    grid on;

    %% Save figure
    outFile = fullfile(fig_path, ['dpvdym-boundaries', '_eMode-', num2str(n_mode), '_m0-', num2str(m0), '.png']);
    saveas(gcf, outFile);
    close(gcf);

end