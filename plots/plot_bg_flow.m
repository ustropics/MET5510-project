%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_bg_flow.m
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

function plot_bg_flow(yy, zz, jj, kk, Ubar, BPVy, m0, n_mode, fig_path)

    params = cfg(); % load model parameters

    % Define output file name
    file_str = ['background-flow', '_eMode-', num2str(n_mode), '_m0-', num2str(m0), '.png'];

    % Extract boundary PV gradients from BPVy
    qy_surf = BPVy(:, 1) / params.HH;  % surface PV gradient
    qy_trop = -BPVy(:, end) / params.HH;  % tropopause PV gradient

    % Compute scaled PV gradients for plotting
    qy_surf_plot = qy_surf * 1e11;  % (jj+1) x 1 vector
    qy_trop_plot = qy_trop * 1e11;  % (jj+1) x 1 vector
    beta_plot = params.beta * 1e11;  % scalar

    % Create the combined figure
    figure('units', 'inch', 'position', [4, 2, 20, 8], ...
        'Name', sprintf('Background Flow'), 'Visible', 'off');

    % Title for the whole figure
    % annotation('textbox', [0.3, 1, 0.4, 0.1], 'String', ...
    %     ['Hoskins-West Modified Eady Model''s background flow' newline '$\Delta T = 60$; $U_0 = 0$'], ...
    %     'FontSize', 14, 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'Interpreter', 'latex');

    %% Subplot 1: Ubar (mean zonal wind)
    subplot(1, 3, 1);
    contourf(yy, zz, Ubar', 20, 'LineStyle', 'none'); 
    colorbar;
    xlabel('Latitude (degrees)');
    ylabel('Height (km)');
    title('$U_{bar}$', 'Interpreter', 'latex');
    caxis([min(Ubar(:)) max(Ubar(:))]); % set color axis limits

    %% Subplot 2: Interior ∂q/∂y = 0
    subplot(1, 3, 2);
    contourf(BPVy(2:jj,2:kk)', 'LineStyle','none');
    hold on;
    contour(BPVy(2:jj,2:kk)',[0 0], 'LineStyle', '-', 'color', 'white');
    colorbar;
    xlabel('Latitude (degrees)');
    ylabel('Height (km)');
    title('$\frac{\partial \bar{q}}{\partial y} = 0$ (Interior)', 'Interpreter', 'latex');


    %% Subplot 3: Boundary d(PVbar)/dy with beta
    subplot(1, 3, 3);
    hold on;

    % Plot PV gradients at boundaries
    plot(yy, qy_surf_plot, 'r-', 'LineWidth', 2);  % surface
    plot(yy, qy_trop_plot, 'b-', 'LineWidth', 2);  % tropopause
    plot(yy, beta_plot * ones(size(yy)), 'k-', 'LineWidth', 2);  % beta (match yy length)

    xlabel('Latitude (degrees)');
    ylabel('$\frac{d(\overline{PV})}{dy} \times 10^{-11} \, (s^{-1} m^{-1})$', 'Interpreter', 'latex');
    title('$\frac{d(\overline{PV})}{dy}$ at surf./trop./beta (red/blue/black)', 'Interpreter', 'latex');

    % Add equations as annotations
    annotation('textbox', [0.67, 0.5, 0.25, 0.1], 'String', ...
        '$(\frac{\partial q}{\partial y})_{\mathrm{Trop.}} = \frac{f_0^2}{N^2 H} (\frac{\partial \bar{u}}{\partial z})_{\mathrm{Trop.}}$', ...
        'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');
    annotation('textbox', [0.67, 0.3, 0.25, 0.1], 'String', ...
        '$(\frac{\partial q}{\partial y})_{\mathrm{Surf.}} = -\frac{f_0^2}{N^2 H} (\frac{\partial \bar{u}}{\partial z})_{\mathrm{Surf.}}$', ...
        'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter', 'latex');

    % Set y-limits using min and max of the plotted data
    ylim_min = min([min(qy_surf_plot), min(qy_trop_plot), beta_plot]) * 1.2;
    ylim_max = max([max(qy_surf_plot), max(qy_trop_plot), beta_plot]) * 1.2;
    ylim([ylim_min, ylim_max]);

    hold off;

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',15);

    %% Save figure
    outFile = fullfile(fig_path, file_str);
    saveas(gcf, outFile);
    close(gcf);
    
end