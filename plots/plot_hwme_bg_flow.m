%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_hwme_bg_flow.m

% DESCRIPTION: Plots the Hoskins-West Modified Eady model's background flow:
% Ubar, interior d(PVbar)/dy = 0, and boundary d(PVbar)/dy with beta.
%
% INPUT:
%   yy - Latitude coordinates (degrees)
%   zz - Height coordinates (km)
%   Ubar - Mean zonal wind field (jj+1 x kk+1 array, m/s)
%   BPVy - Potential vorticity gradient field (jj+1 x kk+1 array, m^-1 s^-1)
%   model - Model name (string)
%   m0 - Wave number
%
% OUTPUT:
%   Saves combined plot to 'output/plots/hwme_background_flow.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_hwme_bg_flow(yy, zz, jj, kk, Ubar, BPVy, model, m0)

    % Clear any variable named 'title' to avoid shadowing
    clear title;

    % Load parameters to access global or config values
    params = hwme_config();

    % Extract boundary PV gradients from BPVy (assuming first and last rows are boundaries)
    qy_surf = BPVy(:, 1) / params.HH;  % Surface PV gradient (scaled by 1/HH for s^-1 m^-1)
    qy_trop = -BPVy(:, end) / params.HH;  % Tropopause PV gradient (flipped sign to match convention)

    % Compute scaled PV gradients for plotting (×10^11 to match typical units)
    qy_surf_plot = qy_surf * 1e11;  % (jj+1) x 1 vector
    qy_trop_plot = qy_trop * 1e11;  % (jj+1) x 1 vector
    beta_plot = params.beta * 1e11;  % Scalar (ensure beta is a scalar)

    % Create the combined figure
    figure('units', 'inch', 'position', [4, 2, 20, 8], ...
        'Name', sprintf('%s Model Background Flow', model), 'Visible', 'off');

    % Title for the whole figure
    % annotation('textbox', [0.3, 1, 0.4, 0.1], 'String', ...
    %     ['Hoskins-West Modified Eady Model''s background flow' newline '$\Delta T = 60$; $U_0 = 0$'], ...
    %     'FontSize', 14, 'HorizontalAlignment', 'center', 'EdgeColor', 'none', 'Interpreter', 'latex');

    %% Subplot 1: Ubar (mean zonal wind)
    subplot(1, 3, 1);
    contourf(yy, zz, Ubar', 20, 'LineStyle', 'none');  % Transpose to match latitude vs height
    colorbar;
    xlabel('Latitude (degrees)');
    ylabel('Height (km)');
    title('$U_{bar}$', 'Interpreter', 'latex');
    caxis([min(Ubar(:)) max(Ubar(:))]);  % Adjust color limits based on data

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
    plot(yy, qy_surf_plot, 'r-', 'LineWidth', 2);  % Surface
    plot(yy, qy_trop_plot, 'b-', 'LineWidth', 2);  % Tropopause
    plot(yy, beta_plot * ones(size(yy)), 'k-', 'LineWidth', 2);  % Beta (broadcast to match yy length)

    xlabel('Latitude (degrees)');
    ylabel('$\frac{d(\overline{PV})}{dy} \times 10^{-11} \, (s^{-1} m^{-1})$', 'Interpreter', 'latex');
    title('$\frac{d(\overline{PV})}{dy}$ at surf./trop./beta (red/blue/black)', 'Interpreter', 'latex');

    % Add equations as annotations with proper TeX syntax
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

    saveas(gcf, ['output', filesep, 'figures', filesep, 'hwme_background_', model, '_m0_', num2str(m0), '.png']);
    close(gcf);

end