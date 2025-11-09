%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_bg_flow.m

% DESCRIPTION: Generates a three-panel figure visualizing the background state
%              of the Hoskins-West Modified Eady (HWME) model:
%              1) mean zonal wind Ubar(y,z),
%              2) interior meridional PV gradient ∂q̄/∂y = 0,
%              3) boundary PV gradients at surface/tropopause with β.
%              Useful for verifying thermal wind balance and PV structure.
%              The figure is saved as a PNG including wavenumber and mode number
%              in title and filename.

% INPUT:
%   yy       - Latitude coordinates (degrees), size (jj+1)
%   zz       - Height coordinates (km), size (kk+1)
%   jj       - Number of meridional grid points (interior-adjacent)
%   kk       - Number of vertical grid points (including boundaries)
%   Ubar     - Background zonal wind (jj+1 x kk+1 array, m/s)
%   BPVy     - Meridional PV gradient field (jj+1 x kk+1 array, s⁻¹)
%   m0       - Zonal wavenumber for title/filename
%   n_mode   - Mode number for title/filename
%   fig_path - Directory path for saving figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function combined_bg_flow(yy, zz, jj, kk, Ubar, BPVy, m0, n_mode, fig_path)
    
    %% --------------------------------------------------------------------
    %% 1. Load parameters and extract boundary PV gradients
    %% --------------------------------------------------------------------
    params = cfg(); % load model parameters

    % Extract boundary PV gradients from BPVy
    qy_surf = BPVy(:, 1) / params.HH;      % surface PV gradient
    qy_trop = -BPVy(:, end) / params.HH;   % tropopause PV gradient (sign flip)

    % Scale for plotting (×10¹¹)
    qy_surf_plot = qy_surf * 1e11;         % (jj+1) vector
    qy_trop_plot = qy_trop * 1e11;         % (jj+1) vector
    beta_plot    = params.beta * 1e11;     % scalar

    % Compute min/max for y-limits
    qy_min = min([min(qy_surf_plot), min(qy_trop_plot), beta_plot]);
    qy_max = max([max(qy_surf_plot), max(qy_trop_plot), beta_plot]);

    fprintf('\nBoundary PV gradient values from BPVy:\n')
    fprintf('Maximum Values: %.2e and Minimum Values: %.2e\n', qy_max, qy_min)

    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------
    figure('units','inch','position',[4,2,20,8],'Visible','off');

    %% --------------------------------------------------------------------
    %% 3.1 Subplot 1: Ubar (mean zonal wind)
    %% --------------------------------------------------------------------
    subplot(1, 3, 1);
    contourf(yy, zz, Ubar', 20, 'LineStyle','none');
    colorbar;
    xlabel('Latitude (degrees)');
    ylabel('Height (km)');
    title('$U_{\mathrm{bar}}$ (m s$^{-1}$)', 'Interpreter','latex');
    caxis([min(Ubar(:)) max(Ubar(:))]);

    %% --------------------------------------------------------------------
    %% 3.2 Subplot 2: Interior ∂q̄/∂y = 0
    %% --------------------------------------------------------------------
    subplot(1, 3, 2);
    contourf(yy(2:jj), zz(2:kk), BPVy(2:jj,2:kk)', 'LineStyle','none');
    hold on;
    contour(yy(2:jj), zz(2:kk), BPVy(2:jj,2:kk)', [0 0], 'k-', 'LineWidth', 0.8);
    colorbar;
    xlabel('Latitude (degrees)');
    ylabel('Height (km)');
    title('$\frac{\partial \bar{q}}{\partial y} = 0$ (Interior)', 'Interpreter','latex');
    hold off;

    %% --------------------------------------------------------------------
    %% 3.3 Subplot 3: Boundary PV gradients with β
    %% --------------------------------------------------------------------
    subplot(1, 3, 3);
    hold on;
    plot(yy, qy_surf_plot, 'r-', 'LineWidth', 2);   % surface
    plot(yy, qy_trop_plot, 'b-', 'LineWidth', 2);   % tropopause
    plot(yy, beta_plot*ones(size(yy)), 'k-', 'LineWidth', 2); % beta

    xlabel('Latitude (degrees)');
    ylabel('$\frac{d(\overline{\mathrm{PV}})}{dy} \times 10^{-11}$ (s$^{-1}$ m$^{-1}$)', 'Interpreter','latex');
    title('$\frac{d(\overline{PV})}{dy}$ at surf./trop./beta (red/blue/black)', 'Interpreter', 'latex');

    % Add equations as annotations
    annotation('textbox', [0.67, 0.55, 0.25, 0.1], 'String', ...
        '$(\frac{\partial q}{\partial y})_{\mathrm{Trop.}} = \frac{f_0^2}{N^2 H} (\frac{\partial \bar{u}}{\partial z})_{\mathrm{Trop.}}$', ...
        'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter','latex');
    annotation('textbox', [0.67, 0.35, 0.25, 0.1], 'String', ...
        '$(\frac{\partial q}{\partial y})_{\mathrm{Surf.}} = -\frac{f_0^2}{N^2 H} (\frac{\partial \bar{u}}{\partial z})_{\mathrm{Surf.}}$', ...
        'FontSize', 10, 'EdgeColor', 'none', 'HorizontalAlignment', 'center', 'Interpreter','latex');

    % Set y-limits with padding
    ylim([qy_min*1.2, qy_max*1.2]);
    hold off;

    %% --------------------------------------------------------------------
    %% 4. Global formatting
    %% --------------------------------------------------------------------
    set(findall(gcf,'-property','FontSize'),'FontSize',15);

    %% --------------------------------------------------------------------
    %% 5. Save figure
    %% --------------------------------------------------------------------
    outFile = fullfile(fig_path, ['background-flow', ...
                 '_eMode-', num2str(n_mode), ...
                 '_m0-', num2str(m0), '.png']);
    
    fprintf('\nSaving three-panel background state plot to: %s\n', outFile);
    saveas(gcf, outFile);
    close(gcf);
end