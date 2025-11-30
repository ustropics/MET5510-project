%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_dpvdym_bnd.m

% DESCRIPTION: Plots meridional PV gradient (∂q̄/∂y) at surface, mid-level,
%              and tropopause (scaled ×10¹² s⁻¹) along with the planetary
%              vorticity gradient β (scaled) as line plots vs. latitude.

% INPUT:
%   yy       - Latitude coordinates (degrees), size (jj+1)
%   BPVy     - Meridional PV gradient field (jj+1 x kk+1 array, s⁻¹)
%   beta     - Planetary vorticity gradient (s⁻¹ m⁻¹)
%   kk       - Number of vertical levels
%   m0       - Zonal wavenumber for title/filename
%   n_mode   - Mode number for title/filename
%   fig_path - Directory path for saving figure

% OUTPUT:
%   Saves: fig_path/dpvdym-boundaries_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_dpvdym_bnd(yy, BPVy, beta, kk, m0, n_mode, fig_path)
    
    %% --------------------------------------------------------------------
    %% 1. Extract and scale boundary PV gradients
    %% --------------------------------------------------------------------
    dpvdy_surface = BPVy(:, 1) * 1e12;           % surface
    dpvdy_mid     = BPVy(:, floor(kk/2)+1) * 1e12; % mid-troposphere
    dpvdy_top     = BPVy(:, end) * 1e12;         % tropopause

    R_earth = 6.371e6;
    beta_scaled = beta * 1e12 * R_earth;         % scale to match units
    
    fprintf('\nBoundary meridional gradient of background potential vorticity:\n')
    fprintf('Boundary PV gradients (×10⁻¹² s⁻¹):\n');
    fprintf('  Surface:  max = %.2f, min = %.2e\n', max(dpvdy_surface), min(dpvdy_surface));
    fprintf('  Mid:      max = %.2f, min = %.2f\n', max(dpvdy_mid), min(dpvdy_mid));
    fprintf('  Top:      max = %.2e, min = %.2e\n', max(dpvdy_top), min(dpvdy_top));

    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------
    figure('units','inch','position',[4,2,18,14],'Visible','off');
    hold on;

    plot(yy, dpvdy_surface, 'b-', 'LineWidth', 2, 'DisplayName', 'Surface');
    plot(yy, dpvdy_mid,     'r-', 'LineWidth', 2, 'DisplayName', 'Mid-Level');
    plot(yy, dpvdy_top,     'g-', 'LineWidth', 2, 'DisplayName', 'Tropopause');
    plot(yy, beta_scaled*ones(size(yy)), 'k--', 'LineWidth', 2, 'DisplayName', 'β (scaled)');

    xlabel('Latitude (degrees)');
    ylabel('$\frac{\partial \bar{q}}{\partial y}$ ($\times 10^{-12}$ s$^{-1}$)', 'Interpreter','latex');
    title_str = ['d(PVbar)/dy at Boundaries with Beta', newline ...
    'zonal wave # = ', num2str(m0), ...
    ', eMode # = ', num2str(n_mode)];

    title(title_str)
    legend('Location','best');
    grid on;

    set(findall(gcf,'-property','FontSize'),'FontSize',20);

    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------
    outFile = fullfile(fig_path, ['dpvdym-boundaries', ...
                 '_eMode-', num2str(n_mode), ...
                 '_m0-', num2str(m0), '.png']);
    
    fprintf('Saving boundary PV gradient plot to: %s\n', outFile);
    saveas(gcf, outFile);
    close(gcf);
end