%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_momentum_flux.m
%
% DESCRIPTION: Generates a meridional-height (y-z) contour plot of the 
% zonally-averaged meridional eddy momentum flux <v'u'>. This diagnostic 
% measures momentum convergence/divergence by waves, driving mean flow 
% acceleration (Eliassen-Palm flux divergence).
%
% INPUT:
% - vg: 3D meridional wind perturbation (m/s), size (ii+1, jj+1, kk+1)
% - ug: 3D zonal wind perturbation (m/s), size (ii+1, jj+1, kk+1)
% - m0: Zonal wavenumber (integer)
% - n_mode: Selected eigenmode index
% - fig_path: Directory path for saving the figure
%
% OUTPUT:
% - Saves figure to: fig_path/mean_momentum_flux_eMode-n_m0-m.png
%
% MATH/FUNCTIONS:
% - <v'u'> = mean(vg .* ug, 1)  →  squeeze → (jj+1 × kk+1)
% - Units: (m/s)² = m²/s²
%
% VARIABLES:
% - zvu_calc: Zonally-averaged momentum flux <v'u'> (m²/s²)
% - vg: Meridional wind anomaly
% - ug: Zonal wind anomaly (from -XVy)
% - contourf: Filled contour, no lines
% - colorbar: Scale in m²/s²
% - Figure size: 16×12 inches, font size 20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_momentum_flux(vg, ug, m0, n_mode, fig_path)

    %% Compute zonal mean of momentum flux: <v'u'> = mean(v'u', x)
    zvu_calc = squeeze(mean(vg.*ug,1));

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(zvu_calc', 'LineStyle', 'none');
    colorbar;
    colormap(cmap_RdYlBl(256))
    xlabel('Latitude')
    ylabel('z-level')
    % set(gca, 'xtick', 0:30:360)

    % Set title
    title_str = ['Mean Zonal Flow', ' (zonal wave # = ', ...
     num2str(m0), ', eMode # = ', num2str(n_mode), ')'];

    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

    %% Save figure
    outFile = fullfile(fig_path, ['momentum-flux', '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);
    fprintf('Saving zonal wind plot to: %s\n', outFile); % Debug output
    saveas(gcf, outFile);
    close(gcf);
    
end
