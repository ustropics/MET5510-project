%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_meridional_hflux.m
%
% DESCRIPTION: Generates a meridional-height (y-z) contour plot of the 
% zonally-averaged meridional eddy heat flux <v'T'> in the quasi-geostrophic 
% model. This diagnostic quantifies poleward heat transport by baroclinic 
% waves, critical for understanding energy conversion and wave maintenance.
%
% INPUT:
% - vg: 3D meridional wind perturbation (m/s), size (ii+1, jj+1, kk+1)
% - temp: 3D temperature perturbation (K), size (ii+1, jj+1, kk+1)
% - m0: Zonal wavenumber (integer)
% - n_mode: Selected eigenmode index
% - fig_path: Directory path for saving the figure
%
% OUTPUT:
% - Saves figure to: fig_path/mean_meridional_heat_flux_eMode-n_m0-m.png
%
% MATH/FUNCTIONS:
% - <v'T'> = mean(vg .* temp, 1)  →  squeeze → (jj+1 × kk+1) matrix
% - Units: (m/s)·K = m·K/s
%
% VARIABLES:
% - zvt_calc: Zonally-averaged meridional heat flux <v'T'> (m·K/s)
% - vg: Meridional geostrophic wind anomaly (from XVx)
% - temp: Temperature anomaly (from XVz)
% - contourf: Filled contour plot with no lines
% - colorbar: Shows scale in m·K/s
% - Figure size: 16×12 inches, font size 20
% - saveas: PNG format, automatically closed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_meridional_hflux(vg, temp, m0, n_mode, fig_path)

    %% Compute zonal mean of meridional heat flux: <v'T'> = mean(v'T', x)
    zvt_calc = squeeze(mean(vg.*temp,1));

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(zvt_calc', 'LineStyle', 'none');
    colorbar;
    colormap('turbo');
    xlabel('Latitude')
    ylabel('z-level')
    % set(gca, 'xtick', 0:30:360)

    % Set title
    title_str = ['Meridional Heat Flux', ' (zonal wave # = ', ...
     num2str(m0), ', eMode # = ', num2str(n_mode), ')'];
    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

        
    %% Save figure
    outFile = fullfile(fig_path, ['meridional-hflux', '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);
    fprintf('Saving zonal wind plot to: %s\n', outFile); % Debug output
    saveas(gcf, outFile);
    close(gcf);
    
end
