%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_vertical_heat_flux.m
%
% DESCRIPTION: Generates a meridional-height (y-z) contour plot of the 
% zonally-averaged vertical eddy heat flux <w'T'>. This term drives diabatic 
% heating/cooling and is key in the QG thermodynamic energy equation and 
% wave-mean flow interaction.
%
% INPUT:
% - wfield: 3D vertical velocity (m/s), size (ii+1, jj+1, kk+1)
% - temp: 3D temperature perturbation (K), size (ii+1, jj+1, kk+1)
% - m0: Zonal wavenumber (integer)
% - n_mode: Selected eigenmode index
% - fig_path: Directory path for saving the figure
%
% OUTPUT:
% - Saves figure to: fig_path/mean_vertical_heat_flux_eMode-n_m0-m.png
%
% MATH/FUNCTIONS:
% - <w'T'> = mean(wfield .* temp, 1)  →  squeeze → (jj+1 × kk+1)
% - Units: (m/s)·K = m·K/s
%
% VARIABLES:
% - zwt_calc: Zonally-averaged vertical heat flux <w'T'> (m·K/s)
% - wfield: Vertical velocity from omega equation
% - temp: Temperature anomaly
% - contourf: Filled contour, no lines
% - colorbar: Scale in m·K/s
% - Figure size: 16×12 inches, font size 20

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_vertical_hflux(wfield, temp, m0, n_mode, fig_path)

    %% Compute zonal mean of vertical heat flux: <w'T'> = mean(w'T', x)
    zwt_calc = squeeze(mean(wfield.*temp,1));

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(zwt_calc', 'LineStyle', 'none');
    colorbar;
    colormap('turbo')
    xlabel('Latitude')
    ylabel('z-level')
    % set(gca, 'xtick', 0:30:360)

    % Set title
    title_str = ['Vertical Heat Flux ', '(zonal wave # = ', ...
     num2str(m0), ', eMode # = ', num2str(n_mode), ')'];

    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

        
    %% Save figure
    outFile = fullfile(fig_path, ['vertical-hflux', '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);
    fprintf('Saving zonal wind plot to: %s\n', outFile); % Debug output
    saveas(gcf, outFile);
    close(gcf);
    
end
