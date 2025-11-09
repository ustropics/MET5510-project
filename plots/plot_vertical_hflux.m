%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_vertical_hflux.m

% DESCRIPTION: Generates a latitude-height (y-z) filled contour plot of the
%              zonally-averaged vertical heat flux <w'T'>. This diagnostic
%              drives diabatic heating/cooling in the QG thermodynamic
%              equation. The figure is saved as a PNG with wavenumber and mode.

% INPUT:
%   wfield   - 3D vertical velocity (ii+1 x jj+1 x kk+1 array, m/s)
%   temp     - 3D temperature perturbation (ii+1 x jj+1 x kk+1 array, K)
%   m0       - zonal wavenumber for title/filename
%   n_mode   - mode number for title/filename
%   fig_path - directory path for saving figure

% OUTPUT:
%   Saves: fig_path/vertical-hflux_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_vertical_hflux(zz, yy, wfield, temp, m0, n_mode, fig_path)

    %% --------------------------------------------------------------------
    %% 1. Extract the data and compute limits
    %% --------------------------------------------------------------------

    %% Compute zonal mean of vertical heat flux: <w'T'> = mean(w'T', x)
    data = squeeze(mean(wfield.*temp,1));

    % get the maximum and minimum values
    vmin = min(data(:));   % absolute minimum
    vmax = max(data(:));   % absolute maximum
    
    % print maximum and minimum values
    fprintf('\nZonal mean of vertical heat flux:\n')
    fprintf('Maximum Value: %.3e and Minimum Value: %.3e\n', vmax, vmin)
    
    % sets the +/- value to add to contourf (0.2 = ~20%)  
    step = 0.2;                              
    vmin = floor(vmin/step)*step;               
    vmax = ceil (vmax/step)*step;

    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------

    % Create the figure and set it's size [left, bottom, width, height]
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(yy, zz, data');
    hold on
    contour(yy, zz, data', 'LineColor', 'k', 'LineStyle', '-');
    hold off

    colormap(flipud(cmap_spectral(256))); % set our colormap choice
    colorbar;

    % caxis([vmin vmax]); % sets limits from section 1
    
    % Set x and y-label text as well as tick spacing for grid labels
    xlabel('Latitude (degrees)')
    ylabel('Height (km)')
    % set(gca, 'xtick', 0:30:360)

    % Create title string with input variables and set it
    title_str = ['Vertical Heat Flux ', newline ...
        'zonal wave # = ', num2str(m0), ...
        ', eMode # = ', num2str(n_mode)];

    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------

    % Create filename to save our figure and set output directory
    outFile = fullfile(fig_path, ['vertical-hflux', ...
        '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);

    fprintf('Saving vertical heat flux plot to: %s\n', outFile);

    % Comment out 'close(gcf)' if you want to see the plot and save it
    saveas(gcf, outFile);
    close(gcf);
    
end
