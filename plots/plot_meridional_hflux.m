%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_meridional_hflux.m

% DESCRIPTION: Generates a latitude-height (y-z) filled contour plot of the
%              zonally-averaged meridional eddy heat flux <v'T'>. This
%              diagnostic quantifies poleward heat transport by baroclinic
%              waves, critical for energy conversion and wave maintenance.

% INPUT:
%   vg       - 3D meridional wind perturbation (ii+1 x jj+1 x kk+1 array, m/s)
%   temp     - 3D temperature perturbation (ii+1 x jj+1 x kk+1 array, K)
%   m0       - zonal wavenumber for title/filename
%   n_mode   - mode number for title/filename
%   fig_path - directory path for saving figure

% OUTPUT:
%   Saves: fig_path/meridional-hflux_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_meridional_hflux(zz, yy, vg, temp, m0, n_mode, fig_path)

    %% --------------------------------------------------------------------
    %% 1. Extract the data and compute limits
    %% --------------------------------------------------------------------

    %% Compute zonal mean of meridional heat flux: <v'T'> = mean(v'T', x)
    data = squeeze(mean(vg.*temp,1));

    % get the maximum and minimum values
    vmin = min(data(:));   % absolute minimum
    vmax = max(data(:));   % absolute maximum
    
    % print maximum and minimum values
    fprintf('\nZonal mean of meridional heat flux:\n')
    fprintf('Maximum Value: %.2e and Minimum Value: %.2e\n', vmax, vmin)
    
    % sets the +/- value to add to contourf (0.2 = ~20%)  
    step = 0.2;                              
    vmin = floor(vmin/step)*step;               
    vmax = ceil (vmax/step)*step;

    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------    
   
    % Create the figure and set it's size [left, bottom, width, height]
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(yy, zz, data', 'LineStyle', 'none');
    hold on
    contour(yy, zz, data', 'LineColor', 'k', 'LineStyle', '-');
    hold off
    
    colormap(flipud(cmap_spectral(256))); % set our colormap choice
    colorbar;

    % Set x and y-label text as well as tick spacing for grid labels    
    xlabel('Latitude (degree)')
    ylabel('Height (km)')
    % set(gca, 'xtick', 0:30:360)

    % Create title string with input variables and set it
    title_str = ['Meridional Heat Flux', newline ...
        'zonal wave # = ', num2str(m0), ... 
        ', eMode # = ', num2str(n_mode)];

    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);
        
    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------

    % Create filename to save our figure and set output directory
    outFile = fullfile(fig_path, ['meridional-hflux', ...
        '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);

    fprintf('Saving meridional heat flux plot to: %s\n', outFile);

    % Comment out 'close(gcf)' if you want to see the plot and save it    
    saveas(gcf, outFile);
    close(gcf);
    
end
