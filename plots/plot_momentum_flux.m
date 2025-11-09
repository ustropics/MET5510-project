%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_momentum_flux.m

% DESCRIPTION: Generates a latitude-height (y-z) filled contour plot of the
%              zonally-averaged meridional eddy momentum flux <v'u'>. This
%              diagnostic drives mean flow acceleration via EP flux divergence.
%              The figure is saved as a PNG with wavenumber and mode number.

% INPUT:
%   vg       - 3D meridional wind perturbation (ii+1 x jj+1 x kk+1 array, m/s)
%   ug       - 3D zonal wind perturbation (ii+1 x jj+1 x kk+1 array, m/s)
%   m0       - zonal wavenumber for title/filename
%   n_mode   - mode number for title/filename
%   fig_path - directory path for saving figure

% OUTPUT:
%   Saves: fig_path/momentum_flux_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_momentum_flux(zz, yy, vg, ug, m0, n_mode, fig_path)

    %% --------------------------------------------------------------------
    %% 1. Extract the data and compute limits
    %% --------------------------------------------------------------------

    % Get 2d data for zonal mean of momentum flux: <v'u'> = mean(v'u', x)
    data = squeeze(mean(vg.*ug,1));

    % get the maximum and minimum values
    vmin = min(data(:));   % absolute minimum
    vmax = max(data(:));   % absolute maximum
    
    % print maximum and minimum values
    fprintf('\nZonal mean of momentum flux:\n')
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
    
    colormap(flipud(cmap_twilight(256))); % set our colormap choice
    colorbar;

    % Set x and y-label text as well as tick spacing for grid labels   
    xlabel('Latitude (degrees)')
    ylabel('Height (km)')
    set(gca, 'xtick', 0:30:360)

    % Create title string with input variables and set it
    title_str = ['Momentum Flux', newline ...
        'zonal wave # = ', num2str(m0), ...
        ', eMode # = ', num2str(n_mode)];

    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------

    % Create filename to save our figure and set output directory
    outFile = fullfile(fig_path, ['momentum_flux', ...
        '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);

    fprintf('Saving momentum flux plot to: %s\n', outFile);

    % Comment out 'close(gcf)' if you want to see the plot and save it
    saveas(gcf, outFile);
    close(gcf);
    
end
