%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_temp.m

% DESCRIPTION: Plots temperature perturbation (filled contour) at a
%              specified vertical level across longitude and latitude.

% INPUT:
%   xx       - longitude coordinates (degrees)
%   yy       - latitude coordinates (degrees)
%   temp     - temperature field (ii+1 x jj+1 x kk+1 array, K)
%   hlevel   - vertical level index (e.g., 1, 25, 51)
%   m0       - zonal wavenumber for title/filename
%   n_mode   - mode number for title/filename
%   fig_path - directory path for saving figure

% OUTPUT:
%   Saves: fig_path/temperature_hlevel-<h>_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_temp(xx, yy, temp, hlevel, m0, n_mode, fig_path)

    %% --------------------------------------------------------------------
    %% 1. Extract the data and compute limits
    %% --------------------------------------------------------------------

    % Get 2d data for our temperature contours at specific hlevel
    data = squeeze(temp(:,:,hlevel));

    % get the maximum and minimum values
    vmin = min(data(:));   % absolute minimum
    vmax = max(data(:));   % absolute maximum
    
    % print maximum and minimum values
    fprintf('\nTemperature contours at hlevel = %d:\n', hlevel)
    fprintf('Maximum Value: %.2f and Minimum Value: %.2f\n', vmax, vmin)
    
    % sets the +/- value to add to contourf (0.2 = ~20%)  
    step = 0.2;                              
    vmin = floor(vmin/step)*step;               
    vmax = ceil (vmax/step)*step;

    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------

    % Create the figure and set it's size [left, bottom, width, height]
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, yy, data', 'LineStyle', 'none');
    hold on
    contour(xx, yy, data', 'LineColor', 'k', 'LineStyle', '-');
    hold off

    colormap(flipud(cmap_twilight(256))); % set our colormap choice
    colorbar;
    caxis([vmin vmax]); % sets limits from section 1

    % Set x and y-label text as well as tick spacing for grid labels
    xlabel('Longitude (degrees)')
    ylabel('Latitude (degrees)')
    set(gca, 'xtick', 0:30:360)
    
    % Create title string with input variables and set it
    title_str = ['Temperature Perturbation', newline ...
        'hlevel = ', num2str(hlevel), ...
        ', zonal wave # = ', num2str(m0), ... 
        ', eMode # = ', num2str(n_mode),];

    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);   

    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------

    % Create filename to save our figure and set output directory
    outFile = fullfile(fig_path, ['temperature', ...
        '_hlevel-', num2str(hlevel), ...
        '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);
    
    fprintf('Saving temperature plot to: %s\n', outFile)

    % Comment out 'close(gcf)' if you want to see the plot and save it    
    saveas(gcf, outFile);
    close(gcf);
    
end