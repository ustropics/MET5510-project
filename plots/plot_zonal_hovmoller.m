%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_zonal_hovmoller.m

% DESCRIPTION: Plots a Hovmoller diagram of zonal wind perturbation at a
%              specified vertical level and latitude index across longitude
%              and time.

% INPUT:
%   xx         - longitude coordinates (degrees)
%   time       - time coordinates (days)
%   ug_hovmoler- Hovmöller data for zonal wind (ii+1 x nt array, m/s)
%   hlat       - latitude index (for title)
%   hlevel     - vertical level index (for title and filename)
%   m0         - zonal wavenumber for title/filename
%   n_mode     - mode number for title/filename
%   fig_path   - directory path for saving figure

% OUTPUT:
%   Saves: fig_path/zonal-ug-hovmoller_hlevel-<h>_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_zonal_hovmoller(xx, time, ug_hovmoler, hlat, hlevel, m0, n_mode, fig_path)

    %% --------------------------------------------------------------------
    %% 1. Extract the data and compute limits
    %% --------------------------------------------------------------------

    % Get data for hovmoller diagram of zonal wind at specified hlevel
    data = ug_hovmoler'*10;

    % get the maximum and minimum values
    vmin = min(data(:));   % absolute minimum
    vmax = max(data(:));   % absolute maximum
    
    % print maximum and minimum values
    fprintf('\nHovmoller diagram of zonal wind at specified hlevel:\n')
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
    contourf(xx, time, data, 'LineStyle', 'none');
    hold on
    contour(xx, time, data, 'LineColor', 'k', 'LineStyle', '-')

    
    colormap(cmap_PuOr(256)); % set our colormap choice
    colorbar;

    caxis([vmin vmax]); % sets limits from section 1

    % Set x and y-label text as well as tick spacing for grid labels
    xlabel('Longitude (degrees)')
    ylabel('Time (days)')

    % Create title string with input variables and set it
    title_str = ['Zonal Wind Hovmoller Diagram', newline ... 
        'hlevel = ', num2str(hlat), ...
        ', zonal wave # = ', num2str(m0), ...
        ', eMode # = ', num2str(n_mode)];
    
    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------

    % Create filename to save our figure and set output directory
    outFile = fullfile(fig_path, ['zonal-ug-hovmoller', ...
        '_hlevel-', num2str(hlevel), ... 
        '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);

    fprintf('Saving Hovmoller diagram for zonal wind plot to: %s\n', outFile)

    % Comment out 'close(gcf)' if you want to see the plot and save it
    saveas(gcf, outFile);
    close(gcf);
    
end