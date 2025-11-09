%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_gph_hovmoller.m

% DESCRIPTION: Plots a Hovmoller diagram of geopotential height perturbation
%              at a specified vertical level across longitude and time.

% INPUT:
%   xx       - longitude coordinates (degrees)
%   time     - time coordinates (days)
%   gpt_h_hovmoler - Hovmöller data matrix (ii+1 x nt)
%   m0       - zonal wavenumber for title/filename
%   n_mode   - mode number for title/filename
%   fig_path - directory path for saving figure
%   hlevel   - vertical level index

% OUTPUT:
%   Saves: fig_path/gph-hovmoller_hlevel-<h>_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_gph_hovmoller(xx, time, gpt_h_hovmoler, m0, n_mode, fig_path, hlevel)

    %% --------------------------------------------------------------------
    %% 1. Extract the data and compute limits
    %% --------------------------------------------------------------------

    % Plot data for Hovmoller diagram for geopotential heights
    data = gpt_h_hovmoler; % get our data to plot

    % get the maximum and minimum values
    vmin = min(data(:));   % absolute minimum
    vmax = max(data(:));   % absolute maximum
    
    % print maximum and minimum values
    fprintf('\nHovmoller diagram for geopotential heights:\n')
    fprintf('Maximum Value: %.2f and Minimum Value: %.2f\n', vmax, vmin)
    
    % sets the +/- value to add to contourf (0.2 = ~20%)  
    step = 0.1;                              
    vmin = floor(vmin/step)*step;               
    vmax = ceil (vmax/step)*step;

    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, time, data', 'linestyle', 'none');
    hold on
    contour(xx, time, data', 'LineColor', 'k', 'LineStyle', '-');
    hold off

    colormap(cmap_coolwarm(256)); % set our colormap choice
    colorbar;

    % caxis([vmin vmax]); % sets limits from section 1

    % Set x and y-label text as well as tick spacing for grid labels    
    xlabel('Longitude (degrees)') 
    ylabel('Time (days)')

    % Set x and y-label text as well as tick spacing for grid labels
    title_str = ['Geopotential Height Hovmoller Diagram', newline ...
        'hlevel = ', num2str(hlevel), ...
        ', zonal wave # = ', num2str(m0), ...
        ', eMode # = ', num2str(n_mode)];

    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);
    
    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------

    % Create filename to save our figure and set output directory
    outFile = fullfile(fig_path, ['gph-hovmoller', ...
        '_hlevel-', num2str(hlevel), ... 
        '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);

    fprintf('Saving hovmoller for geopotential height plot to: %s\n', outFile)

    % Comment out 'close(gcf)' if you want to see the plot and save it
    saveas(gcf, outFile);
    close(gcf);

end