%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_meridional_xsec.m

% DESCRIPTION: Plots meridional wind vertical cross-section at mid-latitude
%              (jj/2+1) across longitude and height using filled contours.

% INPUT:
%   xx       - longitude coordinates (degrees)
%   zz       - height coordinates (km)
%   vg       - meridional wind field (ii+1 x jj+1 x kk+1 array, m/s)
%   jj       - number of latitude grid points
%   m0       - zonal wavenumber for title/filename
%   n_mode   - mode number for title/filename
%   fig_path - directory path for saving figure

% OUTPUT:
%   Saves: fig_path/meridional-vg-cross-section_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_meridional_xsec(xx, zz, vg, jj, m0, n_mode, fig_path)
    
    %% --------------------------------------------------------------------
    %% 1. Extract the data and compute limits
    %% --------------------------------------------------------------------

    lat = floor(jj/2)+1;

    % Get 2d data for our meridional wind cross-section
    data = squeeze(vg(:,lat,:)); % get our data to plot

    % get the maximum and minimum values
    vmin_data = min(data(:));   % absolute minimum
    vmax_data = max(data(:));   % absolute maximum
    
    % print maximum and minimum values
    fprintf('\nMeridional wind cross-section:\n')
    fprintf('Maximum Value: %.2f and Minimum Value: %.2f\n', vmax_data, vmin_data)

    % sets the +/- value to add to contourf (0.2 = ~20%) 
    step = 0.05;                                 
    vmin = floor(vmin_data/step)*step;               
    vmax = ceil (vmax_data/step)*step;

    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------

    % Create the figure and set it's size [left, bottom, width, height]
    figure('units', 'inch', 'position', [4,2,18,14], 'Visible', 'off')
    contourf(xx, zz, data');
    hold on
    contour(xx, zz, data', 'LineColor', 'k', 'LineStyle', '-');
    hold off
    
    colormap(cmap_PRGn(256)); % set our colormap choice
    colorbar;
    caxis([vmin vmax]); % sets limits from section 1

    % Set x and y-label text as well as tick spacing for grid labels
    xlabel('Longitude (degrees)')
    ylabel('Height (km)')
    set(gca, 'xtick', 0:30:360)
    set(gca, 'ytick', 0:2:10)

    % Create title string with input variables and set it
    title_str = ['Meridional Wind Vertical Cross-Section', newline ...
        'latitude = ', num2str(lat), ...
        ', zonal wave # = ', num2str(m0), ...
        ', eMode # = ', num2str(n_mode), ...
        ', max val = ', num2str(vmax_data, '%.3f'), ' m/s',  ...
        ', min val = ', num2str(vmin_data, '%.3f'), ' m/s'];

    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);
    
    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------

    % Create filename to save our figure and set output directory
    outFile = fullfile(fig_path, ['meridional-vg_xsec', ...
        '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);
    
    fprintf('Saving meridional wind vertical cross-section plot to: %s\n', outFile)

    % Comment out 'clouse(gcf)' if you want to see the plot and save it
    saveas(gcf, outFile);
    close(gcf);
    
end