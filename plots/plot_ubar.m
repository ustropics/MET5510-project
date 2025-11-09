%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_ubar.m

% DESCRIPTION: Plots background zonal wind Ubar(y,z) in the latitude-height
%              plane using filled contours. The figure is saved as a PNG
%              including wavenumber and mode number in filename.

% INPUT:
%   yy       - latitude coordinates (degrees), size (jj+1)
%   zz       - height coordinates (km), size (kk+1)
%   Ubar     - background zonal wind (jj+1 x kk+1 array, m/s)
%   m0       - zonal wavenumber for filename
%   n_mode   - mode number for filename
%   fig_path - directory path for saving figure

% OUTPUT:
%   Saves: fig_path/ubar_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_ubar(yy, zz, Ubar, m0, n_mode, fig_path)

    %% --------------------------------------------------------------------
    %% 1. Extract the data and compute limits
    %% --------------------------------------------------------------------

    % Get data to plot background zonal wind (Ubar)
    data = Ubar;

    % get the maximum and minimum values
    vmin = min(data(:));   % absolute minimum
    vmax = max(data(:));   % absolute maximum
    
    % print maximum and minimum values
    fprintf('\nBackground zonal wind (Ubar):\n')
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
    contourf(yy, zz, data', 'LineStyle', 'none');
    hold on
    contour(yy, zz, data', 'LineColor', 'k', 'LineStyle', '-');
    hold off

    colormap(flipud(cmap_spectral(256))); % set our colormap choice
    colorbar;
    caxis([vmin vmax]); % sets limits from section 1

    % Set x and y-label text as well as tick spacing for grid labels    
    xlabel('Latitude (degrees)')
    ylabel('Height (km)')

    % Create title string and set it
    title_str = ['Background Zonal Wind (Ubar)', newline ...
    'zonal wave # = ', num2str(m0), ...
    ', eMode # = ', num2str(n_mode)];
    title(title_str);

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);
    
    %% --------------------------------------------------------------------
    %% 3. Save figure
    %% --------------------------------------------------------------------

    % Create filename to save our figure and set output directory
    outFile = fullfile(fig_path, ['ubar', ...
        '_eMode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '.png']);

    fprintf('Saving background zonal wind (Ubar) to: %s\n', outFile)

    % Comment out 'close(gcf)' if you want to see the plot and save it    
    saveas(gcf, outFile);
    close(gcf);
    
end