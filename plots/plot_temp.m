%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_temperature.m

% DESCRIPTION: Plots the temperature contour at a specified vertical level
% across longitude and latitude, and saves it as an image, including hlevel in title and filename.

% INPUT:
% - xx: Longitude coordinates (degrees)
% - yy: Latitude coordinates (degrees)
% - temp: Temperature field (ii+1 x jj+1 x kk+1 array, K)
% - hlevel: Vertical level index (e.g., 1, 25, 51) for title, filename, and data slice
% - model: Model name for filename
% - m0: Wavenumber for filename
% - n_mode: Mode number for filename
% - fig_path: Directory path for saving figure

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function plot_temp(xx, yy, temp, hlevel, m0, n_mode, fig_path)

    
    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, yy, squeeze(temp(:,:,hlevel))', 'LineStyle', 'none');
    clim([-0.8 0.8]); % sets limit to contourf
    colormap(flipud(cmap_twilight(256)));
    colorbar;
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca, 'xtick', 0:30:360)
    
    title_str = ['Temperature at (hlevel = ', num2str(hlevel), ...
        ', zonal wave # = ', num2str(m0), ', eMode # = ', num2str(n_mode), ')'];

    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);   

    %% Save figure
    outFile = fullfile(fig_path, ['temperature', '_hlevel-', num2str(hlevel), ...
        '_eMode-', num2str(n_mode), '_m0-', num2str(m0), '.png']);
    
    fprintf('Saving temperature plot to: %s\n', outFile)
    saveas(gcf, outFile);
    close(gcf);
    
end