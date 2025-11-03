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

% OUTPUT:
% - Saves plot to 'fig_path/model_temperature_nmode-n_mode_m0-m0_hlevel-hlevel.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_temperature(xx, yy, temp, hlevel, model, m0, n_mode, fig_path)

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, yy, squeeze(temp(:,:,hlevel))', 'LineStyle', 'none');
    colorbar;
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca, 'xtick', 0:30:360)


    
    title(['Temperature at hlevel = ', num2str(hlevel)]);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);   

    %% Save figure
    outFile = fullfile(fig_path, [model, '_temperature_', '_nmode-', num2str(n_mode), ...
        '_m0-', num2str(m0), '_hlevel-', num2str(hlevel), '.png']);
    fprintf('Saving temperature plot to: %s\n', outFile); % Debug output
    saveas(gcf, outFile);
    close(gcf);
    
end