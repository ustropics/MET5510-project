%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_ug_hovmoller.m

% DESCRIPTION: Plots the Hovmoller diagram for zonal wind using precomputed
% time-longitude data, and saves it as an image, including hlevel in filename.

% INPUT:
% - xx: Longitude coordinates (degrees)
% - time: Time coordinates (days)
% - ug_hovmoler: Precomputed Hovmoller data for zonal wind (ii+1 x 51 array, m/s)
% - hlat: Latitude index (for title only)
% - hlevel: Vertical level index (for title and filename)
% - model: Model name for filename
% - m0: Wavenumber for filename
% - n_mode: Mode number for filename
% - fig_path: Directory path for saving figure

% OUTPUT:
% - Saves plot to 'fig_path/model_ug_hovmoller_nmode-n_mode_m0-m0_hlevel-hlevel.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_ug_hovmoller(xx, time, ug_hovmoler, hlat, hlevel, model, m0, n_mode, fig_path)

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, time, ug_hovmoler'*10, 'LineStyle', 'none');
    colorbar;
    xlabel('Longitude')
    ylabel('Time (days)')


    title_str = ['Hovmoller Diagram for Zonal Wind (hlevel = ', num2str(hlat), ...
        ', wave # = ', num2str(m0), ', eMode # = ', num2str(n_mode)];
    title(title_str);

    % Set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

    %% Save figure
    outFile = fullfile(fig_path, ['ug_hovmoller',  '_hlevel-', num2str(hlevel), ... 
        '_nmode-', num2str(n_mode), '_m0-', num2str(m0), '.png']);

    fprintf('Saving temperature plot to: %s\n', outFile)
    saveas(gcf, outFile);
    close(gcf);
    
end