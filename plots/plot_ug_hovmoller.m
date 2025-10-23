%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_ug_hovmoller.m

% DESCRIPTION: Plots the Hovmoller diagram for zonal wind using precomputed
% time-longitude data, and saves it as an image.

% INPUT:
% - xx: Longitude coordinates (degrees)
% - time: Time coordinates (days)
% - ug_hovmoler: Precomputed Hovmoller data for zonal wind (ii+1 x 51 array, m/s)
% - hlat: Latitude index (for title only)
% - hlevel: Vertical level index (for title only)

% OUTPUT:
% - Saves plot to 'output/figures/ug_hovmoller.png'

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
    title(['Hovmoller Diagram for Zonal Wind at lat=', num2str(hlat), ', level=', num2str(hlevel)]);

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);

    %% Save figure
    outFile = fullfile(fig_path, [model, '_ug_hovmoller_', '_nmode-', num2str(n_mode), '_m0-', num2str(m0), '.png']);
    saveas(gcf, outFile);
    close(gcf);
    
end