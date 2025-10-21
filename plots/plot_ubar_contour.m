%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_ubar_contour.m

% DESCRIPTION: Plots the background zonal wind (Ubar) contour with latitude on
% the x-axis and height on the y-axis, and saves it as an image.

% INPUT:
% - yy: Latitude coordinates (degrees)
% - zz: Height coordinates (km)
% - Ubar: Background zonal wind field (jj+1 x kk+1 array, m/s)

% OUTPUT:
% - Saves plot to 'output/plots/ubar_contour.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_ubar_contour(yy, zz, Ubar, model, m0, n_mode, fig_path)

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(yy, zz, Ubar', 'LineStyle', 'none');
    colorbar;
    xlabel('Latitude (degrees)')
    ylabel('Height (km)')
    title('Background Zonal Wind (Ubar)');

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);
    
    %% Save figure
    outFile = fullfile(fig_path, [model, '_ubar_', '_nmode-', num2str(n_mode), '_m0-', num2str(m0), '.png']);
    saveas(gcf, outFile);
    close(gcf);
    
end