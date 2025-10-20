%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: rossby_ubar_contour.m

% Description: Plots the background zonal wind (Ubar) contour with latitude on
% the x-axis and height on the y-axis, and saves it as an image.

% Input:
% - yy: Latitude coordinates (degrees)
% - zz: Height coordinates (km)
% - Ubar: Background zonal wind field (jj+1 x kk+1 array, m/s)

% Output:
% - Saves plot to 'output/plots/ubar_contour.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_ubar_contour(yy, zz, Ubar, model, m0)
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(yy, zz, Ubar', 'LineStyle', 'none');
    colorbar;
    xlabel('Latitude (degrees)')
    ylabel('Height (km)')
    title('Background Zonal Wind (Ubar)');
    
    % Save plot
    saveas(gcf, ['output', filesep, 'figures', filesep, 'ubar_', model, '_m0_', num2str(m0), '.png']);
    close(gcf);
end