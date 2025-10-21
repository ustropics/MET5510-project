%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: rossby_meridional_wind.m

% Description: Plots the meridional wind contour at the surface (k=1) across
% longitude and latitude, and saves it as an image.

% Input:
% - xx: Longitude coordinates (degrees)
% - yy: Latitude coordinates (degrees)
% - vg: Meridional wind field (ii+1 x jj+1 x kk+1 array, m/s)

% Output:
% - Saves plot to 'output/plots/meridional_wind.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_meridional_wind(xx, yy, vg, model, m0)
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, yy, squeeze(vg(:,:,1))', 'LineStyle', 'none');
    colorbar;
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca, 'xtick', 0:30:360)
    title('Meridional Wind at Surface');

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);
    
    % Save plot
    saveas(gcf, ['output', filesep, 'figures', filesep, 'meridional_wind_', model, '_m0_', num2str(m0), '.png']);
    close(gcf);
end