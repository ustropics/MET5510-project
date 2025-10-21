%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: rossby_temperature.m

% Description: Plots the temperature contour at the middle vertical level
% (k=kk/2+1) across longitude and latitude, and saves it as an image.

% Input:
% - xx: Longitude coordinates (degrees)
% - yy: Latitude coordinates (degrees)
% - temp: Temperature field (ii+1 x jj+1 x kk+1 array, K)
% - kk: Number of height grid points

% Output:
% - Saves plot to 'output/plots/temperature.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_temperature(xx, yy, temp, kk, model, m0)
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, yy, squeeze(temp(:,:,floor(kk/2)+1))', 'LineStyle', 'none');
    colorbar;
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca, 'xtick', 0:30:360)
    title('Temperature at Mid-Level');

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);   

    % Save plot
    saveas(gcf, ['output', filesep, 'figures', filesep, 'temperature_', model, '_m0_', num2str(m0), '.png']);
    close(gcf);
end
