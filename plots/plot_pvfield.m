%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_pvfield.m

% DESCRIPTION: Plots the potential vorticity contour at the surface (k=1)
% across longitude and latitude, and saves it as an image.

% INPUT:
% - xx: Longitude coordinates (degrees)
% - yy: Latitude coordinates (degrees)
% - pvfield: Potential vorticity field (ii+1 x jj+1 x kk+1 array, s^-1)

% OUTPUT:
% - Saves plot to 'output/figures/pvfield.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_pvfield(xx, yy, pvfield, model, m0)

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, yy, squeeze(pvfield(:,:,1))', 'LineStyle', 'none');
    colorbar;
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca, 'xtick', 0:30:360)
    title('Potential Vorticity at Surface');

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);

    %% Save plot
    saveas(gcf, ['output', filesep, 'figures', filesep, 'pvfield_', model, '_m0_', num2str(m0), '.png']);
    close(gcf);
end