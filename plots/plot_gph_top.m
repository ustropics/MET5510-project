%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: rossby_gph_top.m

% Description: Plots the geopotential height contour at the top boundary
% (k=kk+1) across longitude and latitude, and saves it as an image.

% Input:
% - xx: Longitude coordinates (degrees)
% - yy: Latitude coordinates (degrees)
% - gpt_h: Geopotential height field (ii+1 x jj+1 x kk+1 array, m)

% Output:
% - Saves plot to 'output/plots/gph_top.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rossby_gph_top(xx, yy, gpt_h, model, m0)
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, yy, squeeze(gpt_h(:,:,end))', 'LineStyle', 'none');
    colorbar;
    xlabel('Longitude')
    ylabel('Latitude')
    set(gca, 'xtick', 0:30:360)
    title('Geopotential Height at Top Boundary');

    % Save plot
    saveas(gcf, ['output', filesep, 'figures', filesep, 'gph_top_', model, '_m0_', num2str(m0), '.png']);
    close(gcf);
end