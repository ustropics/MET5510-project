%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_gph.m

% DESCRIPTION: Plots the geopotential height contour at the middle latitude
% (jj/2+1) across longitude and height, and saves it as an image.

% INPUT:
% - xx: Longitude coordinates (degrees)
% - zz: Height coordinates (km)
% - gpt_h: Geopotential height field (ii+1 x jj+1 x kk+1 array, m)
% - jj: Number of latitude grid points

% OUTPUT:
% - Saves plot to 'output/figures/geopotential_height.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_gph(xx, zz, gpt_h, jj, model, m0, n_mode, fig_path)

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, zz, squeeze(gpt_h(:,floor(jj/2)+1,:))', 'LineStyle', 'none');
    colorbar;
    xlabel('Longitude')
    ylabel('Height (km)')
    set(gca, 'xtick', 0:30:360)
    set(gca, 'ytick', 0:2:10)
    title('Geopotential Height at Middle Latitude');
    
    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);

    %% Save figure
    outFile = fullfile(fig_path, [model, '_gph', '_nmode-', num2str(n_mode), '_m0-', num2str(m0), '.png']);
    saveas(gcf, outFile);
    close(gcf);

end