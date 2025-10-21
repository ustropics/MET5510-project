%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_vg_cross_section.m

% DESCRIPTION: Plots the meridional wind vertical cross-section at the middle
% latitude (jj/2+1) across longitude and height, and saves it as an image.

% INPUT:
% - xx: Longitude coordinates (degrees)
% - zz: Height coordinates (km)
% - vg: Meridional wind field (ii+1 x jj+1 x kk+1 array, m/s)
% - jj: Number of latitude grid points

% OUTPUT:
% - Saves plot to 'output/figures/vg_cross_section.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_vg_cross_section(xx, zz, vg, jj, model, m0, n_mode, fig_path)
    
    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, zz, squeeze(vg(:,floor(jj/2)+1,:))', 'LineStyle', 'none');
    colorbar;
    xlabel('Longitude')
    ylabel('Height (km)')
    set(gca, 'xtick', 0:30:360)
    set(gca, 'ytick', 0:2:10)
    title('Meridional Wind Vertical Cross-Section');

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);

    %% Save figure
    outFile = fullfile(fig_path, [model, '_vg_cross_section_', '_nmode-', num2str(n_mode), '_m0-', num2str(m0), '.png']);
    saveas(gcf, outFile);
    close(gcf);
    
end