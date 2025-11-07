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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_vg_xsec(xx, zz, vg, jj, m0, n_mode, fig_path)
    
    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, zz, squeeze(vg(:,floor(jj/2)+1,:))', 'LineStyle', 'none');
    colorbar;
    xlabel('Longitude')
    ylabel('Height (km)')
    set(gca, 'xtick', 0:30:360)
    set(gca, 'ytick', 0:2:10)

    title_str = ['Meridional Wind Vertical Cross-Section (', ...
        'zonal wave # = ', num2str(m0), ', eMode # = ', num2str(n_mode)];
    title(title_str);

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);

    %% Save figure
    outFile = fullfile(fig_path, ['vg-cross-section', '_eMode-', num2str(n_mode), '_m0-', num2str(m0), '.png']);
    
    fprintf('Saving temperature plot to: %s\n', outFile)
    saveas(gcf, outFile);
    close(gcf);
    
end