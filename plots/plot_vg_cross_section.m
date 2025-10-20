%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: rossby_vg_cross_section.m

% Description: Plots the meridional wind vertical cross-section at the middle
% latitude (jj/2+1) across longitude and height, and saves it as an image.

% Input:
% - xx: Longitude coordinates (degrees)
% - zz: Height coordinates (km)
% - vg: Meridional wind field (ii+1 x jj+1 x kk+1 array, m/s)
% - jj: Number of latitude grid points

% Output:
% - Saves plot to 'output/plots/vg_cross_section.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_vg_cross_section(xx, zz, vg, jj, model, m0)
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    contourf(xx, zz, squeeze(vg(:,floor(jj/2)+1,:))', 'LineStyle', 'none');
    colorbar;
    xlabel('Longitude')
    ylabel('Height (km)')
    set(gca, 'xtick', 0:30:360)
    set(gca, 'ytick', 0:2:10)
    title('Meridional Wind Vertical Cross-Section');

    % Save plot
    saveas(gcf, ['output', filesep, 'figures', filesep, 'vg_cross_section_', model, '_m0_', num2str(m0), '.png']);
    close(gcf);
end