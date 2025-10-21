%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_dpvdym_int.m

% DESCRIPTION: Plots the meridional gradient of background potential vorticity
% (d(PVbar)/dy) in the interior (scaled by 10^-12 s^-1) with latitude on the
% x-axis and height on the y-axis, and saves it as an image.

% INPUT:
% - yy: Latitude coordinates (degrees)
% - zz: Height coordinates (km)
% - BPVy: Meridional gradient of background PV (jj+1 x kk+1 array, s^-1)

% OUTPUT:
% - Saves plot to 'output/figures/dpvdym_int.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_dpvdym_int(yy, zz, BPVy, model, m0)

    %% Create figure
    figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
    dpvdym_scaled = BPVy * 1e12; % Scale to 10^-12 s^-1
    contourf(yy, zz, dpvdym_scaled', 'LineStyle', 'none');
    colorbar;
    xlabel('Latitude (degrees)')
    ylabel('Height (km)')
    title('d(PVbar)/dy (Interior, 10^-12 s^-1)');

    % set global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);

    %% Save figure
    saveas(gcf, ['output', filesep, 'figures', filesep, 'dpvdym_int_', model, '_m0_', num2str(m0), '.png']);
    close(gcf);

end