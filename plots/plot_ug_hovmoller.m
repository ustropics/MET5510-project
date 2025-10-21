%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Filename: rossby_ug_hovmoller.m

% Description: Plots the Hovmoller diagram for zonal wind using precomputed
% time-longitude data, and saves it as an image.

% Input:
% - xx: Longitude coordinates (degrees)
% - time: Time coordinates (days)
% - ug_hovmoler: Precomputed Hovmoller data for zonal wind (ii+1 x 51 array, m/s)
% - hlat: Latitude index (for title only)
% - hlevel: Vertical level index (for title only)

% Output:
% - Saves plot to 'output/plots/ug_hovmoller.png'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_ug_hovmoller(xx, time, ug_hovmoler, hlat, hlevel, model, m0)
figure('units', 'inch', 'position', [4,2,16,12], 'Visible', 'off')
contourf(xx, time, ug_hovmoler', 'LineStyle', 'none');
colorbar;
xlabel('Longitude')
ylabel('Time (days)')
title(['Hovmoller Diagram for Zonal Wind at lat=', num2str(hlat), ', level=', num2str(hlevel)]);

% set global font size
set(findall(gcf, '-property', 'FontSize'), 'FontSize',20);

% Save plot
saveas(gcf, ['output', filesep, 'figures', filesep, 'ug_hovmoller_', model, '_m0_', num2str(m0), '.png']);
close(gcf);
end