function plot_background_flow(yy, zz, Ubar, qy_surf, qy_trop, model, m0)
% PLOT_BACKGROUND_FLOW Plots the Eady model's background flow: Ubar, interior
% d(PVbar)/dy = 0, and boundary d(PVbar)/dy with beta.
%
% Inputs:
%   yy - Latitude coordinates (degrees)
%   zz - Height coordinates (km)
%   Ubar - Mean zonal wind field
%   qy_surf - Surface PV gradient (s^-1 m^-1)
%   qy_trop - Tropopause PV gradient (s^-1 m^-1)
%   model - Model name (string)
%   m0 - Wave number
%
% Output:
%   Saves combined plot to 'output/plots/eady_background_flow.png'

% Clear any variable named 'title' to avoid shadowing
clear title;

% Create the combined figure
figure('units', 'inch', 'position', [4, 2, 16, 6], ...
       'Name', sprintf('%s Model Background Flow', model));

%% Subplot 1: Ubar (mean zonal wind)
subplot(1, 3, 1);
contourf(yy, zz, Ubar', 20, 'LineStyle', 'none');  % Transpose for latitude vs height
colorbar;
xlabel('Latitude (degrees)');
ylabel('Height (km)');
title('$U_{\mathrm{bar}}$', 'Interpreter', 'latex');
caxis([0 max(Ubar(:))]);

%% Subplot 2: Interior ∂q/∂y = 0
subplot(1, 3, 2);
rectangle('Position', [min(yy), min(zz), max(yy)-min(yy), max(zz)-min(zz)], ...
          'FaceColor', [0 0.8 0.6], 'EdgeColor', 'none');
axis([min(yy) max(yy) min(zz) max(zz)]);
colorbar;
xlabel('Latitude (degrees)');
ylabel('Height (km)');
title('$\frac{\partial \bar{q}}{\partial y} = 0$ (Interior)', 'Interpreter', 'latex');
caxis([-0.1 0.1]);

%% Subplot 3: Boundary d(PVbar)/dy with beta
subplot(1, 3, 3);
hold on;

% Define scaled boundary PV gradients (×10^11 for plotting)
qy_surf_plot = qy_surf * 1e11;
qy_trop_plot = qy_trop * 1e11;

% Plot PV gradients at boundaries
plot([min(yy) max(yy)], [qy_surf_plot qy_surf_plot], 'r-', 'LineWidth', 2);  % Surface
plot([min(yy) max(yy)], [qy_trop_plot qy_trop_plot], 'b-', 'LineWidth', 2);  % Tropopause
plot([min(yy) max(yy)], [0 0], 'k-', 'LineWidth', 2);                        % Beta = 0

% Axis labels and title
xlabel('Latitude (degrees)');
ylabel('$\frac{d(\overline{PV})}{dy} \times 10^{-11} \, (s^{-1} m^{-1})$', 'Interpreter', 'latex');
title('$\frac{d(\overline{PV})}{dy}$ at surf./trop./beta (red/blue/black)', 'Interpreter', 'latex');

% Add text annotations
text(mean(yy), qy_trop_plot * 1.15, ...
    '$\left(\frac{\partial q}{\partial y}\right)_{\mathrm{Trop.}} = \frac{f_0^2}{N^2 H} \left(\frac{\partial \bar{u}}{\partial z}\right)_{\mathrm{Trop.}}$', ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 10);

text(mean(yy), qy_surf_plot * 1.15, ...
    '$\left(\frac{\partial q}{\partial y}\right)_{\mathrm{Surf.}} = -\frac{f_0^2}{N^2 H} \left(\frac{\partial \bar{u}}{\partial z}\right)_{\mathrm{Surf.}}$', ...
    'Interpreter', 'latex', 'HorizontalAlignment', 'center', 'FontSize', 10);

% Set y-limits
ylim([qy_surf_plot*1.5, qy_trop_plot*1.5]);
yline(0, 'k-', 'LineWidth', 1);

hold off;

% Save the plot
if ~exist('output/plots', 'dir')
    mkdir('output/plots');
end
saveas(gcf, ['output', filesep, 'figures', filesep, 'eady_background_', model, '_m0_', num2str(m0), '.png']);
close(gcf);
end