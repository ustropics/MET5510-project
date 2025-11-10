%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FILE DESCRIPTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FILENAME: plot_zonal_hovmoller.m

% DESCRIPTION: Plots a Hovmoller diagram of zonal wind perturbation at a
%              specified vertical level and latitude index across longitude
%              and time.

% INPUT:
%   xx         - longitude coordinates (degrees)
%   time       - time coordinates (days)
%   ug_hovmoler- Hovmöller data for zonal wind (ii+1 x nt array, m/s)
%   hlat       - latitude index (for title)
%   hlevel     - vertical level index (for title and filename)
%   m0         - zonal wavenumber for title/filename
%   n_mode     - mode number for title/filename
%   fig_path   - directory path for saving figure

% OUTPUT:
%   Saves: fig_path/zonal-ug-hovmoller_hlevel-<h>_eMode-<n>_m0-<m>.png

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plot_zonal_hovmoller(xx, time, ug_hovmoler, hlat, hlevel, m0, n_mode, Lx, eigVal2, fig_path)
    %% --------------------------------------------------------------------
    %% 1. Extract the data and compute phase speed
    %% --------------------------------------------------------------------
    data = ug_hovmoler' * 10;  % scale by 10 for visibility (common in perturbations)

    % Compute phase speed
    omega = -imag(eigVal2(n_mode)); % angular frequency (rad/s)
    k = 2*pi*m0 / Lx; % wavenumber (rad/m)
    phase_speed = omega / k; % m/s
    phase_speed_day = Lx / (phase_speed * 86400); % days per wavelength

    % Get data limits
    vmin_data = min(data(:));
    vmax_data = max(data(:));
    fprintf('\nHovmoller diagram of zonal wind at hlevel = %d, hlat = %d:\n', hlevel, hlat);
    fprintf('Max: %.2f, Min: %.2f (scaled x10)\n', vmax_data, vmin_data);

    % Round to nearest 0.2 step
    step = 0.2;
    vmin = floor(vmin_data / step) * step;
    vmax = ceil(vmax_data / step) * step;

    %% --------------------------------------------------------------------
    %% 2. Create figure
    %% --------------------------------------------------------------------
    figure('units', 'inch', 'position', [4, 2, 18, 14], 'Visible', 'off');
    contourf(xx, time, data, 'LineStyle', 'none');
    hold on;
    contour(xx, time, data, 'LineColor', 'k', 'LineStyle', '-');

    colormap(cmap_PuOr(256));
    colorbar;
    caxis([vmin vmax]);

    % Labels
    xlabel('Longitude (degrees)');
    ylabel('Time (days)');
    set(gca, 'XTick', 0:60:360);

% Create title string with input variables and set it
    title_str = ['Zonal Wind Hovmoller Diagram', newline ...
        'latitude = ', num2str(hlat), ...
        ', hlevel = ', num2str(hlevel), ...
        ', zonal wave # = ', num2str(m0), ...
        ', eMode # = ', num2str(n_mode), ...
        ', phase speed = ', num2str(phase_speed_day, '%.2f'), ' days'];
    
    title(title_str);

    % Font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

    %% --------------------------------------------------------------------
    %% 3. Add HORIZONTAL LINE at phase_speed_day (on TIME axis)
    %% --------------------------------------------------------------------
    y_line = phase_speed_day;
    x_limits = xlim;  % [xmin, xmax] — NO parentheses!

    % Draw horizontal line across full longitude
    line(x_limits, [y_line y_line], ...
         'Color', 'k', 'LineStyle', '--', 'LineWidth', 4);

    % Add label on the right
    text(x_limits(2)*0.98, y_line, ...
         sprintf(' %.2f days', phase_speed_day), ...
         'Color', 'k', 'FontWeight', 'bold', ...
         'HorizontalAlignment', 'right', ...
         'VerticalAlignment', 'middle', 'FontSize', 18, ...
         'BackgroundColor', 'white', 'EdgeColor', 'k', 'Margin', 2);

    %% --------------------------------------------------------------------
    %% 4. Save figure
    %% --------------------------------------------------------------------
    outFile = fullfile(fig_path, ...
        sprintf('zonal-ug-hovmoller_hlevel-%d_eMode-%d_m0-%d.png', hlevel, n_mode, m0));
    fprintf('Saving Hovmoller diagram to: %s\n', outFile);
    saveas(gcf, outFile);
    close(gcf);
end