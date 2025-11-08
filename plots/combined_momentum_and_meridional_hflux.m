function combined_momentum_and_meridional_hflux(vg, ug, temp, m0, n_mode, fig_path)
% PLOT_MOMENTUM_AND_MERIDIONAL_HFLUX
% Overlays <v'T'> (meridional heat flux) as contours on <v'u'> (momentum flux) filled plot
%
% INPUTS:
%   vg     - 3D meridional wind perturbation (m/s)
%   ug     - 3D zonal wind perturbation (m/s)
%   temp   - 3D temperature perturbation (K)
%   m0     - Zonal wavenumber
%   n_mode - Eigenmode index
%   fig_path - Directory to save figure
%
% OUTPUT:
%   Saves: fig_path/combined_momentum_meridional_hflux_eMode-n_m0-m.png

    %% Compute zonal means
    zvu_calc = squeeze(mean(vg .* ug, 1));        % <v'u'>  (m²/s²)
    zvt_calc = squeeze(mean(vg .* temp, 1));      % <v'T'>  (m·K/s)

    %% Create figure
    figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
    hold on;

    % Filled contour: Momentum flux <v'u'>
    contourf(zvu_calc', 'LineStyle', 'none');
    colormap(cmap_twilight(256));
    c1 = colorbar;
    c1.Label.String = '<v''u''> (m^2/s^2)';
    c1.Label.FontSize = 20;

    % Overlay: Contour lines for meridional heat flux <v'T'>
    [C, h] = contour(zvt_calc', 'LineColor', 'k', 'LineWidth', 1.3);
    clabel(C, h, 'FontSize', 14, 'Color', 'k', 'LabelSpacing', 300);

    % Axis labels
    xlabel('Latitude', 'FontSize', 20);
    ylabel('z-level', 'FontSize', 20);

    % Title
    title_str = sprintf('Momentum Flux (filled) & Meridional Heat Flux (contours)\n(zonal wave # = %d, eMode # = %d)', ...
                        m0, n_mode);
    title(title_str, 'FontSize', 20);

    % Global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

    % Optional: Customize contour levels for <v'T'> based on data range
    % Example:
    % vT_range = max(abs(zvt_calc), [], 'all');
    % levels = linspace(-vT_range, vT_range, 13);
    % [C, h] = contour(zvt_calc', levels, 'LineColor', 'k', 'LineWidth', 1.3);

    %% Save figure
    outFile = fullfile(fig_path, ...
        sprintf('combined_momentum_meridional_hflux_eMode-%d_m0-%d.png', n_mode, m0));
    fprintf('Saving combined momentum + meridional heat flux plot to: %s\n', outFile);
    saveas(gcf, outFile);
    close(gcf);

end