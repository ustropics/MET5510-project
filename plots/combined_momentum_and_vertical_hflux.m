function combined_momentum_and_vertical_hflux(vg, ug, wfield, temp, m0, n_mode, fig_path)
% PLOT_MOMENTUM_AND_VERTICAL_HFLUX
% Overlays <w'T'> contours on <v'u'> filled contour plot
%
% INPUTS:
%   vg, ug     - 3D perturbations for momentum flux <v'u'>
%   wfield     - 3D vertical velocity
%   temp       - 3D temperature perturbation
%   m0         - Zonal wavenumber
%   n_mode     - Eigenmode index
%   fig_path   - Directory to save figure
%
% OUTPUT:
%   Saves: fig_path/combined_flux_eMode-n_m0-m.png

    %% Compute zonal means
    zvu_calc = squeeze(mean(vg .* ug, 1));        % <v'u'>  (jj+1 × kk+1)
    zwt_calc = squeeze(mean(wfield .* temp, 1));  % <w'T'>  (jj+1 × kk+1)

    %% Create figure
    figure('units', 'inch', 'position', [4, 2, 16, 12], 'Visible', 'off');
    hold on;

    % Filled contour: Momentum flux <v'u'>
    contourf(zvu_calc', 'LineStyle', 'none');
    colormap(cmap_twilight(256));
    c1 = colorbar;
    c1.Label.String = '<v''u''> (m²/s²)';
    c1.Label.FontSize = 20;

    % Overlay: Contour lines for vertical heat flux <w'T'>
    [C, h] = contour(zwt_calc', 'LineColor', 'k', 'LineWidth', 1.2);
    clabel(C, h, 'FontSize', 14, 'Color', 'k', 'LabelSpacing', 300);

    % Axis labels
    xlabel('Latitude', 'FontSize', 20);
    ylabel('z-level', 'FontSize', 20);

    % Title
    title_str = sprintf('Momentum Flux (filled) & Vertical Heat Flux (contours)\n(zonal wave # = %d, eMode # = %d)', ...
                        m0, n_mode);
    title(title_str, 'FontSize', 20);

    % Global font size
    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 20);

    % Optional: adjust contour levels for <w'T'> if needed
    % [C, h] = contour(zwt_calc', linspace(min(zwt_calc,[],'all'), max(zwt_calc,[],'all'), 15), ...

    %% Save figure
    outFile = fullfile(fig_path, ...
        sprintf('combined_flux_eMode-%d_m0-%d.png', n_mode, m0));
    fprintf('Saving combined flux plot to: %s\n', outFile);
    saveas(gcf, outFile);
    close(gcf);

end